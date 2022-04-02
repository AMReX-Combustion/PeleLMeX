#include "DiagFramePlane.H"
#include "AMReX_VisMF.H"
#include <AMReX_FPC.H>
#include "AMReX_PlotFileUtil.H"
#include <regex>

#if __GNUC__ > 7
#include <filesystem>
namespace filesystem = std::filesystem;
#else
#include <experimental/filesystem>
namespace filesystem = std::experimental::filesystem;
#endif

void
printLowerDimIntVect(std::ostream &a_File,
                     const amrex::IntVect &a_IntVect,
                     int skipDim)
{
    int doneDim = 0;
    a_File << '(';
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim != skipDim) {
            a_File << a_IntVect[idim];
            doneDim ++;
            if (doneDim < AMREX_SPACEDIM-1) a_File << ",";
        }
    }
    a_File << ')';
}

void
printLowerDimBox(std::ostream &a_File,
                 const amrex::Box &a_box,
                 int skipDim)
{
    a_File << '(';
    printLowerDimIntVect(a_File, a_box.smallEnd(),skipDim);
    a_File << ' ';
    printLowerDimIntVect(a_File, a_box.bigEnd(),skipDim);
    a_File << ' ';
    printLowerDimIntVect(a_File, a_box.type(),skipDim);
    a_File << ')';
}

void
DiagFramePlane::init(const std::string &a_prefix)
{
    amrex::ParmParse pp(a_prefix);

    // Plane normal
    pp.get("normal", m_normal);
    AMREX_ASSERT(m_normal>=0 && m_normal<AMREX_SPACEDIM);

    // Plane center
    amrex::Vector<amrex::Real> center;
    pp.getarr("center",center,0,pp.countval("center"));
    if (center.size() == AMREX_SPACEDIM) {
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            m_center[idim] = center[idim];
        }
    } else if (center.size() == 1) {
        m_center[m_normal] = center[0];
    }
    
    // IO
    pp.query("int", m_interval);
    pp.query("per", m_per);
    m_diagfile = "DiagFramePlane";
    pp.query("file",m_diagfile);
    AMREX_ASSERT(m_interval>0 || m_per>0.0);

    // Interpolation
    std::string intType = "Quadratic";
    pp.query("interpolation", intType);
    if (intType == "Linear") {
        m_interpType = Linear;
    } else if (intType == "Quadratic") {
        m_interpType = Quadratic;
    } else {
        amrex::Abort("Unknown interpolation type for "+a_prefix);
    }
}

void 
DiagFramePlane::prepare(int a_nlevels,
                        const amrex::Vector<amrex::Geometry> &a_geoms,
                        const amrex::Vector<amrex::BoxArray> &a_grids,
                        const amrex::Vector<amrex::DistributionMapping> &a_dmap)
{
    if (first_time) {    
        // Store the level0 geometry
        auto initDomain  = a_geoms[0].Domain();
        auto initRealBox = a_geoms[0].ProbDomain();
        const amrex::Real* dxlcl = a_geoms[0].CellSize();
        initDomain.setRange(2,0,1);
        initRealBox.setLo(2,0.0);
        initRealBox.setHi(2,dxlcl[2]);
        m_geomLev0.define(initDomain,initRealBox,a_geoms[0].Coord(),amrex::Array<int,AMREX_SPACEDIM>({AMREX_D_DECL(0,0,0)}));
    }
    first_time = false;

    // Resize internal vectors
    m_intwgt.resize(a_nlevels);
    m_k0.resize(a_nlevels);

    // On each level, find the k0 where the plane lays
    // and the weight of the directionnal interpolation
    for (int lev = 0; lev < a_nlevels; lev++) {
        const amrex::Real* dx     = a_geoms[lev].CellSize();
        const amrex::Real* problo = a_geoms[lev].ProbLo();
        // How many dx away from the lowest cell-center ?
        amrex::Real dist = (m_center[m_normal] - (problo[m_normal]+0.5*dx[m_normal]))/dx[m_normal];
        int k0 = static_cast<int>(std::round(dist));
        dist -= static_cast<amrex::Real>(k0);
        m_k0[lev] = k0;
        if (m_interpType == Quadratic) {
            // Quadratic interp. weights on k0-1, k0, k0+1
            m_intwgt[lev][0] = 0.5 * (dist - 1.0) * (dist - 2.0);
            m_intwgt[lev][1] = dist * (2.0 - dist);
            m_intwgt[lev][2] = 0.5 * dist * (dist - 1.0);;
        } else if (m_interpType == Linear) {
            // linear interp. weights on k0-1, k0, k0+1
            if ( dist > 0.0 ) {
               m_intwgt[lev][0] = 0.0;
               m_intwgt[lev][1] = 1.0 - dist;
               m_intwgt[lev][2] = dist;
            } else if  ( dist < 0.0 ) {
               m_intwgt[lev][0] = -dist;
               m_intwgt[lev][1] = 1.0 + dist;
               m_intwgt[lev][2] = 0.0;
            } else {
               m_intwgt[lev][0] = 0.0;
               m_intwgt[lev][1] = 1.0;
               m_intwgt[lev][2] = 0.0;
            }
        }
    }

    // Assemble the 2D slice boxArray
    m_sliceBA.resize(a_nlevels);
    m_sliceDM.resize(a_nlevels);
    m_dmConvert.resize(a_nlevels);
    for (int lev = 0; lev < a_nlevels; lev++) {
        amrex::Vector<int> pmap;
        amrex::BoxList bl(a_grids[lev].ixType());
        bl.reserve(a_grids[lev].size());
        amrex::Vector<int> dmConvertLev;
        for (int i = 0; i < a_grids[lev].size(); ++i) {
            auto cBox = a_grids[lev][i];
            amrex::IntVect ploc(AMREX_D_DECL(cBox.smallEnd(0),cBox.smallEnd(1),cBox.smallEnd(2)));
            ploc[m_normal] = m_k0[lev];
            if (cBox.contains(ploc)) {
                amrex::Box zNormalBax;
                int idx = 0;
                for (int idim = 0; idim < AMREX_SPACEDIM;++idim) {
                    if (idim != m_normal) { 
                        zNormalBax.setRange(idx,cBox.smallEnd(idim),cBox.length(idim));
                        idx++;
                    }
                }
                zNormalBax.setRange(AMREX_SPACEDIM-1,0,1);
                bl.push_back(zNormalBax);
                pmap.push_back(a_dmap[lev][i]);
                dmConvertLev.push_back(i);
            }
            m_dmConvert[lev] = dmConvertLev;
        }
        m_sliceBA[lev].define(bl);
        m_sliceDM[lev].define(pmap);
    }
}

void
DiagFramePlane::processDiag(int a_nstep,
                            const amrex::Real &a_time,
                            const amrex::Vector<const amrex::MultiFab*> &a_state,
                            const amrex::Vector<std::string> &a_stateVar)
{
    // Interpolate data to slice
    amrex::Vector<amrex::MultiFab> planeData(a_state.size());
    for (int lev = 0; lev < a_state.size(); ++lev) {
        planeData[lev].define(m_sliceBA[lev], m_sliceDM[lev], a_state[0]->nComp(), 0);
        int p0 = m_k0[lev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(planeData[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto &bx = mfi.tilebox();
            const int state_idx = m_dmConvert[lev][mfi.index()];
            auto const& state = a_state[lev]->const_array(state_idx,0);
            auto const& plane = planeData[lev].array(mfi);
            if (m_normal == 0) {
                amrex::ParallelFor(bx, a_state[0]->nComp(), [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    plane(i,j,k,n) = m_intwgt[lev][0] * state(p0-1,i,j,n) +
                                     m_intwgt[lev][1] * state(p0  ,i,j,n) +
                                     m_intwgt[lev][2] * state(p0+1,i,j,n);
                });
            } else if (m_normal == 1) {
                amrex::ParallelFor(bx, a_state[0]->nComp(), [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    plane(i,j,k,n) = m_intwgt[lev][0] * state(i,p0-1,j,n) +
                                     m_intwgt[lev][1] * state(i,p0  ,j,n) +
                                     m_intwgt[lev][2] * state(i,p0+1,j,n);
                });
            } else if (m_normal == 2) {
                amrex::ParallelFor(bx, a_state[0]->nComp(), [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    plane(i,j,k,n) = m_intwgt[lev][0] * state(i,j,p0-1,n) +
                                     m_intwgt[lev][1] * state(i,j,p0  ,n) +
                                     m_intwgt[lev][2] * state(i,j,p0+1,n);
                });
            }
        }
    }

    // Count the number of level where the cut exists
    int nlevs = 0;
    for (int lev = 0; lev < a_state.size(); lev++) {
        if (!m_sliceBA[lev].empty()) nlevs +=1;
    }

    // Build up a z-normal 2D Geom
    amrex::Vector<amrex::Geometry> pltGeoms(nlevs);
    pltGeoms[0] = m_geomLev0;
    amrex::Vector<amrex::IntVect> ref_ratio;
    amrex::IntVect rref(AMREX_D_DECL(2,2,1));
    for (int lev = 1; lev < nlevs; ++lev) {
        pltGeoms[lev] = amrex::refine(pltGeoms[lev-1], rref);
        ref_ratio.push_back(rref);
    }

    // File name based on tep or time
    std::string diagfile;
    if (m_interval > 0) {
        diagfile = amrex::Concatenate(m_diagfile,a_nstep,6);
    }
    if (m_per > 0.0) {
        diagfile = m_diagfile+std::to_string(a_time);
    }
    amrex::Vector<int> step_array(nlevs,a_nstep);
    Write2DMultiLevelPlotfile(diagfile, nlevs, GetVecOfConstPtrs(planeData), a_stateVar,
                                   pltGeoms, a_time, step_array, ref_ratio);
    
}

void
DiagFramePlane::Write2DMultiLevelPlotfile(const std::string &a_pltfile,
                                          int a_nlevels,
                                          const amrex::Vector<const amrex::MultiFab*> &a_slice,
                                          const amrex::Vector<std::string> &a_varnames,
                                          const amrex::Vector<amrex::Geometry> &a_geoms,
                                          const amrex::Real &a_time,
                                          const amrex::Vector<int> &a_steps,
                                          const amrex::Vector<amrex::IntVect> &a_rref)
{
    const std::string levelPrefix = "Level_";
    const std::string mfPrefix = "Cell";
    const std::string versionName = "HyperCLaw-V1.1";

    bool callBarrier(false);
    amrex::PreBuildDirectorHierarchy(a_pltfile, levelPrefix, a_nlevels, callBarrier);
    amrex::ParallelDescriptor::Barrier();

    if (amrex::ParallelDescriptor::MyProc() == amrex::ParallelDescriptor::NProcs()-1) {
        amrex::Vector<amrex::BoxArray> boxArrays(a_nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = a_slice[level]->boxArray();
        }   

        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::string HeaderFileName(a_pltfile + "/Header");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |   
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good()) amrex::FileOpenFailed(HeaderFileName);
        Write2DPlotfileHeader(HeaderFile, a_nlevels, boxArrays, a_varnames,
                              a_geoms, a_time, a_steps, a_rref, versionName,
                              levelPrefix, mfPrefix);
    }
    for (int level = 0; level < a_nlevels; ++level) {
        VisMF2D(*a_slice[level], amrex::MultiFabFileFullPrefix(level, a_pltfile, levelPrefix, mfPrefix+"2D"));
        amrex::VisMF::Write(*a_slice[level], amrex::MultiFabFileFullPrefix(level, a_pltfile, levelPrefix, mfPrefix));
        ReWriteLevelVisMFHeader(amrex::MultiFabFileFullPrefix(level, a_pltfile, levelPrefix, ""));
    }
}

void
DiagFramePlane::Write2DPlotfileHeader(std::ostream &HeaderFile,
                                      int nlevels,
                                      const amrex::Vector<amrex::BoxArray> &bArray,
                                      const amrex::Vector<std::string> &varnames,
                                      const amrex::Vector<amrex::Geometry> &geom,
                                      const amrex::Real &time,
                                      const amrex::Vector<int> &level_steps,
                                      const amrex::Vector<amrex::IntVect> &ref_ratio,
                                      const std::string &versionName,
                                      const std::string &levelPrefix,
                                      const std::string &mfPrefix)
{
    int finest_level(nlevels - 1);
    HeaderFile.precision(17);

    int lowerSpaceDim = AMREX_SPACEDIM-1;

    HeaderFile << versionName << '\n';
    HeaderFile << varnames.size() << '\n';
    for (int ivar = 0; ivar < varnames.size(); ++ivar) {
        HeaderFile << varnames[ivar] << "\n";
    }
    HeaderFile << lowerSpaceDim << '\n';
    HeaderFile << time << '\n';
    HeaderFile << finest_level << '\n';
    for (int idim = 0; idim < lowerSpaceDim; ++idim) {
        HeaderFile << geom[0].ProbLo(idim) << ' ';
    }
    HeaderFile << '\n';
    for (int idim = 0; idim < lowerSpaceDim; ++idim) {
        HeaderFile << geom[0].ProbHi(idim) << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i < finest_level; ++i) {
        HeaderFile << ref_ratio[i][0] << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        printLowerDimBox(HeaderFile, geom[i].Domain(), 2);
        HeaderFile << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        HeaderFile << level_steps[i] << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        for (int idim = 0; idim < lowerSpaceDim; ++idim) {
            HeaderFile << geom[i].CellSize()[idim] << ' ';
        }
        HeaderFile << '\n';
    }
    HeaderFile << (int) geom[0].Coord() << '\n';
    HeaderFile << "0\n";

    for (int level = 0; level <= finest_level; ++level) {
        HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
        HeaderFile << level_steps[level] << '\n';

        const amrex::IntVect& domain_lo = geom[level].Domain().smallEnd();
        for (int i = 0; i < bArray[level].size(); ++i)
        {   
            // Need to shift because the RealBox ctor we call takes the
            // physical location of index (0,0,0).  This does not affect
            // the usual cases where the domain index starts with 0.
            const amrex::Box& b = amrex::shift(bArray[level][i], -domain_lo);
            amrex::RealBox loc = amrex::RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
            for (int idim = 0; idim < lowerSpaceDim; ++idim) {
                HeaderFile << loc.lo(idim) << ' ' << loc.hi(idim) << '\n';
            }   
        }   
        HeaderFile << amrex::MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
    }   
}

void
DiagFramePlane::ReWriteLevelVisMFHeader(const std::string &a_HeaderPath) {

    std::string OldHeaderFileName(a_HeaderPath+ "Cell_H");
    amrex::Vector<char> oldfileCharPtr;
    amrex::ParallelDescriptor::ReadAndBcastFile(OldHeaderFileName, oldfileCharPtr);
    
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::VisMF::IO_Buffer io_buffer_new(amrex::VisMF::IO_Buffer_Size);
        std::string HeaderFileName(a_HeaderPath+ "Cell_temp_H");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer_new.dataPtr(), io_buffer_new.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |   
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good()) amrex::FileOpenFailed(HeaderFileName);

        std::string fileCharPtrString(oldfileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);
        std::string line, word;
        
        // Version, How, nComp, nGrow
        int version, how, nComp, nGrow;
        is >> version;
        is.ignore(std::streamsize(10000), '\n');
        is >> how;
        is.ignore(std::streamsize(10000), '\n');
        is >> nComp;
        is.ignore(std::streamsize(10000), '\n');
        is >> nGrow;

        HeaderFile << version << "\n";
        HeaderFile << how << "\n";
        HeaderFile << nComp << "\n";
        HeaderFile << nGrow << "\n";

        // BoxArray
        int nbox = -1, dummy = 0;
        is.ignore(std::streamsize(10000), '(') >> nbox >> dummy;
        is.ignore(std::streamsize(10000), '\n');
        HeaderFile << '(' << nbox << ' ' << 0 << '\n';
        for (int i = 0; i < nbox; ++i) {
            amrex::Box fullDimBox;
            is >> fullDimBox;
            printLowerDimBox(HeaderFile, fullDimBox, 2);
            HeaderFile << '\n';
        }
        is.ignore(std::streamsize(10000), '\n');
        std::getline(is, line);
        HeaderFile << line << "\n";

        // Get list of Datafiles
        std::getline(is, line);
        HeaderFile << line << "\n";
        int nFabs = std::stoi(line);
        amrex::Vector<std::string> dataFiles;
        for (int i = 0; i < nFabs; ++i) {
            std::getline(is, line);
            HeaderFile << line << "\n";
            std::istringstream lis(line);
            std::string w1, w2, w3;
            lis >> w1;
            lis >> w2;
            lis >> w3;
            int offset = std::stoi(w3);
            if (offset==0) dataFiles.push_back(w2);
        }

        // Just pass from istream to ostream `till the end
        while (is) {
            std::getline(is, line);
            HeaderFile << line << "\n";
        }
        HeaderFile.close();

        // Replace header file
        filesystem::copy(HeaderFileName, OldHeaderFileName, filesystem::copy_options::overwrite_existing);
        filesystem::remove(HeaderFileName);

        // Replace 3D data file by 2D ones
        for (int i = 0; i < dataFiles.size(); ++i) {
            std::string newname = a_HeaderPath+dataFiles[i];
            newname = std::regex_replace(newname, std::regex("Cell_"), "Cell2D_");
            filesystem::rename(newname, a_HeaderPath+dataFiles[i]);
        }
    }
}

void
DiagFramePlane::VisMF2D(const amrex::MultiFab& a_mf,
                        const std::string& a_mf_name)
{
     auto whichRD = amrex::FArrayBox::getDataDescriptor();

     amrex::Long bytesWritten(0);

     std::string filePrefix(a_mf_name + "_D_");

     bool calcMinMax = false;
     amrex::VisMF::Header::Version currentVersion = amrex::VisMF::Header::Version_v1;
     amrex::VisMF::How how = amrex::VisMF::How::NFiles;
     amrex::VisMF::Header hdr(a_mf, how, currentVersion, calcMinMax);

     int nOutFiles = std::max(1, std::min(amrex::ParallelDescriptor::NProcs(), 256));
     bool groupSets = false;
     bool setBuf = true;

     amrex::NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf);

     // Check if mf has sparse data
     bool useSparseFPP = false;
     const amrex::Vector<int> &pmap = a_mf.DistributionMap().ProcessorMap();
     std::set<int> procsWithData;
     amrex::Vector<int> procsWithDataVector;
     for(int i(0); i < pmap.size(); ++i) {
         procsWithData.insert(pmap[i]);
     }    
     if(static_cast<int>(procsWithData.size()) < nOutFiles) {
         useSparseFPP = true;
         for(std::set<int>::iterator it = procsWithData.begin(); it != procsWithData.end(); ++it) {
             procsWithDataVector.push_back(*it);
         }    
     } 

     if ( useSparseFPP ) {
         nfi.SetSparseFPP(procsWithDataVector);
     } else {
         nfi.SetDynamic();
     }
     for( ; nfi.ReadyToWrite(); ++nfi) {
         int whichRDBytes(whichRD->numBytes()), nFABs(0);
         amrex::Long writeDataItems(0), writeDataSize(0);
         for(amrex::MFIter mfi(a_mf); mfi.isValid(); ++mfi) {
             const amrex::FArrayBox &fab = a_mf[mfi];
             std::stringstream hss;
             write_2D_header(hss, fab, fab.nComp());
             bytesWritten += static_cast<std::streamoff>(hss.tellp());
             bytesWritten += fab.box().numPts() * a_mf.nComp() * whichRDBytes;
             ++nFABs;
         }

         for(amrex::MFIter mfi(a_mf); mfi.isValid(); ++mfi) {
             int hLength(0);
             const amrex::FArrayBox &fab = a_mf[mfi];
             writeDataItems = fab.box().numPts() * a_mf.nComp();
             writeDataSize = writeDataItems * whichRDBytes;
             std::stringstream hss; 
             write_2D_header(hss, fab, fab.nComp());
             hLength = static_cast<std::streamoff>(hss.tellp());
             auto tstr = hss.str();
             nfi.Stream().write(tstr.c_str(), hLength);    // ---- the fab header
             nfi.Stream().flush();
             amrex::Real const* fabdata = fab.dataPtr();
#ifdef AMREX_USE_GPU
             std::unique_ptr<amrex::FArrayBox> hostfab;
             if (fab.arena()->isManaged() || fab.arena()->isDevice()) {
                 hostfab = std::make_unique<amrex::FArrayBox>(fab.box(), fab.nComp(),
                                                              The_Pinned_Arena());
                 amrex::Gpu::dtoh_memcpy_async(hostfab->dataPtr(), fab.dataPtr(),
                                               fab.size()*sizeof(Real));
                 amrex::Gpu::streamSynchronize();
                 fabdata = hostfab->dataPtr();
             }    
#endif
             nfi.Stream().write((char *) fabdata, writeDataSize);
             nfi.Stream().flush();
          }
     }
}

void
DiagFramePlane::write_2D_header(std::ostream&           os,
                                const amrex::FArrayBox& f,
                                int                     nvar)
{
    os << "FAB " << amrex::FPC::NativeRealDescriptor();
    amrex::StreamRetry sr(os, "FABio_write_header", 4);
    while(sr.TryOutput()) {
        printLowerDimBox(os, f.box(), 2);
        os << ' ' << nvar << '\n';
    }
}
