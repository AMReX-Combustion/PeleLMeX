#include <PeleLM.H>
#include <AMReX_ParmParse.H>

// -----------------------------------------------------------
// Read a binary file
// INPUTS/OUTPUTS:
// iname => filename
// nx    => input resolution
// ny    => input resolution
// nz    => input resolution
// data  <= output data
// -----------------------------------------------------------
void
read_binary(
  const std::string& iname,
  const size_t nx, 
  const size_t ny, 
  const size_t nz, 
  const size_t ncol,
  amrex::Vector<double>& data /*needs to be double*/)
{
  std::ifstream infile(iname, std::ios::in | std::ios::binary);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }

  for (size_t i = 0; i < nx * ny * nz * ncol; i++) {
    infile.read(reinterpret_cast<char*>(&data[i]), sizeof(data[i]));
  }
  infile.close();
}

AMREX_FORCE_INLINE
std::string
read_file(std::ifstream& in) 
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

// -----------------------------------------------------------
// Read a csv file
// INPUTS/OUTPUTS:
// iname => filename
// nx    => input resolution
// ny    => input resolution
// nz    => input resolution
// data  <= output data
// -----------------------------------------------------------
void
read_csv(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<amrex::Real>& data)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  std::istringstream iss(memfile);

  // Read the file
  size_t nlines = 0;
  std::string firstline;
  std::string line;
  std::getline(iss, firstline); // skip header
  while (getline(iss, line)) {
    ++nlines;
  }

  // Quick sanity check
  if (nlines != nx * ny * nz) {
    amrex::Abort(
      "Number of lines in the input file (= " + std::to_string(nlines) +
      ") does not match the input resolution (=" + std::to_string(nx) + ")");
  }

  // Read the data from the file
  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline); // skip header
  int cnt = 0;
  while (std::getline(iss, line)) {
    std::istringstream linestream(line);
    std::string value;
    while (getline(linestream, value, ',')) {
      std::istringstream sinput(value);
      sinput >> data[cnt];
      cnt++;
    }
  }
}


void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");

   // Chamber conditions
   pp.query("P_mean", PeleLM::prob_parm->P_mean);
   pp.query("T_mean", PeleLM::prob_parm->T_mean);
   pp.query("Y_CH4_chamber", PeleLM::prob_parm->Y_CH4_chamber);
   pp.query("Y_O2_chamber", PeleLM::prob_parm->Y_O2_chamber);

   // Injection parameters
   pp.query("nholes", PeleLM::prob_parm->nholes);
   pp.query("cone_angle", PeleLM::prob_parm->cone_angle);
   pp.query("centx", PeleLM::prob_parm->centx);
   pp.query("centy", PeleLM::prob_parm->centy);
   pp.query("r_circ", PeleLM::prob_parm->r_circ);
   pp.query("r_hole", PeleLM::prob_parm->r_hole);
   pp.query("T_jet", PeleLM::prob_parm->T_jet);
   pp.query("vel_jet", PeleLM::prob_parm->vel_jet);
   pp.query("injection_start", PeleLM::prob_parm->inj_start);
   pp.query("injection_duration", PeleLM::prob_parm->inj_dur);
   pp.query("tau", PeleLM::prob_parm->tau);
   pp.query("Z", PeleLM::prob_parm->Z);

   // Simplifications
   pp.query("useSymmetry", PeleLM::prob_parm->doQuarterDomain);

   // HIT init
   amrex::ParmParse ppic("ic");
   ppic.query("hitIC", PeleLM::prob_parm->hitIC);
   ppic.query("input_resolution", PeleLM::prob_parm->input_resolution);
   ppic.query("uin_norm", PeleLM::prob_parm->uin_norm);
   ppic.query("lscale", PeleLM::prob_parm->lscale);
   ppic.query("offset", PeleLM::prob_parm->offset);
   ppic.query("urms0", PeleLM::prob_parm->urms0);
   amrex::Vector<amrex::Real> win_lo(
     AMREX_SPACEDIM, std::numeric_limits<amrex::Real>::lowest());
   amrex::Vector<amrex::Real> win_hi(
     AMREX_SPACEDIM, std::numeric_limits<amrex::Real>::max());
   ppic.queryarr("win_lo", win_lo, 0, AMREX_SPACEDIM);
   ppic.queryarr("win_hi", win_hi, 0, AMREX_SPACEDIM);
   for (int i = 0; i < AMREX_SPACEDIM; i++) {
     PeleLM::prob_parm->win_lo[i] = win_lo[i];
     PeleLM::prob_parm->win_hi[i] = win_hi[i];
   }
   ppic.query("win_slope", PeleLM::prob_parm->win_slope);

   if (PeleLM::prob_parm->hitIC) {
      amrex::Print() << "Initializing HIT data \n";

      std::string datafile;
      ppic.query("input_name",datafile);
      int binfmt = 0;             // Default is ASCII format
      ppic.query("input_binaryformat",binfmt);

      // Read initial velocity field
      const size_t nx = PeleLM::prob_parm->input_resolution;
      const size_t ny = PeleLM::prob_parm->input_resolution;
      const size_t nz = PeleLM::prob_parm->input_resolution;
      amrex::Vector<double> data(nx * ny * nz * 6); /* this needs to be double */
      if (binfmt) {
        read_binary(datafile, nx, ny, nz, 6, data);
      } else {
        read_csv(datafile, nx, ny, nz, data);
      }

      // Extract position and velocities
      amrex::Vector<amrex::Real> xinput;
      amrex::Vector<amrex::Real> uinput;
      amrex::Vector<amrex::Real> vinput;
      amrex::Vector<amrex::Real> winput;
      amrex::Vector<amrex::Real> xdiff;
      amrex::Vector<amrex::Real> xarray;
      
      xinput.resize(nx * ny * nz);
      uinput.resize(nx * ny * nz);
      vinput.resize(nx * ny * nz);
      winput.resize(nx * ny * nz);
      
      for (int i = 0; i < xinput.size(); i++) {
        xinput[i] = (data[0 + i * 6] + PeleLM::prob_parm->offset) / 
                    PeleLM::prob_parm->lscale;
        uinput[i] = data[3 + i * 6] * PeleLM::prob_parm->urms0 / PeleLM::prob_parm->uin_norm;
        vinput[i] = data[4 + i * 6] * PeleLM::prob_parm->urms0 / PeleLM::prob_parm->uin_norm;
        winput[i] = data[5 + i * 6] * PeleLM::prob_parm->urms0 / PeleLM::prob_parm->uin_norm;
      }
      
      // Get the xarray table and the differences.
      xarray.resize(nx);
      for (int i = 0; i < xarray.size(); i++) {
        xarray[i] = xinput[i];
      }
      xdiff.resize(nx);
      std::adjacent_difference(
        xarray.begin(),
        xarray.end(),
        xdiff.begin());
      xdiff[0] = xdiff[1];
      
      // Make sure the search array is increasing
      if (not std::is_sorted(
            xarray.begin(),
            xarray.end())) {
        amrex::Abort("Error: non ascending x-coordinate array.");
      }
      PeleLM::prob_parm->Linput = xarray[nx - 1] + 0.5 * xdiff[nx - 1];
      
      // Initialize PeleLM::prob_parm containers
      PeleLM::prob_parm->d_xarray = (amrex::Real*) amrex::The_Arena()->alloc(nx*sizeof(amrex::Real));
      PeleLM::prob_parm->d_xdiff = (amrex::Real*) amrex::The_Arena()->alloc(nx*sizeof(amrex::Real));
      PeleLM::prob_parm->d_uinput = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real)); 
      PeleLM::prob_parm->d_vinput = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
      PeleLM::prob_parm->d_winput = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
      
      // Copy into PeleLM::prob_parm
      amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                xarray.begin(),
                                xarray.end(),
                       PeleLM::prob_parm->d_xarray);
      amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                xdiff.begin(),
                                xdiff.end(),
                       PeleLM::prob_parm->d_xdiff);
      amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                uinput.begin(),
                                uinput.end(),
                       PeleLM::prob_parm->d_uinput);
      amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                vinput.begin(),
                                vinput.end(),
                       PeleLM::prob_parm->d_vinput);
      amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                winput.begin(),
                                winput.end(),
                       PeleLM::prob_parm->d_winput);
   }
   
}
