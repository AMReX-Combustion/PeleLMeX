#ifndef _DIAGCOND_H_
#define _DIAGCOND_H_

#include "DiagBase.H"

class DiagConditional: public DiagBase::Register<DiagConditional>
{
public:
  static std::string identifier() { return "DiagConditional"; }

  enum conditionalType {Average, Integral, Sum};

  void init(const std::string &a_prefix,
            std::string_view a_diagName) override;

  void prepare(int a_nlevels,
               const amrex::Vector<amrex::Geometry> &a_geoms,
               const amrex::Vector<amrex::BoxArray> &a_grids,
               const amrex::Vector<amrex::DistributionMapping> &a_dmap,
               const amrex::Vector<std::string> &a_varNames) override;

  void processDiag(int a_nstep,
                   const amrex::Real &a_time,
                   const amrex::Vector<const amrex::MultiFab*> &a_state,
                   const amrex::Vector<std::string> &a_varNames) override;

  int getFieldIndex(const std::string &a_field,
                    const amrex::Vector<std::string> &a_varNames);
  amrex::Real MFVecMin(const amrex::Vector<const amrex::MultiFab*> &a_state,
                       int comp);  
  amrex::Real MFVecMax(const amrex::Vector<const amrex::MultiFab*> &a_state,
                       int comp);  
  void writeAverageDataToFile(const amrex::Vector<amrex::Real> &a_condAbs,
                              const amrex::Vector<amrex::Real> &a_cond,
                              const amrex::Vector<amrex::Real> &a_condSq);

  void close() override {}

private:
  // List of variables
  conditionalType m_condType;                // Type of conditional calculation
  std::string m_cFieldName{""};              // Condition Field name
  amrex::Vector<std::string> m_fieldNames;   // Processed Field names
  amrex::Gpu::DeviceVector<int> m_fieldIndices_d;
  int m_nBins{-1};                           // Number of bins
  bool m_useFieldMinMax{true};               // Use min/max from field
  amrex::Real m_lowBnd{0.0};                 // Low bound
  amrex::Real m_highBnd{0.0};                // High bound

  // Geometrical data
  amrex::Vector<amrex::Geometry> m_geoms;    // Squirrel away the geoms
  amrex::Vector<amrex::IntVect> m_refRatio;
};

#endif
