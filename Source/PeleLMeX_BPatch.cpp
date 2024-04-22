#include "PeleLMeX_BPatch.H"
#include "PelePhysics.H"

BPatch::BPatch(const std::string& patch_name, const amrex::Geometry& geom)
  : m_patchname(std::move(patch_name))
{

  std::string ppbpatch = "bpatch." + m_patchname;
  amrex::ParmParse ps(ppbpatch);

  const auto dx = geom.CellSizeArray();
  auto prob_lo = geom.ProbLoArray();
  auto prob_hi = geom.ProbHiArray();

  // Reading patch related info
  ps.get("patchtype", m_patchtype);
  ps.get("boundary_direction", m_bpdata_h.m_boundary_dir);
  ps.get("boundary_lo_or_hi", m_bpdata_h.m_boundary_lo_hi);

  // Checking if patch type matches defined ones
#if (AMREX_SPACEDIM == 1)
  if (m_patchtype != "full-boundary") {
    amrex::Abort(
      "\nPatch type for 1-D domains should be specified as full-boundary.\n");
  }
#elif (AMREX_SPACEDIM == 2)
  if (m_patchtype != "line") {
    amrex::Abort("\nPatch type for 2-D domains should be specified as line.\n");
  }
#else
  if (
    m_patchtype != "full-boundary" && m_patchtype != "circle" &&
    m_patchtype != "rectangle" && m_patchtype != "circle-annular" &&
    m_patchtype != "rectangle-annular") {
    amrex::Abort(
      "\nPatch type do not match allowed values "
      "(full-boundary,circle,rectangle,circle-annular,rectangle-annular) \n");
  }
#endif

  if (
    m_bpdata_h.m_boundary_dir != 0 && m_bpdata_h.m_boundary_dir != 1 &&
    m_bpdata_h.m_boundary_dir != 2) {
    amrex::Abort("\nBoundary direction should be 0,1 or 2");
  }

  if (m_bpdata_h.m_boundary_lo_hi != 0 && m_bpdata_h.m_boundary_lo_hi != 1) {
    amrex::Abort("\nBoundary high low should be 0 or 1");
  }

#if (AMREX_SPACEDIM == 2)
  if (m_bpdata_h.m_boundary_dir == 2) {
    amrex::Abort("\nFor 2D problems, patch boundary direction cannot be 2");
  }
#endif

  // Define patch_num
  if (m_patchtype == "full-boundary") {
    m_bpdata_h.m_patchtype_num = 0;
  } else if (m_patchtype == "circle") {
    m_bpdata_h.m_patchtype_num = 1;
  } else if (m_patchtype == "rectangle") {
    m_bpdata_h.m_patchtype_num = 2;
  } else if (m_patchtype == "circle-annular") {
    m_bpdata_h.m_patchtype_num = 3;
  } else if (m_patchtype == "rectangle-annular") {
    m_bpdata_h.m_patchtype_num = 4;
  } else {
    amrex::Abort("\nError! Unknown patch type");
  }

  // Define patch variables

  if (m_patchtype == "full-boundary") {
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
      m_bpdata_h.m_patch_rectangle_lo[n] = prob_lo[n];
      m_bpdata_h.m_patch_rectangle_hi[n] = prob_hi[n];
    }
    m_bpdata_h.m_patch_rectangle_lo[m_bpdata_h.m_boundary_dir] = 0.0;
    m_bpdata_h.m_patch_rectangle_hi[m_bpdata_h.m_boundary_dir] = 0.0;
  } else if (m_patchtype == "circle") {
    for (int n = 0; n < 3; ++n) {
      ps.get("patch_circle_center", m_bpdata_h.m_patch_circle_center[n], n);
    }
    m_bpdata_h.m_patch_circle_center[m_bpdata_h.m_boundary_dir] = 0.0;
    ps.get("patch_circle_radius", m_bpdata_h.m_patch_circle_radius);
    if (m_bpdata_h.m_patch_circle_radius <= 0.0) {
      amrex::Abort("\nPatch radius for circle should be greater than 0");
    }
  } else if (m_patchtype == "rectangle") {
    for (int n = 0; n < 3; ++n) {
      ps.get("patch_rectangle_lo", m_bpdata_h.m_patch_rectangle_lo[n], n);
    }
    for (int n = 0; n < 3; ++n) {
      ps.get("patch_rectangle_hi", m_bpdata_h.m_patch_rectangle_hi[n], n);
    }
    m_bpdata_h.m_patch_rectangle_lo[m_bpdata_h.m_boundary_dir] = 0.0;
    m_bpdata_h.m_patch_rectangle_hi[m_bpdata_h.m_boundary_dir] = 0.0;
  } else if (m_patchtype == "circle-annular") {
    for (int n = 0; n < 3; ++n) {
      ps.get("patch_circ_ann_center", m_bpdata_h.m_patch_circ_ann_center[n], n);
    }
    m_bpdata_h.m_patch_circ_ann_center[m_bpdata_h.m_boundary_dir] = 0.0;
    ps.get(
      "patch_circ_ann_inner_radius", m_bpdata_h.m_patch_circ_ann_inner_radius);
    ps.get(
      "patch_circ_ann_outer_radius", m_bpdata_h.m_patch_circ_ann_outer_radius);
  } else if (m_patchtype == "rectangle-annular") {
    for (int n = 0; n < 3; ++n) {
      ps.get(
        "patch_rect_ann_outer_lo", m_bpdata_h.m_patch_rect_ann_outer_lo[n], n);
    }
    for (int n = 0; n < 3; ++n) {
      ps.get(
        "patch_rect_ann_outer_hi", m_bpdata_h.m_patch_rect_ann_outer_hi[n], n);
    }
    for (int n = 0; n < 3; ++n) {
      ps.get(
        "patch_rect_ann_inner_lo", m_bpdata_h.m_patch_rect_ann_inner_lo[n], n);
    }
    for (int n = 0; n < 3; ++n) {
      ps.get(
        "patch_rect_ann_inner_hi", m_bpdata_h.m_patch_rect_ann_inner_hi[n], n);
    }

    m_bpdata_h.m_patch_rect_ann_inner_hi[m_bpdata_h.m_boundary_dir] = 0.0;
    m_bpdata_h.m_patch_rect_ann_outer_hi[m_bpdata_h.m_boundary_dir] = 0.0;
    m_bpdata_h.m_patch_rect_ann_inner_lo[m_bpdata_h.m_boundary_dir] = 0.0;
    m_bpdata_h.m_patch_rect_ann_outer_lo[m_bpdata_h.m_boundary_dir] = 0.0;
  } else {
    amrex::Abort("\nError! Unknown patch type");
  }

  m_bpdata_h.num_species = ps.countval("species");
  ps.getarr("species", speciesList);

  if (m_bpdata_h.num_species > 0) {
    speciesList.resize(m_bpdata_h.num_species);
    m_bpdata_h.speciesIndex = (int*)amrex::The_Pinned_Arena()->alloc(
      m_bpdata_h.num_species * sizeof(int));
    m_bpdata_h.speciesFlux = (amrex::Real*)amrex::The_Pinned_Arena()->alloc(
      m_bpdata_h.num_species * sizeof(amrex::Real));
  } else {
    amrex::Abort("\nError! No species provided to plot flux at boundary patch");
  }

  amrex::Vector<std::string> names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(names);
  names.resize(names.size());

  for (int n = 0; n < names.size(); n++) {
    m_bpdata_h.speciesIndex[n] = -1;
  }

  for (int m = 0; m < m_bpdata_h.num_species; m++) {
    for (int n = 0; n < names.size(); n++) {
      if (speciesList[m] == names[n]) {
        m_bpdata_h.speciesIndex[m] = n;
      }
    }
  }

  for (int n = 0; n < m_bpdata_h.num_species; n++) {
    if (m_bpdata_h.speciesIndex[n] == -1) {
      std::string msg =
        "\nError! Unable to find species index " + std::to_string(n);
      amrex::Abort(msg);
    }
  }

  allocate();
  amrex::Gpu::streamSynchronize();
}
