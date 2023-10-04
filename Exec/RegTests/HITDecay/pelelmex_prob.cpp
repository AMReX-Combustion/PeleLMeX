#include <PeleLMeX.H>
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

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");

  ProbParm local_prob_parm;

  // Read HIT user-defined options from the input file
  pp.query("input_resolution", local_prob_parm.input_resolution);
  if (local_prob_parm.input_resolution == 0) {
    amrex::Abort("for HIT, the input resolution cannot be 0 !");
  }

  std::string datafile;
  pp.query("input_name", datafile);
  int binfmt = 0; // Default is ASCII format
  pp.query("input_binaryformat", binfmt);
  pp.query("urms0", local_prob_parm.urms0);

  amrex::Real lambda0 = 0.5;
  amrex::Real tau = lambda0 / local_prob_parm.urms0;

  // Output initial conditions
  std::ofstream ofs("initialConditions.txt", std::ofstream::out);
  amrex::Print(ofs) << "lambda0, urms0, tau \n";
  amrex::Print(ofs).SetPrecision(17)
    << lambda0 << "," << local_prob_parm.urms0 << "," << tau << std::endl;
  ofs.close();

  // Read initial velocity field
  const size_t nx = local_prob_parm.input_resolution;
  const size_t ny = local_prob_parm.input_resolution;
  const size_t nz = local_prob_parm.input_resolution;
  amrex::Vector<amrex::Real> data(
    nx * ny * nz * 6); /* this needs to be double */
  if (binfmt != 0) {
    read_binary(datafile, nx, ny, nz, 6, data);
  } else {
    read_csv(datafile, nx, ny, nz, data);
  }

  // Extract position and velocities
  amrex::Vector<amrex::Real> xinput(nx * ny * nz);
  amrex::Vector<amrex::Real> uinput(nx * ny * nz);
  amrex::Vector<amrex::Real> vinput(nx * ny * nz);
  amrex::Vector<amrex::Real> winput(nx * ny * nz);
  amrex::Vector<amrex::Real> xdiff(nx);
  amrex::Vector<amrex::Real> xarray(nx);

  for (long i = 0; i < xinput.size(); i++) {
    xinput[i] = data[0 + i * 6];
    uinput[i] =
      data[3 + i * 6] * local_prob_parm.urms0 / local_prob_parm.uin_norm;
    vinput[i] =
      data[4 + i * 6] * local_prob_parm.urms0 / local_prob_parm.uin_norm;
    winput[i] =
      data[5 + i * 6] * local_prob_parm.urms0 / local_prob_parm.uin_norm;
  }

  // Get the xarray table and the differences.
  for (long i = 0; i < xarray.size(); i++) {
    xarray[i] = xinput[i];
  }
  std::adjacent_difference(xarray.begin(), xarray.end(), xdiff.begin());
  xdiff[0] = xdiff[1];

  // Make sure the search array is increasing
  if (not std::is_sorted(xarray.begin(), xarray.end())) {
    amrex::Abort("Error: non ascending x-coordinate array.");
  }

  // Pass data to the local_prob_parm
  local_prob_parm.Linput = xarray[nx - 1] + 0.5 * xdiff[nx - 1];

  local_prob_parm.d_xarray =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
  local_prob_parm.d_xdiff =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
  local_prob_parm.d_uinput =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
  local_prob_parm.d_vinput =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
  local_prob_parm.d_winput =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));

  for (unsigned long i = 0; i < nx; i++) {
    local_prob_parm.d_xarray[i] = xarray[i];
    local_prob_parm.d_xdiff[i] = xdiff[i];
  }
  for (unsigned long i = 0; i < nx * ny * nz; i++) {
    local_prob_parm.d_uinput[i] = uinput[i];
    local_prob_parm.d_vinput[i] = vinput[i];
    local_prob_parm.d_winput[i] = winput[i];
  }

  // Initialize PeleLM::prob_parm container
  PeleLM::prob_parm->d_xarray =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
  PeleLM::prob_parm->d_xdiff =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
  PeleLM::prob_parm->d_uinput =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
  PeleLM::prob_parm->d_vinput =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));
  PeleLM::prob_parm->d_winput =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));

  // Copy into PeleLM::prob_parm: CPU only for now
  PeleLM::prob_parm->d_xarray =
    (amrex::Real*)amrex::The_Arena()->alloc(nx * sizeof(amrex::Real));
  std::memcpy(
    &PeleLM::prob_parm->d_xarray, &local_prob_parm.d_xarray,
    sizeof(local_prob_parm.d_xarray));
  std::memcpy(
    &PeleLM::prob_parm->d_xdiff, &local_prob_parm.d_xdiff,
    sizeof(local_prob_parm.d_xdiff));
  std::memcpy(
    &PeleLM::prob_parm->d_uinput, &local_prob_parm.d_uinput,
    sizeof(local_prob_parm.d_uinput));
  std::memcpy(
    &PeleLM::prob_parm->d_vinput, &local_prob_parm.d_vinput,
    sizeof(local_prob_parm.d_vinput));
  std::memcpy(
    &PeleLM::prob_parm->d_winput, &local_prob_parm.d_winput,
    sizeof(local_prob_parm.d_winput));
  PeleLM::prob_parm->Linput = local_prob_parm.Linput;
  PeleLM::prob_parm->input_resolution = local_prob_parm.input_resolution;
  PeleLM::prob_parm->urms0 = local_prob_parm.urms0;
  PeleLM::prob_parm->uin_norm = local_prob_parm.uin_norm;
}
