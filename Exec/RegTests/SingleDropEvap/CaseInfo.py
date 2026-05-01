import os
import re


class Droplet:
    def __init__(self, T, dia, fuel_names, Y=None, vel=None, Reyn=None):
        self.T = T
        self.dia = dia
        self.fuel_names = fuel_names
        self.dep_fuel_names = None
        if Y is None:
            self.Y = [1.0, 0.0]
        else:
            self.Y = Y
        self.fixed = True
        if abs(sum(self.Y) - 1.0) > 1.0e-8:
            errorstatement = "Liquid mass fractions must sum to 1"
            raise ValueError(errorstatement)
        if vel is None:
            self.vel = 0.0
        else:
            self.vel = vel
        if Reyn is None:
            self.Reyn = 0.0
        else:
            self.Reyn = Reyn


class GasPhase:
    def __init__(self, T, P, vel=None):
        self.T = T
        self.P = P
        if vel is None:
            self.vel = 0.0
        else:
            self.vel = vel


class CaseInfo:
    def __init__(
        self,
        name,
        dname,
        droplet: Droplet,
        gas: GasPhase,
        LiqPropsType,
        xyunits,
        end_time=None,
        dt=1e-2,
        plot_per=0.1,
        domain=[1.0, 1.0, 1.0],
        cell_num=[32, 32, 32],
        reftype=None,
        PeleMP_PsatModel="Antoine",
        **kwargs,
    ):

        # Model specifics
        self.LiqPropsType = LiqPropsType
        if LiqPropsType.lower() == "gcm":
            self.PeleMP_PsatModel = None
        else:
            self.PeleMP_PsatModel = PeleMP_PsatModel
        self.droplet = droplet
        self.gas = gas
        self.time = end_time
        self.dt = dt
        self.plot_per = plot_per
        self.domain = domain
        self.cell_num = cell_num
        self.num_liq_spec = len(droplet.fuel_names)
        self.use_file_y0 = False

        # File paths, names, etc.
        FILE_PATH = os.path.dirname(os.path.abspath(__file__))
        self.name = name
        self.ref_name = name  # Name to use for reference data lookup
        self.dname = dname
        self.case_dir = f"{LiqPropsType.upper()}_{name}"
        if LiqPropsType.lower() == "mp":
            if PeleMP_PsatModel.lower() == "antoine":
                self.case_dir += "_Antoine"
            else:
                self.case_dir += "_CC"
        if "use_manifold" in kwargs.keys():
            if kwargs["use_manifold"]:
                self.case_dir += "_Manifold"
        self.case_path = os.path.join(FILE_PATH, self.case_dir)
        self.input_file = os.path.join(self.case_path, f"input_{name}.inp")
        self.input_spray = os.path.join(self.case_path, f"input_{name}_spray.inp")

        # If reference is experimental or computational results
        if reftype is None:
            self.reftype = "exp"
        else:
            self.reftype = "comp"

        # Unit conversions for plotting
        if xyunits[0] == "s/mm2":
            self.xconv = 1.0 / (self.droplet.dia * 1.0e3) ** 2
            self.xlabel = "$t/d_0^2$ [s/mm$^2$]"
        elif xyunits[0] == "ms":
            self.xconv = 1.0e3
            self.xlabel = "$t$ [ms]"
        elif xyunits[0] == "runge":
            # This is nu_gas / r_0**2 * 1E-2
            if self.droplet.T == 273.15:
                nu_gas = 1.346452e-05
            else:
                nu_gas = 1.469687e-05
            self.xconv = 10 ** (-2) * nu_gas / (self.droplet.dia / 2.0) ** 2
            self.xlabel = r"$(t \nu / r_0^2) \times 10^{-2}$"
        else:
            self.xconv = 1.0
            self.xlabel = "$t$ [s]"
        if "2" in xyunits[1]:
            self.yexp = 2.0
        else:
            self.yexp = 1.0
        if "dd0" in xyunits[1]:
            self.yconv = 1.0 / self.droplet.dia
            if self.yexp == 2:
                self.ylabel = "$d^2/d_0^2$"
            else:
                self.ylabel = "$d/d_0$"
        elif xyunits[1] == "r2_mm":
            self.yconv = 1.0e3 * 0.5
            self.ylabel = "$r^2$ [mm$^2$]"

        # Check domain parameters
        diff_dxdy = abs(
            self.cell_num[0] / self.domain[0] - self.cell_num[1] / self.domain[1]
        )
        diff_dxdz = abs(
            self.cell_num[0] / self.domain[0] - self.cell_num[2] / self.domain[2]
        )
        if (diff_dxdy > 0.0) or (diff_dxdz > 0.0):
            error = "Uniform grid spacing required"
            raise ValueError(error)

    def set_end_time(self, time):
        if self.time is None:
            a_time = time * 1.05
            self.time = 0.1 * round(a_time / 0.1)

        self.plot_int = round(self.plot_per / self.dt)


def SpecifyCase(case_name, LiqPropsType, PeleMP_PsatModel="Antoine", **kwargs):
    if case_name.lower() == "nomura":
        case = Nomura(LiqPropsType, PeleMP_PsatModel, **kwargs)
    elif case_name.lower() == "wonglin":
        case = WongLin(LiqPropsType, PeleMP_PsatModel, **kwargs)
    elif case_name.lower() == "daif":
        case = Daif(LiqPropsType, PeleMP_PsatModel, **kwargs)
    elif case_name.lower() == "rungehep":
        case = RungeHep(LiqPropsType, PeleMP_PsatModel, **kwargs)
    elif case_name.lower() == "rungedec":
        case = RungeDec(LiqPropsType, PeleMP_PsatModel, **kwargs)
    elif case_name.lower() == "rungemix":
        case = RungeMix(LiqPropsType, PeleMP_PsatModel, **kwargs)
    elif "rungejp8" in case_name.lower():
        if "-h" in case_name.lower():
            kwargs["hychem"] = True
        elif "-d" in case_name.lower():
            kwargs["detailed"] = True
        case = RungeJP8(LiqPropsType, PeleMP_PsatModel, **kwargs)
    elif "burger" in case_name.lower():
        if "-h" in case_name.lower():
            kwargs["hychem"] = True
        if "50bar" in case_name.lower():
            kwargs["50bar"] = True
        elif "10bar" in case_name.lower():
            kwargs["10bar"] = True
        elif "1bar" in case_name.lower():
            kwargs["1bar"] = True
        case = Burger(LiqPropsType, PeleMP_PsatModel, **kwargs)
    else:
        raise ValueError(f"Unknown case name: {case_name}")
    return case


def Nomura(LiqPropsType, PeleMP_PsatModel="Antoine", **kwargs):
    drop = Droplet(298.0, 7.0e-4, ["NC7H16", "NC10H22"], [1.0, 0.0])
    gas = GasPhase(471, 1.0e5, vel=0.0)
    case = CaseInfo(
        f"Nomura",
        "Nomura et al.",
        drop,
        gas,
        LiqPropsType,
        xyunits=["s/mm2", "dd02"],
        end_time=2.94,
        PeleMP_PsatModel=PeleMP_PsatModel,
        **kwargs,
    )
    return case


def WongLin(LiqPropsType, PeleMP_PsatModel="Antoine", **kwargs):
    end_time = 4
    drop = Droplet(315.0, 1.961e-3, ["NC7H16", "NC10H22"], [0.0, 1.0], Reyn=17)
    gas = GasPhase(1000.0, 1.01325e5)
    case = CaseInfo(
        "WongLin",
        "Wong & Lin",
        drop,
        gas,
        LiqPropsType,
        xyunits=["s", "dd0"],
        end_time=end_time,
        PeleMP_PsatModel=PeleMP_PsatModel,
        **kwargs,
    )
    return case


def Daif(LiqPropsType, PeleMP_PsatModel="Antoine", **kwargs):
    drop = Droplet(291.4, 1.334e-3, ["NC7H16", "NC10H22"], [0.7375, 0.2625])
    gas = GasPhase(348.0, 1.01325e5, vel=3.10)
    case = CaseInfo(
        "Daif",
        "Daif et al.",
        drop,
        gas,
        LiqPropsType,
        xyunits=["s", "r2_mm"],
        dt=2e-3,
        PeleMP_PsatModel=PeleMP_PsatModel,
        **kwargs,
    )
    return case


def RungeMix(LiqPropsType, PeleMP_PsatModel="Antoine", **kwargs):
    drop = Droplet(272, 5.94e-4, ["NC7H16", "NC10H22"], [0.5, 0.5])
    gas = GasPhase(272, 1.01325e5, vel=2.5)
    case = CaseInfo(
        "RungeMix",
        "Runge et al.",
        drop,
        gas,
        LiqPropsType,
        xyunits=["runge", "dd02"],
        dt=5e-3,
        plot_per=1,
        PeleMP_PsatModel=PeleMP_PsatModel,
        **kwargs,
    )
    return case


def RungeDec(LiqPropsType, PeleMP_PsatModel="Antoine", **kwargs):
    drop = Droplet(272, 5.88e-4, ["NC7H16", "NC10H22"], [0.0, 1.0])
    gas = GasPhase(272, 1.01325e5, vel=2.5)
    case = CaseInfo(
        "RungeDec",
        "Runge et al.",
        drop,
        gas,
        LiqPropsType,
        xyunits=["runge", "dd02"],
        dt=5e-3,
        plot_per=1,
        PeleMP_PsatModel=PeleMP_PsatModel,
        **kwargs,
    )
    return case


def RungeHep(LiqPropsType, PeleMP_PsatModel="Antoine", **kwargs):
    drop = Droplet(272, 5.7e-4, ["NC7H16", "NC10H22"], [1.0, 0.0])
    gas = GasPhase(272, 1.01325e5, vel=2.5)
    case = CaseInfo(
        "RungeHep",
        "Runge et al.",
        drop,
        gas,
        LiqPropsType,
        xyunits=["runge", "dd02"],
        dt=5e-3,
        plot_per=0.25,
        PeleMP_PsatModel=PeleMP_PsatModel,
        **kwargs,
    )
    return case


def RungeJP8(LiqPropsType, PeleMP_PsatModel="Antoine", **kwargs):
    name = "RungeJP8"
    if "hychem" in kwargs.keys():
        if kwargs["hychem"]:
            name += "_HyChem"
    elif "detailed" in kwargs.keys():
        if kwargs["detailed"]:
            name += "_Detailed"
    drop = Droplet(294.15, 6.36e-4, ["POSF10264"], [1.0])
    gas = GasPhase(294.15, 1.01325e5, vel=3.0)
    case = CaseInfo(
        name, 
        "Runge et al.",
        drop,
        gas,
        LiqPropsType,
        xyunits=["runge", "dd02"],
        dt=5e-3,
        plot_per=1,
        PeleMP_PsatModel=PeleMP_PsatModel,
        **kwargs,
    )
    case.use_file_y0 = True
    # Use same reference data for RungeJP8, RungeJP8-H, and RungeJP8-D
    case.ref_name = "RungeJP8"
    return case

def Burger(LiqPropsType, PeleMP_PsatModel="Antoine", **kwargs):
    name = "Burger1Bar"  # Default to 1 bar
    P = 1e5 
    dt = 1e-3
    end_time = 0.09
    if "50bar" in kwargs.keys():
        if kwargs["50bar"]:
            name = "Burger50Bar"
            P = 50.0 * 1e5
    elif "10bar" in kwargs.keys():
        if kwargs["10bar"]:
            name = "Burger10Bar"
            P = 10.0 * 1e5
    elif "1bar" in kwargs.keys():
        if kwargs["1bar"]:
            name = "Burger1Bar"
    
    # Add HyChem suffix if specified
    if "hychem" in kwargs.keys():
        if kwargs["hychem"]:
            name += "_HyChem"
    
    drop = Droplet(300, 1e-4, ["POSF10325"], [1.0])
    gas = GasPhase(800, P, vel=0.0)
    case = CaseInfo(
        name, 
        "Burger et al.",
        drop,
        gas,
        LiqPropsType,
        xyunits=["s", "dd02"],
        dt=dt,
        end_time=end_time,
        plot_per=0.001,
        cell_num=[128, 128, 128],
        PeleMP_PsatModel=PeleMP_PsatModel,
        **kwargs,
    )
    case.use_file_y0 = True
    # Use base name without _HyChem suffix for reference data lookup
    case.ref_name = name.replace("_HyChem", "")
    return case

def CreateInputFile(case):
    FILE_PATH = os.path.dirname(os.path.abspath(__file__))

    # Boundary conditions depend on particle movement
    fixed_parts = True
    if (case.gas.vel > 0.0) or (case.droplet.Reyn > 0.0):
        lo_bc = "Inflow Interior Interior"
        hi_bc = "Outflow Interior Interior"
        is_periodic = "0 1 1"
    else:
        lo_bc = "Outflow Outflow Outflow"
        hi_bc = "Outflow Outflow Outflow"
        is_periodic = "0 0 0"

    # Read general input file
    gen_input_file = os.path.join(FILE_PATH, case.gen_input_file)

    with open(gen_input_file, "r") as f:
        gen_lines = f.readlines()

    new_lines = []
    num_psat_lines = 0
    for line in gen_lines:
        # Domain definition
        if "geometry.is_periodic" in line:
            new_line = f"geometry.is_periodic = {is_periodic}\n"
        elif "geometry.prob_lo" in line:
            dom_lo = [0.0, 0.0, 0.0]
            new_line = (
                f"geometry.prob_lo = {dom_lo[0]:.1f} {dom_lo[1]:.1f} {dom_lo[2]:.1f}\n"
            )
        elif "geometry.prob_hi" in line:
            dom_hi = case.domain
            new_line = (
                f"geometry.prob_hi = {dom_hi[0]:.1f} {dom_hi[1]:.1f} {dom_hi[2]:.1f}\n"
            )

        # BC Flags
        elif "peleLM.lo_bc" in line:
            new_line = f"peleLM.lo_bc = {lo_bc}\n"
        elif "peleLM.hi_bc" in line:
            new_line = f"peleLM.hi_bc = {hi_bc}\n"

        # AMR Control
        elif "amr.n_cell" in line:
            n_cell = case.cell_num
            new_line = f"amr.n_cell = {n_cell[0]:d} {n_cell[1]:d} {n_cell[2]:d}\n"
        elif "amr.plot_per" in line:
            new_line = f"amr.plot_per = {case.plot_per:d}\n"

        # Problem
        elif "prob.P_mean" in line:
            new_line = f"prob.P_mean = {case.gas.P:.1f}\n"
        elif "prob.T0_gas" in line:
            new_line = f"prob.T0_gas = {case.gas.T:.1f}\n"
        elif "prob.vel_gas" in line:
            new_line = f"prob.vel_gas = {case.gas.vel:f}\n"
        elif "prob.part_temp" in line:
            new_line = f"prob.part_temp = {case.droplet.T:.1f}\n"
        elif "prob.part_dia" in line:
            new_line = f"prob.part_dia = {case.droplet.dia:g}\n"
        elif "prob.part_vel" in line:
            if case.droplet.vel > 0.0:
                new_line = f"prob.part_vel = {case.droplet.vel:f} 0.0 0.0\n"
            else:
                new_line = f"prob.part_vel = {case.droplet.vel:1f} 0.0 0.0\n"
        elif "prob.part_Re" in line:
            new_line = f"prob.part_Re = {case.droplet.Reyn:.1f}\n"

        # Time stepping
        elif "amr.stop_time" in line:
            new_line = f"amr.stop_time = {case.time:g}\n"
        elif "amr.fixed_dt" in line:
            new_line = f"amr.fixed_dt = {case.dt:g}\n"

        # IO Control
        elif "amr.plot_file" in line:
            new_line = f'amr.plot_file = "{case.case_dir}/plt"\n'
        elif "amr.plot_int" in line:
            new_line = f"amr.plot_int = {case.plot_int:d}\n"

        # Spray particle data
        elif "particles.write_ascii_files" in line:
            new_line = "particles.write_ascii_files = 1\n"
        elif "particles.fixed_parts" in line:
            new_line = f"particles.fixed_parts = {fixed_parts:d}\n"
        elif "particles.Y_0" in line:
                # particles.Y_0 will be in sprayProps{case.LiqPropsType}_*.inp file
                new_line = "\n"
        elif "particles.fuel_species" in line:
                # particles.fuel_species will be in sprayProps{case.LiqPropsType}_*.inp file
                new_line = "\n"
        elif "FILE" in line:
            new_line = f"FILE = {case.case_dir}/input_{case.name}_spray.inp\n"
        else:
            new_line = line
        new_lines.append(new_line)

    # Save to output file
    with open(case.input_file, "w") as f:
        f.writelines(new_lines)

    # Edit spray_input_file
    spray_input_file = os.path.join(FILE_PATH, case.spray_input_file)

    with open(spray_input_file, "r") as f:
        spray_lines = f.readlines()

    new_spray_lines = []
    for line in spray_lines:
        if "particles.Y_0" in line:
            if case.use_file_y0:
                new_line = line
            else:
                new_line = f"particles.Y_0 = "
                for y in case.droplet.Y:
                    new_line += f"{y:.2f} "
                new_line += "\n"
        elif "particles.fuel_species" in line:
            if ("jp8" in case.name.lower()) and ("hychem" not in case.name.lower()):
                new_line = line
                # Set fuel_names to list after "particles.fuel_species = "
                case.droplet.fuel_names = line.split("=")[1].strip().split()
            elif ("burger" in case.name.lower()):
                new_line = line
                # Set fuel_names to list after "particles.fuel_species = "
                case.droplet.fuel_names = line.split("=")[1].strip().split()
            else: 
                new_line = "particles.fuel_species = "
                for fuel in case.droplet.fuel_names:
                    new_line += f"{fuel} "
                new_line += "\n"
        elif "particles.dep_fuel_species" in line:
            case.droplet.dep_fuel_names = line.split("=")[1].strip().split()
            case.droplet.unique_dep_fuel_names = list(
                set(case.droplet.dep_fuel_names)
            )
            new_line = line
        elif "# Units" in line:
            new_line = line
            new_line += f"# Notes: Y_0 modified for {case.name} case\n"
        elif re.search(r"particles\S*_psat", line):
            if case.LiqPropsType.lower() == "mp":
                # Only edit for PeleMP case
                if case.PeleMP_PsatModel.lower() == "antoine":
                    new_line = line
                    num_psat_lines += 1
                else:
                    # Clausius-Clapeyron relation, ignore existing line
                    new_line = ""
        else:
            new_line = line
        new_spray_lines.append(new_line)
    
    # Check that Psat lines were found for PeleMP if needed
    if (case.LiqPropsType.lower() == "mp") and (
        case.PeleMP_PsatModel.lower() == "antoine"
    ):
        if num_psat_lines < case.num_liq_spec:
            error = f"Expected {case.num_liq_spec} particles.SP_psat lines, found {num_psat_lines}"
            raise ValueError(error)

    # Save to output file
    with open(case.input_spray, "w") as f:
        f.writelines(new_spray_lines)


def CreateManifoldFiles(case, cmlm_dir):
    if not os.path.exists(cmlm_dir):
        raise RuntimeError(f"CMLM installation not found at specified path: {cmlm_dir}")
    fuels = case.droplet.fuel_names
    if case.droplet.dep_fuel_names is not None:
        dep_fuels = case.droplet.unique_dep_fuel_names
    else:
        dep_fuels = fuels

    # first find latent heat for each fuel from input files
    ifile = case.input_spray
    delta_h_vap = [0.0] * len(dep_fuels)
    found = [False] * len(dep_fuels)
    with open(ifile, "r") as f:
        for line in f.readlines():
            for i, fuel in enumerate(dep_fuels):
                if line.startswith(f"particles.{fuel}_latent"):
                    delta_h_vap[i] = float(line.split("=")[1].split("#")[0])
                    found[i] = True
    for ifound, fuel in zip(found, dep_fuels):
        if not ifound:
            raise RuntimeError(f"Latent heat not found in input files for dep species: {fuel}")

    # full input data
    table_file = os.path.join(case.case_path, "table.ctb")
    table_metadata_file = os.path.join(case.case_path, "table_metadata.txt")
    input_data = {
        "phys": {
            "mechanism": "../../../Submodules/PelePhysics/Mechanisms/liquid_fuels_nonreacting/mechanism.yaml",
            "pressure": case.gas.P,
            "X_ox": "O2:1.0, N2:3.76",
            "T_ox": case.gas.T,
            "X_fuel": [f"{fuel}:1.0" for fuel in dep_fuels],
            "liq_temp_fuel": [case.droplet.T] * len(dep_fuels),
            "delta_h_vap": delta_h_vap,
            "T_min": 100.0,
        },
        "table": {
            "use_fmix": False,
            "grid": [50] * len(dep_fuels),
            "filename": table_file,
            "metadata_file": table_metadata_file,
        },
    }

    # write input data to toml file
    import toml

    table_toml_file = os.path.join(case.case_path, "table_generation.toml")
    with open(table_toml_file, "w") as tomlfile:
        toml.dump(input_data, tomlfile)

    # run table generation scripts
    command = f"python {os.path.join(cmlm_dir,'run_scripts/ctable/create_spray_table_nd.py')} {table_toml_file}"
    print("Running Table Generation Command:")
    print(command)
    error = os.system(command)
    if error:
        raise RuntimeError(f"Table generation failed with error code {error}")

    # append necessary manifold information to input files
    with open(case.input_file, "a") as f:
        f.write("\n\n")
        f.write("# --------------------- # \n")
        f.write("# MANIFOLD MODEL INPUTS # \n")
        f.write("# --------------------- # \n")
        f.write("\n")
        f.write("manifold.model = Table \n")
        f.write(f"manifold.table.filename = {table_file} \n")
        f.write(f"manifold.metadata_filename = {table_metadata_file} \n")
        f.write("manifold.compute_temperature = true \n")
        f.write("manifold.has_species_mw = true \n")
        f.write("manifold.v = 1 \n")
        if "jp8" in case.name.lower() and "hychem" not in case.name.lower():
            f.write(
                "particles.dep_manifold_species = "
                + " ".join([f"ZMIX0" for i in range(len(fuels))])
                + "\n"
            )
        else:
            f.write(
                "particles.dep_manifold_species = "
                + " ".join([f"ZMIX{i}" for i in range(len(fuels))])
                + "\n"
            )
        f.write("peleLM.use_wbar = 0 \n")
