import os
import re


class Droplet:
    def __init__(self, T, dia, fuel_names, Y=None, vel=None, Reyn=None):
        self.T = T
        self.dia = dia
        self.fuel_names = fuel_names
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

        # File paths, names, etc.
        FILE_PATH = os.path.dirname(os.path.abspath(__file__))
        self.name = name
        self.dname = dname
        self.case_dir = f"{LiqPropsType.upper()}_{name}"
        if LiqPropsType.lower() == "mp":
            if PeleMP_PsatModel.lower() == "antoine":
                self.case_dir += "_Antoine"
            else:
                self.case_dir += "_CC"
        self.case_path = os.path.join(FILE_PATH, self.case_dir)
        self.input_file = os.path.join(self.case_path, f"input_{name}.inp")
        if LiqPropsType.lower() == "gcm":
            self.input_gcm = os.path.join(self.case_path, f"input_{name}_gcm.inp")

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


def SpecifyCase(case_name, LiqPropsType, PeleMP_PsatModel="Antoine"):
    if case_name.lower() == "nomura":
        case = Nomura(LiqPropsType, PeleMP_PsatModel)
    elif case_name.lower() == "wonglin":
        case = WongLin(LiqPropsType, PeleMP_PsatModel)
    elif case_name.lower() == "daif":
        case = Daif(LiqPropsType, PeleMP_PsatModel)
    elif case_name.lower() == "rungehep":
        case = RungeHep(LiqPropsType, PeleMP_PsatModel)
    elif case_name.lower() == "rungedec":
        case = RungeDec(LiqPropsType, PeleMP_PsatModel)
    elif case_name.lower() == "rungemix":
        case = RungeMix(LiqPropsType, PeleMP_PsatModel)
    elif case_name.lower() == "rungejp8":
        case = RungeJP8(LiqPropsType, PeleMP_PsatModel)
    else:
        raise ValueError(f"Unknown case name: {case_name}")
    return case


def Nomura(LiqPropsType, PeleMP_PsatModel="Antoine"):
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
    )
    return case


def WongLin(LiqPropsType, PeleMP_PsatModel="Antoine"):
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
    )
    return case


def Daif(LiqPropsType, PeleMP_PsatModel="Antoine"):
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
    )
    return case


def RungeMix(LiqPropsType, PeleMP_PsatModel="Antoine"):
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
    )
    return case


def RungeDec(LiqPropsType, PeleMP_PsatModel="Antoine"):
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
    )
    return case


def RungeHep(LiqPropsType, PeleMP_PsatModel="Antoine"):
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
    )
    return case


def RungeJP8(LiqPropsType, PeleMP_PsatModel="Antoine"):
    drop = Droplet(294.15, 6.36e-4, ["POSF10264"], [1.0])
    gas = GasPhase(294.15, 1.01325e5, vel=3.0)
    case = CaseInfo(
        "RungeJP8",
        "Runge et al.",
        drop,
        gas,
        LiqPropsType,
        xyunits=["runge", "dd02"],
        dt=2e-3,
        plot_per=1,
        PeleMP_PsatModel=PeleMP_PsatModel,
    )
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
            if case.LiqPropsType.lower() == "mp":
                # Only edit gen_input for PeleMP case
                new_line = "particles.Y_0 = "
                for y in case.droplet.Y:
                    new_line += f"{y:.2f} "
                new_line += "\n"
            else:
                # particles.Y_0 is in gcm_input_file
                new_line = "\n"
        elif "particles.fuel_species" in line:
            new_line = "particles.fuel_species = "
            for n in case.droplet.fuel_names:
                new_line += f"{n} "
            new_line += "\n"
        elif re.search(r"particles\S*_psat", line):
            if case.LiqPropsType.lower() == "mp":
                # Only edit gen_input for PeleMP case
                if case.PeleMP_PsatModel.lower() == "antoine":
                    new_line = line
                    num_psat_lines += 1
                else:
                    # Clausius-Clapeyron relation, ignore existing line
                    new_line = ""

        elif "FILE" in line:
            if case.LiqPropsType.lower() == "gcm":
                new_line = f"FILE = {case.case_dir}/input_{case.name}_gcm.inp\n"
            else:
                # Ignore existing FILE line for PeleMP case
                new_line = ""
        else:
            new_line = line
        new_lines.append(new_line)

    # Check that Psat lines were found for PeleMP if needed
    if (case.LiqPropsType.lower() == "mp") and (
        case.PeleMP_PsatModel.lower() == "antoine"
    ):
        if num_psat_lines < case.num_liq_spec:
            error = f"Expected {case.num_liq_spec} particles.SP_psat lines, found {num_psat_lines}"
            raise ValueError(error)

    # Save to output file
    with open(case.input_file, "w") as f:
        f.writelines(new_lines)

    # For GCM cases edit particles.Y_0 in gcm_input_file
    if case.LiqPropsType.lower() == "gcm":

        gcm_input_file = os.path.join(FILE_PATH, case.gcm_input_file)

        with open(gcm_input_file, "r") as f:
            gcm_lines = f.readlines()

        new_gcm_lines = []
        for line in gcm_lines:
            if "particles.Y_0" in line:
                # Edit gcm_input_file for GCM case
                new_line = f"particles.Y_0 = "
                for y in case.droplet.Y:
                    new_line += f"{y:.2f} "
                new_line += "\n"
            elif "# Units" in line:
                new_line = line
                new_line += f"# Notes: Y_0 modified for {case.name} case\n"
            else:
                new_line = line
            new_gcm_lines.append(new_line)

        # Save to output file
        with open(case.input_gcm, "w") as f:
            f.writelines(new_gcm_lines)
