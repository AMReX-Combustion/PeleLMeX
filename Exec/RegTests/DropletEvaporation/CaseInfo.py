
import os

class Droplet:
    def __init__(self, T, dia, fuel_names, Y = None, vel = None, Reyn = None):
        self.T = T
        self.dia = dia
        self.fuel_names = fuel_names
        if (Y is None):
            self.Y = [1., 0.]
        else:
            self.Y = Y
        self.fixed = True
        if (abs(sum(self.Y) - 1.) > 1.E-8):
            errorstatement = "Liquid mass fractions must sum to 1"
            raise ValueError(errorstatement)
        if (vel is None):
            self.vel = 0.
        else:
            self.vel = vel
        if (Reyn is None):
            self.Reyn = 0.
        else:
            self.Reyn = Reyn

class GasPhase:
    def __init__(self, T, P, vel = None):
        self.T = T
        self.P = P
        if (vel is None):
            self.vel = 0.
        else:
            self.vel = vel
class CaseInfo:
    def __init__(self, name, droplet: Droplet, gas: GasPhase, xyunits, end_time=None, dt=1e-2, plot_per=0.1, domain=None, cell_num=None, reftype=None):
        FILE_PATH = os.path.dirname(os.path.abspath(__file__))
        self.name = name
        self.case_dir = os.path.join(FILE_PATH,name)
        self.input_file = os.path.join(self.case_dir,f"input_{name}.inp")
        self.droplet = droplet
        self.gas = gas
        self.time = end_time
        self.dt = dt
        self.plot_per = plot_per

        # If reference is experimental or computational results
        if (reftype is None):
            self.reftype = "exp"
        else:
            self.reftype = "comp"
        
        # Unit conversions for plotting
        if (xyunits[0] == "s/mm2"):
            self.xconv = 1. / (self.droplet.dia * 1.E3)**2
            self.xlabel = "$t/d_0^2$ [s/mm$^2$]"
        elif (xyunits[0] == "ms"):
            self.xconv = 1.E3
            self.xlabel = "$t$ [ms]"
        else:
            self.xconv = 1.
            self.xlabel = "$t$ [s]"
        if ("2" in xyunits[1]):
            self.yexp = 2.
        else:
            self.yexp = 1.
        if ("dd0" in xyunits[1]):
            self.yconv = 1. / self.droplet.dia
            if (self.yexp == 2):
                self.ylabel = "$d^2/d_0^2$"
            else:
                self.ylabel = "$d/d_0$"
        elif(xyunits[1] == "r2_mm"):
            self.yconv = 1.E3 * 0.5
            self.ylabel = "$r^2$ [mm$^2$]"

        # Set domain parameters
        if (domain is None):
            self.domain = [1., 1.]
        else:
            self.domain = domain
        if (cell_num is None):
            self.cell_num = [128, 128]
        else:
            self.cell_num = cell_num
        if (abs(self.cell_num[0] / self.domain[0] - self.cell_num[1] / self.domain[1]) > 0.):
            error = "Uniform grid spacing required"
            raise ValueError(error)
    def set_end_time(self, time):
        if self.time is None:
            a_time = time * 1.05
            self.time = 0.1 * round(a_time/0.1)
        
        self.plot_int = round(self.plot_per/self.dt)

def Nomura(temp):
    if temp == 471:
        end_time = 2.94
    elif temp == 741:
        end_time = 1.225
    else:
        end_time = None
    drop = Droplet(298., 7.0e-4, ["NC7H16","NC10H22"], [1.0,0.0])
    gas = GasPhase(temp, 1.0e5, vel = 0.0)
    case = CaseInfo(f"Nomura_{int(temp)}", drop, gas, xyunits=["s/mm2", "dd02"], end_time=end_time)
    return case

def WongLin():
    end_time = 4
    drop = Droplet(315., 1.961e-3, ["NC7H16","NC10H22"], [0.0,1.0], Reyn=17)
    gas = GasPhase(1000, 1.01325e5)
    case = CaseInfo("WongLin", drop, gas, xyunits=["s", "dd0"], end_time=end_time)
    return case

def CreateInputFile(case):
    FILE_PATH = os.path.dirname(os.path.abspath(__file__))

    # Boundary conditions depend on particle movement
    fixed_parts = 1
    if case.droplet.vel > 0:
        fixed_parts = 0
    if (fixed_parts):
        lo_bc = "Inflow Interior"
        hi_bc = "Outflow Interior"
        is_periodic = "0 1"
    else:
        lo_bc = "Interior Interior"
        hi_bc = "Interior Interior"
        is_periodic = "1 1"

    # Read general input file
    gen_file = os.path.join(FILE_PATH,"input_general.inp")
    with open(gen_file, 'r') as f:
        gen_lines = f.readlines()

    new_lines = []
    for line in gen_lines:
        # Domain definition
        if "geometry.is_periodic" in line:
            new_line = f"geometry.is_periodic = {is_periodic}\n"
        elif "geometry.prob_lo" in line:
            dom_lo = [0.0, 0.0]
            new_line = f"geometry.prob_lo = {dom_lo[0]:.1f} {dom_lo[1]:.1f}\n"
        elif "geometry.prob_hi" in line:
            dom_hi = [case.domain[0], case.domain[1]]
            new_line = f"geometry.prob_hi = {dom_hi[0]:.1f} {dom_hi[1]:.1f}\n"
        
        # BC Flags
        elif "peleLM.lo_bc" in line:
            new_line = f"peleLM.lo_bc = {lo_bc}\n"
        elif "peleLM.hi_bc" in line:
            new_line = f"peleLM.hi_bc = {hi_bc}\n"
        
        # AMR Control
        elif "amr.n_cell" in line:
            new_line = f"amr.n_cell = {case.cell_num[0]:d} {case.cell_num[1]:d}\n"
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
            new_line = f'amr.plot_file = "{case.name}/plt"\n'
        elif "amr.plot_int" in line:
            new_line = f"amr.plot_int = {case.plot_int:d}\n"

        # Spray particle data
        elif "particles.fixed_parts" in line:
            new_line = f"particles.fixed_parts = {fixed_parts:d}\n"
        elif "prob.Y_drop" in line:
            new_line = f"prob.Y_drop = "
            for y in case.droplet.Y:
                new_line += f"{y:.2f} "
            new_line += "\n"
        elif "particles.fuel_species" in line:
            new_line = f"particles.fuel_species = "
            for n in case.droplet.fuel_names:
                new_line += f"{n} "
            new_line += "\n"
        else:
            new_line = line
        new_lines.append(new_line)
        
    # Save to output file
    with open(case.input_file, 'w') as f:
        f.writelines(new_lines)

        
        
