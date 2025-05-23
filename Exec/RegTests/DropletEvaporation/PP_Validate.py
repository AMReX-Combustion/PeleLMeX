import os
from CaseInfo import *
from PP_ExtractData import *
import matplotlib.pyplot as plt

#case = Nomura(471)
case = WongLin(17)

# Plotting parameters
marker_s = 15
line_w = 3
font_s = 16

[refdvals, reftvals, refyvals] = ExtractRefVals(case)
# Set end time based on reference values
time = refdvals[-1, 0] / case.xconv
case.set_end_time(time)

# Create a new directory for plt and spray files
#os.system("mkdir -p {}".format(case.name))
#os.system("rm -r {}/plt*".format(case.name))
params = CreateInputParams(case)

executable = "None"
for f in os.listdir("./"):
    if (f.startswith("Pele") and f.endswith(".ex")):
        executable = f
if (not os.path.exists(executable)):
    error = "Pele executable not found"
    raise ValueError(error)

#input_file = "gen-input.inp"
#run_cmd = "mpiexec -np 6 ./"
#os.system("{}{} {} {}".format(run_cmd, executable, input_file, params))

outfile = case.name + "/pele_vals.csv"
pele_vals = ExtractData(case, outfile)

numplots = 1
if (reftvals is not None):
    numplots += 1
if (refyvals is not None):
    numplots += 1
ylabels = [case.ylabel, "$T$ [K]", "$Y_f$"]

if numplots == 1:
    plt.figure()
    plt.plot(pele_vals[:,0], pele_vals[:,1], label="Pele", color='red', linewidth = line_w)
    refarr = refdvals
    if (refarr is not None):
        if (case.reftype == "exp"):
            plt.scatter(refarr[:,0], refarr[:,1], marker='o', s=80, facecolor='none',
                        edgecolor='black',label="Ref",linewidth=round(line_w/2))
        else:
            plt.plot(refarr[:,0], refarr[:,1], label="Ref", color='black', linewidth = line_w)
    plt.ylabel(ylabels[0], fontsize=font_s)
    plt.xlabel(case.xlabel, fontsize=font_s)
    plt.xlim(min(pele_vals[:,0]),max(pele_vals[:,0]))
    plt.tick_params(labelsize=font_s)
    plt.legend(fontsize=font_s)
    plt.grid()
    
else:
    fig, axs = plt.subplots(1,numplots,figsize=(numplots*6.4,4.8))
    for i in range(numplots):
        axs[i].plot(pele_vals[:,0], pele_vals[:,i+1], label="Pele", color='red', linewidth = line_w)
        if (i == 0):
            refarr = refdvals
        elif (i == 1):
            refarr = reftvals
        elif (i == 2):
            refarr = refyvals
        if (refarr is not None):
            if (case.reftype == "exp"):
                axs[i].scatter(refarr[:,0], refarr[:,1], marker='o', s=80, facecolor='none',
                        edgecolor='black',label="Ref",linewidth=round(line_w/2))
            else:
                axs[i].plot(refarr[:,0], refarr[:,1], label="Ref", color='black', linewidth = line_w)
        axs[i].set_ylabel(ylabels[i], fontsize=font_s)
        axs[i].tick_params(labelsize=font_s)
        axs[i].set_xlim(min(pele_vals[:,0]),max(pele_vals[:,0]))
        axs[i].grid()
        axs[i].set_xlabel(case.xlabel, fontsize=font_s)
    plt.legend(fontsize=font_s)
    

plt.tight_layout()    
plt.savefig(case.name + "/results.png")
plt.show()
