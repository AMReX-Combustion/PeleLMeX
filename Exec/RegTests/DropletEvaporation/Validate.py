import os
from CaseInfo import *
from ExtractData import *
import matplotlib.pyplot as plt

# Select test case:
# Nomura at 471 K -> case = Nomura(471)
# Nomura at 741 K -> case = Nomura(741)
# Wong and Lin with Decane and Re=17 -> case = WongLin()

case = Nomura(471)

# Run new or extract existing simulation data?
run_new = True

# Plotting parameters
marker_s = 15
line_w = 3
font_s = 16

# Get reference values from experiments
[refdvals, reftvals, refyvals] = ExtractRefVals(case)

# Set end time based on reference values
time = refdvals[-1, 0] / case.xconv
case.set_end_time(time)

if run_new:
    # Create a new directory for plt and spray files
    FILE_PATH = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(case.case_dir):
        os.makedirs(case.case_dir)

    # Remove existing plt and .p3d files
    else:
        os.system(
            f"rm -rf {case.name}/plt* {case.name}/*.p3d {case.name}/pele_vals.csv"
        )

    # Create case-specific input file
    CreateInputFile(case)

    # Get the Pele executable
    exe = ""
    for f in os.listdir(FILE_PATH):
        if f.startswith("Pele") and f.endswith(".ex"):
            exe = f
    if not os.path.exists(exe):
        error = "Pele executable not found"
        raise ValueError(error)

    # Run the case
    os.system(f"mpiexec -np 4 ./{exe} {case.input_file}")
else:
    # Check that the case directory exists
    if not os.path.exists(case.case_dir):
        raise ValueError(f"Case directory not found: {case.case_dir}")

outfile = os.path.join(case.case_dir, "pele_vals.csv")
pele_vals = ExtractData(case, outfile)

numplots = 1
if reftvals is not None:
    numplots += 1
if refyvals is not None:
    numplots += 1
ylabels = [case.ylabel, "$T$ [K]", "$Y_f$"]

if numplots == 1:
    plt.figure()
    plt.plot(
        pele_vals[:, 0], pele_vals[:, 1], label="Pele", color="red", linewidth=line_w
    )
    refarr = refdvals
    if refarr is not None:
        if case.reftype == "exp":
            plt.scatter(
                refarr[:, 0],
                refarr[:, 1],
                marker="o",
                s=80,
                facecolor="none",
                edgecolor="black",
                label="Ref",
                linewidth=round(line_w / 2),
            )
        else:
            plt.plot(
                refarr[:, 0], refarr[:, 1], label="Ref", color="black", linewidth=line_w
            )
    plt.ylabel(ylabels[0], fontsize=font_s)
    plt.xlabel(case.xlabel, fontsize=font_s)
    plt.xlim(min(pele_vals[:, 0]), max(pele_vals[:, 0]))
    plt.tick_params(labelsize=font_s)
    plt.legend(fontsize=font_s)
    plt.grid()

else:
    fig, axs = plt.subplots(1, numplots, figsize=(numplots * 6.4, 4.8))
    for i in range(numplots):
        axs[i].plot(
            pele_vals[:, 0],
            pele_vals[:, i + 1],
            label="Pele",
            color="red",
            linewidth=line_w,
        )
        if i == 0:
            refarr = refdvals
        elif i == 1:
            refarr = reftvals
        elif i == 2:
            refarr = refyvals
        if refarr is not None:
            if case.reftype == "exp":
                axs[i].scatter(
                    refarr[:, 0],
                    refarr[:, 1],
                    marker="o",
                    s=80,
                    facecolor="none",
                    edgecolor="black",
                    label="Ref",
                    linewidth=round(line_w / 2),
                )
            else:
                axs[i].plot(
                    refarr[:, 0],
                    refarr[:, 1],
                    label="Ref",
                    color="black",
                    linewidth=line_w,
                )
        axs[i].set_ylabel(ylabels[i], fontsize=font_s)
        axs[i].tick_params(labelsize=font_s)
        axs[i].set_xlim(min(pele_vals[:, 0]), max(pele_vals[:, 0]))
        axs[i].grid()
        axs[i].set_xlabel(case.xlabel, fontsize=font_s)
    plt.legend(fontsize=font_s)


plt.tight_layout()
plt.savefig(os.path.join(case.case_dir, "results.png"))
plt.show()
