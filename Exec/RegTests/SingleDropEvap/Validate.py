import os
from CaseInfo import *
from ExtractData import *
import matplotlib.pyplot as plt

"""
Script for validating PelePhysics spray model
Test cases:
| Case         | Fuel           |
| ------------ | -------------- |
| Nomura(471)  | heptane        |
| Nomura(741)  | heptane        |
| WongLin()    | decane         |
| Daif()       | heptane/decane |
| RungeHep()   | heptane        |
| RungeDec()   | decane         |
| RungeMix()   | heptane/decane |
| RungeJP8()   | POSF10264      |
"""

# Case object
case = WongLin()

# Run new or extract existing simulation data?
run_new = True

# Number of processors to run on
num_proc = 6

# Plotting parameters
marker_s = 40
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
    elif (num_proc > 1) and ("MPI" not in exe):
        error = f"Pele not compiled with MPI and num_proc = {num_proc}"
        raise ValueError(error)

    # Run the case
    if ("MPI" in exe) and (num_proc > 1):
        os.system(f"mpiexec -np {num_proc} ./{exe} {case.input_file}")
    else:
        os.system(f"./{exe} {case.input_file}")

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
                s=marker_s,
                facecolor="none",
                edgecolor="black",
                label=case.dname,
                linewidth=round(line_w / 2),
            )
            # Plot uncertainty if available
            if refarr.shape[1] == 4:
                uncrt = refarr[~np.isnan(refarr).any(axis=1)]
                plt.scatter(
                    uncrt[:, 0],
                    uncrt[:, 2],
                    marker="_",
                    color="black",
                    label=None,
                    linewidth=round(line_w / 2),
                )
                plt.scatter(
                    uncrt[:, 0],
                    uncrt[:, 3],
                    marker="_",
                    color="black",
                    label=None,
                    linewidth=round(line_w / 2),
                )
                for k in range(len(uncrt)):
                    tval = [uncrt[k, 0], uncrt[k, 0]]
                    uline = [uncrt[k, 2], uncrt[k, 3]]
                    plt.plot(tval, uline, "k-", linewidth=round(line_w / 2))
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
                    s=marker_s,
                    facecolor="none",
                    edgecolor="black",
                    label=case.dname,
                    linewidth=round(line_w / 2),
                )
                # Plot uncertainty if available
                if refarr.shape[1] == 4:
                    uncrt = refarr[~np.isnan(refarr).any(axis=1)]
                    axs[i].scatter(
                        uncrt[:, 0],
                        uncrt[:, 2],
                        marker="_",
                        color="black",
                        label=None,
                        linewidth=round(line_w / 2),
                    )
                    axs[i].scatter(
                        uncrt[:, 0],
                        uncrt[:, 3],
                        marker="_",
                        color="black",
                        label=None,
                        linewidth=round(line_w / 2),
                    )
                    for k in range(len(uncrt)):
                        tval = [uncrt[k, 0], uncrt[k, 0]]
                        uline = [uncrt[k, 2], uncrt[k, 3]]
                        axs[i].plot(tval, uline, "k-", linewidth=round(line_w / 2))
            else:
                axs[i].plot(
                    refarr[:, 0],
                    refarr[:, 1],
                    label=case.dname,
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
