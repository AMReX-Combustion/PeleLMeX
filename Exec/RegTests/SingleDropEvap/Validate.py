import os
from CaseInfo import *
from ExtractData import *
import matplotlib.pyplot as plt

"""
Script for validating PelePhysics spray model
Test cases:
| Case Name  | Fuel           | Requirements for SPRAY_FUEL_NUM                |
| ---------- | -------------- | ---------------------------------------------- |
| Nomura     | heptane        | SPRAY_FUEL_NUM = 2                             |
| WongLin    | decane         | SPRAY_FUEL_NUM = 2                             |
| Daif       | heptane/decane | SPRAY_FUEL_NUM = 2                             |
| RungeHep   | heptane        | SPRAY_FUEL_NUM = 2                             |
| RungeDec   | decane         | SPRAY_FUEL_NUM = 2                             |
| RungeMix   | heptane/decane | SPRAY_FUEL_NUM = 2                             |
| RungeJP8   | POSF10264      | SPRAY_FUEL_NUM = 1                             |
| ---------- | -------------- | ---------------------------------------------- |
"""
# Case to run
case_name = "WongLin"

# Liquid properties model: "mp" or "gcm"
LiqPropsType = "gcm"

# Psat model for PeleMP: "Antoine" or "Clasius-Clapeyron"
PeleMP_PsatModel = "Antoine"

# Run new or extract existing simulation data?
run_new = True

# Number of processors to run on
num_proc = 6

# Create case instance
case = SpecifyCase(case_name, LiqPropsType, PeleMP_PsatModel)

# General input file
case.gen_input_file = f"single-drop-evap-{LiqPropsType.lower()}.inp"
if "jp8" in case.name.lower():
    case.gcm_input_file = f"sprayPropsGCM_mixture_jp8.inp"
else:
    case.gcm_input_file = f"sprayPropsGCM_heptane-decane.inp"

# Get reference values from experiments
[refdvals, reftvals, refyvals] = ExtractRefVals(case)

# Set end time based on reference values
time = refdvals[-1, 0] / case.xconv
case.set_end_time(time)

if run_new:
    # Create a new directory for plt and spray files
    FILE_PATH = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(case.case_path):
        os.makedirs(case.case_path)

    # Remove existing plt and .p3d files
    else:
        os.system(
            f"rm -rf {case.case_path}/plt* {case.case_path}/*.p3d {case.case_path}/pele_vals.csv"
        )

    # Create case-specific input file
    CreateInputFile(case)

    # Check GNUmakefile for correct compilation flags
    with open(os.path.join(FILE_PATH, "GNUmakefile"), "r") as f:
        lines = f.readlines()
    gcm_flag = False
    for line in lines:
        if "SPRAY_GCM" in line:
            if "TRUE" in line:
                gcm_flag = True
        elif "SPRAY_FUEL_NUM" in line:
            if "jp8" in case.name.lower():
                if "1" not in line:
                    error = "GNUmakefile SPRAY_FUEL_NUM must be 1 for JP-8"
                    raise ValueError(error)
            else:
                if "2" not in line:
                    error = "GNUmakefile SPRAY_FUEL_NUM must be 2 for heptane/decane"
                    raise ValueError(error)
    if case.LiqPropsType.lower() == "gcm" and not gcm_flag:
        error = "GNUmakefile SPRAY_GCM must be TRUE for GCM liquid properties model"
        raise ValueError(error)
    elif case.LiqPropsType.lower() == "mp" and gcm_flag:
        error = "GNUmakefile SPRAY_GCM must be FALSE for MP liquid properties model"
        raise ValueError(error)

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
    if not os.path.exists(case.case_path):
        raise ValueError(f"Case directory not found: {case.case_path}")

# Extract Pele simulation data
outfile = os.path.join(case.case_path, "pele_vals.csv")
pele_vals = ExtractData(case, outfile)

# Plotting parameters
marker_s = 40
line_w = 3
font_s = 16
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
plt.savefig(os.path.join(case.case_path, "results.png"))
plt.show()
