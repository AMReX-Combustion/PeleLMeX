import os
from CaseInfo import *
from ExtractData import *
import matplotlib.pyplot as plt
import argparse

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

parser = argparse.ArgumentParser(
    description="Run single droplet evaporation cases and compare to experimental data"
)

cases = ["WongLin", "Nomura", "Daif", "RungeHep", "RungeDec", "RungeMix", "RungeJP8"]
parser.add_argument(
    "--case_name",
    "-c",
    type=str,
    default="WongLin",
    choices=cases,
    help="Case name to run, default: WongLin",
)
prop_models = ["mp", "gcm"]
parser.add_argument(
    "--liq_props_type",
    "-l",
    type=str,
    default="gcm",
    choices=prop_models,
    help="Liquid properties model, default: gcm",
)
psat_models = ["Antoine", "Clausius-Clapeyron", "CC"]
parser.add_argument(
    "--mp_psat_model",
    "-p",
    type=str,
    default="Antoine",
    choices=psat_models,
    help="Psat model for PeleMP properties (Antoine, Clausius-Clapeyron, or CC), default: Antoine",
)
parser.add_argument(
    "--use_manifold",
    "-m",
    action="store_true",
    help="Use Manifold chemistry/EOS instead of Detailed chemistry/EOS",
)
parser.add_argument(
    "--cmlm_path",
    type=str,
    default="./cmlm",
    help="Path to CMLM install, required only for Manifold chemistry, default: ./cmlm",
)
parser.add_argument(
    "--dont_run_new",
    "-d",
    action="store_true",
    help="Plot previously computed data instead of running new simulation",
)
parser.add_argument(
    "--build_new", "-b", action="store_true", help="Build executable to run case"
)
parser.add_argument(
    "--num_proc",
    "-n",
    type=int,
    default=6,
    help="number of processors for parallel runs, default: 6",
)
args = parser.parse_args()


# Case to run
case_name = args.case_name

# Liquid properties model: "mp" or "gcm"
LiqPropsType = args.liq_props_type

# Psat model for PeleMP: "Antoine" or "Clausius-Clapeyron"
PeleMP_PsatModel = args.mp_psat_model
# Map CC shortcut to full name
if PeleMP_PsatModel == "CC":
    PeleMP_PsatModel = "Clausius-Clapeyron"

# Use manifold model for EOS (requires CMLM dependency)
use_manifold = args.use_manifold
cmlm_path = args.cmlm_path

# Run new or extract existing simulation data?
run_new = not args.dont_run_new

# Build new executable for case if needed
build_new = args.build_new

# Number of processors to run on
num_proc = args.num_proc

# Create case instance
case = SpecifyCase(case_name, LiqPropsType, PeleMP_PsatModel, use_manifold=use_manifold)

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
    if use_manifold:
        CreateManifoldFiles(case, cmlm_path)

    # Build the executable if needed
    if build_new:
        build_flags = f" -j {num_proc} "
        if case.LiqPropsType.lower() == "gcm":
            build_flags += " SPRAY_GCM=TRUE"
        elif case.LiqPropsType.lower() == "mp":
            build_flags += " SPRAY_GCM=FALSE"
        if "jp8" in case.name.lower():
            build_flags += " SPRAY_FUEL_NUM=1"
        else:
            build_flags += " SPRAY_FUEL_NUM=2"
        if use_manifold:
            build_flags += " USE_MANIFOLD=TRUE"
        else:
            build_flags += " USE_MANIFOLD=FALSE"
        error = os.system(f"make {build_flags}")
        if error:
            raise RuntimeError(f"Compilation of PeleLMeX failed with code {error}")

    # Get the Pele executable
    exe_files = [
        f for f in os.listdir(FILE_PATH) if f.startswith("Pele") and f.endswith(".ex")
    ]
    exe = None
    found = 0
    for f in exe_files:
        # We continue past invalid executables for our configuration
        if case.LiqPropsType.lower() == "gcm" and ".SprayGCM." not in f:
            continue
        if case.LiqPropsType.lower() == "mp" and ".SprayMP." not in f:
            continue
        if "jp8" in case.name.lower():
            if ".1SprayFuel." not in f:
                continue
        else:
            if ".2SprayFuel." not in f:
                continue
        if use_manifold and ".Manifold" not in f:
            continue
        elif not use_manifold and ".Manifold" in f:
            continue
        found += 1
        exe = f

    if found == 0:
        error = "Valid Pele executable for case not found"
        raise ValueError(error)
    elif found == 1:
        print(f"Running with executable: {exe}")
    else:
        print(f"Found {found} valid executables, using the last: {exe}")

    if (num_proc > 1) and ("MPI" not in exe):
        error = f"Pele not compiled with MPI and num_proc = {num_proc}"
        raise ValueError(error)

    # Run the case
    if ("MPI" in exe) and (num_proc > 1):
        error = os.system(f"mpiexec -np {num_proc} ./{exe} {case.input_file}")
    else:
        error = os.system(f"./{exe} {case.input_file}")
    if error:
        raise RuntimeError(f"Pele simulation failed with error code {error}")

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
