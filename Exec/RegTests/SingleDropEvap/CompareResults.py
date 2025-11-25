import os
import numpy as np
from CaseInfo import *
from ExtractData import *
import matplotlib.pyplot as plt
import pandas as pd
import argparse

"""
Script for plotting and comparing PelePhysics liquid properties models
Test cases:
| Case Name  | Fuel                 | Notes                                    |
| ---------- | -------------------- | ---------------------------------------- |
| Nomura     | heptane              |                                          |
| WongLin    | decane               |                                          |
| Daif       | heptane/decane       |                                          |
| Runge      | heptane, decane, mix | Plots RungeHep, RungeDec and RungeMix    |
| RungeJP8   | POSF10264            | Plots JP8 case only                      |
| ---------- | -------------------- | ---------------------------------------- |
"""

parser = argparse.ArgumentParser(
    description="Quantitatively compare results for given case to experimental data"
)

cases = ["WongLin", "Nomura", "Daif", "Runge", "RungeJP8"]
parser.add_argument(
    "--case_name",
    "-c",
    type=str,
    default="WongLing",
    choices=cases,
    help="Case name to run, default: WongLin",
)

args = parser.parse_args()

# Case to run
case_name = args.case_name

def find_case_dirs(case_name):
    # List all directories in file path
    FILE_PATH = os.path.dirname(os.path.abspath(__file__))
    dir_list = [
        d for d in os.listdir(FILE_PATH) if os.path.isdir(os.path.join(FILE_PATH, d))
    ]
    # Find directories that match the case name
    matching_dirs = [d for d in dir_list if case_name.lower() in d.lower()]
    if not matching_dirs:
        raise ValueError(
            f"No matching case directories found for case name: {case_name}"
        )

    # Sort alphabetically
    matching_dirs.sort()

    # If case_name.lower() == "runge", ensure all three sub-cases are present
    if case_name.lower() == "runge":
        sub_cases = ["Hep", "Dec", "Mix"]
        for sub_case in sub_cases:
            if not any(sub_case.lower() in d.lower() for d in matching_dirs):
                raise ValueError(f"Missing sub-case directory for Runge: {sub_case}")
        # Remove any with "jp8" in the name
        matching_dirs = [d for d in matching_dirs if "jp8" not in d.lower()]
        # Sort to ensure consistent order
        matching_dirs.sort(
            key=lambda x: ["hep", "dec", "mix"].index(
                next(sub for sub in ["hep", "dec", "mix"] if sub in x.lower())
            )
        )

    return matching_dirs


def setup(case_name):
    matching_dirs = find_case_dirs(case_name)
    cases = []
    leg_lab = []
    leg_col = []
    line_sty = []
    for k, d in enumerate(matching_dirs):
        if "gcm" in d.lower():
            name = d.split("_")[1]
            leg_lab.append("PeleGCM")
            cases.append(SpecifyCase(name, "gcm"))
            leg_col.append("tab:blue")
            cases[k].case_path = d
            if "_manifold" in d.lower():
                line_sty.append(":")
                leg_lab[k] += " + Manifold"
            else:
                line_sty.append("-")

        elif "mp" in d.lower():
            name = d.split("_")[1]
            leg_lab.append("PeleMP")
            line_sty.append("--")
            if "_manifold" in d.lower():
                leg_lab[k] += " + Manifold"
                line_sty[k] = ":"
            if "antoine" in d.lower():
                cases.append(SpecifyCase(name, "mp", "Antoine"))
                leg_col.append("tab:red")
                leg_lab[k] += ": Antoine"
            elif "cc" in d.lower():
                cases.append(SpecifyCase(name, "mp", "Clasius-Clapeyron"))
                leg_col.append("tab:orange")
                leg_lab[k] += ": Clasius-Clapeyron"
            else:
                raise ValueError(f"Unknown Psat model in directory name: {d}")
            cases[k].case_path = d
        else:
            raise ValueError(f"Unknown liquid properties model in directory name: {d}")

        # Specifics for Runge sub-cases
        if "rungehep" in name.lower():
            leg_lab[k] += ": Heptane"
            leg_col[k] = "tab:red"
        elif "rungedec" in name.lower():
            leg_lab[k] += ": Decane"
            leg_col[k] = "tab:blue"
        elif "rungemix" in name.lower():
            leg_lab[k] += ": Mix"
            leg_col[k] = "tab:orange"
        elif "rungejp8" in name.lower():
            if ":" not in leg_lab[k]:
                leg_lab[k] += ":"
            if "hychem" in d.lower():
                leg_lab[k] += " HyChem"
                line_sty[k] = "-."
            elif "detailed" in d.lower():
                leg_lab[k] += " Detailed"
                line_sty[k] = "--"
            else:
                leg_lab[k] += " Many-to-One"
    return cases, leg_lab, leg_col, line_sty


# Get cases and plotting parameters
cases, leg_lab, leg_col, line_sty = setup(case_name)

# Line and marker styles
marker_s = 40
line_w = 3
font_s = 16


def case_info(case):
    refdvals, reftvals, _ = ExtractRefVals(case)
    # Set end time based on reference values if available
    if refdvals is not None:
        time = refdvals[-1, 0] / case.xconv
        case.set_end_time(time)
    if not os.path.exists(case.case_path):
        raise ValueError(f"Case directory not found: {case.case_path}")
    outfile = os.path.join(case.case_path, "pele_vals.csv")
    if not os.path.isfile(outfile):
        try:
            pele_vals = ExtractData(case, outfile)
        except Exception as e:
            raise RuntimeError(f"Error extracting data for case {case.name}: {e}")
    else:
        df = pd.read_csv(outfile)
        pele_vals = df.to_numpy()
    return refdvals, reftvals, pele_vals


# Determine number of plots based on reference data
refdvals, reftvals, _ = ExtractRefVals(cases[0])
numplots = 1
if reftvals is not None:
    numplots += 1

ylabels = [cases[0].ylabel if hasattr(cases[0], "ylabel") else "$d/d_0$"]
if numplots == 2:
    ylabels.append("$T$ [K]")

# Use wider figure for JP8 case to accommodate legend
figwidth = 9.5 if case_name.lower() == "rungejp8" else 6.4
fig, axs = (
    plt.subplots(1, numplots, figsize=(numplots * figwidth, 4.8), constrained_layout=True)
    if numplots > 1
    else (plt.figure(figsize=(figwidth, 4.8), constrained_layout=True), [plt.gca()])
)

# Plot simulation lines first
for k, case in enumerate(cases):
    refdvals, reftvals, pele_vals = case_info(case)
    # Diameter plot
    i = 0
    axs[i].plot(
        pele_vals[:, 0],
        pele_vals[:, i + 1],
        line_sty[k],
        label=leg_lab[k],
        color=leg_col[k],
        linewidth=line_w,
    )
    axs[i].set_ylabel(ylabels[i], fontsize=font_s)
    axs[i].tick_params(labelsize=font_s)
    axs[i].set_xlim(min(pele_vals[:, 0]), max(pele_vals[:, 0]))
    axs[i].grid()
    axs[i].set_xlabel(
        case.xlabel if hasattr(case, "xlabel") else "Time", fontsize=font_s
    )

    # Temperature plot if available
    if numplots == 2 and reftvals is not None:
        i = 1
        axs[i].plot(
            pele_vals[:, 0],
            pele_vals[:, i + 1],
            line_sty[k],
            label=None,
            color=leg_col[k],
            linewidth=line_w,
        )
        axs[i].set_ylabel(ylabels[i], fontsize=font_s)
        axs[i].tick_params(labelsize=font_s)
        axs[i].set_xlim(min(pele_vals[:, 0]), max(pele_vals[:, 0]))
        axs[i].grid()
        axs[i].set_xlabel(
            cases[k].xlabel if hasattr(cases[k], "xlabel") else "Time", fontsize=font_s
        )

# Plot reference data last so legend entry is last
if case_name.lower() == "runge":
    # Diameter reference for each sub-case
    for idx in range(0, 6, 2):  # 0, 2, 4 (Heptane, Decane, Mix)
        refdvals, _, _ = case_info(cases[idx])
        i = 0
        label = f"{cases[idx].dname}" if idx == 0 else None
        if refdvals is not None:
            axs[i].scatter(
                refdvals[:, 0],
                refdvals[:, 1],
                marker="o",
                s=marker_s,
                facecolor="none",
                edgecolor="black",
                label=label,
                linewidth=round(line_w / 2),
            )
            # Plot uncertainty if available
            if refdvals.shape[1] == 4:
                uncrt = refdvals[~np.isnan(refdvals).any(axis=1)]
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
                for j in range(len(uncrt)):
                    tval = [uncrt[j, 0], uncrt[j, 0]]
                    uline = [uncrt[j, 2], uncrt[j, 3]]
                    axs[i].plot(tval, uline, "k-", linewidth=round(line_w / 2))
    # Temperature reference only for Heptane
    refdvals, reftvals, _ = case_info(cases[0])
    if numplots == 2 and reftvals is not None:
        i = 1
        axs[i].scatter(
            reftvals[:, 0],
            reftvals[:, 1],
            marker="o",
            s=marker_s,
            facecolor="none",
            edgecolor="black",
            label=None,
            linewidth=round(line_w / 2),
        )
        if reftvals.shape[1] == 4:
            uncrt = reftvals[~np.isnan(reftvals).any(axis=1)]
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
            for j in range(len(uncrt)):
                tval = [uncrt[j, 0], uncrt[j, 0]]
                uline = [uncrt[j, 2], uncrt[j, 3]]
                axs[i].plot(tval, uline, "k-", linewidth=round(line_w / 2))
else:
    # Non-Runge cases: reference data only for first case
    refdvals, reftvals, _ = case_info(cases[0])
    i = 0
    if refdvals is not None:
        axs[i].scatter(
            refdvals[:, 0],
            refdvals[:, 1],
            marker="o",
            s=marker_s,
            facecolor="none",
            edgecolor="black",
            label=f"{cases[0].dname}",
            linewidth=round(line_w / 2),
        )
        if refdvals.shape[1] == 4:
            uncrt = refdvals[~np.isnan(refdvals).any(axis=1)]
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
            for j in range(len(uncrt)):
                tval = [uncrt[j, 0], uncrt[j, 0]]
                uline = [uncrt[j, 2], uncrt[j, 3]]
                axs[i].plot(tval, uline, "k-", linewidth=round(line_w / 2))
    if numplots == 2 and reftvals is not None:
        i = 1
        axs[i].scatter(
            reftvals[:, 0],
            reftvals[:, 1],
            marker="o",
            s=marker_s,
            facecolor="none",
            edgecolor="black",
            label=None,
            linewidth=round(line_w / 2),
        )
        if reftvals.shape[1] == 4:
            uncrt = reftvals[~np.isnan(reftvals).any(axis=1)]
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
            for j in range(len(uncrt)):
                tval = [uncrt[j, 0], uncrt[j, 0]]
                uline = [uncrt[j, 2], uncrt[j, 3]]
                axs[i].plot(tval, uline, "k-", linewidth=round(line_w / 2))

# Single legend from axs[0]
handles, labels = axs[0].get_legend_handles_labels()
if case_name.lower() == "rungejp8":
    # Place legend outside for JP8 case
    axs[-1].legend(handles, labels, fontsize=font_s-2, loc="center left", bbox_to_anchor=(1, 0.5))
else:
    # Place legend inside for other cases
    axs[-1].legend(handles, labels, fontsize=font_s-2, loc="best")

plt.show()
