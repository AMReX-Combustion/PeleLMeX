import os
import numpy as np
from CaseInfo import *
from ExtractData import *
import matplotlib.pyplot as plt

case = "Runge"  # "Nomura", "WongLin", "Daif", or "Runge"

# Set up cases and labels automatically
if case.lower() == "nomura":
    cases = [Nomura("gcm"), Nomura("mp")]
    leg_lab = ["PeleGCM", "PeleMP"]
    leg_col = ["red", "red"]
    line_sty = ["-", "--"]
elif case.lower() == "wonglin":
    cases = [WongLin("gcm"), WongLin("mp")]
    leg_lab = ["PeleGCM", "PeleMP"]
    leg_col = ["red", "red"]
    line_sty = ["-", "--"]
elif case.lower() == "daif":
    cases = [Daif("gcm"), Daif("mp")]
    leg_lab = ["PeleGCM", "PeleMP"]
    leg_col = ["red", "red"]
    line_sty = ["-", "--"]
elif case.lower() == "runge":
    cases = [
        RungeHep("gcm"), RungeHep("mp"),
        RungeDec("gcm"), RungeDec("mp"),
        RungeMix("gcm"), RungeMix("mp")
    ]
    leg_lab = [
        "PeleGCM: Heptane", "PeleMP: Heptane",
        "PeleGCM: Decane", "PeleMP: Decane",
        "PeleGCM: Mix", "PeleMP: Mix"
    ]
    leg_col = ["red", "red", "blue", "blue", "orange", "orange"]
    line_sty = ["-", "--", "-", "--", "-", "--"]
else:
    raise ValueError(f"Unknown case: {case}")

marker_s = 40
line_w = 3
font_s = 16

def case_info(case):
    refdvals, reftvals, pele_vals = ExtractRefVals(case)
    # Set end time based on reference values if available
    if refdvals is not None:
        time = refdvals[-1, 0] / case.xconv
        case.set_end_time(time)
    if not os.path.exists(case.case_dir):
        raise ValueError(f"Case directory not found: {case.case_dir}")
    outfile = os.path.join(case.case_dir, "pele_vals.csv")
    pele_vals = ExtractData(case, outfile)
    return refdvals, reftvals, pele_vals

# Determine number of plots based on reference data
refdvals, reftvals, _ = ExtractRefVals(cases[0])
numplots = 1
if reftvals is not None:
    numplots += 1

ylabels = [cases[0].ylabel if hasattr(cases[0], "ylabel") else "$d/d_0$"]
if numplots == 2:
    ylabels.append("$T$ [K]")

fig, axs = plt.subplots(1, numplots, figsize=(numplots * 6.4, 4.8)) if numplots > 1 else (plt.figure(), [plt.gca()])

# Plot simulation lines first
for k in range(len(cases)):
    refdvals, reftvals, pele_vals = case_info(cases[k])
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
    axs[i].set_xlabel(cases[k].xlabel if hasattr(cases[k], "xlabel") else "Time", fontsize=font_s)

    # Temperature plot if available
    if (numplots == 2) and k < 2:
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
        axs[i].set_xlabel(cases[k].xlabel if hasattr(cases[k], "xlabel") else "Time", fontsize=font_s)

# Plot reference data last so legend entry is last
if case.lower() == "runge":
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

# Single legend from axs[0] placed in rightmost subplot
handles, labels = axs[0].get_legend_handles_labels()
axs[-1].legend(handles, labels, fontsize=font_s, loc="best")

plt.tight_layout()
plt.show()