import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

# A script to plot data from the EB_PipeFlow case
# when running a simple Poiseuille flow.
# Flow properties are hardcoded to: G = 82944.0,
# mu = 0.0576, radius = 0.01, umax = 36.0
# Data from LMeX run at increasing resolution must be extracted
# using fextract a-priori and stored in nr* folders.

cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]

def theory_ooa(order, res, orig):
    return orig * (res[0] / res) ** order

def eval_u_exact(r, dr):
    return 82944.0 / (4.0 * 0.0576) * (0.01 ** 2 - r ** 2)

if __name__ == "__main__":
    radius = 0.01
    umax = 36
    r = np.linspace(-radius, radius, 200)
    u_exact = eval_u_exact(r, 0.0001)

    # Plot the exact data (pointwise, not cell-averaged)
    plt.figure("u")
    plt.rc('font', size=12)
    plt.plot(
        r / radius, u_exact/umax, lw=2, color=cmap[-1], label="Exact"
    )

    # Read in LMeX data
    resolution=[8,16,32,64]

    errors = np.zeros((2, len(resolution)))
    for k, r in enumerate(resolution):
        f = open('./nr{}/nr{}prof.dat'.format(r,r), 'r')
        rlabel="res{}".format(r)
        rad = []
        u_x = []
        for line in f:
            line = line.strip()
            columns = line.split()
            rad.append(float(columns[0]))
            u_x.append(float(columns[1]))
        plt.plot(
            np.array(rad)/radius, np.array(u_x)/umax, linestyle='--', lw=1, color=cmap[k], label=rlabel
        )
        errors[0, k] = r
        errors[1, k] = np.sqrt(
            np.sum(
                (
                    np.array(u_x) / umax
                    - eval_u_exact(np.array(rad), 0.0001) / umax
                )
                ** 2
            )
            / r
        )

    plt.xlabel(r"$r / R$",fontsize=14)
    plt.ylabel(r"$u / Umax$",fontsize=14)
    plt.grid(which='both',color='#E6E3E3', linestyle=':', linewidth=1.0)
    plt.legend(bbox_to_anchor=(0.5, 0.4), loc=1, borderaxespad=0.)
    plt.savefig("PoiseuilleVelProf.png")

    plt.figure("error")
    plt.loglog(
        errors[0, :],
        errors[1, :],
        label="PeleLMeX",
    )
    p2 = theory_ooa(2, errors[0, :], 0.95*errors[1, 0])
    plt.loglog(errors[0, :], p2, color='k', linestyle='--', label=f"2nd-order")
    plt.grid(which='both',color='#E6E3E3', linestyle=':', linewidth=1.0)
    plt.xlabel("Resolution")
    plt.ylabel("Error L2-norm")
    plt.legend(bbox_to_anchor=(0.9, 0.9), loc=1, borderaxespad=0.)
    plt.savefig("PoiseuilleConvergence.png")
