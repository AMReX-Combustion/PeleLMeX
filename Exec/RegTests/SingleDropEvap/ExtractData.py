import os
import csv
import numpy as np
import pandas as pd
import re


def ExtractData(case, outfile):
    FILE_PATH = os.path.dirname(os.path.abspath(__file__))
    CASE_PATH = os.path.join(FILE_PATH, case.name)
    pltfiles = []
    sprayfiles = []
    for x in os.listdir(CASE_PATH):
        if x.startswith("plt"):
            pltfiles.append(os.path.join(CASE_PATH,x))
        if x.endswith(".p3d"):
            sprayfiles.append(os.path.join(CASE_PATH,x))

    def get_step(fn):
        match = re.search(r"(?:plt|spray)(\d+)", fn)
        return int(match.group(1)) if match else -1

    pltfiles = sorted(pltfiles, key=get_step)
    sprayfiles = sorted(sprayfiles, key=get_step)
    time = []
    # Open Header in the plt files and extract the solution time
    for cf in pltfiles:
        curfile = cf + "/Header"
        Lines = []
        with open(curfile, "r") as fn:
            Lines = fn.readlines()
            numcomp = int(Lines[1])
            timeline = numcomp + 3
        time.append(float(Lines[timeline]))

    # Column designations in the spray*.p3d files
    numspec = len(case.droplet.Y)  # Liquid fuel components
    dims = 2 # Solution dimensions
    if "z_velocity\n" in Lines:
        dims = 3  
    loccols = dims - 1
    velcols = loccols + 1
    tcol = velcols + dims
    dcol = tcol + 1
    mfcol = dcol + 1
    numcol = mfcol + numspec
    timevals = []
    vals = np.zeros([len(sprayfiles), numcol])
    crow = 5
    for i, cf in enumerate(sprayfiles):
        with open(cf, "r") as fn:
            Lines = fn.readlines()
            if len(Lines) >= crow + 1:
                timevals.append(time[i])
                sline = Lines[crow].split()
                for col in range(numcol):
                    vals[i][col] = float(sline[col])

    # Unit conversions for plotting
    xconv = case.xconv
    yconv = case.yconv
    yexp = case.yexp

    modvals = []
    with open(outfile, "w+") as new_file:
        new_file.write("t, dd0, T, Y1, Y2\n")
        csv_writer = csv.writer(new_file, delimiter=",", lineterminator="\n")
        for k, tv in enumerate(timevals):
            dia = vals[k][dcol]
            dd0 = (dia * yconv) ** yexp
            T = vals[k][tcol]
            Y1 = vals[k][mfcol]
            if len(case.droplet.Y) == 2:
                Y2 = vals[k][mfcol + 1]
                outvals = [tv * xconv, dd0, T, Y1, Y2]
            else:
                outvals = [tv * xconv, dd0, T, Y1]
            modvals.append(outvals)
            csv_writer.writerow(outvals)
    modvals = np.array(modvals)
    return modvals


def ExtractRefVals(case):
    FILE_PATH = os.path.dirname(os.path.abspath(__file__))

    ldir = os.path.join(FILE_PATH, f"ref_files/{case.name}")
    fnames = ["refdvals.csv", "refTvals.csv", "refYvals.csv"]
    reffiles = []

    def getdata(fname):
        df = pd.read_csv(fname)
        return df.to_numpy()

    dvals = None
    tvals = None
    yvals = None
    cname = fnames[0]
    if cname in os.listdir(ldir):
        dvals = getdata(os.path.join(ldir, cname))
    cname = fnames[1]
    if cname in os.listdir(ldir):
        tvals = getdata(os.path.join(ldir, cname))
    cname = fnames[2]
    if cname in os.listdir(ldir):
        yvals = getdata(os.path.join(ldir, cname))
    return [dvals, tvals, yvals]
