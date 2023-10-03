#!/usr/bin/env python3.6

import sys
import os
import shutil
import argparse
import socket
import numpy as np
import subprocess
import fnmatch

USAGE = """
    A script to extract weak scaling data from PeleLMeX logfiles.
"""

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("-n","--test_name", type=str, default="Case1", metavar="test-name",
                        help="Name of the test. Default = Case1")

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    return args

if __name__ == "__main__":
    args = parse_args(arg_string=sys.argv[1:])

    # log file pattern
    logprefix = '{}_*'.format(args.test_name)
    nodeCounts = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192]
    #nodeCounts = [1,2,4,8,16,32,64,128,256,512,8192]

    caselist = []
    for nCount in nodeCounts:
        caselist.append("{:04d}".format(nCount))

    # get list of cases from folders
    cases = [ os.path.relpath(f.path,"./") for f in os.scandir(".") if f.is_dir() ]
    cases.sort()

    # data contqiner
    runTimes = []
    reactionTimes = []
    diffusionTimes = []
    MacProjTimes = []
    NodalProjTimes = []
    VelAdvTimes = []
    ScalAdvTimes = []
    AvgDownTimes = []
    ParCopyTimes = []

    for case in caselist:
        print(case)
        folder = "./{}".format(case)

        logfile = ""
        for root, dirs, files in os.walk(folder):
            for name in files:
                if fnmatch.fnmatch(name,logprefix):
                    logfile = os.path.join(root, name)
                    break

        if (logfile == ""):
            print("WARNING ! Could not find logfile in {}".format(case))
            continue

        # Get runTime
        cmd = "cat {}".format(logfile)+" | grep 'Total Time:' | awk -F: '{print $2}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currRunTime = procRunTime.communicate()[0].decode("utf-8").strip()
        runTimes.append(currRunTime)
        # Get component times: reaction, diffusion, MacProj, NodalProj, ScalAdv, VelAdv, Sync
        cmd = "cat {}".format(logfile)+" | grep 'PeleLMeX::advance::reactions' | awk 'NR%2==0' | awk 'NR%2!=0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currReactTime = procRunTime.communicate()[0].decode("utf-8").strip()
        reactionTimes.append(currReactTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLMeX::advance::diffusion' | awk 'NR%2==0' | awk 'NR%2!=0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currDiffTime = procRunTime.communicate()[0].decode("utf-8").strip()
        diffusionTimes.append(currDiffTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLMeX::advance::mac' | awk 'NR%2==0' | awk 'NR%2!=0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currMACTime = procRunTime.communicate()[0].decode("utf-8").strip()
        MacProjTimes.append(currMACTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLMeX::velocityProjection()' | awk 'NR%2==0' | awk 'NR%2!=0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currProjTime = procRunTime.communicate()[0].decode("utf-8").strip()
        NodalProjTimes.append(currProjTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLMeX::advance::scalars_adv' | awk 'NR%2==0' | awk 'NR%2!=0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currScalAdvTime = procRunTime.communicate()[0].decode("utf-8").strip()
        ScalAdvTimes.append(currScalAdvTime)
        cmd = "cat {}".format(logfile)+" | grep 'PeleLMeX::advance::velocity' | awk 'NR%2==0' | awk 'NR%2!=0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currVelAdvTime = procRunTime.communicate()[0].decode("utf-8").strip()
        VelAdvTimes.append(currVelAdvTime)
        cmd = "cat {}".format(logfile)+" | grep 'amrex::average_down ' | awk 'NR%2==0' | awk 'NR%2!=0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currAvgTime = procRunTime.communicate()[0].decode("utf-8").strip()
        AvgDownTimes.append(currAvgTime)
        cmd = "cat {}".format(logfile)+" | grep 'average_down_faces' | awk 'NR%2==0' | awk 'NR%2!=0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currAvgTime = procRunTime.communicate()[0].decode("utf-8").strip()
        prevAvgTime = AvgDownTimes[-1]
        AvgDownTimes[-1] = float(prevAvgTime) + float(currAvgTime)
        cmd = "cat {}".format(logfile)+" | grep 'FabArray::ParallelCopy()' | awk 'NR%2==0' | awk 'NR%2!=0' | awk '{print $4}'"
        procRunTime = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        currParCopyTime = procRunTime.communicate()[0].decode("utf-8").strip()
        ParCopyTimes.append(currParCopyTime)

    fout = open("LMeX_ScalingData.dat", "w")
    fout.write(" NodeCount RunTime Efficiency Reaction Diffusion MacProj NodalProj ScalAdv VelUpdate AvgDown ParCopy\n")
    for n in range(len(caselist)):
        fout.write("{} {} {} {} {} {} {} {} {} {} {}\n".format(nodeCounts[n], runTimes[n], float(runTimes[0])/float(runTimes[n]), reactionTimes[n], diffusionTimes[n], MacProjTimes[n], NodalProjTimes[n], ScalAdvTimes[n], VelAdvTimes[n], AvgDownTimes[n], ParCopyTimes[n]))
    fout.close()

