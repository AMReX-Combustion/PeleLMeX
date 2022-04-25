#!/usr/bin/env python3

# Template post-processing script for PeleLM convergence analysis
# Must be used after multirun.py script
# Input are limited by the regression framework.

# Usage:
#   ./compareOutput.py --test_name DummyTest

# Input:
#   * --test_name: a TESTNAME that will looked for during the postprocessing
#   * --max_error: the maximum allow able error

# "Internal" user input
#   * vars : a list of the variables of interest (no check is done on whether it exists in plt ...)

# Output:
#  If the error in vars is higher than the max_error, an error statement is returned

# Head's up :
#   - The script will get a copy of the post-processing program (if not already there) in the testing folder. The name of this folder is assumed to be the TESTNAME.
#   - The plt files naming convention is: ${TESTNAME}/64_1_plt***** where 64 is the box size and 1 is the first test
#   - Errors are parsed from the screen output of the standard fcompare. Beware of any change of these programs.

import sys
import os
import fnmatch
import shutil
import argparse
import numpy as np

USAGE = """
    Template post-processing script for PeleLM convergence analysis
"""

def pproc(args):

    # User data
    vars=["density", "rho.Y(NC10H22)", "x_velocity", "temp", "rhoh"]

    # Get a local copy of post-processing executable
    run_dir = os.getcwd()
    pproc_exe = "fcompare.llvm.ex"

    # Check the test name: current folder name is default
    if ( args.test_name == "None" ):
        args.test_name = "testfiles"

    # Run the postprocessing
    pltfiles = ['', '', '', '']
    # Put all plot files in list
    # Note: this assumes only the final plot files are there and plt00000 has been removed
    for f in os.listdir(args.test_name):
        if ( not fnmatch.fnmatch(f, '*old*')):
            if (f.startswith("32_1_plt")):
                pltfiles[0] = f
            elif (f.startswith("32_2_plt")):
                pltfiles[1] = f
            elif (f.startswith("64_1_plt")):
                pltfiles[2] = f
            elif (f.startswith("64_2_plt")):
                pltfiles[3] = f
    # We have 3 comparisons to make,
    # 32_1 to 32_2, 64_1 to 64_2, and 32_1 to 64_1
    numcomps = 3
    comp1 = [0, 2, 0]
    comp2 = [1, 3, 2]
    errors = np.empty([numcomps,len(vars)+1])
    for comp in range(numcomps):
        indx1 = comp1[comp]
        indx2 = comp2[comp]
        plt1 = args.test_name + "/" + pltfiles[indx1]
        plt2 = args.test_name + "/" + pltfiles[indx2]
        print("Comparing " + plt1 + " with " + plt2)
        outfile = "error_{}.analysis.out".format(comp)
        os.system("./{} -n 2 -a {} {} > {}".format(os.path.basename(pproc_exe), plt1, plt2, outfile))
        # Extract errors on each variable
        with open(outfile) as fp:
            for i, line in enumerate(fp):
                if (i >= 5):
                    var = line.split()[0]
                    for v in range(len(vars)):
                        if ( var == vars[v] ):
                            errors[comp,v+1] = line.split()[2]
        os.system("rm {}".format(outfile))
    # Check error
    passed = checkError(errors, args.test_name, vars, args.max_error)

def checkError(data, test_name, vars, maxerror):
    # Evaluate order
    for v in range(len(vars)):
        for i in range(len(data[:,0])):
            if (data[i,v] > maxerror):
                errorStatement="{} did not reach target convergence order: {} > {}".format(vars[v],data[i,v],maxerror)
                raise ValueError(errorStatement)

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("--test_name", type=str, default="None", metavar="test-name",
                        help="name of the test. Default = current folder name")

    parser.add_argument("--max_error", type=float, default=1.E-12,
                        help="max error for tests.")

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    return args

if __name__ == "__main__":
    args = parse_args(arg_string=sys.argv[1:])
    pproc(args)
