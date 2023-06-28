# PeleLMeX weak scaling

This folder contains two scripts used to perform weak scaling on PeleLMeX
FlameSheet case (or any other case periodic in x and y). One need to compile
the code with TINY\_PROFILER = TRUE in order to get the timing data.

Prior to using the scripts, the user first design a case on a single node,
setting up the domain dimension, number of grid cells on the base level,
number of AMR levels and tagging criterion such that the case fits nicely
in the available memory and perform the targeted number of time steps.

Once this single node input file is ready, use the `WeakScaling.py` script
to launch a user-specified number of jobs, increasing domain size/cells count
along with the number of node. The usage is as follows:

```
python3 WeakScaling.py -n TestName -e Pele**.ex -b batch.Script -i input_file_LMeX --extra_file extraFile1 extraFile2 ...
```

A `TestName` folder will be created, in which subfolders for each of the user-provided list of node counts specified
in the `WeakScaling.py` script (variable `nodeCounts` at l. 169) will be created. The PeleLMeX executable,
the batch script, the input file and any extra files are copied into each subfolders and the batch and input files
are adapted for weak scaling. Note that the number of nodes is expected to double between each consecutive entry of 
the `nodeCounts` list, consistently with a doubling of the domain size aternating between x and y in the input file
to maintain a constant load per node. The script will also submit jobs to the queue manager depending using the provided
batch script (subject to update if the queue manager change). A number of HPC machine are recognized by the script 
(Summit, Perlmutter, Frontier, ...) but one need to update the script for new platforms.

The second script, `ExtractScalingData.py` can be used when several or all the job submitted by the first script
have completed. Copy the script file into the top-level `TestName` folder and use: 

```
python3 ExtractScalingData.py -n TestName
```

Once again, the user need to update the `nodeCounts` variable in the script to part or all of the node counts
jobs launched by the first script. The resulting ASCII file will contain overall PeleLMeX timing and efficiency as
well as the timing of targeted pieces of the algorithm using TPROF data.
