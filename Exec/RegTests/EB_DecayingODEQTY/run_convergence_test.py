import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def generate_input_file(final_time, dt, test, temporal_int, template_filename='input.convergence_0'):
    new_input_filename = f'input.convergence_{test}'

    # Read the template input file
    with open(template_filename, 'r') as file:
        lines = file.readlines()

    # Modify amr.fixed_dt
    for i, line in enumerate(lines):
        line_start = line.strip().split()
        if len(line_start) > 0:
            line_start = line_start[0]
            if line_start == 'amr.fixed_dt':
                lines[i] = f'amr.fixed_dt = {dt}\n'
            elif line_start == 'amr.plot_file':
                lines[i] = f'amr.plot_file = "plt_conv{test}_"\n'
            elif line_start == 'peleLM.temporal_int':
                lines[i] = f'peleLM.temporal_int = {temporal_int}\n'
            elif line_start == 'amr.stop_time':
                lines[i] = f'amr.stop_time = {final_time}\n'
    # Write the new input file
    with open(new_input_filename, 'w') as file:
        file.writelines(lines)

    return new_input_filename

def run_PeleLMeX(input_file,num_proc):
    if (num_proc > 1):
        os.system(f'mpirun -np {num_proc} ./PeleLMeX2d.gnu.MPI.ex {input_file}')
    else:
        os.system(f'./PeleLMeX2d.gnu.MPI.ex {input_file}')

def soln(t,IC):
    k = -100.0
    return IC*np.exp(k*t)

def dataAnalysis(vars,file_name):
    data = pd.read_csv(file_name)
    data_vars = data[vars].values
    time = data.time.values
    IC = data_vars[0,:]
    
    approx = data_vars[-1,:]
    exact = soln(time[-1],IC)
    abs_error = abs(exact - approx)

    return abs_error.reshape(1,len(vars))


def removeAllPrevData():
    os.system('rm -r plt* temporals/*')


# Number of runs for convergence test
startWithNewData = 0
nRuns = 4
num_proc = 6
dt_vec = [1e-3, 1e-4, 1e-5, 1e-6]
final_time = 0.05
num_temporal_write = 60
temporals_directory = 'temporals_simpleDecay_SDC0'
temporals_filename = temporals_directory+'/tempExtremas'

NUM_ODE = 3
vars = [f"max_MY_ODE_{i}" for i in range(0,NUM_ODE)]
errors = np.zeros((nRuns,NUM_ODE))

print("\n==============================================================================")
if startWithNewData:
    removeAllPrevData()

for n in range(0, nRuns):
    dt = dt_vec[n]
    Nt = round(final_time / dt)
    temporal_int = round(Nt / num_temporal_write)
    print(f"Running test #{n} with dt = {dt}, Nt =  {Nt}, temporal_int = {temporal_int}")

    # Generate a new input file for the given dt
    input_file = generate_input_file(final_time, dt, n, temporal_int)

    # Run the convergence test in PeleLMeX
    if startWithNewData:
        run_PeleLMeX(input_file,num_proc)

    # Rename the tempExtremas file and remove unused tempState file
    new_temporals_filename = temporals_filename + f'_{n}'
    os.system(f'mv {temporals_filename} {new_temporals_filename}')
    os.system('rm temporals/tempState')

    print(f"Finished running test #{n}\n")
    
    # Read in extrema data from temporals
    errors[n,:] = dataAnalysis(vars,new_temporals_filename)

print("\n==============================================================================\n")

plt.figure()
for j in range(NUM_ODE):
    

    # Fit a line to the log-log data
    log_dt_vec = np.log(dt_vec[:nRuns])
    log_errors = np.log(errors[:,j])
    slope, intercept = np.polyfit(log_dt_vec, log_errors, 1)
    
    # Plot scatter plot of errors and the line of best fit
    best_fit_line = np.exp(intercept) * np.array(dt_vec[:nRuns])**slope
    plt.scatter(dt_vec[:nRuns], errors[:,j], label=f'ODE {j} (rate: {slope:.2f})', marker='o')
    plt.plot(dt_vec[:nRuns], best_fit_line, '--')

    # Reset the color cycle
    #plt.gca().set_prop_cycle(None)

# Set x and y scales to log-log
plt.xscale("log")
plt.yscale("log")

# Add labels and legend
plt.xlabel('$\Delta t$')
plt.ylabel('Error')
plt.legend()

plt.show()