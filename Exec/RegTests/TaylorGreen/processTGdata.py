import sys
sys.path.append('utils')
import numpy as np
import matplotlib.pyplot as plt

tmpState_file = 'temporals/tempState'
ref_file = 'refData'

stateData = np.loadtxt(tmpState_file)
refData = np.loadtxt(ref_file)

# Simulation data
t_norm = 0.00010580074
vel_norm = 30.08579015
rho0 = 1.0
dom_length = 0.02
visc = 0.0005985379607

stateTime = stateData[:,1] / t_norm
stateEnstrophy = stateData[:,4] / (dom_length*dom_length*dom_length) * t_norm * t_norm
stateKinEnergy = stateData[:,3] / (dom_length*dom_length*dom_length * vel_norm * vel_norm)
stateDissipationKin = -np.diff(stateKinEnergy) / np.diff(stateTime)
stateDissipationEns = 2.0 * visc / rho0 * stateEnstrophy

refTime = refData[:,0]
refEnstrophy = refData[:,3]
refKinEnergy = refData[:,1]
refDissipation = refData[:,2]

plt.plot(stateTime,stateEnstrophy,linewidth=2, color='r',label='PeleLMeX')
plt.plot(refTime,refEnstrophy,linewidth=0,marker='o',markevery=40, color='k',label='Spectral 512^3')
plt.xlabel("t* [-]")
plt.ylabel("Enstrophy* [-]")
plt.grid(which='both',color='k', linestyle=':', linewidth=0.4)
plt.legend(bbox_to_anchor=(0.95, 0.9), loc=1, borderaxespad=0.)
plt.savefig("Enstrophy.png")
plt.figure()
plt.plot(stateTime,stateKinEnergy,linewidth=2, color='r',label='PeleLMeX')
plt.plot(refTime,refKinEnergy,linewidth=0,marker='o',markevery=40, color='k',label='Spectral 512^3')
plt.xlabel("t* [-]")
plt.ylabel("Kinetic Energy* [-]")
plt.grid(which='both',color='k', linestyle=':', linewidth=0.4)
plt.legend(bbox_to_anchor=(0.95, 0.9), loc=1, borderaxespad=0.)
plt.savefig("KinEnergy.png")
plt.figure()
plt.plot(stateTime[0:-1],stateDissipationKin,linewidth=2, color='r',label='PeleLMeX - eps_1')
plt.plot(stateTime,stateDissipationEns,linewidth=2, color='b',label='PeleLMeX - eps_2')
plt.plot(refTime,refDissipation,linewidth=0,marker='o',markevery=40, color='k',label='Spectral 512^3')
plt.xlabel("t* [-]")
plt.ylabel("Dissipation* [-]")
plt.grid(which='both',color='k', linestyle=':', linewidth=0.4)
plt.legend(bbox_to_anchor=(0.95, 0.9), loc=1, borderaxespad=0.)
plt.savefig("Dissipation.png")
exit()
