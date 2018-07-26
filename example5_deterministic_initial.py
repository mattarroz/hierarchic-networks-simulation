# coding: utf-8
import stochpy
smod = stochpy.SSA()
smod.Model('Transc-Transl.psc')

smod.DoStochSim(end = 2.0,mode = 'time',trajectories = 100000)
#smod.PlotAverageSpeciesTimeSeries()
#stochpy.plt.title("")
#stochpy.plt.ylabel("Cells")
#stochpy.plt.ylabel("")
#stochpy.plt.xlabel("Time [minutes]")
#stochpy.plt.legend(["mRNA","Protein"],numpoints=1,loc='center right')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
smod.GetRegularGrid()
# species[species_index][trajectory_index][time_step_index]
protarray=np.array(smod.data_stochsim_grid.species[1])
# we transpose to get the time step index as the first index
protarray=protarray.transpose()
# select last time step
#ts=len(protarray)-1
ts=25
plt.hist(protarray[ts],range=[0,max(protarray[ts])],bins=max(protarray[ts]))
plt.title("")
plt.ylabel("Frequency")
plt.xlabel("Number of molecules")
plt.show()
ts=len(protarray)-1


plt.ioff()
with PdfPages('Transc-Transl-poisson.pdf') as pdf:
	fig=plt.hist(protarray[ts],range=[0,max(protarray[ts])],bins=max(protarray[ts]))
	plt.title("")
	plt.ylabel("Frequency")
	plt.xlabel("Number of molecules")
	pdf.savefig(fig)
