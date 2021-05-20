import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import MD_MC_LJ as md
import Parallel_Functions as pf
import time

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ===================
#     Initialize
# ===================
# Set configuration parameters
       # Number of particles
# Set simulation parameters


sigma = 1
T = 0.8
bet = 1/T
N1 = 125
N2 = 125
L1 = 6*sigma
L2 = 6*sigma
vol = L1**3 + L2**3
delV = vol/20
N = N1 + N2
nPart = 3000
nVol = 50
nSwap = 6000
rcut = 2.5*sigma
invdrcut6 = (1/rcut)**6
phicut = invdrcut6*(invdrcut6-1)
cut = [rcut**2,phicut]
maxDr = np.minimum(L1,L2)/10
rho = np.array([N1/L1**3,N2/L2**3])
marr = md.index_alloc(nPart-1,size)
#m = [mbeg,mend]

if rank == 0:
	print("\nGibbs Ensemble Monte Carlo Simulation of the Lennard - Jones Fluid\n")
	print("Temperature of the System = "+str(T))
	print("Volume of the Simulation = "+str(vol))
	print("Number of Particles = "+str(N))
	print("Net Density = "+str(N/vol))


	# Clear files if they already exist.
	if os.path.exists('Data.csv'):
		os.remove('Data.csv')
	'''
	with open('Data.csv', 'a') as f:
		f.write('{0}, {1}\n'.format('Step', 'L1','L2','N1','N2','E1','E2',''))
	'''

	coords1 = md.initCubicGrid(N1,L1)
	coords2 = md.initCubicGrid(N2,L2)

	e1 = pf.LJ_Energy_Master(comm,[coords1,L1,cut,size,marr],1)
	e2 = pf.LJ_Energy_Master(comm,[coords2,L2,cut,size,marr],1)

	data = [coords1,coords2,L1,L2,N1,N2,e1,e2]

	nSteps = 2000
	nTherm = 500
	nFreq = 10

	sum_rho_1 = 0
	sum_rho_2 = 0
	num_rho_1 = 0
	num_rho_2 = 0
	rms_rho_1 = 0
	rms_rho_2 = 0
	t_write = 0
	n_write = 0
	s_write = 0
	t_disp = 0
	n_disp = 0
	s_disp = 0
	t_vol = 0
	n_vol = 0
	s_vol = 0
	t_swap = 0
	n_swap = 0
	s_swap = 0
	for step in range(1,nSteps+1):
		if (data[4]== 0)or(data[5]==0):
			print("Box is Empty")
		if (data[4]+data[5] != N):
			print("STOP")
		ran = np.random.rand()*(nPart + nVol + nSwap)
		#ran = nPart - 3
		if (ran<nPart):

			coords1 = data[0]
			coords2 = data[1]
			L1 = data[2]
			L2 = data[3]
			N1 = data[4]
			N2 = data[5]

			#	Particle Displacement
			
			#	Choose Particle to be displaced
			i = np.random.randint(0,N)

			if (i<N1):
				box = 0
			else:
				box = 1

			#	Calculate old energy
			coords = data[box]
			L = data[box+ 2]
			e1_o = pf.LJ_Energy_Master(comm,[coords,L,cut,size,marr],1)

			#	Propose a Trial
			if i>=N1:
				i = i-N1
			rTrial = coords[i,:] + maxDr*(np.random.rand(3,)-0.5)
			rTrial = md.PBC3D(rTrial,L)
			deltaE = md.LJ_EnergyIncrease(coords,rTrial,i,L,rcut)

			#	Accept or Reject
			if (np.random.rand()<np.exp(-bet*(deltaE))):
				coords[i,:] = rTrial
				data[box] = coords
				data[6+box] = e1_o + deltaE
			


		elif(ran<nPart+nVol):
			#	Volume Change Trial Move

			#	Calculate Energies
			coords1_o = data[0]
			coords2_o = data[1]
			L1_o = data[2]
			L2_o = data[3]
			e1_o = pf.LJ_Energy_Master(comm,[coords1_o,L1_o,cut,size,marr],1)
			e2_o = pf.LJ_Energy_Master(comm,[coords2_o,L2_o,cut,size,marr],1)
			N1 = data[4]
			N2 = data[5]
			
			#	Propose Volume Change
			vol1_o = L1_o**3
			vol2_o = L2_o**3
			#print(vol1_o)
			vol1_n = -1
			while vol1_n<0:
				vol1_n = vol1_o + (2*np.random.rand()-1)*delV

			vol2_n = vol - vol1_n
			L1_n = vol1_n**(1/3)
			L2_n = vol2_n**(1/3)
			coords1_n = coords1_o*(L1_n/L1_o)
			coords2_n = coords2_o*(L2_n/L2_o)
			e1_n = pf.LJ_Energy_Master(comm,[coords1_n,L1_n,cut,size,marr],1)
			e2_n = pf.LJ_Energy_Master(comm,[coords2_n,L2_n,cut,size,marr],1)

			#	Accept or Reject
			p = np.exp(-bet*(e2_n - e1_n))*((vol1_n/vol1_o)**N1)*((vol2_n/vol2_o)**N2)
			if np.random.rand()<p:
				data[0] = coords1_n
				data[1] = coords2_n
				data[2] = L1_n
				data[3] = L2_n
				data[6] = e1_n
				data[7] = e2_n
			
		else:
			#	Particle Swap Trial Move
			#	Select Box
			if np.random.rand()<0.5:
				box = 0
			else:
				box = 1
			#	Propose Trial Move
			Lin  = data[box + 2]
			Lout = data[1-box+2]
			rnew = Lin*np.random.rand(3,)

			#	Calculate Energies
			coords_in_o = data[box]
			coords_in_n = np.vstack((coords_in_o,rnew))
		
			en_in_o = pf.LJ_Energy_Master(comm,[coords_in_o,Lin,cut,size,marr],1)
			en_in_n = pf.LJ_Energy_Master(comm,[coords_in_n,Lin,cut,size,marr],1)

			coords_out_o = data[1-box]
			coords_out_n = np.delete(coords_out_o,np.random.randint(0,data[5-box]),0)
			en_out_o = pf.LJ_Energy_Master(comm,[coords_out_o,Lout,cut,size,marr],1)
			en_out_n = pf.LJ_Energy_Master(comm,[coords_out_n,Lout,cut,size,marr],1)

			Nin  = data[box+4]
			Nout = data[5-box]

			#	Accept or Reject
			p = (Nin/(Nout+1))*((Lout/Lin)**3)*np.exp(-bet*((en_in_n+en_out_n)-(en_in_o+en_out_o)))
			if data[5-box] > 1:
				if np.random.rand()<p:
					data[box+4] = Nin+ 1
					data[5-box] = Nout-1
					data[ box ] = coords_in_n
					data[1-box] = coords_out_n
					data[box+6] = en_in_n
					data[7-box] = en_out_n

		rho_1 = data[4]/(data[2]**3)
		rho_2 = data[5]/(data[3]**3)

		rho = np.vstack((rho,[rho_1,rho_2]))

		if step>nTherm:
			if step%nFreq == 0:
				t1 = time.time()
				sum_rho_1 += rho_1
				sum_rho_2 += rho_2
				num_rho_1 += 1
				num_rho_2 += 1
				rms_rho_1 += rho_1**2
				rms_rho_2 += rho_2**2
				md.Write_to_file(step,data,'Data.csv')
				t2 = time.time()
				t_write += (t2 - t1)
				s_write += (t2 - t1)**2
				n_write += 1

		sys.stdout.write("\rStep Number: {0}\t\t ".format(step))
		sys.stdout.flush()

	plt.plot(rho)
	plt.title("Gibbs Ensemble Simulation of Lennard Jones Fluid")
	#plt.show()
	plt.xlabel("Step Number")
	plt.ylabel("$\\rho^{*}$")
	plt.savefig("Result.png")

	sys.stdout.write("\nDensity 1 = {0}\t Density 2 = {1}".format(sum_rho_1/num_rho_1,sum_rho_2/num_rho_2))
	sys.stdout.write("\nStandard Deviation = {0} Standard Deviation = {1}".format(np.sqrt(rms_rho_1/num_rho_1-(sum_rho_1/num_rho_1)**2),np.sqrt(rms_rho_2/num_rho_2-(sum_rho_2/num_rho_2)**2)))
	sys.stdout.flush()
	pf.Shut_Down(comm)

else:
	pf.LJ_Energy_Slave(comm)