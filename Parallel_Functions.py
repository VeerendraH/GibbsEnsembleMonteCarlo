from mpi4py import MPI
import numpy as np
import sys
import MD_MC_LJ as md

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def Summer(comm,a,program):
	print(a,"h")
	switch = 1
	rank = comm.Get_rank()
	size = comm.Get_size()
	comm.barrier()
	comm.bcast(switch,root=0)

	comm.barrier()
	comm.bcast(a,root=0)

	t = a*rank*np.array([1,2,3,4])
	b = np.array([0,0,0,0])
	comm.Reduce(t,b,op=MPI.SUM,root=0)

	comm.barrier()
	comm.bcast(program,root = 0)

	return b
def Shut_Down(comm):
	program = 0
	switch = 0
	comm.barrier()
	comm.bcast(switch,root=0)
	comm.barrier()
	comm.bcast(program,root=0)
	return
	

def LJ_Energy_Master(comm,a,program):
	
	switch = 1
	rank = comm.Get_rank()
	size = comm.Get_size()
	comm.barrier()
	comm.bcast(switch,root=0)

	comm.barrier()
	comm.bcast(a,root=0)

	t = np.array([Energy_LJ_Rank(a[0],a[1],a[2],rank,a[4])])
	b = comm.reduce(t,op=MPI.SUM,root=0)
	
	comm.barrier()
	comm.bcast(program,root = 0)

	return b[0]
def LJ_Energy_Slave(comm):
	switch = 0
	a = None
	program = 1	#Bulb 1 is on
	while program == 1:	#Day begins as Bulb 1 is on 
		comm.barrier()
		switch = comm.bcast(switch,root = 0) #Look at the bulb 2

		if switch == 1:		#If Bulb 2 is on, start working on the function
			comm.barrier()				#These two sentences equate picking
			a = comm.bcast(a,root=0)	#up the phone for Inputs 
			t = np.array([Energy_LJ_Rank(a[0],a[1],a[2],rank,a[4])])	#Rank-specific Transform
			comm.reduce(t,op=MPI.SUM,root=0)	#Return Outputs
			switch = 0					#Equivalent to turning off Bulb 2
		comm.barrier()							#These two sentences equate to 
		program = comm.bcast(program,root = 0)	#asking if they can leave
		
	return

def Energy_LJ_Rank(coords,L,cut,rank,m):
	energy = 0
	nPart = coords.shape[0]
	rcut2 = cut[0]
	phicut = cut[1]
	mbeg = m[0][rank][0]
	mend = m[1][rank][0]
	for partA in range(mbeg,mend):
		for partB in range(partA+1,nPart):
			dr = coords[partA,:] - coords[partB,:]
			dr = md.distPBC3D(dr,L)
			dr2 = np.sum(dr**2)
			if dr2<rcut2:
				invDr6 = (1/dr2)**3
				energy = energy + (invDr6*(invDr6-1)) - (phicut)

	energy = 4*energy
	return energy


def Wait(comm):
	switch = 0
	a = None
	program = 1	#Bulb 1 is on
	while program == 1:	#Day begins as Bulb 1 is on 
		comm.barrier()
		switch = comm.bcast(switch,root = 0) #Look at the bulb 2

		if switch == 1:		#If Bulb 2 is on, start working on the function
			b = np.array([0,0,0,0])

			comm.barrier()				#These two sentences equate picking
			a = comm.bcast(a,root=0)	#up the phone for Inputs 
			
			t = a*rank*np.array([1,2,3,4])	#Rank-specific Transform

			comm.Reduce(t,b,op=MPI.SUM,root=0)	#Return Outputs
			switch = 0					#Equivalent to turning off Bulb 2

		comm.barrier()							#These two sentences equate to 
		program = comm.bcast(program,root = 0)	#asking if they can leave
		
	return