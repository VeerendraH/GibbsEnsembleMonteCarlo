import numpy as np


def initCubicGrid(N,L):
    coords = np.zeros((N,3))
    nCube = 2

    while (nCube**3 < N):
        nCube = nCube + 1

    index = np.zeros((3,))

    for part in range(0,N):
        coords[part,:] = (index+ 0.5*np.ones((3,)))*(L/nCube)
        index[0] = index[0] + 1
        if (index[0] == nCube):
            index[0] = 0
            index[1] = index[1] + 1
            if (index[1] == nCube):
                index[1] = 0
                index[2] = index[2] + 1
    return coords
def index_alloc(N,size):
	m = (N//size)*np.ones((1,size))
	for i in range(0,N%size):
		m[0,i] = m[0,i] + 1

	mbeg = np.tril(np.ones((size,size)),-1)@m.reshape([-1,1])
	mend = (np.hstack((-1*np.ones((size,1)),np.tril(np.ones((size,size)),0))))@(np.vstack(([1],m.reshape([-1,1])))) + 1
	
	return [(mbeg.astype(int)),(mend.astype(int))]
def LJ_Energy(coords,L,cut,size):
	comm.bcast
	energy = 0
	nPart = coords.shape[0]
	rcut2 = cut[0]
	phicut = cut[1]
	mbeg,mend = index_alloc(nPart-1,size)
	print(mbeg,mend,'Hi')
	for partA in range(mbeg[rank],mend[rank]):
		for partB in range(partA+1,nPart):
			dr = coords[partA,:] - coords[partB,:]
			dr = distPBC3D(dr,L)
			dr2 = np.sum(dr**2)
			if dr2<rcut2:
				invDr6 = (1/dr2)**3
				energy = energy + (invDr6*(invDr6-1)) - (phicut)

	energy = 4*energy
	comm.Barrier()
	comm.Reduce(en,energy,op=MPI.SUM,root = 0)
	return en

def LJ_EnergyIncrease(coords,trialPos,part,L,rcut):
	deltaE = 0
	nPart = coords.shape[0]
	for otherPart in range(0,nPart):
		if (otherPart==part):
			continue

		drNew = coords[otherPart,:] - trialPos
		drOld = coords[otherPart,:] - coords[part,:]

		drNew = distPBC3D(drNew,L)
		drOld = distPBC3D(drOld,L)

		dr2_New = sum(drNew**2)
		dr2_Old = sum(drOld**2)

		if dr2_New < rcut**2:

			invDr6_New = 1/(dr2_New**3)
			invDr6_Old = 1/(dr2_Old**3)

			eNew = invDr6_New*(invDr6_New-1)
			eOld = invDr6_Old*(invDr6_Old-1)

			deltaE = deltaE + eNew - eOld
	return 4*deltaE

def distPBC3D(vec,L):
    hL = L/2
    for dim in range(0,2):
        if (vec[dim] > hL):
            vec[dim] = vec[dim] - L
        elif (vec[dim] < -hL):
            vec[dim] = vec[dim] + L
    return vec

def PBC3D(vec,L):
    for dim in range(0,2):
        if (vec[dim] > L):
            vec[dim] = vec[dim] - L
        elif (vec[dim] < 0):
            vec[dim] = vec[dim] + L
    return vec

def Write_to_file(step,data,filename):
	'''Writes the energy to a file.'''
	'''Order =Step L1 L2 N1 N2 e1 e2 d1 d2 coords1 coords2'''
	coords = np.vstack((data[0],data[1])).flatten() 
	with open(filename, 'a') as f:
		f.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}\n'.format(step,data[2],data[3],data[4],data[5],data[6],data[7],data[4]/(data[2]**3),data[5]/(data[3]**3),','.join([str(num) for num in coords])))
	return