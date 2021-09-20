#python 3
#Carron Leopold 2020
#TOOL that you can need and easy to find on my git account :)

#use
import h5py
import numpy as np
from scipy import sparse
from hmmlearn import hmm

#view
from matplotlib import pyplot as plt
%matplotlib inline
# Seaborn for plotting and styling
import seaborn as sns

#code done for sparse matrix

def bin2d(Data,p,q):    
    """   
    Data = input matrix
    p,q rescaling factors
    Written for sparse
    This implementation is little bit tricky
    """
    n,m=np.shape(Data);
    s=(int(np.ceil(n/p)),int(np.ceil(m/q)))
    i,j,d = sparse.find(Data);
    i=np.int_(np.ceil(i/p))
    j=np.int_(np.ceil(j/q))
    M=sparse.csr_matrix((d,(i,j))) #go to csr to to indexing step
    if p==1:
        M=M[1:,:]
    if q==1:
        M=M[:,1:]
    return M

def SCN(D, max_iter = 10):
	"""
	Input : spare array * int
	Out  : SCN(D)
	Code version from Boost-HiC paper, update for sparse
	"""    
	# Iteration over max_iter    
	for i in range(max_iter):
		D /= np.maximum(1, D.sum(axis = 0))       
		D /= np.maximum(1, D.sum(axis = 1))    
	return (D + D.T)/2 # To make matrix symetric again   


def filteramat(Hicmat,Filterextremum=True,factor=1.5):
	"""
	in : a HiCmat without any transformation, factor of reduction
	out : the HiCmatreduce,thevector of his transformation
	"""
	Hicmatreduce=Hicmat
	#first step : filter empty bin
	sumHicmat=Hicmat.sum(axis = 0)
	segmenter1=sumHicmat>0
	A=np.where(segmenter1)
	Hicmatreduce=Hicmatreduce[A[1],:]
	Hicmatreduce=Hicmatreduce[:,A[1]]
	if Filterextremum:
		#second step : filter lower bin
		sumHicmat=np.sum(Hicmatreduce,0)
		msum=np.mean(sumHicmat)
		mstd=np.std(sumHicmat)
		mini = msum-mstd*factor
		maxi = msum+mstd*factor
		#Make the bolean condition
		newcond=mini < sumHicmat
		newcond2=sumHicmat < maxi
		newcond=np.logical_and(newcond,newcond2)
		B=np.where(newcond)
		#Filter
		Hicmatreduce=Hicmatreduce[B[1],:]
		Hicmatreduce=Hicmatreduce[:,B[1]]
		segmenter1=A[1][B[1]] #Create the binsaved index
	return Hicmatreduce,segmenter1

def reversesegmentation1D(segmentedvec,segmenter,originalsize):
	"""
	in : a matrix, a segmenter, his original size
	type of segmenter : list of coord wich has been saved
	out : the matrix that is reverse segmented
	/!\segmented value will be return as -1 in the vector
	"""	
	reversedvec=np.zeros(originalsize)-1
	i=0
	Li=len(segmenter)
	while i<Li:
		reversedvec[segmenter[i]]=segmentedvec[i]
		i+=1
	return reversedvec

  

def observed_expected(OE):
	""" in : a matrix 
	out : the matrix with point divide by mean of the diagonal
	"""
	i=0
	j=0
	L=len(OE)
	while j<L:
		thediag=np.diag(OE,k=j)
		mtg=np.mean(thediag)
		while i<(L-j):
			v=OE[i,i+j]/mtg
			OE[i,i+j]=v
			OE[i+j,i]=v
			i+=1
		i=0
		j+=1
	return OE

def makecompartimentbyGaussianHMM(chrmat,N=2):
	"""
	generate a simple hmm to have compartiment on correlation map
	out : the model of the hmm, AND the states of the prediction at given N
	"""
	# Run Gaussian HMM
	chrmat[np.isnan(chrmat)]=0
	chrmat[np.isinf(chrmat)]=0
	chrmat=np.nan_to_num(chrmat)
	#print("fitting to HMM and decoding ...")
	# Make an HMM instance and execute fit
	model = hmm.GaussianHMM(n_components=N,n_iter=1000, covariance_type="diag").fit(chrmat) #n_iter has to be >100 after test
	#print(model)
	# Predict the optimal sequence of internal hidden state 
	hidden_states = model.predict(chrmat)
	#print("done")
	return model,hidden_states


  
#Example of command to view


#figure = plt.figure(dpi=100)
#sns.heatmap(ccunseg, vmin=-0.4, vmax=0.4,cmap='seismic')
#plt.plot(comp)



