#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
# %%
S = [1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0]
deltaNu = np.arange(1,40+1,1) #1,40
N = len(S)
V = len(deltaNu)
inverseN = N**(-1)

# %%
#------------------ Defining ACFs: Sum Method

def ACF_Sum():
	ACF = np.zeros(V)
	for j in range(V): #different deltaNu
		for v in range(N):
			ACF[j] += inverseN * S[v]*S[(v+j)%N]

	plt.figure()
	plt.plot(deltaNu , ACF)
	plt.xlim([0,40])
 
	for i in range(3):
		plt.savefig("hi"+str(i))
  
	plt.show(block=True)
	
	print("ACF sum plot finished.")

	return None

#------------------ Defining ACFs: Numpy Correlation Method

def ACF_Numpy():
	corr = np.correlate(S,S, mode="full")
 
	plt.figure()
	plt.plot(np.arange(0,len(corr),1) , corr )
	# plt.xlim([0,40])
	plt.show(block=True)
	print("ACF numpy plot finished.")

	return None

#------------------ Defining ACFs: Fast Fourier Transform Method

def ACF_FFT():
	fft = np.fft.fft(S)
	autocorr = np.fft.ifft(fft * np.conj(fft))
 
	plt.figure()
	plt.plot( np.arange(1, N+1, 1) , inverseN*autocorr )
	plt.xlim([0,40])
	plt.show(block=True)
	print("ACF FFT plot finished.")
 
	return None

ACF_Sum()
ACF_Numpy()
ACF_FFT()

# %%
plt.plot(deltaNu , deltaNu)
for i in range(3):
		plt.savefig("hi"+str(i))

# %%
for i in [[1,2],[3,4]]:
    print("Hi :" + str(i))
# %%
