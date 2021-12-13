#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@authors: pearlman, morais

"""
# ---------------------------------------------------------------------------------------------- Importing Modules

#region
import matplotlib
import matplotlib.pyplot as plt
matplotlib.pyplot.switch_backend('agg')
import getopt, sys, copy, os, h5py

from subprocess import call, check_call, check_output, Popen
import numpy as np
import filterbank
from scipy import signal
from scipy.optimize import curve_fit
from scipy.stats import chisquare
#endregion

print("Done importing modules.")

# ---------------------------------------------------------------------- Defining Filterbank Reader, DM, & Dedispersion Functions

BLOCKSIZE = 1e6

""" Read the filterbank file into memory. Store the data in a dynamically
    accessible h5py file, stored in a binary .hdf5 file. """
def readFilterbank(inputFilename, logFile=""):
    
    if (logFile == ""):
        print("Reading filterbank file (%s)...\n" % inputFilename)
    else:
        logFile.write("Reading filterbank file (%s)...\n\n" % inputFilename)
    
    fb = filterbank.FilterbankFile(inputFilename)
    
    inputHeader = copy.deepcopy(fb.header)
    inputNbits = fb.nbits
    
    totalChans = fb.nchans
    nchans = np.arange(0, fb.nchans-1, 1) # Top of the band is index 0.
    
    freqs = fb.frequencies
    
    startbin = 0
    endbin = fb.nspec
    nspec = np.subtract(endbin, startbin)
    
    nblocks = int(np.divide(nspec, BLOCKSIZE))
    remainder = nspec % BLOCKSIZE
    totalBlocks = nblocks
    
    if (remainder):
        totalBlocks = nblocks + 1
    
    
    
    h5pyFile = h5py.File("%s.hdf5" % inputFilename, "w")
    spectraData = h5pyFile.create_dataset("data", (totalChans, nspec), dtype="float32")
    
    
    iblock = 0
    for iblock in np.arange(0, nblocks, 1):
        
        progress = np.multiply(np.divide(iblock + 1.0, totalBlocks), 100.0)
        
        if (logFile == ""):
            sys.stdout.write("Reading... [%3.2f%%]\r" % progress)
            sys.stdout.flush()
        else:
            logFile.write("Reading... [%3.2f%%]\n" % progress)
        
        
        lobin = int(np.add(np.multiply(iblock, BLOCKSIZE), startbin))
        hibin = int(np.add(lobin, BLOCKSIZE))
        
        spectra = fb.get_spectra(lobin, hibin)
        
        
        for ichan in np.arange(0, totalChans, 1):
            
            spectraData[ichan, lobin:hibin] = spectra[:,ichan]
        
    
    if (remainder):
        
        progress = np.multiply(np.divide(iblock + 2.0, totalBlocks), 100.0)
        
        if (logFile == ""):
            sys.stdout.write("Reading... [%3.2f%%]\r" % progress)
            sys.stdout.flush()
        else:
            logFile.write("Reading... [%3.2f%%]\n" % progress)
        
        
        lobin = int(np.subtract(endbin, remainder))
        hibin = int(endbin)
        
        spectra = fb.get_spectra(lobin, hibin)
        
        for ichan in np.arange(0, totalChans, 1):
            
            spectraData[ichan, lobin:hibin] = spectra[:,ichan]
    
    if (logFile == ""):
        print("\n")
    else:
        logFile.write("\n")
    
    return spectraData, inputHeader, inputNbits, h5pyFile;

print("Done defining filterbank function.")

""" Calculate the DM constant. """
def dm_constant():
    
    dm_constant = 1.0 / 0.000241
    
    return dm_constant

print("Done defining DM function.")

def dedisperse(spectraData, dm=87.757, freq_hi=800.0, tsamp=0.00004096, foff=0.390625):

    # DM correction for each channel (de-disperse with respect to the top
    # of the band).
    # Delay the higher frequency channels until the lowest frequency one
    # has arrived.
    # Crude de-dispersion algorithm.
    freq_lo = freq_hi
    
    # spectraData_dedispersed = copy.deepcopy(spectraData)    Was taking too much RAM ~2.18GiB
    spectraData_dedispersed = np.asarray(spectraData)
           
    for channelIndex in range(0, np.shape(spectraData_dedispersed)[0]):
        tau = np.multiply(dm_constant(), np.multiply(dm, np.subtract(np.divide(1.0, np.power(freq_lo, 2.0)), np.divide(1.0, np.power(freq_hi, 2.0)))))
                
        # Could provide a more sophisticated way of doing this...
        # We're just shifting an integer number of samples, and cycling back to the front.
        # Better way could be to plot each channel time-series and shift in time, one channel at a time.
        num_samples = int(np.round(np.divide(-tau, tsamp)))
        spectraData_dedispersed[channelIndex, :] = np.roll(spectraData_dedispersed[channelIndex, :], num_samples, axis=0)
                
        if (channelIndex < (np.shape(spectraData)[0] - 1)):
            freq_lo = np.subtract(freq_lo, foff)
     
    return spectraData_dedispersed
     
print("Done defining dedispersion function.")


# -------------------------------------------------------------------------------------------- Plotting Filterbank File

#Note: Temporarily working with hdf5 file for faster compiling

#region
#inputFilename = "/psr_scratch/jeffmorais/processed_data/m81r_59327_burst_data/J0957+68_cand_59327_pow_extract_32bit_bp0.5_bavg0.1_burst1.fil"

#spectraData, inputHeader, inputNbits, h5pyFile = readFilterbank(inputFilename)

#print("Input header: ", inputHeader)
#print("Shape of data: ", np.shape(spectraData))
#print("Spectra Data: ", spectraData[:, :])
#endregion

# ----------------------- Opening file

filename = "J0957+68_cand_59327_pow_extract_32bit_bp0.5_bavg0.1_burst1.fil.hdf5"

with h5py.File(filename, "r") as f:

    #List all groups
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]

    #Getting the data
    spectraData = list(f[a_group_key])

#spectraDataPlot = copy.deepcopy(spectraData[:, :]) #takes too much memory
spectraDataPlot = np.array(spectraData)

#Manually Inputting Header
inputHeader = {'telescope_id': 20, 'machine_id': 20, 'fch1': 800.0, 'data_type': 1, 'nchans': 1024, 'ibeam': 1, 'tsamp': 4.096e-05, 'nbeams': 1, 'foff': -0.390625, 'src_raj': 95754.70436096191, 'src_dej': 684901.0125732422, 'tstart': 59327.15853905457, 'az_start': 0.0, 'source_name': 'frb_20200120e', 'nifs': 1, 'rawdatafile': 'J0957+68_cand_59327_pow_extract_32bit_bp0.5_bavg0.1.fil', 'za_start': 0.0, 'nbits': 32}

#Incoherent dedispersion
spectraData_dedispersed = dedisperse(spectraDataPlot, dm=87.757, freq_hi=inputHeader["fch1"], tsamp=inputHeader["tsamp"], foff=np.abs(inputHeader["foff"]))

#h5pyFile.close()

print("Done dedispersing.")

# ------------------------ Plotting

xmin_plot = 14.99
xmax_plot = 15.01

#region
# Dedispersed Time-Series
#fig1 = plt.figure()
#ax1 = fig1.gca()

#ax1.plot(np.linspace(0,30,np.shape(spectraData_dedispersed)[1]),np.mean(spectraData_dedispersed,axis=0))
#ax1.axis(xmin=xmin_plot,xmax=xmax_plot)

#plt.show(block=True)
#plt.savefig('dedispersedTimeSeries.png')
#print("Dedispersed time-series plot finished.")

# Dedispersed Dynamic Spectrum
#fig2 = plt.figure()
#ax2 = fig2.gca()

#ax2.imshow(spectraData_dedispersed, aspect="auto", origin="upper", interpolation="none", extent=(0.0, 30.0, 400.0, 800.0))
#ax2.axis(xmin=xmin_plot,xmax=xmax_plot)

#plt.show(block=True)
#plt.savefig('dedispersedDynamicSpectrum.png')
#print("Dedispersed dynamic spectrum plot finished.")

# Frequency spectrum (REDACTED)
#fig3 = plt.figure()
#ax3 = fig3.gca()

print("Slicing the indices of the dynamic spectra: ", int(np.round(xmin_plot/inputHeader["tsamp"])), int(np.round((xmax_plot/inputHeader["tsamp"] + 1))))
spectraData_dedispersed_zoom = spectraData_dedispersed[:, int(np.round(xmin_plot/inputHeader["tsamp"])) : int(np.round((xmax_plot/inputHeader["tsamp"] + 1)))]
freq_spectrum = np.mean(spectraData_dedispersed_zoom, axis=1)

#ax3.plot(np.arange(0, len(freq_spectrum), 1), freq_spectrum, "b-")

#plt.show(block=True)
#plt.savefig('frequencySpectrum.png')
#print("Frequency spectrum plot finished.")
#endregion

print("This is the dedispersed, zoomed data:", spectraData_dedispersed_zoom)

# ---------------------------------------------------------------------------------------------- Spectral Analysis

# ------------------------- Removing all channels whose frequency spectrum amplitude > 3sigma
from astropy.stats import mad_std, sigma_clipped_stats
stats = sigma_clipped_stats(freq_spectrum, sigma_lower=3, sigma_upper=3) #mean, median, std

print("Generating bad channel list:")

threeSigma = 3*stats[-1]
removedChannels = []
removedChannelAmplitudes = []

for i in range(len(freq_spectrum)):
	if abs(freq_spectrum[i]) > threeSigma:
	#	#print("Channel: " + str(i))
	#	print("Frequency: " + str(freq_spectrum[i]))
		removedChannels.append(i)
		removedChannelAmplitudes.append(freq_spectrum[i])
		freq_spectrum[i] = 0
		spectraData_dedispersed[i-1]=0

#-------------------------- A Note on Variables

#Let freq_spectrum[i] = f(γ), where γ represents one of the channels. To get this in terms of frequency, we use the relation ν = 800 - 0.390625(γ) (ex: 800 - 0.390625*1024 = 400.390625)
#Inverting the relation gives γ = γ(ν), giving: f(γ) = f( (800- ν)/(0.390625) )

#-------------------------- Secondary pass of bad channels: Auto

#region
#Procedure: Zero points along the bright DM vertical bar, then average over each channel and zero channels above a threshold

#tempList = []
#tempSpectra = np.copy(spectraData_dedispersed)

#Zeroing DM bar
#for i in range(len(spectraData_dedispersed.T)):
#	tempList.append(np.mean(spectraData_dedispersed.T[i]))

#maximum = max(tempList) #Finding max of the data points, aka the bright burst vertical line

#for i in range(len(spectraData_dedispersed.T)):
#	if np.mean(spectraData_dedispersed.T[i]) > 0.1*maximum: #Removing data points with >10% amplitude of max
		#print(i)
#		tempSpectra.T[i] = 0 
		#spectraData_dedispersed.T[i] = 0

#Zeroing channels based on average threshold
#for i in range(len(freq_spectrum)):
#	if abs(np.mean(tempSpectra[i])) > 1.0e-1: #Removing channels with average over a threshold
#		spectraData_dedispersed[i-1] = 0
#endregion

#-------------------------------- Secondary pass of bad channels: Manual

def channel(v): # γ = γ(ν)
	return int(round( (800 - v)*(0.390625)**(-1) ))

for v in range(710,751):
	spectraData_dedispersed[channel(v)-1]=0

for v in range(620,641):
        spectraData_dedispersed[channel(v)-1]=0


#---------------------------- Finding bad data channels for manual pass

#region
#--- Dictionary of freq/amplitude

#for i in range(len(spectraData_dedispersed)):
#	channelFrequencyMeans.append(np.mean(spectraData_dedispersed[i]))
	
#for i in range(len(spectraData_dedispersed)):
#	channelDict[i+1] = channelFrequencyMeans[i]


#--- Plotting spectrum & array

#print("Test plot to find bad channels: ")

#plt.plot(indexChannels, channelFrequencyMeans)
#plt.savefig('badChannels.png')
#endregion

#---------------------------- Showing removed channels + plotting clean data

#---- Bad channels

removedChannelArray = np.array([removedChannels,removedChannelAmplitudes])
print("Removed channels & amplitudes:", removedChannelArray)

#---- Clean data plots

print("Cleaned up frequency range:")

plt.plot(np.arange(0, len(freq_spectrum), 1), freq_spectrum, "b-")

plt.show(block=True)
plt.savefig('cleanedUpFreq.png')
print("Cleaned up frequency range plot finished.")

print("Spectra shape: ", np.shape(spectraData_dedispersed))
print("Full spectra array: ", spectraData_dedispersed)
print("First spectra row (frequency channel): ", spectraData_dedispersed[0])
print("First spectra column (time sample):", spectraData_dedispersed.T[0])

plt.imshow(np.asarray(spectraData_dedispersed), aspect="auto", origin="upper", interpolation="none", extent=(0.0, 30.0, 400.0, 800.0)) #testing np.asarray in first argument for memory errors
plt.axis(xmin=xmin_plot,xmax=xmax_plot)

plt.show(block=True)
plt.savefig('cutSpectra.png')
print("Cut dynamic spectra plot finished.")

#------------------------------------------------------------------------------------------------- Auto Correlation Functions

#------------------ Defining ACF fitting functions

def lorentzian(deltaFreq, m, nuDC):
	return (m**2) * ((nuDC)**2 + (deltaFreq)**2)**(-1)

def scintillationWidth(v, A, alpha):
	return A * ( np.asarray(v) * (600)**(-1) )**(alpha)

#------------------ Defining ACFs: Sum Method

def ACF_Sum(fullFreqRange):

	#Preliminary list of sWidth values & frequency values to fit sWidth
	sWidth = []
	freqValues = []

	#Iterating over given frequency pair ranges
	fullFreqRange = np.array(fullFreqRange)
	for i in fullFreqRange: 

		#Defining variables
		low = int(i[0]) #lower frequency cutoff
		high = int(i[1]) #higher frequency cutoff
		S = freq_spectrum[channel(high)-1:channel(low)] #split S(v) array, interchanged slicing order since low > high for channels
		deltaNu = np.arange(1,100+1,1) #different values of Δν
		N = len(S)
		V = len(deltaNu)
		inverseN = N**(-1)

		#Performing ACF for all Δν
		ACF = np.zeros(V)
		for j in range(V): #iterating over Δν
			for v in range(N): #iterating over v
				ACF[j] += inverseN * S[v]*S[(v+j)%N] #the modulo % is to wrap the index around when it gets out of bounds

		#Plotting & saving ACF
		plt.plot(deltaNu , ACF)
		plt.xlim([0,100])

		#Fitting ACFs to a Lorentzian & plotting fits
		popt, pcov = curve_fit(lorentzian, deltaNu, ACF)
		sWidth.append(popt[1])
		freqValues.append(np.mean([low,high]))
		plt.plot(deltaNu, lorentzian(deltaNu, *popt))

	#Fitting sWidth to the scintillation width function
	poptsW, pcovsW = curve_fit(scintillationWidth, freqValues, sWidth)

	plt.savefig('ACF_Sum; Fit.png')

	#Plotting sWidth and fit
	plt.figure()
	plt.scatter(freqValues, sWidth)
	print("FreqValues: ", freqValues) #debugging
	print("sWidth: ", sWidth) #debugging
	plt.plot(freqValues, scintillationWidth(freqValues, *poptsW))
	plt.savefig('ACF_Sum; sWFit.png')

	#Print statements
	print("ACF sum plot finished.")
	print("Alpha : ", poptsW[1])
	print("A : ", poptsW[0])

	return None

#------------------ Defining ACFs: Numpy Correlation Method

def ACF_Numpy(fullFreqRange):

	#Iterating over given frequency pair ranges
	fullFreqRange = np.array(fullFreqRange)
	for i in fullFreqRange: 

		#Defining variables
		low = int(i[0]) #lower frequency cutoff
		high = int(i[1]) #higher frequency cutoff
		S = freq_spectrum[channel(high)-1:channel(low)] #split S(v) array, interchanged slicing order since low > high for channels

		#Performing ACF for all Δν
		corr = np.correlate(S,S, mode="full")

		#Plotting & saving ACF
		plt.plot(np.arange(0,len(corr),1) , corr )
		# plt.xlim([0,100])

	plt.savefig('ACF_Numpy.png')
	print("ACF numpy plot finished.")

	return None

#------------------ Defining ACFs: Fast Fourier Transform Method

def ACF_FFT(fullFreqRange):

	#Iterating over given frequency pair ranges
	fullFreqRange = np.array(fullFreqRange)
	for i in fullFreqRange: 

		#Defining variables
		low = int(i[0]) #lower frequency cutoff
		high = int(i[1]) #higher frequency cutoff
		S = freq_spectrum[channel(high)-1:channel(low)] #split S(v) array, interchanged slicing order since low > high for channels
		N = len(S)
		inverseN = N**(-1)

		#Performing ACF FFT for all Δν
		fft = np.fft.fft(S)
		autocorr = np.fft.ifft(fft * np.conj(fft))

		#Plotting & saving ACF
		plt.plot( np.arange(1, N+1, 1) , inverseN*autocorr )
		plt.xlim([0,100])

	plt.savefig('ACF_FFT.png')
	print("ACF FFT plot finished.")

	return None

#------------------ Experimental Method: Matplotlib Acorr Method

def ACF_MAM(fullFreqRange):

	#Iterating over given frequency pair ranges
	fullFreqRange = np.array(fullFreqRange)
	for i in fullFreqRange:

		#Defining variables
		low = int(i[0]) #lower frequency cutoff
		high = int(i[1]) #higher frequency cutoff
		S = freq_spectrum[channel(high)-1:channel(low)] #split S(v) array, interchanged slicing order since low > high for channels
		
		#Plotting & saving ACF
		plt.acorr(S)
		plt.ylim([0,0.2])

	# plt.savefig('TEST-MAM.png')
	print("ACF MAM plot finished.")

	return None

#------------------ Experimental Method: Scipy Correlate Method

def ACF_SCM(fullFreqRange):

	#Iterating over given frequency pair ranges
	fullFreqRange = np.array(fullFreqRange)
	for i in fullFreqRange:

		#Defining variables
		low = int(i[0]) #lower frequency cutoff
		high = int(i[1]) #higher frequency cutoff
		S = freq_spectrum[channel(high)-1:channel(low)] #split S(v) array, interchanged slicing order since low > high for channels
		N = len(S)
		inverseN = N**(-1)

		#Performing ACF for all Δν
		corr = signal.correlate(S,S)

		#Plotting & saving ACF
		plt.plot(np.arange(0,len(corr),1) , inverseN*corr )
		# plt.xlim([0,100])

	# plt.savefig('TEST_SCM.png')
	print("ACF SCM plot finished.")

	return None

#------------------- Computing, plotting, & fitting ACFs

plt.figure()
ACF_Sum(
	[[760,800-1],[680,720+1],[640,680+1],[600,640+1],[560,600+1],[520,560+1],[480,520+1],[440,580+1]]
)
plt.figure()
ACF_Numpy(
	[[760,800-1],[680,720+1],[640,680+1],[600,640+1],[560,600+1],[520,560+1],[480,520+1],[440,580+1]]
)
plt.figure()
ACF_FFT(
	[[760,800-1],[680,720+1],[640,680+1],[600,640+1],[560,600+1],[520,560+1],[480,520+1],[440,580+1]]
)
plt.figure()
ACF_MAM(
	[[760,800-1],[680,720+1],[640,680+1],[600,640+1],[560,600+1],[520,560+1],[480,520+1],[440,580+1]]
)
plt.figure()
ACF_SCM(
	[[760,800-1],[680,720+1],[640,680+1],[600,640+1],[560,600+1],[520,560+1],[480,520+1],[440,580+1]]
)

#------------------------------------------------------------------------------------------------------------ Objectives

#region
# Item 1:

# Use data in spectraData_dedispersed_zoom to determine which channel indices are bad. Compare to moving average/median, determine +/- 3 sigma values, and identify channel indices that exceed these threshol.
# Zero channels that exceed the threshold.
# Replot the data.
# Write it as a function.
# Repeat over and over until you're satisfied.
# If bad channels still exist, then manually remove them.
# IMPORTANT: Save the list of channels that you zero (both from the algorithm, and also from manual removal).


# Item 2:

# Make clean plots of the other bursts. Remember to save the mask, and +/- 50 ms around the burst. Use np.savez(...) to save it as a .npz file. Then, you can load the cleaned, shortened data using np.load(...)


# Item 3 (HIGH PRIORITY):

# Move to scintillation analysis using clean burst data. Read about decorrelation bandwidth in FRBs. This will be calculated using an ACF (auto-correlation function). Start experimenting yourself with simple code examples.
# Eventually, we'll write algorithms/functions that will produce reliable results that we can run on many different bursts. 
#endregion


#---------------------------------------------------------------------------------------------------------- Code Graveyard

#region
#---------------------------------- Deprecated ACF Functions

#------------------ Defining ACFs: Sum Method

# def ACF_Sum_Single(lowerFreqRange, upperFreqRange):

# 	#Defining variables
# 	low = int(lowerFreqRange) #lower frequency cutoff
# 	high = int(upperFreqRange) #higher frequency cutoff
# 	S = freq_spectrum[channel(high)-1:channel(low)] #split S(v) array, interchanged slicing order since low > high for channels
# 	deltaNu = np.arange(1,100+1,1) #different values of Δν
# 	N = len(S)
# 	V = len(deltaNu)
# 	inverseN = N**(-1)

# 	#Performing ACF for all Δν
# 	ACF = np.zeros(V)
# 	for j in range(V): #iterating over Δν
# 		for v in range(N): #iterating over v
# 			ACF[j] += inverseN * S[v]*S[(v+j)%N] #the modulo % is to wrap the index around if it gets out of bounds

# 	#Plotting & saving ACF
# 	plt.figure()
# 	plt.plot(deltaNu , ACF)
# 	plt.xlim([0,100])
# 	plt.savefig('ACF_Sum.png')
# 	plt.show(block=True)
# 	print("ACF sum plot finished.")

# 	return None

# def ACF_Sum(deltaFreq, lowerFreqRange, upperFreqRange):
# 	low = int(lowerFreqRange) #lower frequency cutoff
# 	high = int(upperFreqRange) #higher frequency cutoff
# 	spectrumArray = freq_spectrum[channel(high)-1:channel(low)] #split S(nu) array, interchanged slicing order since low > high for channels
# 	N = len(spectrumArray)
# 	inverseN = N**(-1)

# 	sumAmplitude = np.zeros(len(deltaFreq)) #r(deltaNu) with empty entries
# 	for j in range(len(deltaFreq)):	#iterating over values of deltaFreq Δν
# 		for i in range(N): #iterating over values of frequency ν
# 			sumAmplitude += spectrumArray[i]*spectrumArray[(i+deltaFreq[j])%N] #the modulo % is to wrap the index around if it gets out of bounds

# 	return (inverseN*sumAmplitude)

#------------------ Defining ACFs: Numpy Correlation Method

# def ACF_Numpy(lowerFreqRange, upperFreqRange):

# 	#Defining variables
# 	low = int(lowerFreqRange) #lower frequency cutoff
# 	high = int(upperFreqRange) #higher frequency cutoff
# 	S = freq_spectrum[channel(high)-1:channel(low)] #split S(v) array, interchanged slicing order since low > high for channels

# 	#Performing ACF for all Δν
# 	corr = np.correlate(S,S, mode="full")

# 	#Plotting & saving ACF
# 	plt.figure()
# 	plt.plot(np.arange(0,len(corr),1) , corr )
# 	# plt.xlim([0,100])
# 	plt.savefig('ACF_Numpy.png')
# 	plt.show(block=True)
# 	print("ACF numpy plot finished.")

# 	return None

# #------------------ Defining ACFs: Fast Fourier Transform Method

# def ACF_FFT(lowerFreqRange, upperFreqRange):

# 	#Defining variables
# 	low = int(lowerFreqRange) #lower frequency cutoff
# 	high = int(upperFreqRange) #higher frequency cutoff
# 	S = freq_spectrum[channel(high)-1:channel(low)] #split S(v) array, interchanged slicing order since low > high for channels

# 	#Performing ACF FFT for all Δν
# 	fft = np.fft.fft(S)
# 	autocorr = np.fft.ifft(fft * np.conj(fft))

# 	#Plotting & saving ACF
# 	plt.figure()
# 	plt.plot( np.arange(1, N+1, 1) , inverseN*autocorr )
# 	plt.xlim([0,100])
# 	plt.savefig('ACF_FFT.png')
# 	plt.show(block=True)
# 	print("ACF FFT plot finished.")

# 	return None

# #------------------ Experimental Method: Matplotlib Acorr Method

# def ACF_MAM(lowerFreqRange, upperFreqRange):

# 	#Defining variables
# 	llow = int(760) #lower frequency cutoff
# 	hhigh = int(800-1) #higher frequency cutoff
# 	SS = freq_spectrum[channel(hhigh)-1:channel(llow)] #split S(v) array, interchanged slicing order since low > high for channels
	
# 	#Plotting & saving ACF
# 	plt.figure()
# 	plt.acorr(SS)
# 	plt.ylim([0,0.2])
# 	plt.savefig('TEST-MAM.png')
# 	plt.show(block=True)
# 	print("ACF MAM plot finished.")

# 	return None

# #------------------ Experimental Method: Scipy Correlate Method

# def ACF_SCM(lowerFreqRange, upperFreqRange):

# 	#Defining variables
# 	llow = int(760) #lower frequency cutoff
# 	hhigh = int(800-1) #higher frequency cutoff
# 	SS = freq_spectrum[channel(hhigh)-1:channel(llow)] #split S(v) array, interchanged slicing order since low > high for channels
# 	NN = len(SS)
# 	inverseN = NN**(-1)

# 	#Performing ACF for all Δν
# 	corr = signal.correlate(S,S)

# 	#Performing ACF for all Δν
# 	corr = np.correlate(S,S, mode="full")

# 	#Plotting & saving ACF
# 	plt.figure()
# 	plt.plot(np.arange(0,len(corr),1) , corr )
# 	# plt.xlim([0,100])
# 	plt.savefig('TEST_SCM.png')
# 	plt.show(block=True)
# 	print("ACF SCM plot finished.")

	# return None

#------------------ Defining ACFs: Numpy Correlation Method

# def ACF_Numpy(deltaFreq, lowerFreqRange, upperFreqRange):
# 	low = int(lowerFreqRange) #lower frequency cutoff
# 	high = int(upperFreqRange) #higher frequency cutoff
# 	spectrumArray = freq_spectrum[channel(high)-1:channel(low)] #split S(nu) array, interchanged slicing order since low > high for channels
# 	N = len(spectrumArray)
# 	inverseN = N**(-1)

# 	sumAmplitude = np.correlate(spectrumArray,spectrumArray, mode="full")[deltaFreq%len(np.correlate(spectrumArray,spectrumArray))]

# 	return (inverseN*sumAmplitude)

#------------------ Defining ACFs: Fast Fourier Transform Method

# def ACF_FFT(lowerFreqRange, upperFreqRange):
# 	low = int(lowerFreqRange) #lower frequency cutoff
# 	high = int(upperFreqRange) #higher frequency cutoff
# 	spectrumArray = freq_spectrum[channel(high)-1:channel(low)] #split S(nu) array, interchanged slicing order since low > high for channels
# 	N = len(spectrumArray)

# 	fft = np.fft.fft(spectrumArray)
# 	autocorr = np.fft.ifft(fft * np.conj(fft))
# 	return (autocorr, N)


# #----------------- Plotting ACFs

# #Sum method
# plt.figure()
# frequencyRange = np.arange(1,100+1,1)
# plt.plot( frequencyRange , ACF_Sum( frequencyRange ,760, 800-1) ) #adding -1 to avoid divergence
# # plt.plot( frequencyRange , ACF_Sum( frequencyRange ,720, 760+1) ) #Eve didn't include this range
# plt.plot( frequencyRange , ACF_Sum( frequencyRange ,680, 720+1) )
# plt.plot( frequencyRange , ACF_Sum( frequencyRange ,640, 680+1) )
# plt.plot( frequencyRange , ACF_Sum( frequencyRange ,600, 640+1) )
# plt.plot( frequencyRange , ACF_Sum( frequencyRange ,560, 600+1) )
# plt.plot( frequencyRange , ACF_Sum( frequencyRange ,520, 560+1) )
# plt.plot( frequencyRange , ACF_Sum( frequencyRange ,480, 520+1) )
# plt.plot( frequencyRange , ACF_Sum( frequencyRange ,440, 580+1) )
# plt.xlim([0,100])
# plt.show(block=True)
# plt.savefig('ACF_Sum.png')
# print("ACF sum plot finished.")

# #Numpy method
# plt.figure()
# frequencyRange = np.arange(1,100+1,1)
# plt.plot( frequencyRange , ACF_Numpy( frequencyRange ,760, 800-1) ) #adding -1 to avoid divergence
# # plt.plot( frequencyRange , ACF_Sum( frequencyRange ,720, 760+1) ) #Eve didn't include this range
# plt.plot( frequencyRange , ACF_Numpy( frequencyRange ,680, 720+1) )
# plt.plot( frequencyRange , ACF_Numpy( frequencyRange ,640, 680+1) )
# plt.plot( frequencyRange , ACF_Numpy( frequencyRange ,600, 640+1) )
# plt.plot( frequencyRange , ACF_Numpy( frequencyRange ,560, 600+1) )
# plt.plot( frequencyRange , ACF_Numpy( frequencyRange ,520, 560+1) )
# plt.plot( frequencyRange , ACF_Numpy( frequencyRange ,480, 520+1) )
# plt.plot( frequencyRange , ACF_Numpy( frequencyRange ,440, 580+1) )
# plt.xlim([0,100])
# plt.show(block=True)
# plt.savefig('ACF_Numpy.png')
# print("ACF numpy plot finished.")

# #FFT method
# plt.figure()
# plt.plot( np.arange(1, ACF_FFT(760,800-1)[1]+1, 1) , ACF_FFT(760, 800-1)[0] ) #adding -1 to avoid divergence
# # plt.plot( np.arange(1, ACF_FFT(720,760)[1]+1, 1) , ACF_FFT(720, 760)[0] ) #Eve didn't include this range
# plt.plot( np.arange(1, ACF_FFT(680,720)[1]+1, 1) , ACF_FFT(680, 720)[0] )
# plt.plot( np.arange(1, ACF_FFT(640,680)[1]+1, 1) , ACF_FFT(640, 680)[0] )
# plt.plot( np.arange(1, ACF_FFT(600,640)[1]+1, 1) , ACF_FFT(600, 640)[0] )
# plt.plot( np.arange(1, ACF_FFT(560,600)[1]+1, 1) , ACF_FFT(560, 600)[0] )
# plt.plot( np.arange(1, ACF_FFT(520,560)[1]+1, 1) , ACF_FFT(520, 560)[0] )
# plt.plot( np.arange(1, ACF_FFT(480,520)[1]+1, 1) , ACF_FFT(480, 520)[0] )
# plt.plot( np.arange(1, ACF_FFT(440,480)[1]+1, 1) , ACF_FFT(440, 480)[0] )
# plt.xlim([0,100])
# plt.show(block=True)
# plt.savefig('ACF_FFT.png')
# print("ACF FFT plot finished.")


# #Experimental method: matplotlib acorr
# llow = int(760) #lower frequency cutoff
# hhigh = int(800-1) #higher frequency cutoff
# SS = freq_spectrum[channel(hhigh)-1:channel(llow)] #split S(nu) array, interchanged slicing order since low > high for channels
# plt.figure()
# plt.acorr(SS)
# plt.ylim([0,0.2])
# plt.show()
# plt.savefig('TEST-Matplot.png')

# #Experimental method: scipy correlate
# NN = len(SS)
# inverseN = NN**(-1)
# def corr(a,deltaNu):
# 	return ( (inverseN*signal.correlate(a,a))[deltaNu] )
# plt.figure()
# plt.plot(frequencyRange,corr(SS,frequencyRange))
# plt.ylim([0,0.2])
# plt.show()
# plt.savefig('TEST-Scipy.png')


# #----------------- Fitting ACFs

# #Note: Currently no measure of chisquared implemented

# scintillationWidthArray = []

# def lorentzian(deltaFreq, m, nuDC):
# 	return (m**2) * ((nuDC)**2 + (deltaFreq)**2)**(-1)

# def scintillationWidth(v, A, alpha):
# 	return A * ( v * (600)**(-1) )**(alpha)

# poptOne, pcovOne = curve_fit(lorentzian, np.arange(1, ACF_FFT(760,800-1)[1]+1, 1), ACF_FFT(760, 800-1)[0])
# scintillationWidthArray.append(poptOne[1])

# poptTwo, pcovTwo = curve_fit(lorentzian, np.arange(1, ACF_FFT(680,720)[1]+1, 1), ACF_FFT(680,720)[0])
# scintillationWidthArray.append(poptTwo[1])

# poptThree, pcovThree = curve_fit(lorentzian, np.arange(1, ACF_FFT(640,680)[1]+1, 1), ACF_FFT(640,680)[0])
# scintillationWidthArray.append(poptThree[1])

# poptFour, pcovFour = curve_fit(lorentzian, np.arange(1, ACF_FFT(600,640)[1]+1, 1), ACF_FFT(600,640)[0])
# scintillationWidthArray.append(poptFour[1])

# poptFive, pcovFive = curve_fit(lorentzian, np.arange(1, ACF_FFT(560,600)[1]+1, 1), ACF_FFT(560,600)[0])
# scintillationWidthArray.append(poptFive[1])

# poptSix, pcovSix = curve_fit(lorentzian, np.arange(1, ACF_FFT(520,560)[1]+1, 1), ACF_FFT(520,560)[0])
# scintillationWidthArray.append(poptSix[1])

# poptSeven, pcovSeven = curve_fit(lorentzian, np.arange(1, ACF_FFT(480,520)[1]+1, 1), ACF_FFT(480,520)[0])
# scintillationWidthArray.append(poptSeven[1])

# poptEight, pcovEight = curve_fit(lorentzian, np.arange(1, ACF_FFT(440,480)[1]+1, 1), ACF_FFT(440,480)[0])
# scintillationWidthArray.append(poptEight[1])

# print("scintillationWidthArray: ", scintillationWidthArray)

# #ACF sum/numpy method fitting (doesn't converge for constant functions)
# # poptOne, pcovOne = curve_fit(lorentzian, frequencyRange, ACF_Sum( frequencyRange ,760, 800-1))
# # scintillationWidthArray.append(poptOne[1])

# # poptTwo, pcovTwo = curve_fit(lorentzian, frequencyRange, ACF_Sum( frequencyRange ,680,720,))
# # scintillationWidthArray.append(poptTwo[1])

# # poptThree, pcovThree = curve_fit(lorentzian, frequencyRange, ACF_Sum( frequencyRange ,640, 680+1))
# # scintillationWidthArray.append(poptThree[1])

# # poptFour, pcovFour = curve_fit(lorentzian, frequencyRange, ACF_Sum( frequencyRange ,600, 640+1))
# # scintillationWidthArray.append(poptFour[1])

# # poptFive, pcovFive = curve_fit(lorentzian, frequencyRange, ACF_Sum( frequencyRange ,560, 600+1))
# # scintillationWidthArray.append(poptFive[1])

# # poptSix, pcovSix = curve_fit(lorentzian, frequencyRange, ACF_Sum( frequencyRange ,520, 560+1))
# # scintillationWidthArray.append(poptSix[1])

# # poptSeven, pcovSeven = curve_fit(lorentzian, frequencyRange, ACF_Sum( frequencyRange ,480, 520+1))
# # scintillationWidthArray.append(poptSeven[1])

# # poptEight, pcovEight = curve_fit(lorentzian, frequencyRange, ACF_Sum( frequencyRange ,440, 480+1))
# # scintillationWidthArray.append(poptEight[1])

# # print("scintillationWidthArray: ", scintillationWidthArray)

# #Now we will fit for the scintillationWidth (sW) at different frequencies. We computed sW for different frequency ranges, so for now
# #we will average over those ranges to give approximate frequency values and fit for sW using those values.

# avg1 = np.average(np.mean(np.arange(760, 800-1, 1)))
# avg2 = np.average(np.mean(np.arange(680, 720+1, 1)))
# avg3 = np.average(np.mean(np.arange(640, 680+1, 1)))
# avg4 = np.average(np.mean(np.arange(600, 640+1, 1)))
# avg5 = np.average(np.mean(np.arange(560, 600+1, 1)))
# avg6 = np.average(np.mean(np.arange(520, 560+1, 1)))
# avg7 = np.average(np.mean(np.arange(480, 520+1, 1)))
# avg8 = np.average(np.mean(np.arange(440, 480+1, 1)))
 
# frequencyValues = [avg1, avg2, avg3, avg4, avg5, avg6, avg7, avg8]

# #vectorizedScintillationWidth = np.vectorize(scintillationWidth) #To be able to pass an array in the argument of a function

# poptSW, pcovSW = curve_fit(scintillationWidth, frequencyValues, scintillationWidthArray, maxfev = 10000)

# print("Alpha : ", poptSW[1])

# #--------------- Test Plot to See Goodness of Fit

# plt.figure()
# plt.plot( np.arange(1, ACF_FFT(760,800-1)[1]+1, 1) , ACF_FFT(760, 800-1)[0] ) #adding -1 to avoid divergence
# plt.plot( np.arange(1, ACF_FFT(760,800-1)[1]+1, 1) , lorentzian(np.arange(1, ACF_FFT(760,800-1)[1]+1, 1), *poptOne) )
# plt.xlim([0,100])
# plt.show(block=True)
# plt.savefig('ACF_curveFit.png')
# print("ACF FFT test fit plot finished.")
#endregion