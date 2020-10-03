###################################################
# Code for ETH course "Numerical Methods for CSE"
# Author: H. Gratten
# Date: Oct 01 2020
# Pyhton script demonstrating the use of soundcompression by means of FFT
# Demo code for Chapter 4 of the course
# https://www.sam.math.ethz.ch/~grsam/NCSE20/NumCSE_Lecture_Document.pdf
##################################################
import wave
import numpy as np
from numpy.fft import fft, ifft
# pip install sounddevice
import sounddevice as sd 
from matplotlib import pyplot as plt

# Uses a stylesheet for plotting.
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": [ "Helvetica" ],
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "lines.linewidth": 1 })

def audioread( filename ):

	"""MATLAB-style reads audio into real-numbered sample matrix y (nsamples x nchannels) and 
	sampling rate Fs.
	"""
	
	wave_obj = wave.open( filename, 'rb' )	
	( nchannels, w, Fs, n, _, _ ) = wave_obj.getparams( )
	y = np.zeros(( n, nchannels ))
	for s in range( n ):
		sample = wave_obj.readframes( 1 )
		for ch in range( nchannels ):
			f = int.from_bytes( 
				sample[ ch * w : ch * w + w ], 
				byteorder = 'little',
				signed = True )
			f /= ( 2 **( 8 * w - 1 ))
			y[ s ][ ch ] = f

	wave_obj.close( )
	return (y, Fs)

def audiowrite( filename, y, Fs ):

	"""MATLAB-style writes audio with 16-bit precision (wave does not support any compression).
	"""

	if len( y.shape ) < 2:
		y = np.expand_dims( y, axis = 1 )

	wave_obj = wave.open( filename, 'wb' )
	n = y.shape[ 0 ]
	nchannels = y.shape[ 1 ]
	w = 2
	wave_obj.setparams(( nchannels, w, Fs, n, 'NONE', 'not compressed' ))
	for s in range( n ):
		for ch in range( nchannels ):
			b = int( round( float( y.real[ s ][ ch ] ) * ( 2 **( 8 * w - 1 ))))
			wave_obj.writeframesraw( b.to_bytes( w, byteorder = 'little', signed = True ))

	wave_obj.close( )

def sound( y, Fs, s = '' ):

	"""Plays sound matrix y (nsamples x nchannels) with playback rate Fs.
	"""

	name = 'sound \'{}\''.format( s )
	print( 'Playing {}...'.format( name ), end = '\r', flush = True )
	sd.play( y.real, Fs, blocking = True )
	print( 'Played {}.{}'.format( name, ' ' * 3 ))

# Reads real sound data.
y, Fs = audioread( 'hello.wav' )
# Converts sound matrix to mono.
y = y[ :, 1 ] 
n = len( y )
print( 'Read wav File: {} samples, rate = {}/s'.format( n, Fs ))
sound( y, Fs, 'hello.wav' )

# Plots sound signal.
plt.figure( )
plt.plot( np.array( range( n )) / Fs, y, 'r-' )
plt.title( 'sampled sound signal' )
plt.xlabel( 'time[s]' )
plt.ylabel( 'sound pressure' )
plt.grid( )
plt.savefig( './soundsignal.eps' )

# Retrieves Fourier coefficients with Fast Fourier transform.
coeff = fft( y )

# Plots all Fourier coefficients.
plt.figure( )
plt.plot( np.array( range( n )), abs( coeff ) ** 2, 'm-' )
plt.title( 'power spectrum of sound signal' )
plt.xlabel( 'index $k$ of Fourier coefficients $c_k$' )
plt.ylabel( '$|c_k|^2$' )
plt.grid( )
plt.axvline( n/2, c = 'k', ls = '--')
plt.text( n/2 + 1000, 31000, 'Nyquist frequency', ha = 'left' )
plt.savefig( './soundpower.eps' )

# Plots Fourier coefficients well below Nyquist frequency (one-sided).
plt.figure( )
plt.plot( np.array( range( 3000 )), abs( coeff[ :3000 ] ) ** 2, 'b-' )
plt.title( 'low frequency power spectrum' )
plt.xlabel( 'index $k$ of Fourier coefficients $c_k$' )
plt.ylabel( '$|c_k|^2$' )
plt.grid( )
plt.savefig( './soundlowpower.eps' )

# Plots filtered signals on a small subset of the time axis.
plt.figure( )
plt.title( 'sound filtering' )
plt.xlabel( 'time[s]' )
plt.ylabel( 'sound pressure' )

for cutoff in [ 5000, 3000, 1000 ]:
	
	# Applies low pass filtering (also pass symmetric counterpart of low frequencies).
	coeff_filter = np.zeros( n, dtype = np.complex64 )
	coeff_filter[ :cutoff ] = coeff[ :cutoff ]
	coeff_filter[ n-cutoff: ] = coeff[ n-cutoff: ]
	
	# Constructs a new sound signal from filtered Fourier coefficients.
	y_filter = ifft( coeff_filter )
	audiowrite( 'hello_cutoff{}.wav'.format( cutoff ), y_filter, Fs )
	sound( y_filter, Fs, 'hello.wav with cut-off frequency {}'.format( cutoff ))
	plt.plot( 
		np.array( np.arange( 30000, 32000 )) / Fs, 
		y_filter.real[ 30000:32000 ],
		label = 'cut-off = {}'.format( cutoff ))

plt.legend( )
plt.savefig( './soundfiltered.eps' )

def compress( c ):

	"""Applies a naive, straightforward amplitude-based compression to a Fourier coefficient c
	"""

	global nnz
	global power_cutoff
	if abs( c ) ** 2 > power_cutoff:
		nnz += 1 # Counts all nonzero entries.
		return c
	else:
	 return complex( 0 )

for pc in [ 50, 200, 1000 ]:

	# Resets globals to be used by compress function.
	power_cutoff = pc 
	nnz = 0

	# Applies the compress function element-wise.
	coeff_filter = np.vectorize( compress )( coeff ) 

	# Constructs a new sound signal from filtered Fourier coefficients.
	y_filter = ifft( coeff_filter )
	
	# The wave library does not support any compression for writing, so one just has to imagine how
	# a simple sparse representation based on Fourier coefficients would look like.
	audiowrite( 'hello_compressed{}.wav'.format( nnz ), y_filter, Fs )
	sound( y_filter, Fs, 'hello.wav with only the largest {}% of coefficients'.format( int( 100*nnz/n )))
