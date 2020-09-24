import wave
import numpy as np
from numpy.fft import fft, ifft
import sounddevice as sd
from matplotlib import pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": [ "Helvetica" ],
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "lines.linewidth": 1 })

def audioread( filename ):
	
	wave_obj = wave.open( filename, 'rb' )	
	( nchannels, w, Fs, n, _, _ ) = wave_obj.getparams( )
	y = np.zeros(( n, nchannels ))
	for s in range( n ):
		sample = wave_obj.readframes( 1 )
		for ch in range( nchannels ):
			f = int.from_bytes( sample[ ch * w : ch * w + w ], byteorder = 'little', signed = True )/( 2 **( 8 * w - 1 ))
			y[ s ][ ch ] = f

	wave_obj.close( )
	return (y, Fs)

#write with 16-bit precision
def audiowrite( filename, y, Fs ):

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

	print( 'Playing sound {}...'.format( s ), end = '\r', flush = True )
	sd.play( y.real, Fs, blocking = True )
	print( ' ' * 100, end = '\r', flush = True )

#Real sound data
y, Fs = audioread( 'hello.wav' )
y = y[ :, 1 ] #Convert to mono
n = len( y )
print( 'Read wav File: {} samples, rate = {}/s'.format( n, Fs ))
sound( y, Fs, 'hello.wav' )

plt.figure( )
plt.plot( np.array( range( n )) / Fs, y, 'r-' )
plt.title( 'sampled sound signal' )
plt.xlabel( 'time[s]' )
plt.ylabel( 'sound pressure' )
plt.grid( )
plt.savefig( './soundsignal.eps' )

#Fast fourier transform
coeff = fft( y )

#Plot all Fourier coefficients
plt.figure( )
plt.plot( np.array( range( n )), abs( coeff ) ** 2, 'm-' )
plt.title( 'power spectrum of sound signal' )
plt.xlabel( 'index $k$ of Fourier coefficients $c_k$' )
plt.ylabel( '$|c_k|^2$' )
plt.grid( )
plt.axvline( n/2, c = 'k', ls = '--')
plt.text( n/2 + 1000, 31000, 'Nyquist frequency', ha = 'left' )
plt.savefig( './soundpower.eps' )

#Plot Fourier coefficients well below Nyquist frequency (one-sided)
plt.figure( )
plt.plot( np.array( range( 3000 )), abs( coeff[ :3000 ] ) ** 2, 'b-' )
plt.title( 'low frequency power spectrum' )
plt.xlabel( 'index $k$ of Fourier coefficients $c_k$' )
plt.ylabel( '$|c_k|^2$' )
plt.grid( )
plt.savefig( './soundlowpower.eps' )

plt.figure( )
plt.title( 'sound filtering' )
plt.xlabel( 'time[s]' )
plt.ylabel( 'sound pressure' )

for cutoff in [ 5000, 3000, 1000 ]:
	
	#Low pass filtering (also pass symmetric counterpart of low frequencies)
	coeff_filter = np.zeros( n, dtype = np.complex64 )
	coeff_filter[ :cutoff ] = coeff[ :cutoff ]
	coeff_filter[ n-cutoff: ] = coeff[ n-cutoff: ]
	y_filter = ifft( coeff_filter )
	audiowrite( 'hello_cutoff{}.wav'.format( cutoff ), y_filter, Fs )
	sound( y_filter, Fs, 'hello.wav with cut-off frequency {}'.format( cutoff ))
	plt.plot( 
		np.array( np.arange( 30000, 32000 )) / Fs, 
		y_filter.real[ 30000:32000 ],
		label = 'cut-off = {}'.format( cutoff ))

plt.legend( )
plt.savefig( './soundfiltered.eps' )