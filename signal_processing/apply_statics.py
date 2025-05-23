import numpy as np
import RSSPython as rs
from scipy import io
from scipy import interpolate
import os
import sys
from numpy.fft import rfft, irfft, rfftfreq


f1 = sys.argv[1]
f2 = sys.argv[2]
vel = int(sys.argv[3])
stat_type = int(sys.argv[4])


rss = rs.RSSdata()
rss.read(f1)
data_all = rss.data
n_samples, n_traces = data_all.shape
original_shape = data_all.shape
sx = rss.srcX
sz = rss.srcZ
rx = rss.GroupX
rz = rss.GroupZ
offset = abs(sx-rx)
print("OFFSET examples",offset[0],offset[1])
nmo = (offset / vel) 

# apply linear moveout
if stat_type==0:
    nmo *=-1

# remove linear moveout
elif stat_type==1:
    nmo *=1

#rss_ = rs.RSSdata()
#rss_.read(f2)
#recmap = rss_.data.squeeze()
#a = recmap.T[:,0]
#a = np.tile(a,2).reshape(2,-1).T
#a[:,0] = a[:,0] + 1
#a[0,0] = 0
#a[:,1] = np.roll(a[:,1],-1)
#a[-1,1] = data_all.shape[1]
#traces = a.astype('int')
#print (traces.size)

########################
########################
def apply_stat(trin, t, dt):
    """
    Applies a static shift to a seismic trace or gather of traces.
    trin : ndarray
        Input seismic data. Can be a 1D array (single trace) or a 2D array 
        (gather of traces). Shape: (n_samples, n_traces).
    t : ndarray
        Time vector corresponding to the input trace(s). Shape: (n_samples,).
    dt : float or ndarray
        Static shift(s) to apply in seconds. Can be a scalar (applied to all traces)
        or an array (one static shift per trace). Shape: (n_traces,).
    --------
    trout : ndarray
        The statically shifted trace(s). Shape: (n_samples,) for a single trace
        or (n_samples, n_traces) for a gather.
    Notes:
    ------
    - The static shift is applied in the frequency domain u
    """
    print("QC Time vector",t[1],t[0])
    # Ensure `trin` is at least 2D for consistency (1D for single trace, 2D for gather).
    trin = np.atleast_2d(trin)  # Shape: (n_samples, n_traces)
    n_samples, n_traces = trin.shape
    # Handle the static shift(s) `dt`: ensure it matches the number of traces.
    dt = np.atleast_1d(dt)
    if dt.size == 1:
        dt = dt * np.ones(n_traces)  # Broadcast single value to all traces.
    if dt.size != n_traces:
        raise ValueError('dt must have one static shift per trace.')
    # Calculate sample interval `dt_t` and convert static shifts to sample indices.
    dt_t = t[1] - t[0]  # Assumes uniform time sampling.
    dt_sample = np.abs(np.ceil(dt / dt_t).astype(int))  # Integer sample shifts.
    max_static_shift = np.max(dt_sample)
    # Determine FFT size with padding to handle shifts.
    NFFT = int(2 ** np.ceil(np.log2(n_samples + max_static_shift)))  # Power of 2 for speed.
    pad_length = NFFT - n_samples
    trin_padded = np.pad(trin, ((0, pad_length), (0, 0)), mode='constant')  # Zero padding.
    
    # Perform forward FFT on the input data.
    Trin = rfft(trin_padded, axis=0)  # Real-to-complex FFT.
    freqs = rfftfreq(NFFT, d=dt_t)  # Frequencies for FFT bins.
    # Apply phase shift in the frequency domain.
    phase_shift = np.exp(-1j * 2 * np.pi * freqs[:, None] * dt)  # Phase shift operator.
    shifted_spectrum = Trin * phase_shift
    # Perform inverse FFT to get back to the time domain.
    trout = irfft(shifted_spectrum, n=NFFT, axis=0)
    # Adjust output based on `flag` and whether input is a single trace or gather.
    if n_traces == 1:  # Single trace.
        trout = trout.flatten()  # Ensure 1D output for single trace.
        trout = trout[:n_samples]  # Return original length.
    else:  # Multi-trace gather.
        trout = trout[:n_samples, :]  # Trim to original number of samples.
    # Return the shifted trace(s)
    return trout.squeeze()
########################
########################
#rss = rs.RSSdata()
#rss.read(f1)
#r1,r2 = traces[s]
#data = rss.data[:,r1:r2]
#sx = rss.srcX[r1:r2]
#sz = rss.srcZ[r1:r2]
#rx = rss.GroupX[r1:r2]
#rz = rss.GroupZ[r1:r2]
    
#rss_geom = rs.RSSdata()
#rss_geom.Ndims=2
#rss_geom.Nheader=4
#rss_geom.data_format=4
#rss_geom.type=2
#rss_geom.geomN[0] = data.shape[0]
#rss_geom.geomN[1] = data.shape[1]
#rss_geom.geomD[0] = 0.01
#rss_geom.geomD[1] = 1
#rss_geom.geomO[0] = 3
#rss_geom.data   = data
#rss_geom.fullsize = data.shape[1]*data.shape[0]
#rss_geom.srcX   = sx
#rss_geom.srcZ   = sz
#rss_geom.GroupX = rx
#rss_geom.GroupZ = rz

twt = np.linspace(0,40,n_samples)
print("TWI examples",twt[1],twt[0])
print("data_all.shape,twt.shape,nmo.shape,nmo.squeeze().shape")
print(data_all.shape,twt.shape,nmo.shape,nmo.squeeze().shape)
lmo_data = apply_stat(data_all,twt,nmo.squeeze())
lmo_data = np.reshape(lmo_data,original_shape)
rss.data = lmo_data
rss.write(f2)
