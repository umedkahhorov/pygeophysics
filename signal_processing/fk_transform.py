##1
def fk_transform(data,dt,dx,norm='std',kind='np',ishift=1):
    """
    Uses built in np.fft (type='np')  or scipy.fftpack(type='scipy') to perform
    a 2-d f-k transform on a real valued (seismic time-space) matrix.
    Only the positive frequencies are calculated while all kxs are.
    Input traces should be regularly spaced.  Output traces are real-valued and
    represent the magnitude of the complex F-K spectrum.
    
    input:
    data[time,traces]    -seismic data
    dt,dx                -sampling rate alon time and trace axes
    norm['std','minmax',None] -scale
    
    output:
    spec -2D spectrum,freq-vector, kx-vector
    """
    nt,ntr = data.shape[0],data.shape[1]

    if kind == 'np':
        # returns only the positive freq, the t-f transform each
        # each column (or each traces) of the matrix is transformed
        spec_fx  = np.fft.rfft(data,axis=0)
        freq     = np.fft.rfftfreq(nt,dt)          # vector of frequency coordinates for the rows of spec
        spec     = np.fft.ifft(spec_fx.T,axis=0).T # 1D fft along time axis or along rows
    elif kind =='scipy':
        spec_fx  = fftpack.rfft(data,axis=0)
        freq     = np.fft.rfftfreq(nt,dt)
        spec     = fftpack.ifft(spec_fx.T,axis=0).T


    kxnyg    = 1 / (2*dx)
    dkx    = 2 * kxnyg/spec.shape[1]

    kx1 = np.arange(0,kxnyg,dkx)         # wavenumbers posit
    kx2 = np.arange(-kxnyg,-dkx+dkx,dkx) # wavenumbers negat
    kx1 = np.hstack((kx1,kx2))
    ikx = np.argsort(kx1)
    kx  = np.sort(kx1)

    if ishift == 1:
        spec      = spec[:,ikx]

    if norm == None:
        spec = abs(spec)
    elif norm == 'std':
        spec = abs(spec)
        spec = (spec - np.mean(spec)) / np.std(spec)
    elif norm =='minmax':
        a,b = 0,1
        spec = abs(spec)
        spec = (b-a) * ( (spec-np.min(spec)) / (np.max(spec) - np.min(spec)) ) + a

    return spec,freq, kx
