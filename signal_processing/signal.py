import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import fftpack
from scipy import signal

##### Signal processing
##1
def ricker(center_freq,dt,tmax,t0=None,plot=False):
    """
    tmax            time lenght
    nt              number of samples
    dt              sampling interval
    center_freq     dominant freq.
    t0              time for center sample (peak) in sec. 
                    default is center output trace: t_0 = (nt-1)/2 * dt
    Output:
    w               wavelet
    t               time
    """
    nt = tmax/dt - 1

    if t0 is None:
        t0 = (nt-1)/2 * dt

    t = np.arange(0,tmax+dt,dt)
    t = t - t0
    arg = (np.pi*center_freq*t)**2
    w = (1-2*arg)*np.exp(-arg)
    t = t + t0
    if plot==True:
        plt.plot(t,w)
        plt.xlabel('Time (sec.)')
        plt.ylabel('Amplitude')
        plt.grid()
    return w,t


##2
def spectra_1d(w,dt,tmax,fmax=0,output=False):
    """
    Frequency representation of a time signal using the Fourier transform
    """
    from scipy.fftpack import fft
    
    spec = fft(w)               # FFT forward from -Fmax:+Fmax double-sided
    fs  =  1/dt                 # samplimg freq
    fnq =  1/(2*dt)             # Nyquist freq
    if fmax ==0.:
        fmax = fnq
    N   = int(fs * tmax)+1      # sample numbers
    m = int(N/2) + 1            # location of Nyquist freq, single-sided f[0:fnq]

    f  = np.linspace(0,fnq,m)   # freq array one-sided,the positive frequency vector

    P2 = np.abs(spec/max(spec)) # Compute the two-sided spectrum P2,
    #P2 = np.abs(spec/N)         # norm to N and 2*P1[1:m-1] or max(spec)
                                # the normalization of the output for the length of the input signal
    P1 = P2[:m]                 # Single-sided spectrum P1 from P2
    P1[1:m-1] = P1[1:m-1]       # Double spectrum except for zero and Nyquist frequencies
    #P1[1:m-1] = 2*P1[1:m-1]

    plt.plot(f,P1,'r') # plotting the spectrum
    plt.xlim([0,fmax])
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Amplitude')
    plt.grid()
    if output == True:
        return P1,f
##3
# zero phase wavelet
def get_zerophasewav(data,dt,fft_type='2D'):
    """
    create a zero phase wavelet usind 2D shot gather
    data [traces,time] - seismic shot gather
    fft                - fft2(data)--'2D',rfft(data,axis=1)--'1D'
    
    Outputs:
    wav                - wavelet (np.read(inverse fft))
    time               - time array
    freq               - freq array [0:Fnq], Fnq = 1/(2*dt)
    plot               - plot resulted wavelet -- plot=True
    """
    Nt   = data.shape[1]
    time = np.arange(0,Nt*dt,dt)  
    # forward fft
    if fft_type == '2D':
        fft_forward = np.fft.fft2(data)
    elif fft_type == '1D':
        fft_forward = np.fft.rfft(data,axis=1)
        spectrum = np.abs(fft_forward/(fft_forward.shape[0])) * np.exp(-1j*2*np.pi)
        
    # take abs values and multiples to exp(-1j*2*np.pi)
    # r abs value of z(rcosØ,rsinØ), z = r(cosØ+jsinØ) <==> z = r * exp(1jØ)
    # it has the effect of rotating z counterclockwise by an angle of Ø on the complex plane
    spectrum =  np.abs(fft_forward/(fft_forward.shape[0])) * np.exp(-1j*2*np.pi)
    # sum over all traces
    spectrum = np.sum(spectrum,axis=0)
    n = spectrum.size
    freq = np.fft.fftfreq(n, d=dt)
    
    fft_inverse = np.fft.ifft(spectrum) 
    # ff W is a vector, then ifftshift swaps the left and right halves of W.
    # Inverse zero-frequency shift (ifftshift undoes the result of fftshift)
    fft_inverse = np.fft.ifftshift(fft_inverse)   
    zerophasewav = np.real(fft_inverse)
        
    fnq = 1/(2*dt) # # Nyquist freq
    freq = np.linspace(0,fnq,int(Nt))
        
    if fft_type == '1D':
        time = np.linspace(0,Nt*dt,int(Nt/2)+1)
        freq = np.linspace(0,fnq,int(Nt/2)+1)
        
    d = np.array([zerophasewav,time,freq])
    df   = pd.DataFrame(data=d.T,columns=['wavelet','time','frequency'])


    return df

##4 
# Statistical Wavelet Extraction
def get_statwav(data,dt,phase=2*np.pi):
    """
    Statistical Wavelet Extraction. The frequency spectra of the autocorrelation
    of all input traces are stacked to produce the wavelet.

    1. Calculate the autocorrelation of each trace in the analysis window.
    2. Calculate the frequency spectrum of each autocorrelation.
       Take the square root of the modulus of each frequency spectrum,
       the zero Hertz component is muted to zero. This step approximates the wavelet spectrum
    3. Stack the spectra
    4. Add the wavelet phase
    5. Take the inverse Fast Fourier Transform (FFT) to extract the wavelet

    Input
    data                   [traces:times]
    dt                     sampling rate of data in sec
    phase=2*pi (positive)  [np.pi,...]

    Output -- Dataframe pandas
    wavelet, time, frequencies
    """
    traces = data.shape[0]
    autocorrelations = np.zeros_like(data)
    spectrum_autocor = np.zeros_like(data,dtype='complex64')
    for trace in range(traces):
        autocorrelations[trace,:] = np.correlate(data[trace,:],data[trace,:],mode='same')
        auto = autocorrelations[trace,:]

        spec = np.fft.fft(auto)
        spec = np.sqrt(np.abs(spec)) 
        #spec = spec * np.exp(-1j*phase)

        spectrum_autocor[trace,:] = spec

    spec_stacked = np.sum(spectrum_autocor,axis=0)
    spec_stacked = spec_stacked * np.exp(-1j*phase)
    wav = np.fft.ifft(spec_stacked)
    wav = np.fft.ifftshift(wav)

    out  = np.real(wav)
    Nt   = data.shape[1]
    fnq  = 1 / (2*dt)
    freq = np.linspace(0,fnq,int(Nt))
    time = np.linspace(0,Nt*dt,int(Nt))

    data = np.array([out,time,freq])
    df   = pd.DataFrame(data=data.T,columns=['wavelet','time','frequency'])

    return df     #out,time,freq
##5
def trace_length(data,dt,tmax):
    """
    rows = time
    cols = traces
    Trace length the data
    data[time,traces] >> newdata[:tmax,traces]
    include tmax
    """
    ind = int(tmax//dt)
    data = data[:ind+1,:]
    return data
##6
def resample(data,olddt,newdt):
    """
    data > rows - times, cols - traces
    olddt - example-.008s
    newdt - example-.004s
    
    returns resampled data
    """
    from scipy import signal
    tmax = data.shape[0]*olddt - olddt
    Nt = int(tmax/newdt) + 1
    otime = np.arange(0,tmax+olddt,olddt)
    resam_d,_ = signal.resample(data,Nt,t=otime,axis=0)
    return resam_d

##7
def spectra_2d(data,dt):
    """
    return spectrum of data
    inputs:
    data[traces,time]
    dt in seconds
    
    return:
    spec,freq
    """
    ft2d = np.fft.rfft2(data)
    ft1d = np.sum(ft2d,axis=0)
    spec = np.abs(ft1d/ft1d.shape[0])
    spec *=2
    freq = np.fft.rfftfreq(data.shape[1],dt)
    return spec,freq
##8
def standart_scale(data):
    """
    Standardize features by removing the mean and scaling to unit variance
    The standard score of a sample x is calculated as:
    z = (x - u) / s
    
    input  - numpy array 
    output - scaled numpy array
    """
    scaled = (data - np.mean(data)) / np.std(data)
    return scaled

def normalize(data,between=(0,1)):
    """
    normalization 
    default min-max 
    or      [a,b]
    
    input  - numpy array 
    output - scaled numpy array
    """
    a,b = between[0],between[1]
    scaled = (b-a) * ( (data-np.min(data)) / (np.max(data) - np.min(data)) ) + a
    return scaled
