def source_estimation(data,model,model_wavelet,eps=None):
    """
    source estimation using Pratt approach, transfer function estimation
                          d * m'
    transfer func=     ------------
                        |m|² + eps²
        
    """
    ntrace,nt = data.shape[1],data.shape[0]
    transfer_cmp = np.zeros((nt,ntrace),dtype=np.complex)
    
    
    for j in range(ntrace): 
        # compute Fourier transform
        df = np.fft.fft(data[:,j])   # for jth trace
        mf = np.fft.fft(model[:,j])  # for jth trace
        # Computing transfer function
        amp = np.real(mf * np.conj(mf))
        transfer_cmp[:,j] = ((df) * np.conj(mf))  / (amp + eps)
    # Scaling result with number of traces
    scale = 2.0 / ntrace
    transfer_cmp *=scale
    transfer_cmp = np.sum(transfer_cmp,axis=1)
    # Convolving modeled wavelet with transfer function
    wf = np.fft.fft(model_wavelet)
    wavelet_cmp = transfer_cmp * wf
    
    # Inverse Fourier to obtain wavelet
    wavelet = np.fft.ifft(wavelet_cmp)
    # wavelet = np.fft.ifftshift(wavelet) 
    
    # to QC -- > Convolve modeled data with transfer function
    
    return wavelet 
