def apply_phase_shift_3d(data, shift_angle):
    """
    Applies a constant phase shift to a 3D data array in the frequency domain.
    The function pads the input along the time axis to the next power of 2,
    computes the real-to-complex FFT (rfft) to obtain only the nonnegative frequencies,
    applies a constant phase multiplier:
        phase_multiplier = exp(-1j * shift_angle_in_radians)
    and then uses the inverse FFT (irfft) to recover the time-domain data.
    
    Parameters
    ----------
    data : np.ndarray
        A 3D array with shape (50, 50, 400), where the last axis represents the time domain.
    shift_angle : float
        The phase shift angle in degrees to apply uniformly across all frequencies.
    Returns
    -------
    shifted_data : np.ndarray
        The phase-shifted data in the time domain with the same shape as the input.
    """
    # Get original dimensions and number of time samples.
    n1, n2, n_samples = data.shape
    # Determine FFT length as the next power of 2 (for computational efficiency).
    NFFT = int(2 ** np.ceil(np.log2(n_samples)))
    # Zero-pad the data along the time axis (pad at the end).
    pad_width = ((0, 0), (0, 0), (0, NFFT - n_samples))
    data_padded = np.pad(data, pad_width, mode='constant')
    # Compute the real-to-complex FFT along the time axis.
    spectrum = np.fft.rfft(data_padded, n=NFFT, axis=-1)
    # Convert the phase shift from degrees to radians.
    rad = np.deg2rad(shift_angle)
    # Create a constant phase multiplier (applied uniformly to all frequency components).
    phase_multiplier = np.exp(-1j * rad)
    # Apply the constant phase shift to the spectrum.
    shifted_spectrum = spectrum * phase_multiplier
    # Convert back to the time domain using the inverse FFT.
    shifted_data_padded = np.fft.irfft(shifted_spectrum, n=NFFT, axis=-1)
    # Remove the padded portion to recover the original number of time samples.
    shifted_data = shifted_data_padded[..., :n_samples]
    return shifted_data
def stat(trin, t, dt):
    """
    Applies a static shift to a seismic trace or gather of traces.
    Parameters:
    -----------
    trin : ndarray
        Input seismic data. Can be a 1D array (single trace) or a 2D array 
        (gather of traces). Shape: (n_samples, n_traces).
    t : ndarray
        Time vector corresponding to the input trace(s). Shape: (n_samples,).
    dt : float or ndarray
        Static shift(s) to apply in seconds. Can be a scalar (applied to all traces)
        or an array (one static shift per trace). Shape: (n_traces,).
    Returns:
    --------
    trout : ndarray
        The statically shifted trace(s). Shape: (n_samples,) for a single trace
        or (n_samples, n_traces) for a gather.
    Notes:
    ------
    - The static shift is applied in the frequency domain using the FFT.
    """
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
