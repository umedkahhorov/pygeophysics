# example of usage
#d2 = dataH1[:,r1:r2] - extract 2D array gather - shot gather
#refpicks2 =  calculate_refpicks(d2, sta_window=20, lta_window=400, dt=0.01)

# define feature functions
def classic_sta_lta_py(a, nsta, nlta):
    """
    (copied from obspy.signal)
    short-term-average long-term-average method STA/LTA
    Computes the standard STA/LTA from a array a. 
    STA and LTA are defined in samples.

    :type a: NumPy :class: ~numpy.ndarray
    :param a: Seismic Trace
    :type nsta: int
    :param nsta: Length of short time average window in samples
    :type nlta: int
    :param nlta: Length of long time average window in samples
    :rtype: NumPy :class:`~numpy.ndarray`
    :return: STA/LTA relation
    """
    # The cumulative sum can be exploited to calculate a moving average (the
    # cumsum function is quite efficient)
    sta = np.cumsum(a ** 2, dtype=np.float64)

    # Copy for LTA
    lta = sta.copy()

    # Compute the STA and the LTA
    sta[nsta:] = sta[nsta:] - sta[:-nsta]
    sta /= nsta
    lta[nlta:] = lta[nlta:] - lta[:-nlta]
    lta /= nlta

    # Pad zeros
    sta[:nlta - 1] = 0

    # Avoid division by zero by setting zero values to tiny float
    dtiny = np.finfo(0.0).tiny
    idx = lta < dtiny
    lta[idx] = dtiny

    return sta / lta

def mean(a,w):
    m = np.cumsum(a**2)
    m[w:] = m[w:] - m[:-w]
    m /= w
    return m

def power(a,w):
    m = np.cumsum(a**2)
    m[w:] = m[w:] - m[:-w]
    return m
def calculate_refpicks(data, sta_window=10, lta_window=400, dt=0.01):
    """Vectorized STA/LTA picking with timing"""
    print("Starting NumPy vectorized processing...")
    t1 = time.time()
    
    # STA/LTA calculation
    refpicks = np.apply_along_axis(
        lambda trace: classic_sta_lta_py(trace, sta_window, lta_window),
        axis=0,
        arr=data
    )
    
    # Pick time extraction (vectorized)
    refpicks_idx = np.argmax(refpicks, axis=0)  # More efficient than apply_along_axis
    refpicks_time = refpicks_idx * dt  # Convert samples to time
    
    t2 = time.time()
    print(f"NumPy result - Shape: {refpicks_idx.shape}, Time: {t2-t1:.2f} sec")
    
    return refpicks_time
