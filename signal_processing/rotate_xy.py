
import numpy as np

def rotate_xy(xy,radian):
    """Rotate x,y coordinates given angle [radian]
       to x start from O x = x + min(n)
       y = y - std(y) 
    
    | x' | = | cos A -sin A | | x |    
    | y' | = | sin A  cos A | | y |
    Params:
        xy     = list[x,y]
        radian = np.deg2rad(degree)
    Output:
        x,y    - list of corrdinates
    
    
    """
    c, s = np.cos(radian), np.sin(radian)
    R = np.array([[c, s], [-s, c]])
    xyhat = np.matmul(R,xy)
    x = xyhat.T[:,0]
    y = xyhat.T[:,1]
    return x,y
    # make X start from O x = x + min(n)
    # set Y=0, std(Y)~0.3 
