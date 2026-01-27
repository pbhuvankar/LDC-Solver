import numpy as np
def get_ghia_data_re1000():
    """
    Returns Ghia et al. (1982) benchmark data for Re=1000
    U velocity along vertical centerline (x=0.5)
    V velocity along horizontal centerline (y=0.5)
    """
    # U velocity along vertical centerline (x = 0.5)
    y_ghia = np.array([0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 
                       0.2813, 0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 
                       0.9531, 0.9609, 0.9688, 0.9766, 1.0000])
    
    u_ghia = np.array([
            0.00000, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289,
            -0.27805, -0.10648, -0.06080, 0.05702, 0.18719, 0.33304,
            0.46604, 0.51117, 0.57492, 0.65928, 1.00000
        ])
    
    # V velocity along horizontal centerline (y = 0.5)
    x_ghia = np.array([0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563,
                       0.2266, 0.2344, 0.5000, 0.8047, 0.9453, 0.9531,
                       0.9609, 0.9688, 1.0000])
    
    v_ghia = np.array([0.00000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077,
                       0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.22194,
                       -0.21952, -0.21766, 0.00000])
    
    return y_ghia, u_ghia, x_ghia, v_ghia