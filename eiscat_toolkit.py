# -*- coding: utf-8 -*-

import numpy as np

# make colormaps
#from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
esrBlueRed = {'red':   ((0.0, 0.0, 0.0),
                        (0.5, 0.0, 0.0),
                        (1.0, 1.0, 1.0)),

              'green': ((0.0, 0.0, 0.0),
                        (0.5, 0.0, 0.0),
                        (1.0, 0.0, 0.0)),

              'blue':  ((0.0, 1.0, 1.0),
                        (0.5, 0.0, 0.0),
                        (1.0, 0.0, 0.0))
              }

isomorphicTest = {'red':   ((0.0, 0.0, 0.0),
                            (0.33, 0.0, 0.0),
                            (0.67, 1.0, 1.0),
                            (1, 1.0, 1.0)),

              'green': ((0.0, 0.0, 0.0),
                            (0.33, 0.0, 0.0),
                            (0.67, 0.0, 0.0),
                            (1, 1.0, 1.0)),

              'blue':  ((0.0, 0.0, 0.0),
                            (0.33, 1.0, 1.0),
                            (0.67, 0.0, 0.0),
                            (1, 1.0, 1.0)),
              }


esrJet = {'red':       ((0./12, 0.0, 0.0),
                        (1./12, 0.0, 0.0),
                        (2./12, 0.0, 0.0),
                        (3./12, 0.0, 0.0),
                        (4./12, 0.0, 0.0),
                        (4.7/12, 0.5, 0.5),
                        (6./12, 1.0, 1.0),
                        (7.3/12, 1.0, 1.0),
                        (8./12, 1.0, 1.0),
                        (9./12, 1.0, 1.0),
                        (10./12, 1.0, 1.0),
                        (11./12, 1.0, 1.0),
                        (12./12, 1.0, 1.0)),

          'green':     ((0./12, 0.0, 0.0),
                        (1./12, 0.0, 0.0),
                        (2./12, 0.0, 0.0),
                        (3./12, 0.5, 0.5),
                        (4.7/12, 1.0, 1.0),
                        (5./12, 1.0, 1.0),
                        (6./12, 1.0, 1.0),
                        (7.3/12, 0.5, 0.5),
                        (8./12, 0.0, 0.0),
                        (9./12, 0.0, 0.0),
                        (10./12, 0.0, 0.0),
                        (11./12, 0.5, 0.5),
                        (12./12, 1.0, 1.0)),

          'blue':      ((0./12, 0.0, 0.0),
                        (1./12, 0.5, 0.5),
                        (2./12, 1.0, 1.0),
                        (3./12, 0.5, 0.5),
                        (4.7/12, 0.0, 0.0),
                        (5./12, 0.0, 0.0),
                        (6./12, 0.0, 0.0),
                        (7.3/12, 0.0, 0.0),
                        (8./12, 0.0, 0.0),
                        (9./12, 0.5, 0.5),
                        (10./12, 1.0, 1.0),
                        (11./12, 1.0, 1.0),
                        (12./12, 1.0, 1.0)),
              }

# no pink at top
esrJet2 = {'red':       ((0, 0.0, 0.0),  # black
                        (0.25, 0.0, 0.0),  # blue
                        (0.5, 0.0, 0.0),  # green
                        (0.75, 1.0, 1.0),  # yellow
                        (1.0, 1.0, 1.0)),  # red

          'green':     ((0, 0.0, 0.0),  # black
                        (0.25, 0.0, 0.0),  # blue
                        (0.5, 1.0, 1.0),  # green
                        (0.75, 1.0, 1.0),  # yellow
                        (1.0, 0.0, 0.0)),  # red

          'blue':      ((0, 0.0, 0.0),  # black
                        (0.25, 1.0, 1.0),  # blue
                        (0.5, 0.0, 0.0),  # green
                        (0.75, 0.0, 0.0),  # yellow
                        (1.0, 0.0, 0.0)),  # red
              }

test =  {'red':      ((0.0, 0.0, 0.0),
                        (0.5, 0.5, 0.5),
                        (1.0, 1.0, 1.0)),

          'green':     ((0.0, 0.0, 0.0),
                        (0.5, 0.25, 0.25),
                        (1.0, 0.0, 0.0)),

          'blue':      ((0.0, 0.0, 0.0),
                        (0.5, 0.25, 0.25),
                        (1.0, 0.0, 0.0)),
              }

plt.register_cmap(name='isomorphicTest', data=isomorphicTest)

#blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)
plt.register_cmap(name='esrBlueRed', data=esrBlueRed)
plt.register_cmap(name='esrJet', data=esrJet)
plt.register_cmap(name='esrJet2', data=esrJet2)


#@profile
def load_param_single_simple(fn, status=[0, np.inf], trueAzEl=False):
    """
    Loads parameters from a single GUISDAP result file. NOT meant to be a
    replacement for GUISDAPs own load_param.m. Loads physical parameters and
    time, az/el etc. and little else from monostatic experiments.

    Parameters
    ----------

    fn : string, required
        name of file to load

    status : list of length 2
        [max_status, max_residual]
        Status: 0 = OK, 1 = max number of iterations exceeded, 2 = No fit done
        because data too noisy

    trueAzEl : boolean
        If False, azimuth and elevation will be cast into 0-360 and 0-90 degrees

    """

    import datetime as dt
    from scipy.io import loadmat

    params = loadmat(fn, mat_dtype=True)

    # needed for calculations and other computations
    Te_Ti = params['r_param'][:, 2]  # Te/Ti
    errTe_Ti = params['r_error'][:, 2]
    r_status = params['r_status'][:, 0]

    # don't really know what this is
    if 'r_Offsetppd' in params:
        rOff = params['r_Offsetppd']
    elif 'r_phasepush' in params:
        rOff = params['r_phasepush']
    else:
        rOff = np.nan

    c1 = r_status > status[0]
    c2 = params['r_res'][:, 0] > status[1]
    c3 = params['r_error'][:, :8] > 0
    c12 = np.array([np.bitwise_or(c1, c2), ]*8).T

    params['r_error'][:, :8][c12*c3] = np.nan
    params['r_param'][c12*c3] = np.nan

    # Time
    # XXX: Time is a python datetime object, not MATLAB datenum!
    t = params['r_time']
    tStart = dt.datetime(int(t[0, 0]), int(t[0, 1]), int(t[0, 2]), int(t[0, 3]), int(t[0, 4]), int(t[0, 5]), int(round(np.mod(t[0, 5], 1)*1e6, 0)))
    tEnd = dt.datetime(int(t[1, 0]), int(t[1, 1]), int(t[1, 2]), int(t[1, 3]), int(t[1, 4]), int(t[1, 5]), int(round(np.mod(t[1, 5], 1)*1e6, 0)))

    # 1D params
    Az = params['r_az'][0][0]  # azimuth
    El = params['r_el'][0][0]  # elevation
    Pt = params['r_Pt'][0][0]/10000  # transmitter power
    Tsys = np.median(params['r_Tsys'])
    Oppd_Php = rOff

    # 2D params
    Ran = params['r_range'][:, 0]  # range of scattering volume
    Alt = params['r_h'][:, 0]  # altitude of scattering volume
    Ne = params['r_param'][:, 0]
    Ti = params['r_param'][:, 1]
    Te = Te_Ti * Ti
    Vi = -params['r_param'][:, 4]
    Coll = params['r_param'][:, 3]  # collision frequency
    Comp = params['r_dp'][:, 0]  # Ion composition ([O+]/Ne)
    Res = params['r_res'][:, 0]  # residual of the fit (or standard deviation?)

    # 2D errors
    errNe = params['r_error'][:, 0]
    errTi = params['r_error'][:, 1]
    errTe = (errTi/Ti + errTe_Ti/Te_Ti)*Te
    errVi = params['r_error'][:, 4]
    errColl = params['r_error'][:, 3]

    Time = np.array([[tStart, tEnd], ]).T
    par1D = np.array([[Az, El, Pt, Tsys], ])
    if not np.isnan(Oppd_Php):
        par1D = np.append(par1D, Oppd_Php, axis=1)
    par2D = np.column_stack((Ran, Alt, Ne, Te, Ti, Vi, Coll, Comp, Res))
    par2D = np.expand_dims(par2D, 1)
    err2D = np.column_stack((errNe, errTe, errTi, errVi, errColl))
    err2D = np.expand_dims(err2D, 1)
    rpar2D = np.array([])  # XXX currently not implemented

    # cast azimuth and elevation into range 0-360, 0-90 degrees
    if not trueAzEl:
        d = np.where(par1D[:, 1] > 90)
        par1D[d, 1] = 180 - par1D[d, 1]
        par1D[d, 0] = par1D[d, 0] + 180
        par1D[:, 0] = np.mod(par1D[:, 0]+360, 360)

    return Time, par2D, par1D, rpar2D, err2D


def load_param_simple(path, trueAzEl=False):
    """
    Loads parameters from a directory of GUISDAP result files. NOT meant to be
    a replacement for GUISDAPs own load_param.m. Loads physical parameters and
    time, az/el etc. and little else from monostatic experiments.

    Parameters
    ----------

    path: string, required
        name of file to load

    trueAzEl : boolean
        If False, azimuth and elevation will be cast into 0-360 and 0-90 degrees

    """

    import os
    import fnmatch

    mat_files = fnmatch.filter(os.listdir(path), '*.mat')
    n_ip = len(mat_files)  # number of integration periods

    try:
        import frogress
        iterator = frogress.bar(enumerate(mat_files), steps=len(mat_files))
    except:
        iterator = enumerate(mat_files)

    for i, mat_file in iterator:
        s_Time, s_par2D, s_par1D, s_rpar2D, s_err2D = load_param_single_simple(os.path.join(path, mat_file), trueAzEl=trueAzEl)

        if i == 0:  # initialize data structures
            n_ran = len(s_par2D[:, 0, 1])  # number of range gates

            Time = np.empty((2, n_ip), dtype=object)
            par2D = np.empty((n_ran, n_ip, 9))*np.nan
            par1D = np.zeros((n_ip, s_par1D.shape[1]))
            rpar2D = np.array([])
            err2D = np.empty((n_ran, n_ip, 5))*np.nan

        # correct dimensions if number of range gates have changed
        # XXX: Double-check how this is done in the original matlab code
        if par2D[:, 0, 0].shape[0] < s_par2D[:, 0, 0].shape[0]:
            par2D = np.append(par2D, np.ones((s_par2D.shape[0]-par2D.shape[0], par2D.shape[1], par2D.shape[2]))*np.nan, axis=0)
            err2D = np.append(err2D, np.ones((s_err2D.shape[0]-err2D.shape[0], err2D.shape[1], err2D.shape[2]))*np.nan, axis=0)
        elif par2D[:, 0, 0].shape[0] > s_par2D[:, 0, 0].shape[0]:
            s_par2D = np.append(s_par2D, np.ones((par2D.shape[0]-s_par2D.shape[0], s_par2D.shape[1], s_par2D.shape[2]))*np.nan, axis=0)
            s_err2D = np.append(s_err2D, np.ones((err2D.shape[0]-s_err2D.shape[0], s_err2D.shape[1], s_err2D.shape[2]))*np.nan, axis=0)

        # somewhat the same for par1D
        if s_par1D.shape[1] > par1D.shape[1]:
            par1D = np.vstack((par1D, np.empty((par1D.shape[1], 1))*np.nan))
        elif s_par1D.shape[1] < par1D.shape[1]:
            s_par1D = np.append(s_par1D, np.nan)

        # add current data to data structures
        Time[:, i] = s_Time[:, 0]
        par1D[i, :] = s_par1D[0, :]
        par2D[:, i, :] = s_par2D[:, 0, :]
        err2D[:, i, :] = s_err2D[:, 0, :]

    return Time, par2D, par1D, rpar2D, err2D


def gg2gc(gg):
    """transforms coordinates from geographic (lat, lon, h) to geocentric"""

    import math

    factor = math.pi/180  # conversion factor from degrees to radians
    r_earth = 6378.135  # earth radius (km) and flatness factor
    g = 1.00673944  # earth flatness factor

    lat = gg[0]*factor
    lon = gg[1]*factor
    h = gg[2]

    hor = (r_earth/math.sqrt(1+math.tan(lat)**2/g)+h*math.cos(lat))
    gc = [hor*math.cos(lon), hor*math.sin(lon), r_earth/math.sqrt(g+g**2/math.tan(lat)**2)+h*math.sin(lat)]

    return gc


def gc2gg(gc):
    """transforms coordinates from geocentric to geographic (lat, lon, h)"""

    import math

    factor = math.pi/180  # conversion factor from degrees to radians
    r_earth = 6378.135  # earth radius (km) and flatness factor
    g = 1.00673944  # earth flatness factor

    gg = [0, 0, 0]  # initialize gg, not needed in original MATLAB code...

    if gc[0] == 0 and gc[1] == 0:
        print('Beware of the spinning earth axis!')
        gg = [90, 0, gc[2] - r_earth/g]
    else:
        gg[1] = math.atan2(gc[1], gc[0]) / factor
        r0 = math.sqrt(sum([gc[0]*gc[0], gc[1]*gc[1]]))
        xi0 = gc[2] / (r0 * math.sqrt(g))
        xi_iter = r_earth*(g-1)/(g*r0)
        tanxi = xi0
        tanxi = xi0 + xi_iter * tanxi / math.sqrt(1+tanxi**2)
        tanxi = xi0 + xi_iter * tanxi / math.sqrt(1+tanxi**2)
        gg[0] = math.atan(math.sqrt(g)*tanxi)/factor
        gg[2] = math.sqrt(1+g*tanxi**2)*(r0-r_earth/math.sqrt(1+tanxi**2))

    return gg


def loc2gg(site1, loc):
    """transforms the scattering point location given in local coordinates
    loc   [elevation, azimuth, range] at location
    site1 [latitude, longitude, height]  to geographic coordinates
    """

    import math
    import numpy as np

    factor = math.pi/180

    #  first calculate the  transformation matrices
    lat1 = site1[0]*factor
    lon1 = site1[1]*factor
    sinlat = math.sin(lat1)
    coslat = math.cos(lat1)
    sinlon = math.sin(lon1)
    coslon = math.cos(lon1)
    rlocgc = np.array([[sinlat*coslon, -sinlon,  coslat*coslon],
                       [sinlat*sinlon,  coslon,  coslat*sinlon],
                       [-coslat,             0,  sinlat]])

    s1 = loc[0]*factor
    s2 = loc[1]*factor
    s3 = loc[2]
    loc = s3*np.array([-math.cos(s2)*math.cos(s1), math.sin(s2)*math.cos(s1), math.sin(s1)]).T
    gc_site1 = gg2gc(site1)  # Site1 to geogentric
    gc_sp = gc_site1 + np.dot(rlocgc, loc).T  # Add scattering distance in geocentric
    gg_sp = gc2gg(gc_sp)  # Transform back to geographic

    return gg_sp
