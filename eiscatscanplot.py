# -*- coding: utf-8 -*-
"""
Module to plot EISCAT radar scans. Scan detection based on azimuth (primary) or
elevation (if stationary in azimuth)

@author: Christer van der Meeren

Initially tested at the ESR 2014-11-26

"""

from __future__ import print_function
import os
import shutil
from os.path import isfile, join
import re
import logging
import datetime as dt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.ticker import LogLocator, MultipleLocator
from eiscat_toolkit import loc2gg, load_param_single_simple
from helper_functions import find_coord_alt, find_coord_alt_2, bm_drawCmappedPoly, cbar_right
from mpl_toolkits.basemap import Basemap
from time import sleep

plt.ioff()  # performance impact is negligible, but when plotting many plots, this allows to plot and save without showing figure
np.seterr(invalid='ignore', divide='ignore')
plt.rcParams['font.size'] = 12

#==============================================================================
# Some defaults
#==============================================================================

# behaviour
RT = True  # realtime plotting. continuously checks data folder for new files
RT_replotAfterScan = True  # after a realtime scan/plot is finished, close and plot everything at once and save this plot instead. This corrects the flat-projection plots and doesn't take long.
onlyDoScanNo = None  # the only scan number to plot. Takes precedence over which scan to start at (which is given in input when running the code)

# plot-related
radarLoc = [78.153, 16.029, 0.438]  # latitude, longitude, alt [km] of radar (and center of map)
mapWidth = 1.8e6  # map width (and height) in map projection coordinates
figSize = 60  # DPI of figure (makes EVERYTHING larger/smaller)

# debug
debugRT = False  # forces update of plot after each file is read

#logging.basicConfig(level=logging.INFO)  # uncomment to turn on info messages (not debug)
#logging.basicConfig(level=logging.DEBUG)  # uncomment to turn on debug messages

#==============================================================================
# End switchboard
#==============================================================================

# get most recent 32m data directory
baseDataDir = '/analysis/results'  # will look for 32m data folders in this directory
if os.path.isdir(baseDataDir):
    subdirs_32m = [join(baseDataDir, d) for d in os.listdir(baseDataDir) if os.path.isdir(os.path.join(baseDataDir, d)) and '@32m' in d]
    dataFolder = join(baseDataDir, max(subdirs_32m, key=os.path.getmtime))
else:
    dataFolder = None


class ScanDetectionError(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Scan(object):

    def __init__(self, init_Time, init_par1D, init_par2D, init_err2D, scDir, scNo):

        # plot-related attributes
        self.fig = None
        self.axes = []
        self.map = None
        self.plotAlts = []  # contains the altitude lines to be plotted
        self.plottedAlts = {}
        self.plottedAltsFlat = {}
        self.altLats = {}
        self.altLons = {}
        self.altLatsFlat = {}
        self.altLonsFlat = {}
        self.cbar = [None]*4
        self.cbarAxis = [None]*4
        self.plottedBeams = []
        self.mainTitle = None
        self.byline = None
        self.elScanDirectionPlotted = False
        self.scDirElev = None

        # data
        self.Time = init_Time
        self.par1D = init_par1D
        self.par2D = init_par2D
        self.err2D = init_err2D

        # scan direction and number
        self.scDir = scDir
        self.scNo = scNo

        # is scan finished? (Used for determining method for rotating data in flat-projection plots)
        self.finished = False

        self.midTimes = self.Time[0, :] + (self.Time[1, :] - self.Time[0, :]) / 2
        self.scanStart = self.Time[0, 0]
        self.scanEnd = self.Time[-1, -1]

        # convenience shortcuts
        self.tStart = self.Time[0, :]
        self.tEnd = self.Time[1, :]
        self.Az = self.par1D[:, 0]
        self.El = self.par1D[:, 1]
        self.Pt = self.par1D[:, 2]
        self.Tsys = self.par1D[:, 3]
        try:
            self.Oppd_Php = self.par1D[:, 4]
        except IndexError:
            self.Oppd_Php = None
        self.Ran = self.par2D[:, :, 0]
        self.Alt = self.par2D[:, :, 1]
        self.Ne = self.par2D[:, :, 2]
        self.Te = self.par2D[:, :, 3]
        self.Ti = self.par2D[:, :, 4]
        self.Vi = self.par2D[:, :, 5]
        self.Coll = self.par2D[:, :, 6]
        self.Comp = self.par2D[:, :, 7]
        self.Res = self.par2D[:, :, 8]
        self.errNe = self.err2D[:, :, 0]
        self.errTe = self.err2D[:, :, 1]
        self.errTi = self.err2D[:, :, 2]
        self.errVi = self.err2D[:, :, 3]
        self.errColl = self.err2D[:, :, 4]

    def add_data(self, new_Time, new_par1D, new_par2D, new_err2D):

        # correct dimensions if number of range gates have changed
        self.par2D, new_par2D, self.err2D, new_err2D = fix_dimensions(self.par2D, new_par2D, self.err2D, new_err2D)

        self.Time = np.append(self.Time, new_Time, axis=1)
        self.par1D = np.append(self.par1D, new_par1D, axis=0)
        self.par2D = np.append(self.par2D, new_par2D, axis=1)
        self.err2D = np.append(self.err2D, new_err2D, axis=1)

        self.midTimes = self.Time[0, :] + (self.Time[1, :] - self.Time[0, :]) / 2
        self.scanEnd = self.Time[-1, -1]

        self.tStart = self.Time[0, :]
        self.tEnd = self.Time[1, :]
        self.Az = self.par1D[:, 0]
        self.El = self.par1D[:, 1]
        self.Pt = self.par1D[:, 2]
        self.Tsys = self.par1D[:, 3]
        try:
            self.Oppd_Php = self.par1D[:, 4]
        except IndexError:
            self.Oppd_Php = None
        self.Ran = self.par2D[:, :, 0]
        self.Alt = self.par2D[:, :, 1]
        self.Ne = self.par2D[:, :, 2]
        self.Te = self.par2D[:, :, 3]
        self.Ti = self.par2D[:, :, 4]
        self.Vi = self.par2D[:, :, 5]
        self.Coll = self.par2D[:, :, 6]
        self.Comp = self.par2D[:, :, 7]
        self.Res = self.par2D[:, :, 8]
        self.errNe = self.err2D[:, :, 0]
        self.errTe = self.err2D[:, :, 1]
        self.errTi = self.err2D[:, :, 2]
        self.errVi = self.err2D[:, :, 3]
        self.errColl = self.err2D[:, :, 4]

    def _find_az_speed(self, ind1, ind2):
        """Returns the azimuth speed in deg/sec between the beams corresponding
        to the two indices (scan direction from ind1 to ind2)"""
        azimDist = self.Az[ind2]-self.Az[ind1]
        timeDist = self.tStart[ind2]-self.tStart[ind1]
        return azimDist/timeDist.total_seconds()

    def _find_el_speed(self, ind1, ind2):
        """Returns the azimuth speed in deg/sec between the beams corresponding
        to the two indices (scan direction from ind1 to ind2)"""
        elevDist = self.El[ind2]-self.El[ind1]
        timeDist = self.tStart[ind2]-self.tStart[ind1]
        return elevDist/timeDist.total_seconds()

    def _data_skip_before(self, ind):
        """Returns True if there is a data skip before the beam corresponding
        to ind (as determined by the difference between end time of ind-1 and
        the start time of ind)"""
        if ind == 0 or ind == len(self.Az):
            # data skip before first beam and after last beam
            return True
        else:
            # currently the limit is set at 2 seconds between end[i-1] and start[i]
            return abs(self.tStart[ind]-self.tEnd[ind-1]).total_seconds() >= 2

    def _data_skip_after(self, ind):
        """Same as _data_skip_before, but checks for data skip AFTER scan"""
        return self._data_skip_before(ind+1)

    def _vertex_array(self, alts, radarLoc, doAll=False):
        """Returns an array of vertices to be used with :function:`~cvdm.bm_drawCmappedPoly`"""

        # skip already plotted beams
        allBeams = np.array(range(0, len(self.Az)))
        plottedBeams = np.array(self.plottedBeams)
        if doAll:
            self.newBeamNos = allBeams
            self.altLats = {}
            self.altLons = {}
        else:
            self.newBeamNos = np.setdiff1d(allBeams, plottedBeams)

        # initiate data structures
        vertexLats = np.ones((self.Ran.shape[0] * len(self.newBeamNos), 4))
        vertexLons = np.ones((self.Ran.shape[0] * len(self.newBeamNos), 4))
        vertexLats_flat = np.ones((self.Ran.shape[0] * len(self.newBeamNos), 4))
        vertexLons_flat = np.ones((self.Ran.shape[0] * len(self.newBeamNos), 4))
        vertexX_elev = np.ones((self.Ran.shape[0] * len(self.newBeamNos), 4))
        vertexY_elev = np.ones((self.Ran.shape[0] * len(self.newBeamNos), 4))

        p = 0  # polygon counter

        for beamNo in self.newBeamNos:

            # find beam width (az1 and az2)
            if not self._data_skip_before(beamNo) and not self._data_skip_after(beamNo):  # not including first/last beams
                az1 = (self.Az[beamNo] + self.Az[beamNo-1]) / 2
                az2 = (self.Az[beamNo] + self.Az[beamNo+1]) / 2
            elif self._data_skip_before(beamNo) and not self._data_skip_after(beamNo):  # including first beam
                az2 = (self.Az[beamNo] + self.Az[beamNo+1]) / 2
                az1 = self.Az[beamNo] - (az2 - self.Az[beamNo])
            elif not self._data_skip_before(beamNo) and self._data_skip_after(beamNo):  # including last beam
                az1 = (self.Az[beamNo] + self.Az[beamNo-1]) / 2
                az2 = self.Az[beamNo] + (self.Az[beamNo] - az1)
            else:  # data skip both before and after
                azSpeed = self._find_az_speed(beamNo-1, beamNo)
                beamWidth = abs(azSpeed * (self.tEnd[beamNo] - self.tStart[beamNo]).total_seconds())
                az1 = self.Az[beamNo] - beamWidth/2
                az2 = self.Az[beamNo] + beamWidth/2

            # do the same for elevation
            if not self._data_skip_before(beamNo) and not self._data_skip_after(beamNo):  # not including first/last beams
                el1 = (self.El[beamNo] + self.El[beamNo-1]) / 2
                el2 = (self.El[beamNo] + self.El[beamNo+1]) / 2
            elif self._data_skip_before(beamNo) and not self._data_skip_after(beamNo):  # including first beam
                el2 = (self.El[beamNo] + self.El[beamNo+1]) / 2
                el1 = self.El[beamNo] - (el2 - self.El[beamNo])
            elif not self._data_skip_before(beamNo) and self._data_skip_after(beamNo):  # including last beam
                el1 = (self.El[beamNo] + self.El[beamNo-1]) / 2
                el2 = self.El[beamNo] + (self.El[beamNo] - el1)
            else:  # data skip both before and after
                elSpeed = self._find_el_speed(beamNo-1, beamNo)
                beamWidth = abs(elSpeed * (self.tEnd[beamNo] - self.tStart[beamNo]).total_seconds())
                el1 = self.El[beamNo] - beamWidth/2
                el2 = self.El[beamNo] + beamWidth/2

            # get coordinates for altitude lines
            # TODO: Find better solution than below for avoiding "stair lines"?
            if beamNo == 0:
                azList = [az1, az2]
            else:
                azList = [az2]
            for az in azList:
                for alt in alts:
                    lat, lon, h = find_coord_alt(radarLoc, self.El[beamNo], az, alt)
                    azFlat = az if self.El[beamNo] <= 90 else az + 180
                    latFlat, lonFlat, hFlat = find_coord_alt_2(radarLoc, self.El[beamNo], 30, azFlat, alt)
                    try:  # append coordinates to this altitude
                        self.altLats[alt] = np.append(self.altLats[alt], lat)
                        self.altLons[alt] = np.append(self.altLons[alt], lon)
                        self.altLatsFlat[alt] = np.append(self.altLatsFlat[alt], latFlat)
                        self.altLonsFlat[alt] = np.append(self.altLonsFlat[alt], lonFlat)
                    except KeyError:  # first values for this altitude
                        self.altLats[alt] = np.array(lat)
                        self.altLons[alt] = np.array(lon)
                        self.altLatsFlat[alt] = np.array(latFlat)
                        self.altLonsFlat[alt] = np.array(lonFlat)

            # loop over all range gates within this beam
            beamRans = self.Ran[:, beamNo]  # ranges for this beam
            for ranNo, ran in enumerate(beamRans):

                if ranNo == 0:  # first range gate
                    ran2 = (ran + beamRans[ranNo+1]) / 2
                    ran1 = ran - (ran2 - ran)
                elif ranNo == len(beamRans)-1:  # last range gate
                    ran1 = (ran + beamRans[ranNo-1]) / 2
                    ran2 = ran + (ran - ran1)
                else:  # any beam not first or last
                    ran1 = (ran + beamRans[ranNo-1]) / 2
                    ran2 = (ran + beamRans[ranNo+1]) / 2

                # compute lats and lons of vertices and add to structure
                if not np.isnan(ran1) and not np.isnan(ran2):
                    v1 = loc2gg(radarLoc, [self.El[beamNo], az1, ran1])
                    v2 = loc2gg(radarLoc, [self.El[beamNo], az1, ran2])
                    v3 = loc2gg(radarLoc, [self.El[beamNo], az2, ran2])
                    v4 = loc2gg(radarLoc, [self.El[beamNo], az2, ran1])
                    vertexLats[p, :] = [v1[0], v2[0], v3[0], v4[0]]
                    vertexLons[p, :] = [v1[1], v2[1], v3[1], v4[1]]

                    azFlat1 = az1 if self.El[beamNo] <= 90 else az1 + 180
                    azFlat2 = az2 if self.El[beamNo] <= 90 else az2 + 180
                    v1_flat = loc2gg(radarLoc, [30, azFlat1, ran1])
                    v2_flat = loc2gg(radarLoc, [30, azFlat1, ran2])
                    v3_flat = loc2gg(radarLoc, [30, azFlat2, ran2])
                    v4_flat = loc2gg(radarLoc, [30, azFlat2, ran1])
                    vertexLats_flat[p, :] = [v1_flat[0], v2_flat[0], v3_flat[0], v4_flat[0]]
                    vertexLons_flat[p, :] = [v1_flat[1], v2_flat[1], v3_flat[1], v4_flat[1]]

                    v1_elev = (np.cos(np.deg2rad(el1))*ran1, np.sin(np.deg2rad(el1))*ran1)
                    v2_elev = (np.cos(np.deg2rad(el1))*ran2, np.sin(np.deg2rad(el1))*ran2)
                    v3_elev = (np.cos(np.deg2rad(el2))*ran2, np.sin(np.deg2rad(el2))*ran2)
                    v4_elev = (np.cos(np.deg2rad(el2))*ran1, np.sin(np.deg2rad(el2))*ran1)
                    vertexX_elev[p, :] = [v1_elev[0], v2_elev[0], v3_elev[0], v4_elev[0]]
                    vertexY_elev[p, :] = [v1_elev[1], v2_elev[1], v3_elev[1], v4_elev[1]]
                else:
                    vertexLats[p, :] = [np.nan]*4
                    vertexLons[p, :] = [np.nan]*4
                    vertexLats_flat[p, :] = [np.nan]*4
                    vertexLons_flat[p, :] = [np.nan]*4
                    vertexX_elev[p, :] = [np.nan]*4
                    vertexY_elev[p, :] = [np.nan]*4

                p = p + 1

            self.plottedBeams.append(beamNo)

        return vertexLats, vertexLons, vertexLats_flat, vertexLons_flat, vertexX_elev, vertexY_elev

    def plot_alts(self, plotFlat=False, rotFlat=None, mapObj=None, ax=None):
        """Plots altitude lines for all altitude in the list `alts`"""

        mapObj = self.map if mapObj is None else mapObj

        for alt in self.altLats:

            # transform to map coordinates
            if hasattr(mapObj, 'coords'):
                x, y = mapObj(self.altLons[alt], self.altLats[alt], coords='geo')
                xFlat, yFlat = mapObj(self.altLonsFlat[alt], self.altLatsFlat[alt], coords='geo')
            else:
                x, y = mapObj(self.altLons[alt], self.altLats[alt])
                xFlat, yFlat = mapObj(self.altLonsFlat[alt], self.altLatsFlat[alt])

            # perform rotation for flat-projection plots if given
            if rotFlat is not None:
                angle = np.deg2rad(rotFlat[0])
                rotMatrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
                # shift rotation point to origin
                if hasattr(mapObj, 'coords'):
                    x0, y0 = mapObj(rotFlat[2], rotFlat[1], coords='geo')
                else:
                    x0, y0 = mapObj(rotFlat[2], rotFlat[1])
                xFlat = xFlat - x0
                yFlat = yFlat - y0

                # rotate
                for i, xy in enumerate(zip(xFlat, yFlat)):
                    xFlat[i] = np.dot(xy, rotMatrix)[0]
                    yFlat[i] = np.dot(xy, rotMatrix)[1]

                # shift back
                xFlat = xFlat + x0
                yFlat = yFlat + y0

            # text alignment for main plots
            xMid = (mapObj.urcrnrx - mapObj.llcrnrx) / 2
            yMid = (mapObj.urcrnry - mapObj.llcrnry) / 2

            if (x[0] >= xMid and y[0] <= yMid and self.scDir == 'cw') or (x[0] <= xMid and y[0] >= yMid and self.scDir == 'ccw'):
                # lower right quadrant and clockwise scan, or upper left quadrant and ccw
                hor = 'left'
                vert = 'bottom'
#                text_alignment = (0, 0)
            elif (x[0] <= xMid and y[0] <= yMid and self.scDir == 'cw') or (x[0] >= xMid and y[0] >= yMid and self.scDir == 'ccw'):
                # lower left quadrant and clockwise scan, or upper right quadrant and ccw scan
                hor = 'left'
                vert = 'top'
#                text_alignment = (0, 1)
            elif (x[0] <= xMid and y[0] >= yMid and self.scDir == 'cw') or (x[0] >= xMid and y[0] <= yMid and self.scDir == 'ccw'):
                # upper left quadrant and clockwise scan, or lower right quadrant and ccw scan
                hor = 'right'
                vert = 'top'
#                text_alignment = (1, 1)
            else:  # upper right quadrant
                hor = 'right'
                vert = 'bottom'
#                text_alignment = (1, 0)

            # text alignment for flat-projection plots
            if self.scDir == 'cw':
                horFlat = 'center'
                vertFlat = 'top'
#                text_alignment_Flat = (0.5, 1)
                angleFlat = -30
                sFlat = str(alt)
            else:
                horFlat = 'center'
                vertFlat = 'top'
#                text_alignment_Flat = (0.5, 1)
                angleFlat = 30
                sFlat = str(alt)

            # initiate dict item for this alt
            try:
                self.plottedAlts[alt]
            except KeyError:
                self.plottedAlts[alt] = []
                self.plottedAltsFlat[alt] = []

            if ax is None:  # plotting by self.plot()

                for _ax in self.axes[0:4]:
                    self.plottedAlts[alt].append(mapObj.plot(x, y, linewidth=2.5, color='w', ax=_ax)[0])
                    self.plottedAlts[alt].append(mapObj.plot(x, y, linewidth=1.5, color='k', ax=_ax)[0])
                    self.plottedAlts[alt].append(_ax.text(x[0], y[0], s=' ESR scan at ' + str(alt) + '$\,$km ', ha='left', va=vert, path_effects=[pe.withStroke(linewidth=3, foreground='w')]))

                for _ax in self.axes[4:8]:
                    if self.scDir not in ['elev incr', 'elev decr']:  # don't plot in flat-projection for elevation scans
                        self.plottedAltsFlat[alt].append(mapObj.plot(xFlat, yFlat, linewidth=2.5, color='w', ax=_ax)[0])
                        self.plottedAltsFlat[alt].append(mapObj.plot(xFlat, yFlat, linewidth=1.5, color='k', ax=_ax)[0])
                        self.plottedAltsFlat[alt].append(_ax.text(xFlat[0], yFlat[0], s=' ' + str(alt) + '$\,$km ', ha=hor, va=vert, path_effects=[pe.withStroke(linewidth=3, foreground='w')]))

            else:  # plotting by self.plot_overlay() or something else
                if plotFlat:
                    self.plottedAltsFlat[alt].append(mapObj.plot(xFlat, yFlat, linewidth=2.5, color='w', ax=ax)[0])
                    self.plottedAltsFlat[alt].append(mapObj.plot(xFlat, yFlat, linewidth=1.5, color='k', ax=ax)[0])
                    self.plottedAltsFlat[alt].append(ax.text(xFlat[0], yFlat[0], s=sFlat, horizontalalignment=horFlat, verticalalignment=vertFlat, rotation=angleFlat, path_effects=[pe.withStroke(linewidth=3, foreground='w')]))
                else:
                    self.plottedAlts[alt].append(mapObj.plot(x, y, linewidth=2.5, color='w', ax=ax)[0])
                    self.plottedAlts[alt].append(mapObj.plot(x, y, linewidth=1.5, color='k', ax=ax)[0])
                    self.plottedAlts[alt].append(ax.text(x[0], y[0], s=str(alt) + '$\,$km', horizontalalignment=hor, verticalalignment=vert, path_effects=[pe.withStroke(linewidth=3, foreground='w')]))

    def update_alts(self, rotFlat=None):

        # lines in large subplots
        for alt, lines in self.plottedAlts.items():
            for line in lines:
                try:
                    line.set_data(*self.map(self.altLons[alt], self.altLats[alt]))
                except AttributeError:  # Probably a Text object
                    pass

        # lines in small subplots
        for alt, lines in self.plottedAltsFlat.items():

            xFlat, yFlat = self.map(self.altLonsFlat[alt], self.altLatsFlat[alt])

            # perform rotation for flat-projection plots if given
            if rotFlat is not None:
                angle = np.deg2rad(rotFlat[0])
                rotMatrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
                # shift rotation point to origin
                x0, y0 = self.map(rotFlat[2], rotFlat[1])
                xFlat = xFlat - x0
                yFlat = yFlat - y0

                # rotate
                for i, xy in enumerate(zip(xFlat, yFlat)):
                    xFlat[i] = np.dot(xy, rotMatrix)[0]
                    yFlat[i] = np.dot(xy, rotMatrix)[1]

                # shift back
                xFlat = xFlat + x0
                yFlat = yFlat + y0

            for line in lines:
                try:
                    line.set_data(xFlat, yFlat)
                except AttributeError:  # Probably a Text object
                    pass

    def plot(self):

        # for some reason, plotting only two beams doesn't work...
        if len(self.Az) < 3:
            return

        # initiate figure/axes
        if self.fig is None:

            self.fig = plt.figure(figsize=(22, 14), dpi=figSize)

            gs_main = mpl.gridspec.GridSpec(3, 1, height_ratios=[2, 2, 0.8])
            gs_map = mpl.gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs_main[0])
            gs_flat = mpl.gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs_main[1])
            gs_elev = mpl.gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs_main[2])

            self.axes.append(plt.subplot(gs_map[0], aspect='equal', adjustable='box-forced'))
            self.axes.append(plt.subplot(gs_map[1], aspect='equal', adjustable='box-forced', sharex=self.axes[0], sharey=self.axes[0]))
            self.axes.append(plt.subplot(gs_map[2], aspect='equal', adjustable='box-forced', sharex=self.axes[0], sharey=self.axes[0]))
            self.axes.append(plt.subplot(gs_map[3], aspect='equal', adjustable='box-forced', sharex=self.axes[0], sharey=self.axes[0]))
            self.axes.append(plt.subplot(gs_flat[0], aspect='equal', adjustable='box-forced'))
            self.axes.append(plt.subplot(gs_flat[1], aspect='equal', adjustable='box-forced', sharex=self.axes[4], sharey=self.axes[4]))
            self.axes.append(plt.subplot(gs_flat[2], aspect='equal', adjustable='box-forced', sharex=self.axes[4], sharey=self.axes[4]))
            self.axes.append(plt.subplot(gs_flat[3], aspect='equal', adjustable='box-forced', sharex=self.axes[4], sharey=self.axes[4]))
            self.axes.append(plt.subplot(gs_elev[0], aspect='equal', adjustable='box-forced'))
            self.axes.append(plt.subplot(gs_elev[1], aspect='equal', adjustable='box-forced', sharex=self.axes[8], sharey=self.axes[8]))
            self.axes.append(plt.subplot(gs_elev[2], aspect='equal', adjustable='box-forced', sharex=self.axes[8], sharey=self.axes[8]))
            self.axes.append(plt.subplot(gs_elev[3], aspect='equal', adjustable='box-forced', sharex=self.axes[8], sharey=self.axes[8]))

            gs_main.update(left=0.035, right=0.96, bottom=0.06, top=0.95, wspace=0.1, hspace=0)

            rocket_track = np.loadtxt('/home/kstdev/users/rocket/CAPER/CAPERtrack.txt', skiprows=8)
#            rocket_track = np.loadtxt('CAPERtrack.txt', skiprows=8)

            # initiate maps
            self.map = Basemap(width=mapWidth,
                               height=mapWidth,
                               projection='aeqd',
                               lat_0=radarLoc[0],
                               lon_0=radarLoc[1],
                               resolution='l')

            # draw coastlines
            for ax in self.axes[0:4]:
                self.map.drawcoastlines(linewidth=0.5, color="k", ax=ax)

                # draw rocket tracks
                self.map.plot(rocket_track[:, 4], rocket_track[:, 3], latlon=True, ax=ax, color='r', linewidth=2, path_effects=[pe.withStroke(foreground='w', linewidth=4)])
                for timeafterlaunch in range(0, 1001, 100):
                    row = rocket_track[rocket_track[:, 0] == timeafterlaunch, :]
                    row = row.flatten()
                    x, y = self.map(row[4], row[3])
                    self.map.plot(x, y, 'r.', ax=ax, markersize=15)
                    ax.annotate(s=str(int(row[0])), xy=(x, y), xycoords='data', xytext=(-8, 0), textcoords='offset points', ha='right', va='center', path_effects=[pe.withStroke(foreground='w', linewidth=3)], color='r')
                    if timeafterlaunch == 200:
                        ax.annotate(s='CAPER\ntrajectory', xy=(x, y), xycoords='data', xytext=(-50, 0), textcoords='offset points', ha='right', path_effects=[pe.withStroke(foreground='w', linewidth=3)], color='r')

            # draw inivible coastlines in flat plots
            for ax in self.axes[4:8]:
                self.map.drawcoastlines(linewidth=0, color="w", zorder=-100, ax=ax)

            # small plot "titles"
            for ax in self.axes[4:8]:
                ax.set_title(u'Flattened to 30\u00B0 elev (for mixed az/el scans)')

        # get coordinates of things to plot
        vertexLats, vertexLons, vertexLats_flat, vertexLons_flat, vertexX_elev, vertexY_elev = self._vertex_array(self.plotAlts, radarLoc)

        # plotting parameters
        #          data   clims          ticks                           colormap       logarithmic colormap
#        toPlot = [('Ne',  (0, 1e12),     MultipleLocator(base=2e11),     'jet',         False),
        toPlot = [('Ne',  (1e10, 1e12),  LogLocator(subs=range(1, 10)),  'nipy_spectral_pinktop',         True),
                  ('Vi',  (-500, 500),   [-500, -250, 0, 250, 500],      'coolwarm',  False),
                  ('Te',  (0, 4000),     [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000],    'nipy_spectral_pinktop',         False),
                  ('Ti',  (0, 3000),     [0, 500, 1000, 1500, 2000, 2500, 3000],          'nipy_spectral_pinktop',         False)]

        for i, e in enumerate(toPlot):

            data = np.ma.masked_invalid(eval("self." + e[0])[:, self.newBeamNos])  # mask nans (and infs)
            if removeLargeErrs:
                data = np.ma.masked_where(eval("self.err" + e[0])[:, self.newBeamNos] > abs(eval("self." + e[0])[:, self.newBeamNos]), data)  # mask where error > value
            data = data.flatten(order='F')

            # main plots
            pCol = bm_drawCmappedPoly(self.map,
                                      vertexLons,
                                      vertexLats,
                                      data,
                                      clims=e[1],
                                      cmap=e[3],
                                      logMap=e[4],
                                      linewidths=0)
            self.axes[i].add_collection(pCol)
            if None in self.cbar:
                self.cbar[i], self.cbarAxis[i] = cbar_right(self.axes[i], leftOffset=0.003, width=0.01, relHeight=1, mappable=pCol, ticks=e[2])

            # flat-projection plots
            if self.scDir in ['cw', 'ccw']:
                pCol = bm_drawCmappedPoly(self.map,
                                          vertexLons_flat,
                                          vertexLats_flat,
                                          data,
                                          clims=e[1],
                                          cmap=e[3],
                                          logMap=e[4],
                                          linewidths=0)
                self.axes[i+4].add_collection(pCol)

            # elev plots
            self.axes[i+8].set_xlim((-800, 800))
            self.axes[i+8].set_ylim((100, 700))
            self.axes[i+8].set_yticks([200, 400, 600])
            self.axes[i+8].set_xlabel('Flat ground distance [km]')
            self.axes[i+8].set_ylabel('Altitude [km]')
            self.axes[i+8].grid('on')
            self.axes[i+8].set_title('Elevation-only scan')

            if self.scDir in ['elev incr', 'elev decr']:
                verts = map(zip, vertexX_elev, vertexY_elev)
                if e[4]:
                    norm = mpl.colors.LogNorm(vmin=e[1][0], vmax=e[1][1])
                else:
                    norm = mpl.colors.Normalize(vmin=e[1][0], vmax=e[1][1])
                pCol = mpl.collections.PolyCollection(verts, array=data, cmap=e[3], norm=norm, linewidths=0)
                self.axes[i+8].add_collection(pCol)

                # show cardinal directions on elevation scans
                if not self.elScanDirectionPlotted:
                    azDir = self.Az[0]
                    if self.scDir == 'elev decr':
                        azDir += 180

                # function to find out if azDir is within 22.5 degrees of some direction
                az_near = lambda x: np.cos(np.deg2rad(azDir-x)) >= np.cos(np.deg2rad(22.5))

                # dictionary of direction string for positive/negative elevation at a given azimuth
                dirDict = {0: ['N', 'S'],
                           45: ['NE', 'SW'],
                           90: ['E', 'W'],
                           135: ['SE', 'NW'],
                           180: ['S', 'N'],
                           225: ['SW', 'NE'],
                           270: ['W', 'E'],
                           315: ['NW', 'SE']}

                # find the correct direction string
                for direction in dirDict:
                    if az_near(direction):
                        dirs = dirDict[direction]
                        self.scDirElev = '-'.join(dirs)
                        break

                left, right = dirs if self.scDir == 'elev decr' else dirs[::-1]
                self.axes[i+8].text(0, -0.25, left, ha='left', fontweight='bold', transform=self.axes[i+8].transAxes)
                self.axes[i+8].text(1, -0.25, right, ha='right', fontweight='bold', transform=self.axes[i+8].transAxes)

        self.elScanDirectionPlotted = True

        # big plot titles
        self.axes[0].set_title(u'Electron density [m$\mathregular{^{−3}}$]')
        self.axes[1].set_title('Ion velocity [m/s, red = away from radar]')
        self.axes[2].set_title('Electron temperature [K]')
        self.axes[3].set_title('Ion temperature [K]')

        # plot altitude levels
        if self.plottedAlts:
            self.update_alts()
        else:
            self.plot_alts()

        # set main title of plot
        mainTitleStr = u'Scan #{0} ({1})   {2}\u2212{3}'.format(self.scNo, self.scDirElev or self.scDir, self.tStart[0].strftime('%d %b %Y %H:%M:%S'), self.tEnd[-1].strftime('%H:%M:%S'))
        if self.mainTitle is None:
            self.mainTitle = self.fig.text(0.5, 0.98, mainTitleStr, horizontalalignment='center', verticalalignment='top', fontsize='xx-large')
        else:
            self.mainTitle.set_text(mainTitleStr)

        if self.byline is None:
            self.byline = self.fig.text(0.98, 0.98, 'Plotting software by Christer van der Meeren (cmeeren@gmail.com)\nCode available at https://github.com/cmeeren/eiscatscanplot', ha='right', va='top')

        # draw/update plot if realtime or only plotting a single scan
        if ((RT or debugRT) and not self.finished) or onlyDoScanNo is not None:
            plt.pause(0.001)

    def plot_overlay(self, mapObj, ax=None, data='Ne', altLines=None, radarLoc=None, removeLargeErrs=False, clims=None, cmap=None, logMap=None, flatProjection=False, defScanSpeedPerIP=1.0*3.2, **kwargs):

        altLines = altLines or []
        radarLoc = radarLoc or [78.153, 16.029, 0.438]

        # input handling: data
        if data is 'Ne':
            clims = [1e10, 1e12] if clims is None else clims
            cmap = 'jet' if cmap is None else cmap
            logMap = True if logMap is None else logMap
        elif data is 'Vi':
            clims = [-500, 500] if clims is None else clims
            cmap = 'esrBlueRed' if cmap is None else cmap
            logMap = False if logMap is None else logMap
        elif data is 'Te':
            clims = [0, 4000] if clims is None else clims
            cmap = 'jet' if cmap is None else cmap
            logMap = False if logMap is None else logMap
        elif data is 'Ti':
            clims = [0, 3000] if clims is None else clims
            cmap = 'jet' if cmap is None else cmap
            logMap = False if logMap is None else logMap
        else:
            print('Invalid data specification:', data)
            return

        # input handling: ax
        ax = plt.gca() if ax is None else ax

        # get (ccw) rotation angle for flat-projection plots (ideally, zenith should be vertical)
        ind_maxEl = np.argmin(np.abs(self.El-90))
        # ind_maxEl corresponds to the last data point of the upleg or the first data point of the downleg
        if np.abs(self.El[ind_maxEl]-90) < 30:
            if np.abs(self.El[ind_maxEl]-90) < defScanSpeedPerIP or self._data_skip_before(ind_maxEl) or self._data_skip_after(ind_maxEl):
                # zenith dump is close enough to 90 deg (or the data skips), plot should be centered at 90 degrees
                # use congruence of scan to find azimuth where el=90
                # TODO: What if ind_maxEl corresponds to one of first two indices?
                az1 = self.Az[ind_maxEl-2]
                az2 = self.Az[ind_maxEl-1]
                el1 = self.El[ind_maxEl-2]
                el2 = self.El[ind_maxEl-1]
                smallPlotRotAngle = az1 + ((90 - el1) * (az2 - az1)) / (el2 - el1)
            else:
                # plot should be centered not at 90 degrees, but at the highest elevation
                smallPlotRotAngle = self.Az[ind_maxEl]
#            if self.El[ind_maxEl-1] < self.El[ind_maxEl]:  # upleg
#                smallPlotRotAngle = self.Az[ind_maxEl]+(self.El[ind_maxEl]-90)
#            else:
#                smallPlotRotAngle = self.Az[ind_maxEl-1]+(self.El[ind_maxEl-1]-90)
        else:
            # max elevation of scan not close to 90, make symmetric instead
            smallPlotRotAngle = (self.Az[0] + self.Az[-1]) / 2

        if flatProjection and hasattr(mapObj, 'coords') and mapObj.coords != 'geo':
            print('WARNING: Flat-projected ESR plot will only work on geographic axes')

        # get coordinates of patch vertices
        vertexLats, vertexLons, vertexLats_flat, vertexLons_flat, vertexX_elev, vertexY_elev = self._vertex_array(altLines, radarLoc, doAll=True)

        if flatProjection:
            lats = vertexLats_flat
            lons = vertexLons_flat
            rot = [-smallPlotRotAngle, radarLoc[0], radarLoc[1]]
        else:
            lats = vertexLats
            lons = vertexLons
            rot = None

        # get data to plot
        data = np.ma.masked_invalid(eval("self." + data)[:, self.newBeamNos])  # mask nans (and infs)

        # mask where error > value, if specified
        if removeLargeErrs:
            data = np.ma.masked_where(eval("self.err" + data)[:, self.newBeamNos] > abs(eval("self." + data)[:, self.newBeamNos]), data)

        # flatten data to be used with custom function below
        data = data.flatten(order='F')

        # plot
        pCol = bm_drawCmappedPoly(mapObj,
                                  lons,
                                  lats,
                                  data,
                                  clims=clims,
                                  cmap=cmap,
                                  logMap=logMap,
                                  rot=rot,
                                  **kwargs)
        ax.add_collection(pCol)

        if flatProjection:
            self.plot_alts(mapObj=mapObj, ax=ax, plotFlat=True, rotFlat=rot)
        else:
            self.plot_alts(mapObj=mapObj, ax=ax)

        return pCol

    def saveFig(self, webAccessFolder='', webAccessFolderExternal=''):

        if self.fig is not None:
            fn = '{0}-{1}_{2}_scan{3:03d}.png'.format(self.tStart[0].strftime('%Y-%m-%d_%H%M'), self.tEnd[-1].strftime('%H%M'), self.scDirElev or self.scDir, self.scNo)
            saveTo = os.path.abspath(os.path.expanduser(savePath))
            if not os.path.isdir(saveTo):
                os.makedirs(saveTo)
            plt.savefig(os.path.join(saveTo, fn), dpi=figSize)
            if RT or debugRT:
                print('Saved file ' + os.path.join(saveTo, fn))
#                logging.info('Saved file ' + os.path.join(saveTo, fn))

            # make web page for external access
            if webAccessFolder:
                make_html(webAccessFolder, webAccessFolderExternal, saveTo, lastFile=fn)

        else:
            logging.warn('Scan #{}: No figure to save (perhaps because number of beams in scan is < 3)'.format(self.scNo))

    def closeFig(self):

        if self.fig is not None:
            logging.info('Closing figure')
            plt.close(self.fig)
        else:
            logging.info('Scan #{}: No figure to close (perhaps because number of beams in scan is < 3)'.format(self.scNo))

        # plot-related attributes
        self.fig = None
        self.axes = []
        self.map = None
        self.plottedAlts = {}
        self.plottedAltsFlat = {}
        self.altLats = {}
        self.altLons = {}
        self.altLatsFlat = {}
        self.altLonsFlat = {}
        self.cbar = [None]*4
        self.cbarAxis = [None]*4
        self.plottedBeams = []
        self.mainTitle = None
        self.byline = None
        self.elScanDirectionPlotted = False


def scan_dir(currentAzim, lastAzim):
    '''Returns the direction of the scan (cw, ccw or stat) given two azimuths.'''
    if currentAzim is None or lastAzim is None:
        return None
    elif currentAzim > lastAzim:  # azimuth increases (not crossing 360/0 degrees)
        return 'cw'
    elif currentAzim < lastAzim:  # azimuth decreases (not crossing 360/0 degrees)
        return 'ccw'
    elif currentAzim == lastAzim:  # stationary
        return 'stat'


def scan_dir_elev(currentElev, lastElev):
    '''Returns the direction of the scan (cw, ccw or stat) given two azimuths.'''
    if currentElev is None or lastElev is None:
        return None
    elif currentElev > lastElev:  # azimuth increases (not crossing 360/0 degrees)
        return 'elev incr'
    elif currentElev < lastElev:  # azimuth decreases (not crossing 360/0 degrees)
        return 'elev decr'
    elif currentElev == lastElev:  # stationary
        return 'stat'


def fix_dimensions(par2D, new_par2D, err2D, new_err2D):
    """ correct dimensions if number of range gates have changed"""

    if par2D[:, 0, 0].shape[0] < new_par2D[:, 0, 0].shape[0]:
        par2D = np.append(par2D, np.ones((new_par2D.shape[0]-par2D.shape[0], par2D.shape[1], par2D.shape[2]))*np.nan, axis=0)
        err2D = np.append(err2D, np.ones((new_err2D.shape[0]-err2D.shape[0], err2D.shape[1], err2D.shape[2]))*np.nan, axis=0)
    if par2D[:, 0, 0].shape[0] > new_par2D[:, 0, 0].shape[0]:
        new_par2D = np.append(new_par2D, np.ones((err2D.shape[0]-new_par2D.shape[0], new_par2D.shape[1], new_par2D.shape[2]))*np.nan, axis=0)
        new_err2D = np.append(new_err2D, np.ones((err2D.shape[0]-new_err2D.shape[0], new_err2D.shape[1], new_err2D.shape[2]))*np.nan, axis=0)

    return par2D, new_par2D, err2D, new_err2D


def scan_parse(dataFolder, savePath,
               doPlot=False, onlyDoScanNo=None, startAt=None, removeLargeErrs=False, RT=False, RT_replotAfterScan=True,
               scanWidth=120, defScanSpeedPerIP=0.62*3, alts=None, radarLoc=None, mapWidth=1.8e6, figSize=72,
               debugRT=False, webAccessFolder='', webAccessFolderExternal=''):
    '''docstring'''

    alts = [250, 500] if alts is None else alts
    radarLoc = radarLoc or [78.153, 16.029, 0.438]
    startAt = '1' if onlyDoScanNo else startAt

    oldfiles = set([])  # already loaded data files
    currentScDir = None   # current scan direction
    lastScDir = None  # last scan direction
    currentAzim = None  # current azimuth
    lastAzim = None  # last azimuth
    currentScDirElev = None   # current scan direction in elevation
    lastScDirElev = None  # last scan direction in elevation
    currentElev = None  # current elevation
    lastElev = None  # last elevation

    # not really needed, but included to make the code editor stop complaining about using uninitiated variables
    Time_prev = None
    par1D_prev = None
    par2D_prev = None
    err2D_prev = None

    scanNoThis = 0

    # list to save scans to if not plotting
    allScans = []

    # find out which scan no / time we should start at
    if startAt is not None and ':' in startAt:
        startHour, startMinute = map(int, startAt.split(':'))
        startTime = dt.time(startHour, startMinute)
        startAtScanNo = 1
    elif startAt is not None:
        startAtScanNo = int(startAt)
        startTime = dt.time(0, 0)
    else:
        startAtScanNo = 1
        startTime = dt.time(0, 0)

    try:
        while True:
            allfiles = set([fn for fn in os.listdir(dataFolder) if (isfile(join(dataFolder, fn)) and re.match('\d*\.mat', fn))])
            newfiles = allfiles - oldfiles

            if len(allfiles) < 2:
                if RT:
                    print('Less than 2 files in folder, waiting 3 sec...')
                    sleep(3)
                else:
                    print('Less than 2 files in folder, exiting')
                    return
            elif len(newfiles) > 0:
                if RT:
                    print('\n{} new files discovered, analyzing...'.format(len(newfiles)))
                else:
                    print('Analyzing {} files...'.format(len(newfiles)))

                newfiles_iterator = sorted(list(newfiles))

                # try to do a progress bar if not realtime
                if not RT and not debugRT:
                    try:
                        import frogress
                        newfiles_iterator = frogress.bar(newfiles_iterator)
                    except:
                        pass

                ###########################################
                # LOOP OVER ALL FILES CURRENTLY IN FOLDER #
                ###########################################

                for fn in newfiles_iterator:

                    # load params
                    Time, par2D, par1D, rpar2D, err2D = load_param_single_simple(join(dataFolder, fn), trueAzEl=True)

                    # avoid crashing if a new beam somehow starts before the last beam
                    if Time_prev is not None and Time[0, 0] <= Time_prev[0, 0]:
                        logging.warning(fn + ': Timestamp of file is before timestamp previous file (analysis error?), skipping file')
                        continue

                    # get current azim and scan direction
                    currentAzim = par1D[0, 0]
                    currentScDir = scan_dir(currentAzim, lastAzim)
                    currentElev = par1D[0, 1]
                    currentScDirElev = scan_dir_elev(currentElev, lastElev)
                    logging.debug('currentAzim = {} | lastAzim = {} | currentScDir = {} | lastScDir = {} | currentElev = {} | lastElev = {} | currentScDirElev = {} | lastScDirElev = {} | scanNoThis = {}'.format(currentAzim, lastAzim, currentScDir, lastScDir, currentElev, lastElev, currentScDirElev, lastScDirElev, scanNoThis))

                    if Time[0, 0].time() < startTime:
                        startAtScanNo = scanNoThis + 1

                    # logic to determine if this scan should be plotted
                    skipThisScan = (onlyDoScanNo is not None and scanNoThis != onlyDoScanNo) or (startAtScanNo is not None and scanNoThis < startAtScanNo)
                    skipNextScan = (onlyDoScanNo is not None and scanNoThis+1 != onlyDoScanNo) or (startAtScanNo is not None and scanNoThis+1 < startAtScanNo)
                    doOnlyThisScan = onlyDoScanNo is not None and scanNoThis == onlyDoScanNo

                    # CONTROL logics to determine where we are in relation to scans
                    noScanYet = currentScDir is None and lastScDir is None and currentScDirElev is None and lastScDirElev is None
                    staticPeriod = currentScDir == 'stat' and lastScDir in ['stat', None] and currentScDirElev == 'stat' and lastScDirElev in ['stat', None]
                    startOfVeryFirstScan = (currentScDir in ['cw', 'ccw'] and lastScDir is None) or (currentScDirElev in ['elev incr', 'elev decr'] and lastScDirElev is None)
                    endOfStaticPeriod = lastScDir == 'stat' and lastScDirElev == 'stat' and (currentScDir in ['cw', 'ccw'] or currentScDirElev in ['elev incr', 'elev decr'])
                    sameScanAsBefore = (lastScDir == currentScDir and lastScDir in ['cw', 'ccw']) or (lastScDir == 'stat' and currentScDir == 'stat' and lastScDirElev == currentScDirElev and lastScDirElev in ['elev incr', 'elev decr'])
                    endOfScan_newScan = (currentScDir in ['cw', 'ccw'] or currentScDirElev in ['elev incr', 'elev decr']) and (lastScDir in ['cw', 'ccw'] or lastScDirElev in ['elev incr', 'elev decr']) and (currentScDir != lastScDir or (currentScDirElev != lastScDirElev and currentScDir == 'stat'))
                    endOfScan_static = currentScDir == 'stat' and currentScDirElev == 'stat' and ((lastScDir in ['cw', 'ccw']) ^ (lastScDirElev in ['elev incr', 'elev decr']))

#                    # old control logics commented out below, will only detect azimuthal scans
#                    noScanYet = currentScDir is None and lastScDir is None
#                    staticPeriod = currentScDir == 'stat' and lastScDir in ['stat', None]
#                    startOfVeryFirstScan = currentScDir in ['cw', 'ccw'] and lastScDir is None
#                    endOfStaticPeriod = currentScDir in ['cw', 'ccw'] and lastScDir == 'stat'
#                    sameScanAsBefore = lastScDir == currentScDir and lastScDir in ['cw', 'ccw']
#                    endOfScan_newScan = lastScDir != currentScDir and lastScDir in ['cw', 'ccw'] and currentScDir in ['cw', 'ccw']
#                    endOfScan_static = lastScDir != currentScDir and lastScDir in ['cw', 'ccw'] and currentScDir == 'stat'

                    if (endOfScan_newScan or startOfVeryFirstScan or endOfStaticPeriod) and Time[-1, -1].time() < startTime:
                        skipNextScan = True

                    # check that only one of the above is True
                    controlList = [noScanYet, staticPeriod, startOfVeryFirstScan, endOfStaticPeriod, sameScanAsBefore, endOfScan_newScan, endOfScan_static]
                    if sum(controlList) != 1:
                        raise ScanDetectionError('none or >1 of control logics are True: {}'.format(controlList))

                    # remove altitude lines if elevation-only scan
                    if startOfVeryFirstScan or endOfStaticPeriod or endOfScan_newScan:
                        scDir = currentScDir if currentScDir != 'stat' else currentScDirElev
                        if scDir in ['elev incr', 'elev decr']:
                            useAlts = [400]
                        else:
                            useAlts = alts

                    #################################################
                    # DETERMINE WHAT GETS GONE WITH THE LOADED DATA #
                    #################################################

                    if noScanYet:
                        logging.info('{}: Skipping data file (very first integration period)'.format(fn))
                        pass

                    if staticPeriod:
                        logging.info('{}: Skipping data file (radar stationary)'.format(fn))
                        pass

                    elif startOfVeryFirstScan or endOfStaticPeriod:
                        logging.info('{}: New scan detected (very first scan or end of stationary period)'.format(fn))
                        scanNoThis += 1
                        if not skipNextScan:
                            if RT or debugRT:
                                print('\nInitializing scan #{} (scan start: {})'.format(scanNoThis, Time[0, 0]))
                            par2D_prev, par2D, err2D_prev, err2D = fix_dimensions(par2D_prev, par2D, err2D_prev, err2D)
                            thisScan = Scan(np.append(Time_prev, Time, axis=1), np.append(par1D_prev, par1D, axis=0), np.append(par2D_prev, par2D, axis=1), np.append(err2D_prev, err2D, axis=1), scDir, scanNoThis)
                            thisScan.plotAlts = useAlts

                    elif sameScanAsBefore:
                        logging.info('{}: Adding data to same scan'.format(fn))
                        if not skipThisScan:
                            thisScan.add_data(Time, par1D, par2D, err2D)
                            if debugRT and doPlot:
                                thisScan.plot()

                    elif endOfScan_newScan or endOfScan_static:
                        if endOfScan_newScan:
                            logging.info('{}: End of scan (new scan detected)'.format(fn))
                        else:
                            logging.info('{}: End of scan (entering static period)'.format(fn))
                        if not skipThisScan:
                            logging.info('   Plotting scan')
                            thisScan.finished = True
                            if (RT or debugRT) and RT_replotAfterScan:
                                print('Closing realtime scan and replotting')
                                thisScan.closeFig()
                            if doPlot:
                                thisScan.plot()
                                thisScan.saveFig(webAccessFolder, webAccessFolderExternal)
                                if not doOnlyThisScan:
                                    thisScan.closeFig()
                                else:  # return scan in addition to plotting if only plotting this scan
                                    return thisScan
                            else:  # no plotting, append scan to list instead
                                if not doOnlyThisScan:  # append to list if parsing many scans
                                    allScans.append(thisScan)
                                else:  # return only this scan if doing one scan
                                    return thisScan
                        if endOfScan_newScan:
                            scanNoThis += 1
                            if not skipNextScan:
                                if RT or debugRT:
                                    print('\nInitializing scan #{} (scan start: {})'.format(scanNoThis, Time[0, 0]))
                                par2D_prev, par2D, err2D_prev, err2D = fix_dimensions(par2D_prev, par2D, err2D_prev, err2D)
                                thisScan = Scan(np.append(Time_prev, Time, axis=1), np.append(par1D_prev, par1D, axis=0), np.append(par2D_prev, par2D, axis=1), np.append(err2D_prev, err2D, axis=1), scDir, scanNoThis)
                                thisScan.plotAlts = useAlts

                    # save current data for next loop
                    lastAzim = currentAzim
                    lastScDir = currentScDir
                    lastElev = currentElev
                    lastScDirElev = currentScDirElev
                    Time_prev = Time
                    par1D_prev = par1D
                    par2D_prev = par2D
                    err2D_prev = err2D

                ###########################################################
                # ALL THE PREVIOUSLY DISCOVERED DATA FILES ARE NOW LOADED #
                ###########################################################

                # update list of loaded files
                oldfiles = allfiles

                # if realtime, plot what we currently have
                if (RT or doOnlyThisScan) and not skipThisScan:
                    logging.info('End of list with new files. Plotting current data...')
                    thisScan.plot()

                # finished processing previously loaded files, look for more (if RT) or quit
                if RT or debugRT:
                    print('Loading done! Looking for more files...')
                else:
                    print('All files loaded, finalizing last scan...')
                    if not skipThisScan:
                        thisScan.finished = True
                        if doPlot:
                            thisScan.plot()
                            thisScan.saveFig(webAccessFolder, webAccessFolderExternal)
                            thisScan.closeFig()
                        else:
                            if doOnlyThisScan:
                                return thisScan
                            else:
                                allScans.append(thisScan)
                                return allScans
                    break
            else:  # no new files found
                print('No new files, waiting 3 sec... (Ctrl-C to quit and finalize current scan)')
                sleep(3)
    except KeyboardInterrupt:  # user has pressed Ctrl-C, finalize current scan and quit
        print('Aborting, finalizing current scan...'),
        if (RT or debugRT) and RT_replotAfterScan:
            print('Closing realtime scan and replotting')
            thisScan.finished = True
            thisScan.closeFig()
        if doPlot:
            thisScan.plot()
            thisScan.saveFig(webAccessFolder, webAccessFolderExternal)
            thisScan.closeFig()
        else:
            if doOnlyThisScan:
                return thisScan
            else:
                allScans.append(thisScan)
        print('Done!')
        return allScans


def make_html(webAccessFolder, webAccessFolderExternal, imageFolder, lastFile):
    '''htmlInternal and htmlExternal point to the same folder, e.g.
          internal: /www_kstdev/display/
          external: http://158.39.70.130/~kstdev/display/
    '''

    webAccessFolder = os.path.abspath(os.path.expanduser(webAccessFolder))
    imageFolder = os.path.abspath(os.path.expanduser(imageFolder))

    # local and remote folders where plots are found
    plotFoldersIn = os.path.join(webAccessFolder, 'eiscatscanplot')
    plotFoldersIn_external = os.path.join(webAccessFolderExternal, 'eiscatscanplot')
    plotsIn = os.path.join(plotFoldersIn, os.path.basename(imageFolder))

    # create plot directory
    if not os.path.isdir(plotsIn):
        os.makedirs(plotsIn)

    # copy last file
    shutil.copyfile(os.path.join(imageFolder, lastFile), os.path.join(plotsIn, lastFile))

    # read html template code from source file
    with open('rt_src.html', 'r') as f:
        s = f.read()

    # make list of files
    files = [fn for fn in os.listdir(imageFolder) if '.png' in fn]

    # copy latest file to permanent name
    shutil.copyfile(os.path.join(imageFolder, files[-1]), os.path.join(webAccessFolder, 'ESRlatest.png'))

    files_external = [os.path.join(plotFoldersIn_external, os.path.basename(imageFolder), fn) for fn in files]

    # insert file list and update latest image
    s = s.replace('[!filenames]', 'filenames = ' + str(files_external))
    s = s.replace('[!latestImg]', files_external[-1])
    s = s.replace('[!latestImgPermalink]', os.path.join(webAccessFolderExternal, 'ESRlatest.png'))
    s = s.replace('[!picFolder]', plotFoldersIn_external)

    # write html page
    with open(os.path.join(webAccessFolder, 'scans.html'), 'w') as f:
        f.write(s)

if __name__ == "__main__":

    # datafolder input
    if dataFolder is None:
        dataFolderPrompt = '\nData folder could not be auto-detected, please specify: >> '
    else:
        dataFolderPrompt = '\nData folder [default: {}] >> '.format(dataFolder)
    while True:
        dataFolderOverride = raw_input(dataFolderPrompt)
        if not dataFolderOverride and dataFolder is not None:
            break  # use default
        elif not dataFolderOverride and dataFolder is None:
            print('Please specify a data folder.')
        elif os.path.isdir(dataFolderOverride):
            dataFolder = dataFolderOverride
            break
        else:
            print('Folder {} doesn\'t exist, please try again.'.format(dataFolderOverride))

    savePath = '~/users/eiscatscanplot/plotted/esr_scans_{}'.format(os.path.basename(os.path.normpath(dataFolder)))  # save plotted figures to this path.

    webAccessFolder = '/www_kstdev/display/'
    webAccessFolderExternal = 'http://158.39.70.130/~kstdev/display/'
#    webAccessFolderExternal = 'file:///C:/www_kstdev/display/'
    webAccessFolderMsg = '{} --> {}'.format(os.path.join(webAccessFolder, 'scan.html'), os.path.join(webAccessFolderExternal, 'scan.html')) if webAccessFolder else 'disabled'

    # guess which scan to start at
    startAt = '1'
    startAtMsg = ''
    # guess from files in savePath
    if os.path.isdir(os.path.abspath(os.path.expanduser(savePath))):
        files = [fn for fn in os.listdir(os.path.abspath(os.path.expanduser(savePath))) if '.png' in fn]
        if not len(files) == 0:
            try:
                scanNos = [fn[-7:-4] for fn in files]
                startAt = str(int(max(scanNos)))  # will overwrite the last plot, which might be incomplete
                startAtMsg = ' (auto-detected)'
            except:
                startAtMsg = ' (unable to parse files in plot folder)'

    scanSpeedDegPerSec = 0.63  # scan speed per second
    IPsec = 6.4  # integration period in seconds
    removeLargeErrs = False   # remove data where error > |value|
    alts = []  # altitude lines to plot [km]. Set to empty list [] to disable

    # additional settings
    while True:
        additionalSettings = raw_input('\nPress Enter to plot or type a number to change the settings below. The following additional settings may be edited.\nScan speed and integration period are only used for beam width\nin realtime plots (saved plots will be correct anyway).\n\n' +
                                       '   1. Save figures to: {}\n'.format(savePath) +
                                       '   2. Web access: {}\n'.format(webAccessFolderMsg) +
                                       '   3. Start at (scan no. x or time HH:MM): {}\n'.format(startAt) +
                                       '   4. Scan speed: {} deg/s\n'.format(scanSpeedDegPerSec) +
                                       '   5. Integration period (GUISDAP): {} s\n'.format(IPsec) +
                                       '   6. Altitude lines: {}\n'.format(' '.join(map(str, alts))) +
                                       '   7. Remove data where error > |value|: {}\n'.format(removeLargeErrs) +
                                       '   8. Realtime plotting: {}\n'.format(RT) +
                                       '\nPlease select a number or press Enter to start plotting >> ')
        if not additionalSettings:
            break
        elif additionalSettings == '1':
            # plot save path
            savePathOverride = raw_input('\nSave plots in folder [default: {}] >> '.format(savePath))
            savePath = savePathOverride or savePath
            savePath = os.path.abspath(os.path.expanduser(savePath))
        elif additionalSettings == '2':
            # 'remote web access to scans
            webAccessFolderOverride = raw_input('\nLocal folder in which to put the html file\n[default: {}, single space to disable] >> '.format(webAccessFolder))
            if webAccessFolderOverride == ' ':
                webAccessFolder = ''
            else:
                latestImagePath = webAccessFolderOverride
                webAccessFolderExternalOverride = raw_input('\nRemote (web) address to the same folder\n[default: {}] >> '.format(webAccessFolderExternal))
                webAccessFolderExternal = webAccessFolderExternalOverride or webAccessFolderExternal
            webAccessFolderMsg = '{} --> {}'.format(os.path.join(webAccessFolder, 'scan.html'), os.path.join(webAccessFolderExternal, 'scan.html')) if webAccessFolder else 'disabled'
        elif additionalSettings == '3':
            # start scan no. or start time
            startAtOverride = raw_input('\nEnter scan number or HH:MM from which to start [default: scan no. {}{}] >> '.format(startAt, startAtMsg))
            if startAtOverride == '0':
                startAt = '1'
            elif startAtOverride:
                startAt = startAtOverride
        elif additionalSettings == '4':
            # scan speed per sec
            scanSpeedDegPerSecOverride = raw_input('\nScan speed in degrees per second [default {}] >> '.format(scanSpeedDegPerSec))
            if scanSpeedDegPerSecOverride:
                scanSpeedDegPerSec = float(scanSpeedDegPerSecOverride)
        elif additionalSettings == '5':
            # analysis integration period
            IPsecOverride = raw_input('\nEffective integration period (of GUISDAP analysis) in seconds. If set to \'1\' in GUISDAP, please enter the experiment integration period. [default: {}] >> '.format(IPsec))
            if IPsecOverride:
                IPsec = float(IPsecOverride)
        elif additionalSettings == '6':
            # altitude lines
            altsOverride = raw_input('\nDraw lines at which altitudes? Space-separated list of numbers, blank for default, single space to disable [default: {}] >> '.format(' '.join(map(str, alts))))
            if altsOverride:
                alts = map(int, altsOverride.split())
        elif additionalSettings == '7':
            # switch error filter
            removeLargeErrs = not removeLargeErrs
            print('Error filter turned {}'.format({True: 'ON', False: 'OFF'}[removeLargeErrs]))
        elif additionalSettings == '8':
            # switch realtime plotting
            RT = not RT
            print('Realtime turned {}'.format({True: 'ON', False: 'OFF'}[RT]))
        else:
            print('Invalid choice')

    defScanSpeedPerIP = scanSpeedDegPerSec*IPsec  # default scan speed in degrees per integration period. Used to assist in rotating flat-projection plots

    print('Starting scan parsing')

    scan_parse(dataFolder=dataFolder, savePath=savePath, doPlot=True, RT=RT,
               onlyDoScanNo=onlyDoScanNo, startAt=startAt, removeLargeErrs=removeLargeErrs, RT_replotAfterScan=RT_replotAfterScan,
               defScanSpeedPerIP=defScanSpeedPerIP, alts=alts, radarLoc=radarLoc, mapWidth=mapWidth, figSize=figSize,
               debugRT=debugRT, webAccessFolder=webAccessFolder, webAccessFolderExternal=webAccessFolderExternal)

    raw_input('\nPlotting finished, figures saved to {}. Press Enter to close >> '.format(savePath))
