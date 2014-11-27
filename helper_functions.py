# -*- coding: utf-8 -*-

def cbar_right(ax,
               valign='top',
               leftOffset=0.02,
               vertOffset=0.,
               width=0.03,
               relHeight=1.,
               use_ColorbarBase=False,
               draw=True,
               rasterize=False,
               **kwargs):
    """
    Creates separate colorbar :class:`~matplotlib.axes.Axes` for full control
    of colorbar placement, and draws a :class:`~matplotlib.pyplot.colorbar`
    on the created :class:`~matplotlib.axes.Axes` to the right of *ax*.
    Default vertical alignment is top.


    Parameters
    ----------

    ax : :class:`~matplotlib.axes.Axes` instance
        colorbar will be positioned relative to (and to the right of) *ax*

    valign : 'top' | 'bottom'
        align the colorbar to either the top or bottom of *ax*

    leftOffset : float
        offset, in figure coordinates, between *ax* and the colorbar

    vertOffset : float
        offset, in figure coordinates, between *ax* and the colorbar (positive
        values offsets toward the centre of *ax*, i.e. downwards if
        valign='top', upwards if valign='bottom')

    width : float
        width of the colorbar in figure coordinates

    relHeight : float
        height of colorbar relative to height of *ax*

    use_ColorbarBase : boolean
        if True, uses :func:`~mpl.colorbar.ColorbarBase` instead of
        :func:`~matplotlib.pyplot.colorbar` for creating the colorbar. This
        allows e.g. creating a colorbar by supplying ``norm`` and ``cmap`` as
        kwargs instead of ``mappable``.

    draw : boolean
        If False, do not force fig.canvas.draw(). Drawing may be needed in
        order to find the correct axes positions.

    rasterize : boolean
        if True, rasterize the colorbar. Can prevent the infamous pdf lines

    kwargs
        keyword arguments sent to :func:`~matplotlib.pyplot.colorbar`


    Returns
    -------

    cb : :class:`~matplotlib.pyplot.colorbar` instance
        the drawn colorbar
    cax : :class:`~matplotlib.axes.Axes` instance
        the new axes onto which the colorbar is drawn


    Note
    ----
    The colorbar `mappable` argument must be supplied as a keyword (in kwargs).

    """

    import matplotlib.pyplot as plt
    import matplotlib as mpl

    # we need to force the figure to draw in order to get the correct positions...
    if draw:
        ax.get_figure().canvas.draw()

    [llx, lly], [urx, ury] = ax.get_position().get_points()

    left = urx + leftOffset
    height = relHeight * (ury - lly)
    if valign == 'top':
        bottom = ury - vertOffset - height
    elif valign == 'bottom':
        bottom = lly + vertOffset

    rect = left, bottom, width, height

    cax = ax.get_figure().add_axes(rect)

    # draw colorbar
    if use_ColorbarBase:
        cb = mpl.colorbar.ColorbarBase(ax=cax, **kwargs)
    else:
        cb = plt.colorbar(cax=cax, **kwargs)

    # rasterize
    if rasterize:
        cb.solids.set_rasterized(True)

    return cb, cax

def bm_drawCmappedPoly(Map,
                       lons,
                       lats,
                       values,
                       clims,
                       cmap=None,
                       logMap=False,
                       resolution=2,
                       rot=None,
                       coordsIn='geo',
                       **kwargs):
    """
    Draws filled polygons on a map with optional interpolation, returns :class:`~matplotlib.collections.PolyCollection`

    Parameters
    ----------

    Map : Basemap instance
        Poly will be drawn on this map

    lons : array of size n_polys x n_verts
        longitudes of polygons, one row for each polygon, one column for each vertice

    lats : array of size n_polys x n_verts
        latitudes of polygons, one row for each polygon, one column for each vertice

    value : list-like of length n_polys
        values used for color mapping the polygon

    clims : tuple of integers
        (cmin, cmax)

    cmap : color map object
        default: jet (matplotlib.cm.jet)

    logMap : boolean
        if True, colormap is logarithmic

    resolution : integer >= 2
        number of points to interpolate each line segment to

    kwargs : keyword arguments sent to :class:`matplotlib.collections.PolyCollection` (most notably ``linewidths`` and ``edgecolors``)
        ..

    Returns
    -------
    p : patch object
        plot with ax.add_patch(p)

    Note
    ----
    Note goes here

    """

    # TODO: Make a PolyCollection instead (don't need to return mappable?)

    import matplotlib as mpl
    import numpy as np

    if cmap is None:
        cmap = mpl.cm.jet

    if lons.shape != lats.shape:
        raise Exception('Shape of lons and lats do not match')
    elif lons.shape[0] != values.shape[0]:
        raise Exception('Shape of values and lons/lats do not match')

# TODO: Interpolate
#    interpLons, interpLats = [lons[0]], [lats[0]]
#    for i in xrange(1, len(lons)):
#        interpLons.extend(list(np.linspace(lons[i-1], lons[i], resolution)[1:]))
#        interpLats.extend(list(np.linspace(lats[i-1], lats[i], resolution)[1:]))
#    # also interpolate the line segment between the last and first lats/lons
#    interpLons.extend(list(np.linspace(lons[-1], lons[0], resolution)[1:]))
#    interpLats.extend(list(np.linspace(lats[-1], lats[0], resolution)[1:]))

    # transform to map coordinates
    if hasattr(Map, 'coords'):
        shape = lons.shape
        x, y = Map(lons.flatten(), lats.flatten(), coords=coordsIn)
        x = x.reshape(shape)
        y = y.reshape(shape)
    else:
        x, y = Map(lons, lats)

    # perform rotation if given
    if rot is not None:
        angle = np.deg2rad(rot[0])
        rotMatrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        # shift rotation point to origin
        if hasattr(Map, 'coords'):
            x0, y0 = Map(rot[2], rot[1], coords=coordsIn)
        else:
            x0, y0 = Map(rot[2], rot[1])
        x = x - x0
        y = y - y0

        # rotate
        datashape = x.shape
        x = x.flatten()
        y = y.flatten()
        for i, xy in enumerate(zip(x, y)):
            x[i] = np.dot(xy, rotMatrix)[0]
            y[i] = np.dot(xy, rotMatrix)[1]
        x = x.reshape(datashape)
        y = y.reshape(datashape)

        # shift back
        x = x + x0
        y = y + y0

    # create array of vertices for PolyCollection
    verts = map(zip, x, y)

    if logMap:
        norm = mpl.colors.LogNorm(vmin=clims[0], vmax=clims[1])
    else:
        norm = mpl.colors.Normalize(vmin=clims[0], vmax=clims[1])

    # create polygon collection
    pcol = mpl.collections.PolyCollection(verts, array=values, cmap=cmap, norm=norm, **kwargs)

    return pcol


def find_coord_alt(loc, el, az, alt, range_start=None, accuracy=1):
    """Given a radar at loc [lat, lon, alt] (altitude in km) pointing with
    elevation el and azimuth az, find the latitude and longitude along the beam
    where the altitude is alt (within accuracy given in km).

    Returns: lat, lon, actualAlt"""

    import math
    from eiscat_toolkit import loc2gg

    # find starting guess for range using flat Earth
    if range_start is None:
        range_start = alt/math.sin(math.radians(el))

    lat, lon, h = loc2gg(loc, [el, az, range_start])

    if abs(h-alt) < accuracy:
        return lat, lon, h
    else:
        range_increment = (alt-h)/math.sin(math.radians(el))
        next_range = range_start + range_increment/2  # have to divide by two, otherwise range will jump back and forth for low elevations
        lat, lon, h = find_coord_alt(loc, el, az, alt, next_range)

    return lat, lon, h


def find_coord_alt_2(loc, elRange, elProject, az, alt, range_start=None, accuracy=1):
    """same as find_coord_alt, but lets you find the proper range using one
    elevation elRange and then project to another elevation elProject."""

    import math
    from eiscat_toolkit import loc2gg

    # find starting guess for range using flat Earth
    if range_start is None:
        range_start = alt/math.sin(math.radians(elRange))

    lat, lon, h = loc2gg(loc, [elRange, az, range_start])

    if abs(h-alt) < accuracy:
        return loc2gg(loc, [elProject, az, range_start])
    else:
        range_increment = (alt-h)/math.sin(math.radians(elRange))
        next_range = range_start + range_increment/2  # have to divide by two, otherwise range will jump back and forth for low elevations
        lat, lon, h = find_coord_alt(loc, elRange, az, alt, next_range)

    return loc2gg(loc, [elProject, az, range_start])
