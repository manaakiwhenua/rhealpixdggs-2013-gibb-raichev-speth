"""
This Python 3.3 module implements the rHEALPix discrete global grid system.

CHANGELOG:

- Alexander Raichev (AR), 2012-11-12: Initial version based upon grids.py.
- AR, 2012-12-10: Corrected centroid() and moved some methods from graphics.py to here.
- AR, 2012-12-19: Tested all the methods and added examples.
- AR, 2013-01-01: Added ellipsoidal functionality to neighbor() and neighbors().
- AR, 2013-01-14: Added intersects_meridian(), cell_latitudes(), cells_from_meridian(), cells_from_parallel(), cells_from_region().
- AR, 2013-01-16: Changed the string keyword 'surface' to a boolean keyword 'plane'.
- AR, 2013-03-11: Added minimal_cover(), boundary(), interior(), triangle(), nw_vertex().
- AR, 2013-03-14: Fixed bug in nw_vertex().
- AR, 2013-07-23: Ported to Python 3.3.
- Robert Gibb (RG), 2020-07-13: Issue #1 Multiple tests fail due to rounding errors
- RG, 2020-07-31: Issue #1 Seperated Sage code out of rhealpixdggs-py

NOTES:

All lengths are measured in meters and all angles are measured in radians
unless indicated otherwise.

By 'ellipsoid' throughout, i mean an ellipsoid of revolution and *not* a general (triaxial) ellipsoid.

Points lying on the plane are given in rectangular (horizontal, vertical) coordinates, and points lying on the ellipsoid are given in geodetic (longitude, latitude) coordinates unless indicated otherwise.

DGGS abbreviates 'discrete global grid system'.

Except when manipulating positive integers, I avoid the modulo function '%'
and insted write everything in terms of 'floor()'.
This is because Python interprets the sign of '%' differently than
Java or C, and I don't want to confuse people who are translating this code
to those languages.

EXAMPLES:

Create the (1, 2)-rHEALPix DGGS with N_side = 3 that is based on the WGS84 ellipsoid. Use degrees instead of the default radians for angular measurements ::

    >>> from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
    >>> E = WGS84_ELLIPSOID
    >>> rdggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=2, N_side=3)
    >>> print(rdggs)
    rHEALPix DGGS:
        N_side = 3
        north_square = 1
        south_square = 2
        max_areal_resolution = 1
        max_resolution = 15
        ellipsoid:
            R_A = 6374581.467096525
            a = 6378137.0
            b = 6356752.314140356
            e = 0.0578063088401125
            f = 0.003352810681182319
            lat_0 = 0
            lon_0 = 0
            radians = False
            sphere = False

NOTES::  .. Issue #1 was ..
        ellipsoid:
            R_A = 6374581.4671 *
            a = 6378137.0
            b = 6356752.314140356 *
            e = 0.0578063088401
            f = 0.003352810681182319

Pick a (longitude-latitude) point on the ellipsoid and find the resolution 1 cell that contains it ::

    >>> p = (0, 45)
    >>> c = rdggs.cell_from_point(1, p, plane=False); print(c)
    N8

Find the ellipsoidal (edge) neighbors of this cell ::

    >>> for (direction, cell) in sorted(c.neighbors(plane=False).items()):
    ...     print(direction, cell)
    east N5
    south_east Q0
    south_west P2
    west N7

Find the planar (edge) neighbors of this cell ::

    >>> for (direction, cell) in sorted(c.neighbors('plane').items()):
    ...     print(direction, cell)
    down P2
    left N7
    right Q0
    up N5

Find all the resolution 1 cells intersecting the longitude-latitude aligned ellipsoidal quadrangle with given northwest and southeast corners ::

    >>> nw = (0, 45)
    >>> se = (90, 0)
    >>> cells = rdggs.cells_from_region(1, nw, se, plane=False)
    >>> for row in cells:
    ...     print([str(cell) for cell in row])
    ['N8', 'N5', 'N2']
    ['Q0', 'Q1', 'Q2', 'R0']
    ['Q3', 'Q4', 'Q5', 'R3']

Compute the ellipsoidal nuclei of these cells ::

    >>> for row in cells:
    ...     for cell in row:
    ...         print(cell, cell.nucleus(plane=False))
    N8 (0.0, 58.47067782962736)
    N5 (45.000000000000036, 58.47067782962736)
    N2 (90.00000000000003, 58.470677829627355)
    Q0 (14.999999999999998, 26.438744923100096)
    Q1 (45.0, 26.438744923100096)
    Q2 (74.99999999999999, 26.438744923100096)
    R0 (105.00000000000001, 26.438744923100096)
    Q3 (14.999999999999998, 3.560649871414923e-15)
    Q4 (45.0, 3.560649871414923e-15)
    Q5 (74.99999999999999, 3.560649871414923e-15)
    R3 (105.00000000000001, 3.560649871414923e-15)

NOTES::  .. Issue #1 was ..
    N8 (0.0, 58.470677829627363) *
    N5 (45.000000000000036, 58.470677829627363) *
    N2 (90.000000000000028, 58.470677829627355) *
    Q0 (14.999999999999998, 26.438744923100096)
    Q1 (45.0, 26.438744923100096)
    Q2 (74.999999999999986, 26.438744923100096)
    R0 (105.00000000000001, 26.438744923100096)
    Q3 (14.999999999999998, 3.560649871414923e-15)
    Q4 (45.0, 3.560649871414923e-15)
    Q5 (74.999999999999986, 3.560649871414923e-15) *
    R3 (105.00000000000001, 3.560649871414923e-15)

Create a (0, 0)-rHEALPix DGGS with N_side = 3 based on the WGS84 ellipsoid.
Use degrees instead of the default radians for angular measurements and
orient the DGGS so that the planar origin (0, 0) is on Auckland, New Zealand ::

    >>> p = (174, -37)  # Approximate Auckland lon-lat coordinates
    >>> from rhealpixdggs.ellipsoids import *
    >>> E = Ellipsoid(a=WGS84_A, f=WGS84_F, radians=False, lon_0=p[0], lat_0=p[1])
    >>> rdggs = RHEALPixDGGS(E, N_side=3, north_square=0, south_square=0)
    >>> print(rdggs)
    rHEALPix DGGS:
        N_side = 3
        north_square = 0
        south_square = 0
        max_areal_resolution = 1
        max_resolution = 15
        ellipsoid:
            R_A = 6374581.467096525
            a = 6378137.0
            b = 6356752.314140356
            e = 0.0578063088401125
            f = 0.003352810681182319
            lat_0 = -37
            lon_0 = 174
            radians = False
            sphere = False

NOTES::  .. Issue #1 was ..
        ellipsoid:
            R_A = 6374581.4671 *
            a = 6378137.0
            b = 6356752.314140356
            e = 0.0578063088401 *
            f = 0.003352810681182319

    >>> print(rdggs.cell_from_point(1, p, plane=False))
    Q3

"""
# *****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http: //www.gnu.org/licenses/
# *****************************************************************************
# Import third-party modules.
"""
from numpy import array, base_repr, ceil, log, pi, sign
from scipy import integrate
"""

# Import standard modules.
"""
from itertools import product
from random import uniform, randint
from colorsys import hsv_to_rgb
"""

# Import my modules.
"""
import rhealpixdggs.pj_rhealpix as pjr
import rhealpixdggs.projection_wrapper as pw
from rhealpixdggs.ellipsoids import (
    WGS84_ELLIPSOID,
    WGS84_ELLIPSOID_RADIANS,
    UNIT_SPHERE,
    UNIT_SPHERE_RADIANS,
)
from rhealpixdggs.utils import my_round
"""

# class RHEALPixDGGS(object):
# plot_cells has been extracted from dggs.py,
# where it was defined within the RHEALPixDGGS class.
# No changes have been made to the code,
# nor has not been tested in the new context outside the class.

def plot_cells(
    self, cells, surface="plane", label=True, fontsize=15, saturation=0.5
):
    """
    Plot the given list of cells on the given surface.
    The cells should all come from the same rHEALPix DGGS.
    Inessential graphics method.
    Requires Sage graphics methods.

    INPUT:

    - `cells` - A list of cells from a common rHEALPix DGGS.
    - `surface` - (Optional; default = 'plane').
      One of the strings 'plane', 'plane_lonlat', 'cube', or 'ellipsoid'.
      Surface to draw cells on.
    - `label` - (Optional; default = True). If True, then label cells
      with their names. If False, then don't.
    - `saturation` - (Optional) Number between 0 and 1 indicating the
      saturation value of the cell color.
    """
    from numpy import array, pi, sign

    from sage.all import (
        Graphics,
        text,
        text3d,
        line,
        polygon,
        parametric_plot3d,
        RealNumber,
        Integer,
    )

    # Make Sage types compatible with Numpy.
    RealNumber = float
    Integer = int

    P = Graphics()
    if not cells:
        return P
    # Draw cells.
    if surface == "plane":

        def texty(s, p):
            return text(s, p, color="black", fontsize=fontsize)

        for cell in cells:
            outline = cell.vertices(plane=True)
            # Draw cell boundary.
            P += line(outline + [outline[0]], color="black")
            # Draw cell interior in color.
            P += polygon(outline, rgbcolor=cell.color(saturation=saturation))
            if label:
                # Label cell.
                anchor = cell.nucleus(plane=True)
                P += texty(str(cell), anchor)
    elif surface == "plane_lonlat":

        def texty(s, p):
            return text(s, p, color="black", fontsize=fontsize)

        PI = self.ellipsoid.pi()
        for cell in cells:
            shape = cell.ellipsoidal_shape()
            if shape == "quad":
                outline = cell.boundary(n=2, plane=False, interior=True)
                # Draw cell boundary.
                P += line(outline + [outline[0]], color="black")
                # Draw cell interior.
                P += polygon(outline, rgbcolor=cell.color(saturation=saturation))
            elif shape == "cap":
                phi = cell.vertices(plane=False)[0][1]
                s = sign(phi)
                outline = [
                    (-PI, phi),
                    (-PI, s * PI / 2),
                    (PI, s * PI / 2),
                    (PI, phi),
                ]
                # Draw cell boundary.
                P += line(outline + [outline[0]], color="black")
                # Draw cell interior.
                P += polygon(outline, rgbcolor=cell.color(saturation=saturation))
            elif shape == "skew_quad" or (
                shape == "dart"
                and abs(abs(cell.nucleus(plane=False)[0]) - PI) > PI / 8
            ):
                i = cell.resolution
                n = max(45 // 3 ** i, 3)
                outline = cell.boundary(n=n, plane=False, interior=True)
                # Draw cell boundary.
                P += line(outline + [outline[0]], color="black")
                # Draw cell interior.
                P += polygon(outline, rgbcolor=cell.color(saturation=saturation))
            if label:
                # Label cell.
                if shape == "cap":
                    anchor = (0, phi / 2 + s * PI / 4)
                else:
                    anchor = cell.nucleus(plane=False)
                P += texty(str(cell), anchor)
    elif surface == "ellipsoid":

        def transform(x, y):
            return self.xyz(x, y)

        def texty(s, p):
            return text3d(s, 1.1 * array(p))

        f = (
            lambda x, y: transform(x, y)[0],
            lambda x, y: transform(x, y)[1],
            lambda x, y: transform(x, y)[2],
        )
        for cell in cells:
            i = cell.resolution
            # Draw cell boundary.
            # Number of points on cell edges to interpolate between:
            n = max(20 // 3 ** i, 3)
            outline = [transform(*p) for p in cell.boundary(n=n, plane=True)]
            P += line(outline + [outline[0]], color="black")
            # Draw cell interior.
            # Number of points in cell interior to interpolate between:
            m = max(30 // 3 ** i, 3)
            xr, yr = cell.xy_range()
            P += parametric_plot3d(
                f,
                xr,
                yr,
                color=cell.color(saturation=saturation),
                plot_points=[m, m],
            )
            if label:
                # Label cell.
                anchor = transform(*cell.nucleus(plane=True))
                P += texty(str(cell), anchor)
    else:
        # Draw cells on cube.
        def transform(x, y):
            return self.xyz_cube(x, y)

        def texty(s, p):
            return text3d(s, 1.1 * array(p))

        for cell in cells:
            outline = [transform(*p) for p in cell.vertices(plane=True)]
            # Draw cell boundary.
            P += line(outline + [outline[0]], color="black")
            # Draw cell interior.
            P += polygon(outline, rgbcolor=cell.color(saturation=saturation))
            if label:
                # Label cell.
                anchor = transform(*cell.nucleus(plane=True))
                P += texty(str(cell), anchor)
    return P
