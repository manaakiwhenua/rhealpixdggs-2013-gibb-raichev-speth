"""
This Python 3.3 module implements the rHEALPix map projection.

CHANGELOG:

- Alexander Raichev (AR), 2013-01-26: Refactored code from release 0.3.
- AR, 2013-07-23: Ported to Python 3.3.
- Robert Gibb (RG), 2020-07-13: Issue #1 Multiple tests fail due to rounding errors
- RG, 2020-07-31: Issue #1 Seperated Sage code out of rhealpixdggs-py

NOTE:

All lengths are measured in meters and all angles are measured in radians
unless indicated otherwise.
By 'ellipsoid' below, i mean an oblate ellipsoid of revolution.
"""
# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
from numpy import pi, sign, array, identity, dot, deg2rad, rad2deg

# Import my modules.
from rhealpixdggs.pj_healpix import (
    healpix_sphere,
    healpix_sphere_inverse,
    healpix_ellipsoid,
    healpix_ellipsoid_inverse,
)
from rhealpixdggs.utils import my_round, auth_rad

def rhealpix_diagram(a=1, e=0, north_square=0, south_square=0, shade_polar_region=True):
    """
    Return a Sage Graphics object diagramming the rHEALPix projection
    boundary and polar triangles for the ellipsoid with major radius `a`
    and eccentricity `e`.
    Inessential graphics method.
    Requires Sage graphics methods.
    """
    from sage.all import Graphics, line2d, point, polygon, text, RealNumber, Integer

    # Make Sage types compatible with Numpy.
    RealNumber = float
    Integer = int

    R = auth_rad(a, e)
    g = Graphics()
    color = "black"  # Boundary color.
    shade_color = "blue"  # Polar triangles color.
    north = north_square
    south = south_square
    south_sq = [
        (-R * pi + R * south * pi / 2, -R * pi / 4),
        (-R * pi + R * south * pi / 2, -R * 3 * pi / 4),
        (-R * pi + R * (south + 1) * pi / 2, -R * 3 * pi / 4),
        (-R * pi + R * (south + 1) * pi / 2, -R * pi / 4),
    ]
    north_sq = [
        (-R * pi + R * north * pi / 2, R * pi / 4),
        (-R * pi + R * north * pi / 2, R * 3 * pi / 4),
        (-R * pi + R * (north + 1) * pi / 2, R * 3 * pi / 4),
        (-R * pi + R * (north + 1) * pi / 2, R * pi / 4),
    ]
    # Outline.
    g += line2d(south_sq, linestyle="--", color=color)
    g += line2d(north_sq, linestyle="--", color=color)
    g += line2d(
        [(R * pi, -R * pi / 4), (R * pi, R * pi / 4)], linestyle="--", color=color
    )
    g += line2d(
        [north_sq[0], (-R * pi, R * pi / 4), (-R * pi, -R * pi / 4), south_sq[0]],
        color=color,
    )
    g += line2d([south_sq[3], (R * pi, -R * pi / 4)], color=color)
    g += line2d([north_sq[3], (R * pi, R * pi / 4)], color=color)
    g += point([south_sq[0], south_sq[3]], size=20, zorder=3, color=color)
    g += point([north_sq[0], north_sq[3]], size=20, zorder=3, color=color)
    g += point(
        [(R * pi, -R * pi / 4), (R * pi, R * pi / 4)], size=20, zorder=3, color=color
    )
    g += point(
        [(R * pi, -R * pi / 4), (R * pi, R * pi / 4)], size=10, color="white", zorder=3
    )
    if shade_polar_region:
        # Shade.
        g += polygon(south_sq, alpha=0.1, color=shade_color)
        g += polygon(north_sq, alpha=0.1, color=shade_color)

    # Slice square into polar triangles.
    g += line2d([south_sq[0], south_sq[2]], color="lightgray")
    g += line2d([south_sq[1], south_sq[3]], color="lightgray")
    g += line2d([north_sq[0], north_sq[2]], color="lightgray")
    g += line2d([north_sq[1], north_sq[3]], color="lightgray")

    # Label polar triangles.
    sp = south_sq[0] + R * array((pi / 4, -pi / 4))
    np = north_sq[0] + R * array((pi / 4, pi / 4))
    shift = R * 3 * pi / 16
    g += text(str(south), sp + array((0, shift)), color="red", fontsize=20)
    g += text(
        str((south + 1) % 4),
        sp + array((shift, 0)),
        color="red",
        rotation=90,
        fontsize=20,
    )
    g += text(
        str((south + 2) % 4),
        sp + array((0, -shift)),
        color="red",
        rotation=180,
        fontsize=20,
    )
    g += text(
        str((south + 3) % 4),
        sp + array((-shift, 0)),
        color="red",
        rotation=270,
        fontsize=20,
    )
    g += text(str(north), np + array((0, -shift)), color="red", fontsize=20)
    g += text(
        str((north + 1) % 4),
        np + array((shift, 0)),
        color="red",
        rotation=90,
        fontsize=20,
    )
    g += text(
        str((north + 2) % 4),
        np + array((0, shift)),
        color="red",
        rotation=180,
        fontsize=20,
    )
    g += text(
        str((north + 3) % 4),
        np + array((-shift, 0)),
        color="red",
        rotation=270,
        fontsize=20,
    )
    return g
