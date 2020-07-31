"""
This Python 3.3 module implements the HEALPix map projection as described in [CaRo2007]_.

.. [CaRo2007] Mark R. Calabretta and Boudewijn F. Roukema, Mapping on the healpix grid, Monthly Notices of the Royal
Astronomical Society 381 (2007), no. 2, 865--872.

CHANGELOG:

- Alexander Raichev (AR), 2013-01-26: Refactored code from release 0.3.
- AR, 2013-03-05: In in_healpix_image() increased eps to 1e-10 to decrease out-of-bounds errors i was getting when
drawing figures.
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
from numpy import pi, floor, sqrt, sin, arcsin, sign, array, deg2rad, rad2deg

# Import my modules.
from rhealpixdggs.utils import my_round, auth_lat, auth_rad


def healpix_diagram(a=1, e=0, shade_polar_region=True):
    """
    Return a Sage Graphics object diagramming the HEALPix projection
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
    dl = array((R * pi / 2, 0))
    lu = [(-R * pi, R * pi / 4), (-R * 3 * pi / 4, R * pi / 2)]
    ld = [(-R * 3 * pi / 4, R * pi / 2), (-R * pi / 2, R * pi / 4)]
    g += line2d([(-R * pi, -R * pi / 4), (-R * pi, R * pi / 4)], color=color)
    g += line2d(
        [(R * pi, R * pi / 4), (R * pi, -R * pi / 4)], linestyle="--", color=color
    )
    for k in range(4):
        g += line2d([array(p) + k * dl for p in lu], color=color)
        g += line2d([array(p) + k * dl for p in ld], linestyle="--", color=color)
        g += line2d(
            [
                array(p) + array((k * R * pi / 2 - R * pi / 4, -R * 3 * pi / 4))
                for p in ld
            ],
            color=color,
        )
        g += line2d(
            [
                array(p) + array((k * R * pi / 2 + R * pi / 4, -R * 3 * pi / 4))
                for p in lu
            ],
            linestyle="--",
            color=color,
        )
    pn = array((-R * 3 * pi / 4, R * pi / 2))
    ps = array((-R * 3 * pi / 4, -R * pi / 2))
    g += point(
        [pn + k * dl for k in range(4)] + [ps + k * dl for k in range(4)],
        size=20,
        color=color,
    )
    g += point(
        [pn + k * dl for k in range(1, 4)] + [ps + k * dl for k in range(1, 4)],
        color="white",
        size=10,
        zorder=3,
    )
    npp = [
        (-R * pi, R * pi / 4),
        (-R * 3 * pi / 4, R * pi / 2),
        (-R * pi / 2, R * pi / 4),
    ]
    spp = [
        (-R * pi, -R * pi / 4),
        (-R * 3 * pi / 4, -R * pi / 2),
        (-R * pi / 2, -R * pi / 4),
    ]

    if shade_polar_region:
        for k in range(4):
            g += polygon([array(p) + k * dl for p in npp], alpha=0.1, color=shade_color)
            g += text(
                str(k),
                array((-R * 3 * pi / 4, R * 5 * pi / 16)) + k * dl,
                color="red",
                fontsize=20,
            )
            g += polygon([array(p) + k * dl for p in spp], alpha=0.1, color=shade_color)
            g += text(
                str(k),
                array((-R * 3 * pi / 4, -R * 5 * pi / 16)) + k * dl,
                color="red",
                fontsize=20,
            )
    return g
