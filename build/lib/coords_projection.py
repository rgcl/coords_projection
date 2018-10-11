#!/usr/bin/env python

import numpy as np
from astropy.io import fits
from typing import Tuple, NamedTuple

# shorthand variables
π = np.pi
sin = np.sin
cos = np.cos

# Tuple for deprojected data
Deprojected = NamedTuple('Deprojected', [
    ('rg', float), ('xg', float), ('yg', float), ('psig', float)
])

# Tuple for projected data
Projected = NamedTuple('Projected', [
    ('xsky', float), ('ysky', float)
])


def radius1(
        size: Tuple[int, int],
        inc: float,
        pa: float,
        centre: Tuple[float, float],
        asc: float,
        rmax: float) -> Deprojected:
    """
    Generate deprojected coordinates, radius and azimuthal angle
    in the plane of an artificial disk, knowing the projected sky-plane
    parameters (inclination, centre of the disk, position angle of the major
    axis of the disk).

    :param size:Tuple[int, int]
        a tuple giving the size of the 2D sky image (nx, ny)
    :param inc: float
        inclination in degrees
    :param pa: float
        position angle of the major axis in degrees
    :param centre: Tuple[float, float]
        a tuple (x0, y0) giving the sky position of the centre of the disk
    :param asc: float
        angular scale in arcsec/pixel of a sky-plane pixel
    :param rmax: float
        maximum galactocentric radius to reach (in arcsec)
    :return: Deprojected
    """

    (x0, y0) = centre
    pa0 = pa * π / 180.
    inc0 = inc * π / 180.

    j = np.outer(np.ones(size[0]), np.arange(size[1]))
    i = np.outer(np.arange(size[0]), np.ones(size[1]))

    xg = -(i - x0) * asc * sin(pa0) + (j - y0) * cos(pa0) * asc
    yg = -((i - x0) * asc * cos(pa0) + (j - y0) * asc * sin(pa0)) / cos(inc0)

    rg = np.sqrt((xg ** 2) + (yg ** 2))  # galactocentric radius in arcsec
    tanpsi = yg / xg

    psig = np.arctan(tanpsi)

    psig[np.where(xg < 0.)] += π
    psig[np.where((xg > 0.) & (yg < 0.))] += 2. * π
    uu = np.where((xg == 0.) & (yg > 0.))

    if uu[0].size != 0:
        psig[uu] = π / 2.
    uu2 = np.where((xg == 0.) & (yg < 0.))
    if uu2[0].size != 0:
        psig[uu2] = 3 * π / 2.

    uu = np.where(rg > rmax)

    if uu[0].size != 0:
        psig[uu2] = 3 * π / 2.
        rg[uu] = np.nan
        yg[uu] = np.nan
        xg[uu] = np.nan
        psig[uu] = np.nan

    return Deprojected(rg, xg, yg, psig)


# Legacy compatibility
rayon1 = radius1


def galactocentric_to_sky(
        rg: float,
        theta: float,
        incli: float,
        pa: float,
        centre: Tuple[float, float],
        asc: float) -> Projected:
    """
    Does the opposite of rayon1, as it generates sky-plane positions (xsky, ysky)
    knowing the disk position (radius and azimuth angle inside the disk).
    Projection on the sky also implies knowledge of the inclination and
    position angle of the disk major axis.

    :param rg: float
        azimuthal angle in rad
    :param theta: float
    :param incli: float
    :param pa: float
    :param centre: float
    :param asc: Tuple[float, float]
    :return: Projected
    """

    (x0, y0) = centre

    xsky = - rg / asc * (cos(theta) * sin(pa / 180. * π)) + sin(theta) * cos(pa * π / 180.) * cos(incli * π/180) + x0
    ysky = rg / asc * (cos(theta) * cos(pa / 180. * π) - sin(theta) * sin(pa * π / 180) * cos(incli * π / 1809.)) + y0

    return Projected(xsky, ysky)


# # # # # # # # # # # Example Usage # # # # # # # # # # #

if __name__ == '__main__':

    print('Testing deprojected...')
    deprojected = rayon1(
        size=(256, 256),  # size of the 2D sky image
        inc=60.,  # in degrees
        pa=30.,  # in degrees, counterclockwise from Y sky (North) axis
        centre=(128., 128.),  # coordinates of the disk centre on the sky
        asc=1.,  # pixel size of the sky image in arcsec
        rmax=100.  # arcsec
    )

    hdu = fits.PrimaryHDU(deprojected.rg)
    hdul = fits.HDUList([hdu])
    hdul.writeto('radiusgal.fits', overwrite=True)
    print('> Saving deprojected image to radiusgal.fits\n\n')

    print('Testing galactocentric_to_sky...')
    result = galactocentric_to_sky(
        rg=10.,  # galactocentric radius in arcsec
        theta=π,  # azimuthal angle in rad
        incli=60.,  # in degrees
        pa=30.,  # in degres, conuterclockswise from Y sky (North) axis
        centre=(128., 128) , # coordinates of the disk centre on the sky
        asc=1.  # pixel size of the sky image in arcsec
    )

    print('> Results:')
    print(result)
