"""Blah."""

import astroquery
from astropy import coordinates
import astropy.units as u

def altaz_to_radec(alt_az, pos, time):
    """Blah.

    alt_az ():
    pos ():
    time: LST
    """

    # This seems useful:
    # http://www.brayebrookobservatory.org/BrayObsWebSite/BOOKS/matrix_method_rev_d.pdf
    lat, long = pos
    alt, az = alt_az
    
    HA = time - RA

    sinD = np.sin(alt) * np.sin(lat) + np.cos(alt) * np.cos(lat) * np.cos(AZ)

    cosH = (np.sin(alt) - np.sin(lat) * np.sin(D))/(np.cos(lat) * np.cos(D))

    dec = np.arcsin(sinD)
    RA = time - np.arccos(cosH)
    return [ra, dec]


def find_source(time, pos, alt_az):
    """Blah.

    time (something)
    pos (tuple): (lat, long)
    alt_az (tuple): (alt, az)
    """

    ra, dec = altaz_to_radec(alt_az)

    # c = coordinates.SkyCoord("05h35m17.3s -05d23m28s", frame='icrs')
    altaz_str = str()
    c = coordinates.SkyCoord("05h35m17.3s -05d23m28s", frame='icrs')

    r = 5 * u.arcminute
    result_table = Simbad.query_region(c, radius=r)
