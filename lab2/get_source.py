"""ASTR522 Lab 2: Do some fun stuff.

Should learn to use novas pip package - does some great stuff.
http://aa.usno.navy.mil/software/novas/novas_py/novaspy_intro.php
"""

import astroquery
# import julian
import astropy.units as u
import numpy as np
import pandas as pd
from astropy import coordinates
from astropy.time import Time
from astroquery.simbad import Simbad
from astropy.coordinates import Angle
from datetime import datetime
from astropy.table import Column


# For Mirphak
# alt1 = Angle((25, 46, 35.13), unit=u.deg)
# az1 = Angle((357, 23, 44.27), unit=u.deg)
alt_az_init = (Angle((77, 29, 18.63), unit=u.deg).value,
               Angle((43, 27, 37.09), unit=u.deg).value)
lat_lon_init = (41.552, -72.65)
time_init = datetime(2018, 11, 5, 23, 24, 42, 878762)

# Get right now in MJD (JD is also an option)
# https://astropy.readthedocs.io/en/v0.3/time/index.html
# time = Time(datetime.now()).mjd


def altaz_to_radec_by_hand(alt_az, pos, time):
    """Blah.

    alt_az ():
    pos ():
    time: LST
    date
    """
    # This seems useful:
    # http://www.brayebrookobservatory.org/BrayObsWebSite/BOOKS/matrix_method_rev_d.pdf
    lat, long = pos
    alt, az = alt_az

    # Still need to convert time to the right format

    sin_dec = np.sin(alt) * np.sin(lat) + np.cos(alt) * np.cos(lat) * np.cos(az)
    dec = np.arcsin(sin_dec)

    cosHA = (np.sin(alt) - np.sin(lat) * np.sin(dec))/(np.cos(lat) * np.cos(dec))
    HA = np.arccos(cosHA)
    RA = time - HA

    return [RA, dec]


def find_source(alt_az, lat_lon, time, return_all_sources=True):
    """Blah.

    Args:
        alt_az: Tuple of altitude and azimuth (decimal degrees)
        lat_lon: Tuple of latitude and longitude (decimal degrees)
        time: Clock time (?)
    """

    time = Time(time)
    alt, az = alt_az[0], alt_az[1]
    lat, lon = lat_lon[0], lat_lon[1]

    my_loc = coordinates.EarthLocation.from_geodetic(lat, lon)
    alt, az = coordinates.Angle(alt, unit=u.deg), coordinates.Angle(az, unit=u.deg)

    alt_az = coordinates.SkyCoord(frame='altaz', alt=alt, az=az,
                                  obstime=time, location=my_loc)
    ra_dec = alt_az.transform_to('icrs')



    # c = coordinates.SkyCoord("05h35m17.3s -05d23m28s", frame='icrs')
    """
    # An example
    c = coordinates.SkyCoord("05h35m17.3s -05d23m28s", frame='icrs')
    r = 5 * u.arcminute
    result_table = simbad.query_region(c, radius=r)
    """

    # Get the actual results
    r = 5 * u.arcminute
    results = Simbad.query_region(ra_dec, radius=r)
    # results_pdfdf = results.to_pandas()

    # Can do a bunch of stuff with this now. Mostly look at
    # result_table['MAIN_ID'], result_table.colnames

    # Calculate a goodness score for each captured entry
    # Currently just using Pythagorean distance, but would like to add in a
    # magnitude element as well.
    sources = []
    obs_ra = [ra_dec.ra.hms.h, ra_dec.ra.hms.m, ra_dec.ra.hms.s]
    obs_dec = [ra_dec.dec.dms.d, ra_dec.ra.dms.m, ra_dec.ra.dms.s]
    for i in range(len(results['RA'])):
        # Put things in a format to compare them
        ra = [float(c) for c in str(results['RA'][i]).split()]
        dec = [float(c) for c in str(results['DEC'][i]).split()]
        dRA = np.array(obs_ra) - np.array(ra)
        dDEC = np.array(obs_dec) - np.array(dec)
        dRA_angle = Angle(tuple(dRA), unit='hourangle').value
        dDEC_angle = Angle(tuple(dDEC), unit=u.deg).value
        # This might be preferential to one of the directions bc
        # a degree in hms != a degree in dms?
        d = {'name': results['MAIN_ID'][i],
             'score': np.sqrt(dRA_angle**2 + dDEC_angle**2)
             }
        sources.append(d)
    sources_df = pd.DataFrame(sources)
    best_source = sources_df[sources_df['score'] == sources_df['score'].min()]

    if return_all_sources is True:
        return ((ra_dec.ra.hms, ra_dec.dec.dms), sources_df)
    else:
        return ((ra_dec.ra.hms, ra_dec.dec.dms), best_source['name'])


def find_location(source_name, source_ra_dec, time):
    r = 5 * u.arcminute
    results = Simbad.query_objectids(source_name)











# The End
