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
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from datetime import datetime
from astropy.table import Column
import matplotlib.pyplot as plt


# Some sample things at a sample time
hour, minute = 13, 36
day, month, year = 6, 11, 2018
antares_altaz = (21.7, 185.5)   # 16h31m, -27deg
altair_altaz = (38.4, 117.1)    # 19h51m, 8deg

local_latlong = (41.552, -72.65)

# Used to determine pixel/degree scale in grid search.
res = 1.5

# HELPER FUNCTIONS

def localtime_to_gmst(minute=minute, hour=hour, day=day, month=month, year=year, tz_offset=5):
    # Implement an option to just use current time.
    today = datetime.now()

    utc = datetime(year, month, day, hour + tz_offset, minute, 0)
    T = Time(utc, scale='utc')
    gmst = T.sidereal_time('mean', 'greenwich').degree
    return gmst


def plot_results(mat, name):
    # This stuff isn't right yet.
    lat_coord = (90 + local_latlong[0]) * res
    long_coord = (180 + local_latlong[1]) * res

    plt.matshow(mat, cmap='magma')
    plt.contour(mat)
    plt.plot([lat_coord], [long_coord], 'or')
    plt.savefig(name + '-' + str(res) + '.png', dpi=200)
    plt.show(block=False)


# PROBLEM 1
def altaz_to_radec(alt_az, pos=local_latlong,
                   minute=minute, hour=hour, day=day,
                   month=month, year=year, tz_offset=5):
    """Convert Alt/Az to RA/Dec.

    Args:
        alt_az (tuple of floats): alt, az in decimal degrees
        pos (tuple of floats): local lat, long in decimal degrees
        hour, minute: ints of the desired time.

    Returns:
        ra_dec (tuple): RA and Dec in decimal degrees

    Currently assumes all observations are happening today
    """
    # Retrieve the coordinates and convert them to rads for some trig.
    lat, long = pos[0] * (np.pi/180), pos[1] * (np.pi/180)
    alt, az = alt_az[0] * (np.pi/180), alt_az[1] * (np.pi/180)

    gmst = localtime_to_gmst(minute=minute, hour=hour,
                             day=day, month=month, year=year, tz_offset=5)

    sin_dec = np.sin(alt) * np.sin(lat) + np.cos(alt) * np.cos(lat) * np.cos(az)
    dec = np.arcsin(sin_dec)

    cosHA = (np.sin(alt) - np.sin(lat) * np.sin(dec))/(np.cos(lat) * np.cos(dec))
    HA = np.arccos(cosHA) * (180/np.pi)

    dec *= (180/np.pi)
    ra = gmst + HA + (long * 180/np.pi) if az < np.pi else gmst - HA + (long * 180/np.pi)

    ra_dec = (round(ra, 4), round(dec, 4))
    return ra_dec



def find_source(alt_az, lat_lon=local_latlong,
                minute=minute, hour=hour,
                day=day, month=month, year=year, tz_offset=5,
                return_all_sources=True):
    """Find a source given some coordinates.

    Args:
        alt_az: Tuple of altitude and azimuth (decimal degrees)
        lat_lon: Tuple of latitude and longitude (decimal degrees)
        time: Clock time (?)

    Returns:
        source_name (str):
    """

    alt, az = alt_az[0], alt_az[1]
    lat, lon = lat_lon[0], lat_lon[1]

    # my_loc = coordinates.EarthLocation.from_geodetic(lat, lon)
    ra_dec = altaz_to_radec(alt_az, pos=local_latlong,
                            minute=minute, hour=hour, day=day,
                            month=month, year=year, tz_offset=5)

    # c = coordinates.SkyCoord("05h35m17.3s -05d23m28s", frame='icrs')

    coords = coordinates.SkyCoord(ra=ra_dec[0], dec=ra_dec[1],
                                  unit=(u.deg, u.deg), frame='icrs')
    # Get the actual results
    r = 100 * u.arcminute
    results = Vizier.query_region(coords, radius=r, catalog='V/50')[0]
    df = results.to_pandas()

    candidate_sources = filter(None, [n for n in df['Name']])
    sources = []
    dmax, vmax = 0, 0
    for s in candidate_sources:
        print s, '\n\n'
        source_info = df.loc[df['Name'] == s]
        print source_info['Vmag']
        mag = round(float(source_info['Vmag']), 2)
        source_ra_hms = tuple(map(float, source_info['RAJ2000'][0].split()))
        source_dec_dms = tuple(map(float, source_info['DEJ2000'][0].split()))

        source_ra = Angle(source_ra_hms, unit='hourangle').degree
        source_dec = Angle(source_dec_dms, unit=u.deg).degree

        dist_from_center = np.sqrt((source_ra - ra_dec[0])**2 +
                                   (source_dec - ra_dec[1])**2)

        c1, c2 = 0.5, 0.5
        score = c1 * mag + c2 * dist_from_center
        source_dict = {'Name': source_info['Name'],
                       'RA': source_ra,
                       'DEC': source_dec,
                       'Distance': dist_from_center,
                       'Vmag': source_info['Vmag'],
                       'Score': score}

        sources.append(source_dict)
        dmax = dist_from_center if dist_from_center > dmax else dmax
        vmax = mag if mag > vmax else mag

    sources_df = pd.DataFrame(sources)
    best_source = sources_df[sources_df['score'] == sources_df['score'].min()]

    return sources_df




def find_source_astropy(alt_az, lat_lon, time, return_all_sources=False):
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





def find_location(source_name, source_alt_az,
                  minute, hour, day, month, year,
                  plot_grids=True):
    """Find out where we are on Earth.

    Args:
        source_name (str):
        source_ra_dec (tuple of floats):
        minute, hour, day, month, year (ints):

    Returns:
        lat_long (tuple of floats): your location.
    """

    alt, az = source_alt_az
    source_obj = Vizier.query_object(source_name, catalog='V/50')[0]
    source_ra_dec = (source_obj['RAJ2000'][0], source_obj['DEJ2000'][0])

    source_ra_hms = tuple(map(float, source_ra_dec[0].split()))
    source_dec_dms = tuple(map(float, source_ra_dec[1].split()))

    source_ra = Angle(source_ra_hms, unit='hourangle').degree
    source_dec = Angle(source_dec_dms, unit=u.deg).degree

    lats = np.arange(-90., 90, res)
    longs = np.arange(-180, 180, res)

    ra_grid = np.zeros((len(lats), len(longs)))
    dec_grid = np.zeros((len(lats), len(longs)))
    score_grid = np.zeros((len(lats), len(longs)))

    lat_counter, long_counter = 0, 0
    for i in range(len(lats)):
        for j in range(len(longs)):
            # Need to sort out angular units
            lat, long = lats[i], longs[j]

            ra, dec = altaz_to_radec((alt, az), pos=(lat, long),
                                     minute=minute, hour=hour, day=day,
                                     month=month, year=year, tz_offset=5)

            # pos_grid[i, j] = {'RA': ra, 'DEC': dec}
            ra_grid[i, j] = ra
            dec_grid[i, j] = dec

            # Bad - planar:
            score = np.sqrt((ra - source_ra)**2 + (dec - source_dec)**2)

            # Good - spherical:
            # score = np.arccos(np.sin(dec) * np.sin(source_dec) + np.cos(dec) * np.cos(source_dec) * np.cos(abs(ra - source_ra)))

            score_grid[i, j] = score

            verbose = False
            if verbose is True:
                print('RA, Source RA:', ra, source_ra)
                print('DEC, Source DEC:', dec, source_dec)
                print('Score:', score)
                print('\n')
            else:
                step = long_counter + lat_counter * len(lats)
                print (str(step) + '/' + str(len(lats) * len(longs)))
            long_counter += 1

    outname = 'latlong-gridsearch-results_' + str(res)
    score_df = pd.DataFrame(score_grid)
    score_df.to_csv(outname + '.csv')

    if plot_grids is True:
        lat_coord = (90 + local_latlong[0]) * res
        long_coord = (180 + local_latlong[1]) * res

        plt.matshow(score_grid, cmap='magma')
        plt.contour(score_grid)
        plt.plot([lat_coord], [long_coord], 'or')
        plt.savefig(outname + '.png', dpi=200)
        plt.show(block=False)


    return {'RA': ra_grid, 'DEC': dec_grid, 'SCORE': score_grid}












# The End
