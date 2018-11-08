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
res = 0.1
c1 = 0.5
c2 = 1 - c1

# HELPER FUNCTIONS

def build_source_object():
    # raw_input("Any of the following prompts may be skipped, although some functionality may be lost by doing so. Press Enter to continue.")
    name = raw_input("What is the source's name?\n")
    alt = raw_input("Enter source's altitude (in decimal degrees)\n")
    az = raw_input("Enter source's azimuth (in decimal degrees)\n")
    t_prompt = raw_input("Is the observation being made right now or in the past? (enter 'now' or anything else)\n")

    if t_prompt.lower() == 'now':
        now = datetime.now()
        second, minute, hour = now.second, now.minute, now.hour_now
        day, month, year = now.day, now.month, now.year
    else:
        second = 0
        obs_time = list(raw_input("Enter the date and time of the observation (yyyy, mm, dd, hh, mm)\n"))
        minute, hour = obs_time[4], obs_time[3]
        day, month, year = obs_time[2], obs_time[1], obs_time[0]

    lat = float(raw_input("Enter the latitude of your observation (in decimal degrees)\n"))
    long = float(raw_input("Enter the longitude of your observation (in decimal degrees)\n"))

    # This is tunable, just not yet implemented as such.
    tz_offset = 5

    obs = {'Name': name,
           'Alt': float(alt),
           'Az': float(az),
           'Year': int(year),
           'Month': int(month),
           'Day': int(day),
           'Hour': int(hour),
           'Minute': int(minute),
           'Second': int(second),
           'Lat': float(lat),
           'Long': float(long),
           'Timezone Offset': 5}
    return obs

# Some pre-built source objects
antares = {'Name': 'Antares',
           'Alt': float(21.7),
           'Az': float(a185.5z),
           'Year': int(2018),
           'Month': int(11),
           'Day': int(6),
           'Hour': int(13),
           'Minute': int(36),
           'Second': int(0),
           'Lat': float(41.55),
           'Long': float(-72.65),
           'Timezone Offset': 5}




def localtime_to_gmst(minute, hour, day, month, year, tz_offset=5):
    utc = datetime(year, month, day, hour + tz_offset, minute, 0)
    T = Time(utc, scale='utc')
    gmst = T.sidereal_time('mean', 'greenwich').degree
    return gmst


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
    ra_dec = altaz_to_radec(alt_az, pos=lat_lon,
                            minute=minute, hour=hour, day=day,
                            month=month, year=year, tz_offset=5)

    coords = coordinates.SkyCoord(ra=ra_dec[0], dec=ra_dec[1],
                                  unit=(u.deg, u.deg), frame='icrs')
    # Get the actual results
    # For some reason, if this goes too big it stops seeing the actual source.
    r = 500 * u.arcminute
    results = Vizier.query_region(coords, radius=r, catalog='V/50')[0]
    df = results.to_pandas()

    candidate_sources = filter(None, [n for n in df['HD']])
    sources = []
    dmax, vmax = 0, 0
    for s in candidate_sources:
        source_info = df.loc[df['HD'] == s]
        name = source_info['Name']
        mag = round(float(source_info['Vmag']), 2)

        temp_ra = source_info['RAJ2000'].tolist()[0]
        temp_dec = source_info['DEJ2000'].tolist()[0]
        source_ra_hms = tuple(map(float, temp_ra.split()))
        source_dec_dms = tuple(map(float, temp_dec.split()))
        source_ra = Angle(source_ra_hms, unit='hourangle').degree
        source_dec = Angle(source_dec_dms, unit=u.deg).degree

        dist_from_center = np.sqrt((source_ra - ra_dec[0])**2 +
                                   (source_dec - ra_dec[1])**2)

        score = float(c1 * mag + c2 * dist_from_center)
        source_dict = {'HD': source_info['HD'].values[0],
                       'Name': source_info['Name'].values[0],
                       'RA': source_ra,
                       'DEC': source_dec,
                       'Distance': dist_from_center,
                       'Vmag': source_info['Vmag'],
                       'Score': score}

        sources.append(source_dict)

        dmax = dist_from_center if dist_from_center > dmax else dmax
        vmax = mag if mag > vmax else mag

    for s in range(len(sources)):
        d = sources[s]['Distance']/dmax
        mag = sources[s]['Vmag'].values[0]/vmax
        score = c1 * mag + c2 * d
        sources[s]['Score'] = score
        sources[s]['Scaled-Distance'] = d
        sources[s]['Scaled-Mag'] = mag

    sources_df = pd.DataFrame(sources)


    # Note that this loop is supremely janky, but df.loc'ing wasn't working.
    # best_source = sources_df.loc[sources_df['Score'] == sources_df['Score'].min]
    best_source_idx = 0
    # best_score = np.array([])
    best_score = 10000
    for i in range(len(sources)):
        score = sources[i]['Score']
        if score < best_score:
            best_source_idx = i
            best_score = score

    name = sources_df['Name'].values[0]
    out = {'Coords': ra_dec,
           'HD-Name': 'HD' + str(int(sources[best_source_idx]['HD'])),
           'Name': sources[best_source_idx]['Name'],
           'Score': sources[best_source_idx]['Score'],
           'Scaled-Distance': sources[best_source_idx]['Scaled-Distance'],
           'Scaled-Mag': sources[best_source_idx]['Scaled-Mag']
           }
    return out




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










now = datetime.now()

second_now, minute_now, hour_now = now.second, now.minute, now.hour
day_now, month_now, year_now = now.day, now.month, now.year

def find_location(source_name, source_alt_az,
                  minute=minute_now, hour=hour_now, day=day_now,
                  month=month_now, year=year_now, plot_grids=True):
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


    lat_crit = 0.01
    step = 20
    lats = np.linspace(-90, 90, step)
    longs = np.linspace(-180., 180., 2*step)
    lat_min, lat_max = 0, len(lats) - 1
    long_min, long_max = 0, len(longs) - 1

    real_lat, real_long = 41.55, -72.65

    counter = 0
    while abs(lats[1] - lats[0]) > lat_crit:
        lats = np.linspace(lats[lat_min], lats[lat_max], step)
        longs = np.linspace(longs[long_min], longs[long_max], 2*step)

        score_grid = np.zeros((len(lats), len(longs)))

        # Run the grid
        for i in range(len(lats)):
            for j in range(len(longs)):
                lat, long = lats[i], longs[j]

                ra, dec = altaz_to_radec((alt, az), pos=(lat, long),
                                         minute=minute, hour=hour, day=day,
                                         month=month, year=year, tz_offset=5)

                score = np.sqrt((ra - source_ra)**2 + (dec - source_dec)**2)
                score_grid[i, j] = score

        idx = np.where(score_grid == np.nanmin(score_grid))
        lat_min = idx[0][0] - 2 if idx[0][0] > 1 else 0
        lat_max = idx[0][0] + 2 if idx[0][0] < len(lats) - 2 else len(lats) - 1

        long_min = idx[1][0] - 2 if idx[1][0] > 1 else 0
        long_max = idx[1][0] + 2 if idx[1][0] < len(longs) - 2 else len(longs) - 1

        plt.matshow(score_grid, cmap='magma_r')
        plt.contour(score_grid, cmap='magma')
        xtick_locs = np.arange(0, len(longs), len(longs)/6)
        xtick_labs = [int(longs[i]) for i in xtick_locs]
        plt.xticks(xtick_locs, xtick_labs)

        ytick_locs = np.arange(0, len(lats), len(lats)/6)
        ytick_labs = [int(lats[i]) for i in ytick_locs]
        plt.yticks(ytick_locs, ytick_labs)

        i, j = 0, 0
        while longs[i] < real_long and i < len(longs) - 1:
            i += 1
        while lats[j] < real_lat and j < len(lats) - 1:
            j += 1

        print i, j
        if i != 0 and j != 0 and i != len(longs) and i != len(lats):
            plt.plot([i], [j], 'or')
        plt.savefig('window-evolution' + str(counter) + '.png')
        plt.close

        counter += 1

    return (lats[idx[0][0]], longs[idx[1][0]])
















    outname = 'latlong-gridsearch-results_' + str(res)
    score_df = pd.DataFrame(score_grid)
    score_df.to_csv(outname + '.csv')

    if plot_grids is True:
        lat_coord = (90 + local_latlong[0]) * res
        long_coord = (180 + local_latlong[1]) * res

        plt.contour(score_grid)
        plt.plot([lat_coord], [long_coord], 'or')
        plt.matshow(score_grid, cmap='magma')

        xtick_locs = np.arange(0, len(longs), len(longs)/6)
        xtick_labs = [int(longs[i]) for i in xtick_locs]
        plt.xticks(xtick_locs, xtick_labs)

        # plt.ylim(max(lats), min(lats))
        ytick_locs = np.arange(0, len(lats), len(lats)/10)
        ytick_labs = [int(lats[i]) for i in ytick_locs]
        plt.yticks(ytick_locs, ytick_labs)

        plt.savefig(outname + '.png', dpi=200)
        plt.show(block=False)


    return {'RA': ra_grid, 'DEC': dec_grid, 'SCORE': score_grid}











def find_location_gs(source_name, source_alt_az,
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

    # Run the grid
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

        plt.contour(score_grid)
        plt.plot([lat_coord], [long_coord], 'or')
        plt.matshow(score_grid, cmap='magma')

        xtick_locs = np.arange(0, len(longs), len(longs)/6)
        xtick_labs = [int(longs[i]) for i in xtick_locs]
        plt.xticks(xtick_locs, xtick_labs)

        # plt.ylim(max(lats), min(lats))
        ytick_locs = np.arange(0, len(lats), len(lats)/10)
        ytick_labs = [int(lats[i]) for i in ytick_locs]
        plt.yticks(ytick_locs, ytick_labs)

        plt.savefig(outname + '.png', dpi=200)
        plt.show(block=False)


    return {'RA': ra_grid, 'DEC': dec_grid, 'SCORE': score_grid}












# The End
