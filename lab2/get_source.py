"""ASTR522 Lab 2: Do some fun stuff.

These functions are pretty catastrophically under-tested, and several known
problems exist, but both find* functions are functional under certain
conditions. Written for Python 2

Future work:
- Fix the broken stuff.
- Integrate everything into Source classes with built-in
  alt/az -> RA/Dec conversions and queried catalog info (magnitude, HD #, etc)
  integrated in. That'd be pretty cool.
- Run cost analysis on tuning the tunables. Would be pretty simple but cool.
-
"""
__author__ = "Jonas Powell"
__email__ = "jmpowell@wesleyan.edu"
__status__ = "Prototype"


import copy
import numpy as np
import pandas as pd
import subprocess as sp
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
from astroquery.vizier import Vizier
from astropy.time import Time
from datetime import datetime
from astropy import coordinates


# TUNABLE CONSTANTS
# Set the weights to be used in source selection (230)
c1 = 0.5
c2 = 1 - c1

# Set the latitude accuracy threshold (line 314)
lat_crit = 1

# Set how fine each grid should be (line 300)
step = 20


# Some pre-built source objects:
# Antares: Bad
antares = {'Name': 'Antares',
           'Alt': float(21.7),
           'Az': float(185.5),
           'Year': int(2018),
           'Month': int(11),
           'Day': int(6),
           'Hour': int(13),
           'Minute': int(36),
           'Second': int(0),
           'Lat': float(41.55),
           'Long': float(-72.65),
           'Timezone Offset': 5}

# Fomalhaut: Almost perfect
fomalhaut = {'Alt': 18.2,
             'Az': 170.1,
             'Day': 7,
             'Hour': 18,
             'Lat': 41.22,
             'Long': -72.65,
             'Minute': 57,
             'Month': 11,
             'Name': 'Fomalhaut',
             'Second': 0,
             'Timezone Offset': 5,
             'Year': 2018}

# Kochab: Extremely bad
kochab = {'Alt': 33.5,
          'Az': 342.4,
          'Day': 7,
          'Hour': 19,
          'Lat': 41.55,
          'Long': -72.65,
          'Minute': 5,
          'Month': 11,
          'Name': 'Kochab',
          'Second': 35,
          'Timezone Offset': 5,
          'Year': 2018}

# Ascella: Bad but in a weird way
ascella = {'Alt': 3.7,
           'Az': 224.2,
           'Day': 7,
           'Hour': 19,
           'Lat': 41.55,
           'Long': -72.65,
           'Minute': 18,
           'Month': 11,
           'Name': 'Ascella',
           'Second': 27,
           'Timezone Offset': 5,
           'Year': 2018}


# HELPER FUNCTIONS


def build_source_object():
    """Build an object to contain the observation info.

    This should really be a class, but time only allows this right now.
    """
    name = raw_input("If known, enter the source's name?\n")
    alt = raw_input("Enter source's altitude (in decimal degrees)\n")
    az = raw_input("Enter source's azimuth (in decimal degrees)\n")
    t_prompt = raw_input("Is the observation being made right now \
                         or in the past? (enter 'now' or anything else)\n")

    if t_prompt.lower() == 'now':
        now = datetime.now()
        second, minute, hour = now.second, now.minute, now.hour
        day, month, year = now.day, now.month, now.year
    else:
        second = 0
        obs_time = list(raw_input("Enter the date and time of the observation \
                                  (yyyy, mm, dd, hh, mm)\n"))
        minute, hour = obs_time[4], obs_time[3]
        day, month, year = obs_time[2], obs_time[1], obs_time[0]

    lat = float(raw_input("Enter the latitude of your observation \
                          (in decimal degrees)\n"))
    long = float(raw_input("Enter the longitude of your observation \
                           (in decimal degrees)\n"))
    # Should really just calculate this but lazy so too bad
    tz_offset = raw_input("Enter how many hours the observation's timezone is \
                          away from Greenwich")

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
           'Timezone Offset': tz_offset}
    return obs


def altaz_to_radec(obs):
    """Convert Alt/Az to RA/Dec.

    Args:
        obs (dict): an Observation object, created in build_source_object().

    Returns:
        ra_dec (tuple): RA and Dec in decimal degrees
    """
    # Retrieve the coordinates and convert them to rads for some trig.
    alt, az = obs['Alt'] * (np.pi/180), obs['Az'] * (np.pi/180)
    minute, hour, day = obs['Minute'], obs['Hour'], obs['Day']
    month, year, tz_offset = obs['Month'], obs['Year'], obs['Timezone Offset']
    lat, long = obs['Lat'] * (np.pi/180), obs['Long'] * (np.pi/180)

    # Do some time conversions
    if hour + tz_offset >= 24:
        hour = hour + tz_offset - 24
        day += 1
    utc = datetime(year, month, day, hour + tz_offset, minute, 0)
    T = Time(utc, scale='utc')
    gmst = T.sidereal_time('mean', 'greenwich').degree

    # The trig! Yuck.
    sin_dec = np.sin(alt) * np.sin(lat) + np.cos(alt) * np.cos(lat) * np.cos(az)
    dec = np.arcsin(sin_dec)

    cosHA = (np.sin(alt) - np.sin(lat) * np.sin(dec))/(np.cos(lat) * np.cos(dec))
    HA = np.arccos(cosHA) * (180/np.pi)

    dec *= (180/np.pi)
    if az < np.pi:
        ra = gmst + HA + (long * 180/np.pi)
    else:
        ra = gmst - HA + (long * 180/np.pi)

    if ra < 0:
        ra = 360 + ra

    ra_dec = (round(ra, 4), round(dec, 4))
    return ra_dec


# REAL FUNCTIONS

def find_source(obs, return_all_sources=True):
    """Find a source given some coordinates.

    Args:
        obs (dict): an Observation object, created in build_source_object().

    Returns:
        out (dict): a dictionary of the RA/dec coords we are observing and the
                    Henry Draper number of the star we think is most likely.
    """
    ra_dec = altaz_to_radec(obs)

    coords = coordinates.SkyCoord(ra=ra_dec[0], dec=ra_dec[1],
                                  unit=(u.deg, u.deg), frame='icrs')
    # Get the actual results
    # For some reason, if this goes too big it stops seeing the actual source?!
    r = 100 * u.arcminute
    results = Vizier.query_region(coords, radius=r, catalog='V/50')[0]
    df = results.to_pandas()

    candidate_sources = filter(None, [n for n in df['HD']])
    sources = []
    dmax, vmax = 0, 0
    for s in candidate_sources:
        source_info = df.loc[df['HD'] == s]
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

    # This loop is supremely janky, but df.loc'ing (below) wasn't working.
    # best_source = sources_df.loc[sources_df['Score'] == sources_df['Score'].min]
    best_source_idx = 0
    best_score = 10000
    for i in range(len(sources)):
        score = sources[i]['Score']
        if score < best_score:
            best_source_idx = i
            best_score = score

    print score
    out = {'Coords': ra_dec,
           'HD-Name': 'HD' + str(int(sources[best_source_idx]['HD'])),
           }
    return out


def find_location(obs, step=step, lat_crit=lat_crit, plot_steps=True):
    """Find out where we are on Earth using nested grid searches.

    Args:
        obs (dict): an Observation object, created in build_source_object().

    Returns:
        lat_long (tuple of floats): your location.

    Note that this is still extremely loose and has not been tested with:
        - Different geographic locations
        - Many different stars
        - Edge cases (at all)

    Known problems:
        - Badness in the altaz_to_radec() function that obviously affects this
        - Choosing the wrong gridding scale affects the results significantly
            (for Fomalhaut, steps=100 gives a terrible solution,
                            steps=10 gives an almost perfect solution,
                            steps=6 gives a close but not great solution.)
    """
    source_name = obs['Name']

    source_obj = Vizier.query_object(source_name, catalog='V/50')[0]
    source_ra_dec = (source_obj['RAJ2000'][0], source_obj['DEJ2000'][0])

    source_ra_hms = tuple(map(float, source_ra_dec[0].split()))
    source_dec_dms = tuple(map(float, source_ra_dec[1].split()))

    source_ra = Angle(source_ra_hms, unit='hourangle').degree
    source_dec = Angle(source_dec_dms, unit=u.deg).degree

    lats = np.linspace(-90, 90, step)
    longs = np.linspace(-180., 180., 2*step)
    lat_min, lat_max = 0, len(lats) - 1
    long_min, long_max = 0, len(longs) - 1

    # Set up a new directory to hold the output images in:
    if plot_steps is True:
        new_dir = 'LatLong-Evo-Plots_' + obs['Name'] + '_res' + str(step)
        sp.call(['rm', '-rf', '{}'.format(new_dir)])
        sp.call(['mkdir', '{}'.format(new_dir)])

    counter = 1
    # We want this loop to run as long as we haven't reached out critical
    # latitude resolution and we don't have degenerate solutions (line 376)
    while abs(lats[1] - lats[0]) > lat_crit:
        print "\nCurrently taking step number", counter

        lats = np.linspace(lats[lat_min], lats[lat_max], step)
        longs = np.linspace(longs[long_min], longs[long_max], 2*step)
        print "Latitude range:", str([round(min(lats), 2), round(max(lats), 2)])
        print "Longitude range:", str([round(min(longs), 2), round(max(longs), 2)])

        score_grid = np.zeros((len(lats), len(longs)))

        # Run the grid
        for i in range(len(lats)):
            for j in range(len(longs)):
                lat, long = lats[i], longs[j]
                new_obs = copy.deepcopy(obs)
                new_obs['Lat'] = lat
                new_obs['Long'] = long
                ra, dec = altaz_to_radec(new_obs)

                # Check how far away the sim'ed coords are from the real ones,
                # and use it as a goodness-of-fit metric.
                score = np.sqrt((ra - source_ra)**2 + (dec - source_dec)**2)
                score_grid[i, j] = score

        # Set the new lat/long window for the next step.
        idx = np.where(score_grid == np.nanmin(score_grid))
        lat_min = idx[0][0] - 2 if idx[0][0] > 1 else 0
        lat_max = idx[0][0] + 2 if idx[0][0] < len(lats) - 2 else len(lats) - 1

        long_min = idx[1][0] - 2 if idx[1][0] > 1 else 0
        long_max = idx[1][0] + 2 if idx[1][0] < len(longs) - 2 else len(longs) - 1

        # Plot this step out
        if plot_steps is True:
            plt.matshow(score_grid, cmap='magma_r')
            plt.contour(score_grid, cmap='magma')

            if obs['Lat'] != '' and obs['Long'] != '':
                i, j = 0, 0
                while longs[i] < obs['Long'] and i < len(longs) - 1:
                    i += 1
                while lats[j] < obs['Lat'] and j < len(lats) - 1:
                    j += 1
                # if i != 0 and j != 0 and i != len(longs) and i != len(lats):
                plt.plot([i], [j], 'or')

            # Give the ticks coordinate values
            xtick_locs = np.arange(0, len(longs), len(longs)/6)
            xtick_labs = [int(longs[x]) for x in xtick_locs]
            plt.xticks(xtick_locs, xtick_labs)
            ytick_locs = np.arange(0, len(lats), len(lats)/6)
            ytick_labs = [int(lats[y]) for y in ytick_locs]
            plt.yticks(ytick_locs, ytick_labs)

            plt.xlabel('Longitude', weight='bold')
            plt.ylabel('Latitude', weight='bold')

            plt.title("Grid Searching for "+obs['Name']+"; Step "+str(step),
                      weight='bold')
            plt.savefig(new_dir + '/window-evolution' + str(counter) + '.png')
            plt.close

        # If we start getting degenerate solutions, break the loop.
        # Further progress is meaningless.
        if len(idx[0]) > 1:
            print "Breaking loop"
            break
        counter += 1

    print "\n\nCompleted!"
    if plot_steps is True:
        print "See plots for visualization of progress in", new_dir
    return (lats[idx[0][0]], longs[idx[1][0]])


# The End
