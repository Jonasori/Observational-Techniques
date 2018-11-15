"""
Observational Techniques, Lab 3

Making an HR diagram and stuff.
"""

# Packages
import csv
import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
import astropy.constants as const
from sklearn.linear_model import LinearRegression
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia

# Constants
H0 = 72 * 1e-6                   # km s-1 pc-1
c = const.c.to('km/s').value     # km s-1


# PART 1
def plot_faber_jackson(save=False):
    """Plot the Faber-Jackson relationship.

    The Faber-Jackson relationship defines the relationship between luminosity
    and velocity dispersion for elliptical and S0 galaxies.
    """
    """
    The SQL call:
        SELECT TOP 10000
           p.petroMag_r, p.petroMagErr_r, p.extinction_r,
           p.petroMag_g, p.extinction_g,
           p.petroR50_r, p.petroR90_r, s.velDisp,
           p.dered_r, p.dered_g,
           p.objid, p.ra, p.dec, p.u, p.g, p.r, p.i, p.z,
           p.run, p.rerun, p.camcol, p.field,
           s.specobjid, s.class, s.z as redshift,
           s.plate, s.mjd, s.fiberid
        FROM PhotoObj AS p
           JOIN SpecObj AS s ON s.bestobjid = p.objid
        WHERE
           p.petroR90_r/p.petroR50_r > 2.6
           AND class='galaxy'
           AND p.dered_g - p.dered_r> 1.0
           AND petroMag_r BETWEEN 0 and 19
           AND petroMagErr_r < 0.05
           AND s.z < 0.35
           AND veldisp > 30
           AND veldispErr/veldisp < 0.2
    """
    plt.close()

    # Read in the data
    path = '/Volumes/1TB Storage Drive/Desktop/masters/Observational-Techniques/lab3/'
    raw_query = path + 'SDSS_query.csv'
    new_path = path + 'SDSS_corrected.csv'
    with open(raw_query, 'rb') as inp, open(new_path, 'wb') as outp:
        writer = csv.writer(outp)
        for row in csv.reader(inp):
            if row != ['#Table1']:
                writer.writerow(row)
    df = pd.read_csv(new_path, sep=',')
    df['Abs. Magnitude'] = (df['r'] - df['extinction_r']) - 5 * \
        (np.log10(c * df['redshift'] / H0) - 1)
    # df.columns
    # len(df['petroR50_r'])

    len_test = len(df['Abs. Magnitude'])/2
    len_train = len(df['Abs. Magnitude']) - len_test
    train_x = df['Abs. Magnitude'][:len_train]
    test_x = df['Abs. Magnitude'][len_test:]
    train_y = df['velDisp'][:len_train]
    test_y = df['velDisp'][:len_test]

    train_x = train_x.reshape(len(train_x), 1)
    test_x = test_x.reshape(len(test_x), 1)
    train_y = train_y.reshape(len(train_y), 1)
    test_y = test_y.reshape(len(test_y), 1)

    regr = LinearRegression()
    regr.fit(train_x, train_y)

    predicted_y = regr.predict(df['Abs. Magnitude'].reshape(len(df['Abs. Magnitude']), 1))

    df.plot.scatter('Abs. Magnitude', 'velDisp',
                    c='petroMag_r', cmap='jet',
                    logy=True, alpha=0.2)

    plt.plot(df['Abs. Magnitude'], predicted_y, '.k')

    plt.xlabel('Absolute Magnitude', weight='bold')
    plt.ylabel('Velocity Dispersion', weight='bold')
    plt.title('Faber-Jackson Relation (SDSS Data)', weight='bold')
    plt.gca().invert_xaxis()
    if save is True:
        plt.savefig('faber-jackson.png', dpi=200)
    else:
        plt.show(block=False)
# It'd be nice to have this as a density plot


plot_faber_jackson()


# PART 2
# Parameters by which to ID Pleiades objects.
ra_lims, dec_lims, plx_lims = (10, 30), (-50, -35), (5.5, 10)
# ra_lims, dec_lims, plx_lims = (-30, 80), (-100, -45), (5.5, 10)
# ra_lims, dec_lims, plx_lims = (-60, 80), (-100, -45), (5.5, 10)

# Compile the Hipparcos dataframe


def get_hipparcos_data():
    """Build a data structure for the Hipparcos data."""
    print "Building Hipparcos dataframe."

    coords = SkyCoord(ra='3h47m24s', dec='+24d7m0s',
                      unit=(u.deg, u.deg), frame='icrs')
    r = 2 * u.degree
    hip_results_raw = Vizier.query_region(coords, radius=r,
                                          catalog='I/239/hip_main')[0]
    hip_results = hip_results_raw.to_pandas()

    # Get a sorted list of the columns for reference
    hip_cols = sorted(hip_results.columns)

    hip_df = pd.DataFrame()
    hip_df['Proper Motion (RA)'] = hip_results['pmRA']
    hip_df['Proper Motion (Dec)'] = hip_results['pmDE']
    hip_df['Distance'] = (hip_results['Plx'] * 1e-3)**(-1)
    hip_df['Parallax'] = hip_results['Plx']
    hip_df['Plx. Error'] = hip_results['e_Plx']
    hip_df['mag'] = hip_results['Vmag']
    hip_df['Color'] = hip_results['B-V']
    hip_df['Absolute Magnitude'] = hip_df['mag'] - \
        5 * (np.log10(hip_df['Distance']) - 1)
    hip_df['T Effective'] = [0] * len(hip_df['Color'])
    hip_df['Parallax'] = hip_results['Plx']
    hip_df['Plx. Error'] = hip_results['e_Plx']
    hip_df['Confidence'] = 1 - hip_results['e_Plx']/max(hip_results['e_Plx'])


    # Subset the data based on proper motion and Parallax cuts
    ra_cond1 = hip_results['pmRA'] > ra_lims[0]
    ra_cond2 = hip_results['pmRA'] < ra_lims[1]
    dec_cond1 = hip_results['pmDE'] > dec_lims[0]
    dec_cond2 = hip_results['pmDE'] < dec_lims[1]
    plx_cond1 = hip_results['Plx'] > plx_lims[0]
    plx_cond2 = hip_results['Plx'] < plx_lims[1]

    hip_df = hip_df[ra_cond1 & ra_cond2 & dec_cond1 &
                    dec_cond2 & plx_cond1 & plx_cond2]
    # hip_df.sort_values('Distance')
    pleiades_hip = {'Survey': 'Hipparcos',
                    'Mean Distance': round(np.mean(hip_df['Distance']), 1),
                    'Number of Stars': len(hip_df['Distance']),
                    'text_loc1': (0.25, -1.8),
                    'text_loc2': (0.25, -0.8),
                    'Data': hip_df}
    return pleiades_hip


pleiades_hip = get_hipparcos_data()


# Compile the GAIA dataframe
def get_gaia_data(n_sources=5000):
    """Build a data structure for the Gaia data.

    Really need a better/additional way to identify stars in the pleiades.

    Documentation on the GAIA columns:
    http://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/
    chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html

    'For relations using V-RC, RC-IC, or B-V or for relations linking Sloan
    magnitudes (g or r) and colours (g-r, g-i, or r-i) to Gaia G magnitudes,
    please see arxiv.org/pdf/1008.0815.pdf'

    Basically, since GAIA uses longer wavelength light, their Colors don't map
    to UBVRI bandpasses (p. 4, Table 1).

    Age between 80 and 100 Myr
    https://www.universetoday.com/34994/messier-45-pleiades/
    """
    print "Building GAIA dataframe for {} sources.".format(n_sources)

    gaia_str = "SELECT top {} * FROM gaiadr2.gaia_source \
                WHERE pmra between {} and {} \
                AND pmdec between {} and {} \
                AND parallax between {} and {} \
                AND parallax_error < 2 \
                AND parallax_over_error > 5 \
                ".format(n_sources,
                         ra_lims[0], ra_lims[1],
                         dec_lims[0], dec_lims[1],
                         plx_lims[0], plx_lims[1])

    job = Gaia.launch_job(gaia_str, dump_to_file=True)
    # job = Gaia.launch_job(gaia_str)
    gaia_results_raw = job.get_results()
    # gaia_results_raw['phot_rp_mean_mag'].description

    gaia_results = gaia_results_raw.to_pandas()
    # gaia_cols = sorted(gaia_results.columns)

    print "Acquired data; now building dataframe..."
    gaia_df = pd.DataFrame()
    gaia_df['Distance'] = (gaia_results['parallax'] * 1e-3)**(-1)
    gaia_df['Proper Motion (RA)'] = gaia_results['pmra']
    gaia_df['Proper Motion (Dec)'] = gaia_results['pmdec']
    gaia_df['mag'] = gaia_results['phot_rp_mean_mag']
    gaia_df['Color'] = gaia_results['bp_rp']
    gaia_df['Absolute Magnitude'] = gaia_df['mag'] - \
        5 * (np.log10(gaia_df['Distance']) - 1)
    gaia_df['T Effective'] = gaia_results['teff_val']
    gaia_df['Parallax'] = gaia_results['parallax']
    gaia_df['Plx. Error'] = gaia_results['parallax_error']
    gaia_df['Confidence'] = 1 - gaia_results['parallax_error']/max(gaia_results['parallax_error'])

    pleiades_gaia = {'Survey': 'Gaia',
                     'Mean Distance': round(np.mean(gaia_df['Distance']), 1),
                     'Number of Stars': len(gaia_df['Distance']),
                     'text_loc1': (1.1, -2.2),
                     'text_loc2': (2, -0.2),
                     'Data': gaia_df,
                     'Full Results': gaia_results,
                     'Full Table': gaia_results_raw}

    return pleiades_gaia


pleiades_gaia = get_gaia_data()

cols = sorted(df.columns)
cols


# A generic plotter for the two datasets
def make_plots(survey, save=False):
    """Make the necessary plots for each dataset."""
    plt.close()
    print "Beginning plot-making process for " + survey['Survey'] + "..."

    df = survey['Data']
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

    df.plot.scatter('Proper Motion (RA)', 'Proper Motion (Dec)',
                    ax=axes[0],
                    c='Confidence', cmap='jet', colorbar=False, alpha=0.5
                    )

    # ylims = (max(df['Absolute Magnitude']), min(df['Absolute Magnitude']))
    df.plot.scatter('Color', 'Absolute Magnitude',
                    ax=axes[1],
                    c='Confidence', cmap='jet', alpha=0.5
                    )

    text_str1 = 'Number of stars: ' + str(survey['Number of Stars'])
    text_str2 = 'Mean distance:\n' + str(survey['Mean Distance']) + ' parsecs'

    axes[1].annotate(text_str1, survey['text_loc1'], weight='bold')
    axes[1].annotate(text_str2, survey['text_loc2'], weight='bold')
    axes[1].set_ylim(axes[1].get_ylim()[::-1])

    axes[0].set_title('Proper Motions', weight='bold')
    axes[1].set_title('HR Diagram', weight='bold')
    fig.suptitle('Pleiades Cluster Characteristics, from ' + survey['Survey'] + '\n\n',
                 weight='bold', fontsize=16)

    if save is True:
        outname = './pleiades_' + survey['Survey'] + '.png'
        plt.savefig(outname, dpi=200)
        print "Saved to ", outname
    else:
        print "Showing:"
        plt.show()


# make_plots(pleiades_gaia, save=False)


# The End
