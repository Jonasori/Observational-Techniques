"""Observational Techniques, PS 1.
https://moodle1819.wesleyan.edu/pluginfile.php/268955/mod_resource/content/0/hw1.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import G, M_earth, R_earth, c
r_earth = R_earth.value
c = c.value

pos1 = [19.80159, 155.45581]
pos2 = [17.75652, 64.58376]


### PART A ###
def find_distance(pos1, pos2, to_return='arc'):
    """Blah.
    pos1, pos2 (tuple): Lat, Long in decimal degrees.
    """
    # https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    phi1, lam1 = pos1
    phi2, lam2 = pos2

    phi1 *= np.pi/180
    phi2 *= np.pi/180
    lam1 *= np.pi/180
    lam2 *= np.pi/180
    d_lam = abs(lam2 - lam1)

    d_sig = np.arccos(np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(d_lam))
    sig_error = 0.1 * np.pi/180

    error_angle, error_radius = 1e-5, 0
    error_coeff = np.sqrt((error_angle/d_sig)**2 + (error_radius/r_earth)**2)

    if to_return == 'cord':
        cord = 2 * r_earth * np.sin(d_sig/2)
        return [cord, cord * error_coeff]
    elif to_return == 'arc':
        distance = r_earth * d_sig
        return [distance, distance * error_coeff]
    elif to_return == 'angle':
        return [d_sig, sig_error]
    else:
        return "Invalid return request; choose 'arc', 'cord', or 'angle'."

find_distance(pos1, pos2, to_return='arc')
find_distance(pos1, pos2, to_return='cord')
find_distance(pos1, pos2, to_return='angle')


def find_angres_error(lam, sig_lam, d, sig_d):
    theta = lam/d
    sig_theta = theta * np.sqrt((sig_lam/lam)**2 + (sig_d/d)**2)
    return sig_theta

find_angres_error(0.1, 0, r_earth, find_distance(pos1, pos2)[1])



### PART B ###
def get_photon_delay(sig, declination):
    """Blah.

    sig (float): angle of separation between coordinates
    """
    cord = find_distance(pos1, pos2, return_cord)[0]
    dt = d/c



### PROBLEM 2 ###
def get_relativistic_stuff(v_sat, h_sat):
    special_stuff = v_sat**2 / (2 * c**2)
    general_stuff = G.value * M_earth.value * c**(-2) * (1/R_earth.value - 1/(R_earth.value + h_sat))
    all_stuff = special_stuff + general_stuff
    all_stuff = all_stuff * np.pi * 1e7
    return all_stuff

h_sat = 540 * 1e3
period = 95 * 60 + 28
v_sat = h_sat/period
sat_circumference = 2 * np.pi * (r_earth + h_sat)
get_relativistic_stuff(v_sat, h_sat)


def get_final_time_delay(t=3600):
    R = R_earth.value
    H = 5.4e5
    alt = 4205

    theta_MK = 2 * np.pi * t/86400
    theta_HST_dt = 2 * np.pi * t/5728
    theta_HST = 2*np.pi - ((np.pi/2) + np.arccos(R_earth.value/(R_earth.value + H)) + theta_HST_dt)
    theta_i = theta_MK + theta_HST

    alpha = np.arctan((R + alt) * np.sin(theta_i) / (R + H - (R + alt) * np.cos(theta_i)))
    theta = alpha + theta_HST - np.pi/2
    l = (R + H - (R + alt) * np.cos(theta_i))/(np.cos(alpha))

    d = l * np.sin(theta)
    dt = d/c
    return [d/R_earth.value, dt]


get_final_time_delay()









# The End
