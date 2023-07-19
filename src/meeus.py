# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 16:49:32 2022

@author: Brian

Handles the Meeus ephemeris for the planetary bodies as a dictionary of the form:
{Planet name:
    {L (Mean Longitude): [vals in deg],
     a (semi-major axis): [vals in AU],
     e (eccentricity): [vals],
     i (inclination): [vals in deg],
     LAN aka OMEGA (Longitude of Ascending Node): [vals in deg],
         aka RAAN (Right Angle of Ascending Node)
     LoP aka PI (Longitude of Perihelion): [vals in deg]}
    ...}
    
Notes:
    longitude of Perihelion = Longitude of Ascending Node + argument of perihelion
    Mean Anomaly (aka M) = Mean Longitude - Longitude of Perihelion
    True Anomaly (aka v or f) = Mean Anomaly + C_cen
     
"""

import numpy as np

from .util import *

def JDE_to_T(JDE):
    T = (JDE - J2000) / J_cent
    return T

def sind(deg_in):
    rad_in = np.deg2rad(deg_in)
    rad_out = np.sin(rad_in)
    deg_out = np.rad2deg(rad_out)
    return deg_out
    
def C_cen(e, M):
    return ((2*e - 1/4*e**3 + 5/96*e**5)*sind(M) + (5/4*e**2 - 11/24*e**4)*sind(2*M)
            + (13/12*e**3 - 43/64*e**5)*sind(3*M)
            + 103/96*e**4*sind(4*M) + 1097/960*e**5*sind(5*M))

def get_TA(L, e, LoP, **kwargs):
    M = L - LoP
    TA = M + C_cen(e, M)
    return TA
    
def get_elements(planet, T, angle_units='rad'):
    """Gets the orbital elements

    Parameters
    ----------
    planet : string
        Which planet to get orbital elements for
    T : numeric
        Julian time in centuries from epoch J2000.0
    angle_units : string, optional, either "rad" or "deg"
        What units to use for the returned angle values
        "rad" for radians (default)
        "deg" for degrees

    Raises
    ------
    ValueError
        Invalid planet name

    Returns
    -------
    values : dict
        The orbital element values, keyed by element name for passing into Orbit:
            e, a, i, LAN, AoP, f, tau
    """
    planet = planet.capitalize()
    if planet not in meeus:
        raise ValueError(f"Not a known planet: {planet}")
    coeffs = meeus[planet]
    values = {}
    # Direct application of coefficients
    for k,v in coeffs.items():
        values[k] = v[0] + v[1]*T + v[2]*T**2 + v[3]*T**3
    # Derived value
    values['f'] = get_TA(**values)
    values['AoP'] = values['LoP'] - values['LAN']
    # Unit conversion
    values['a'] = values['a'] * AU
    
    # Cleanup unused
    del values['L'], values['LoP']
    
    # Angle conversion
    if angle_units == 'rad':
        for k in ['i', 'LAN', 'AoP', 'f']:
            values[k] = np.deg2rad(values[k])
    
    # {e, a, i, LAN, AoP, f} -- Note: tau/t0 not included
    return values

meeus = {
 'Mercury': {'L': [252.250906, 149472.6746358, -5.35e-06, 2e-09],
             'a': [0.38709831, 0.0, 0.0, 0.0],
             'e': [0.20563175, 2.0406e-05, -2.84e-08, -1.7e-10],
             'i': [7.004986, -0.0059516, 8.1e-07, 4.1e-08],
             'LAN': [48.330893, -0.1254229, -8.833e-05, -1.96e-07],
             'LoP': [77.456119, 0.1588643, -1.343e-05, 3.9e-08]},
 'Venus': {'L': [181.979801, 58517.815676, 1.65e-06, -2e-09],
           'a': [0.72332982, 0.0, 0.0, 0.0],
           'e': [0.00677188, -4.7766e-05, 9.75e-08, 4.4e-10],
           'i': [3.394662, -0.0008568, -3.244e-05, 1e-08],
           'LAN': [76.67992, -0.278008, -0.00014256, -1.98e-07],
           'LoP': [131.563707, 0.0048646, -0.00138232, -5.332e-06]},
 'Earth': {'L': [100.466449, 35999.3728519, -5.68e-06, 0.0],
           'a': [1.000001018, 0.0, 0.0, 0.0],
           'e': [0.01670862, -4.2037e-05, -1.236e-07, 4e-11],
           'i': [0.0, 0.0130546, -9.31e-06, -3.4e-08],
           'LAN': [174.873174, -0.2410908, 4.067e-05, -1.327e-06],
           'LoP': [102.937348, 0.3225557, 0.00015026, 4.78e-07]},
 'Mars': {'L': [355.433275, 19140.2993313, 2.61e-06, -3e-09],
          'a': [1.523679342, 0.0, 0.0, 0.0],
          'e': [0.09340062, 9.0483e-05, -8.06e-08, -3.5e-10],
          'i': [1.849726, -0.0081479, -2.255e-05, -2.7e-08],
          'LAN': [49.558093, -0.2949846, -0.00063993, -2.143e-06],
          'LoP': [336.060234, 0.4438898, -0.00017321, 3e-07]},
 'Jupiter': {'L': [34.351484, 3034.9056746, -8.501e-05, 4e-09],
             'a': [5.202603191, 1.913e-07, 0.0, 0.0],
             'e': [0.04849485, 0.000163244, -4.719e-07, -1.97e-09],
             'i': [1.30327, -0.0019872, 3.318e-05, 9.2e-08],
             'LAN': [100.464441, 0.1766828, 0.00090387, -7.032e-06],
             'LoP': [14.331309, 0.2155525, 0.00072252, -4.59e-06]},
 'Saturn': {'L': [50.077471, 1222.1137943, 0.00021004, -1.9e-08],
            'a': [9.554909596, -2.1389e-06, 0.0, 0.0],
            'e': [0.05550862, -0.000346818, -6.456e-07, 3.38e-09],
            'i': [2.488878, 0.0025515, -4.903e-05, 1.8e-08],
            'LAN': [113.665524, -0.2566649, -0.00018345, 3.57e-07],
            'LoP': [93.056787, 0.5665496, 0.00052809, 4.882e-06]},
 'Uranus': {'L': [314.055005, 429.8640561, 0.00030434, 2.6e-08],
            'a': [19.218446062, -3.72e-08, 9.8e-10, 0.0],
            'e': [0.0462959, -2.7337e-05, 7.9e-08, 2.5e-10],
            'i': [0.773196, 0.0007744, 3.749e-05, -9.2e-08],
            'LAN': [74.005947, 0.5211258, 0.00133982, 1.8516e-05],
            'LoP': [173.005159, 1.4863784, 0.002145, 4.33e-07]},
 'Neptune': {'L': [304.348665, 219.8833092, 0.00030926, 1.8e-08],
             'a': [30.110386869, -1.663e-07, 6.9e-10, 0.0],
             'e': [0.00898809, 6.408e-06, -8e-10, -5e-11],
             'i': [1.769952, -0.0093082, -7.08e-06, 2.8e-08],
             'LAN': [131.784057, 1.1022057, 0.00026006, -6.36e-07],
             'LoP': [48.123691, 1.4262677, 0.00037918, -3e-09]},
 'Pluto': {'L': [238.92903833, 145.20780515, 0.0, 0.0],
           'a': [39.48211675, -0.00031596, 0.0, 0.0],
           'e': [0.2488273, 5.17e-05, 0.0, 0.0],
           'i': [17.14001206, 4.818e-05, 0.0, 0.0],
           'LAN': [110.30393684, -0.01183482, 0.0, 0.0],
           'LoP': [224.06891629, -0.04062942, 0.0, 0.0]}
 }

if __name__ == "__main__":
    for days in [0, 1, 2]:
        ele = get_elements('earth', days/J_cent, angle_units='deg')
        if days == 0: print(ele)
        print(f"At day {days}:\t{round(ele['f'], 3)}")
    start = get_elements('earth', 0, angle_units='deg')
    end = get_elements('earth', 365.25/J_cent, angle_units='deg')
    ans = end['f'] - start['f']
    print(f"Movement over one year:\t{round(ans, 3)}")
    