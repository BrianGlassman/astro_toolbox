# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 21:57:33 2022

@author: Brian
"""

import numpy as np
from collections import namedtuple, UserDict
import datetime
from matplotlib.patches import ConnectionPatch

deg_sign = u'\N{DEGREE SIGN}'

J_cent = 36525 # Number of days in a Julian century
J2000 = 2451545.0
# AU = 1.49597870700e8 # km / AU - per PlanetaryGravitationalCoefficientsRadii-1 on 6008 Canvas
AU = 1.49597870691e8 # km / AU - per HW9 instructions

class CaseInsensitiveDict(UserDict):
    '''A dictionary that is case-insensitive for keys'''
    def __init__(self, baseDict):
        assert isinstance(baseDict, dict)
        assert all(isinstance(key, str) for key in baseDict.keys())
        assert len(set(k.lower() for k in baseDict.keys())) == len(baseDict.keys()), \
            "Duplicate case-insensitive keys"
        super().__init__(baseDict)
        
    def __getitem__(self, key):
        assert isinstance(key, str)
        for k in super().keys():
            if key.lower() == k.lower():
                return super().__getitem__(k)
        raise KeyError(key)
        
mu = {# km^3/s^2
      'Sun'     : 1.32712440018e11,
      'Venus'   : 3.24858599e5,
      'Earth'   : 3.98600433e5,
      'Mars'    : 4.28283100e4,
      'Jupiter' : 1.266865361e8,
      'Saturn'  : 3.7931208e7,
      'Uranus'  : 5.7939513e6,
      'Neptune' : 6.835100e6,
      'Pluto'   : 8.71e2,
      }
mu = CaseInsensitiveDict(mu)

planetary_radius = {# km
                    'Venus'   : 6051.8,
                    'Earth'   : 6378.14,
                    'Mars'    : 3396.19,
                    'Jupiter' : 71492,
                    'Saturn'  : 60268,
                    'Uranus'  : 25559,
                    'Neptune' : 24764,
                    'Pluto'   : 1188.3,
                    }
planetary_radius = CaseInsensitiveDict(planetary_radius)

# Generated using Orbit.from_meeus April 19, 2022
planetary_period = {# seconds
                    'venus': 19414075.93340966,
                    'earth': 31558244.20761866,
                    'mars': 59354347.98146568,
                    'jupiter': 374492590.3146652,
                    'saturn': 932077469.8750662,
                    'uranus': 2658823893.3606086,
                    'neptune': 5214187832.621648,
                    'pluto': 7829096563.614984,
                    }
planetary_period = CaseInsensitiveDict(planetary_period)

# Minimum distance when designing flybys
# NOTE: defined as altitude, but then immediately converted to total radius
planetary_exclusion = {# km
    'Earth': 300, # HW 6 instructions
    'Mars': 200 # Lab 5 instructions, although it's target PE not necessarily exclusion zone
    }
planetary_exclusion = CaseInsensitiveDict(planetary_exclusion)
planetary_exclusion = {p: v + planetary_radius[p] for p,v in planetary_exclusion.items()}

planetary_approx_orbit = {# AU - from HW9 instructions
    'Mercury': 0.387,
    'Venus': 0.723,
    'Earth': 1,
    'Mars': 1.524,
    'Jupiter': 5.203,
    'Saturn': 9.537,
    'Uranus': 19.191,
    'Neptune': 30.069,
    'Pluto': 39.482,
    }
planetary_approx_orbit = CaseInsensitiveDict(planetary_approx_orbit)

# Default planet colors
colors = { # https://matplotlib.org/stable/gallery/color/named_colors.html
    'Sun': 'yellow',
    'Mercury': 'grey',
    'Venus': 'tan',
    'Earth': 'royalblue',
    'Luna': 'lightsteelblue',
    'Mars': 'orangered',
    'Jupiter': 'sandybrown',
    'Saturn': 'gold',
    'Uranus': 'cyan',
    'Neptune': 'darkblue',
    'Pluto': 'grey',
    }
colors = CaseInsensitiveDict(colors)

deg = np.degrees
rad = np.radians

def days_to_secs(days):
    return days*60*60*24
def secs_to_days(secs):
    return secs/60/60/24

def unit(vector):
    """
    Normalize a vector (i.e. give the unit vector). Returns zero if the vector
    has zero length

    Parameters
    ----------
    vector : 3-vector of numerics
        The input vector

    Returns
    -------
    unit_vector : 3-vector of numerics, magnitude=1
        The unit vector along the input vector's direction
    """
    norm = np.linalg.norm(vector)
    if norm == 0:
        return np.array([0, 0, 0])
    else:
        return vector / norm
norm = np.linalg.norm # Sugar syntax for getting a vector's length

# Note: tau is included even though it's not a true orbital element
Orbital_Elements = namedtuple("Orbital_Elements", ['e','a','i','LAN','AoP','f','tau'])

def fig_3d(*args, **kwargs):
    from matplotlib import pyplot as plt
    # 3D Plotting methodology from:
    # https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html
    fig = plt.figure(*args, **kwargs)
    ax = plt.axes(projection='3d')
    return fig, ax

def square_3d(ax):
    '''Make a 2D or 3D plot have "square" aspect ratio'''
    threeD = ax.name == "3d" # https://stackoverflow.com/a/43563348/14501840
    if threeD:
        # Adapted from https://stackoverflow.com/a/64453375/14501840
        x = np.ptp(ax.get_xlim3d())
        y = np.ptp(ax.get_ylim3d())
        z = np.ptp(ax.get_zlim3d())
        v = max(x, y, z)
        ax.set_box_aspect([v, v, v])
        
        # Remove Z axis ticks if it's too small
        if z / min(x, y) < 1e-3:
            ax.set_zticks([])
    else:
        ax.set_aspect('equal')
        
def plot_point(threeDvector, ax, threeD, *args, **kwargs):
    if threeD:
        ax.plot(threeDvector[0], threeDvector[1], threeDvector[2], '.', *args, **kwargs)
    else:
        ax.plot(threeDvector[0], threeDvector[1], '.', *args, **kwargs)

def plot_line(vector: np.ndarray, ax, threeD: bool, **kwargs):
    """Draw a line
    vector - numpy array with positions in rows.
        ex.
           3D            2D
        [[1 2 3]   or  [[1 2]  <-- point 1
         [4 5 6]]       [3 4]  <-- point 2
    """
    if threeD:
        assert vector.shape == (2, 3)
        ax.plot(vector[:,0], vector[:,1], vector[:,2], **kwargs)
    else:
        assert vector.shape == (2, 2)
        ax.plot(vector[:,0], vector[:,1], **kwargs)

def plot_vector(start, threeDvector, ax, threeD, *args, **kwargs):
    if isinstance(start, (int, float)) and start == 0:
        start = [0, 0, 0]
    if isinstance(threeDvector, (int, float)) and threeDvector == 0:
        threeDvector = [0, 0, 0]
        
    threeDvector = start + threeDvector # Don't use +=, it overwrites the source in calling scope
        
    # FIXME handle args/kwargs
    if args:
        raise NotImplementedError("Can't handle args, only kwargs")
        
    if threeD:
        raise NotImplementedError()
        # TODO
        # use this: https://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-matplotlibs-3d-plot
        # to draw 3D vectors with arrows
    else:
        arrow = ConnectionPatch(start[0:2], threeDvector[0:2], "data", "data",
                                arrowstyle='-|>', *args, **kwargs)
    ax.add_patch(arrow)

def plot_sphere(radius, ax, **plot_kwargs):
    # Adapted from https://github.com/Elucidation/OrbitalElements/blob/master/graphics.py

    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones_like(u), np.cos(v))

    # Plot:
    ax.plot_wireframe(x, y, z, rstride=4, cstride=4, **plot_kwargs)

def normalize_angle(angle):
    """Normalizes an angle to the range [0, 2*pi]"""
    while angle < 0:
        angle += 2*np.pi
    while angle > 2*np.pi:
        angle -= 2*np.pi
    return angle

def _so_calendar_to_jde(date):
    # https://stackoverflow.com/a/39019116/14501840 adapted to include partial days to match jde_to_calendar
    ordinal = date.toordinal()
    jde = ordinal + 1721424.5
    jde += (date.hour + (date.minute + (date.second/60))/60)/24
    return jde
def calendar_to_jde(date):
    # From L2 slides
    assert isinstance(date, datetime.datetime)
    
    # Get date components
    y, m, d = date.year, date.month, date.day
    assert 1901 <= y <= 2099
    assert 1 <= m <= 12
    assert 1 <= d <= 31
    
    # Calculate integer portion of JDE
    int_jde = 367*y - int(7*(y+int((m+9)/12)) / 4) + int(275*m/9) + d + 1721013.5
    
    # Get time components
    h, m, s, ms = date.hour, date.minute, date.second, date.microsecond
    # Convert everything to hours
    s = s + ms*1e-6
    m = m + s/60
    h = h + m/60
    assert 0 <= h <= 24
    
    # Calculate decimal portion of JDE
    dec_jde = h / 24
    
    # Combine
    jde = int_jde + dec_jde
    
    # Compare to Stack Overflow version
    assert jde == _so_calendar_to_jde(date)
    
    return jde
juliandate = calendar_to_jde

def jde_to_calendar(jde):
    # Verified exactly correct using the Calendar Dates and JDE values from HW6 instructions
    remainder = (jde + 0.5) % 1
    ordinal = int(jde - 1721424.5)
    dt = datetime.datetime.fromordinal(ordinal)
    dt = dt + datetime.timedelta(days=remainder)
    return dt
# TODO figure out how to reverse the lecture slide calculations

def jde_to_gmat(jde):
    # from http://gmat.sourceforge.net/docs/nightly/html/SpacecraftEpoch.html
    return jde - 2430000.0

def gmat_to_jde(gmat):
    # from http://gmat.sourceforge.net/docs/nightly/html/SpacecraftEpoch.html
    return gmat + 2430000.0
