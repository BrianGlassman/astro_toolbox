# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 15:51:38 2022

@author: Brian
"""

import numpy as np
from matplotlib import pyplot as plt
from datetime import date

from astro_toolbox import util, Orbit

mu_sun = util.mu['sun']

def _calc(depart_date, arrive_date, depart_planet, arrive_planet):
    if isinstance(depart_date, date):
        depart_date = util.juliandate(depart_date)
    if isinstance(arrive_date, date):
        arrive_date = util.juliandate(arrive_date)
    
    depart_orbit = Orbit.from_meeus(depart_planet, depart_date)
    r_planet_depart, v_planet_depart = depart_orbit.position_velocity()
    arrive_orbit = Orbit.from_meeus(arrive_planet, arrive_date)
    r_planet_arrive, v_planet_arrive = arrive_orbit.position_velocity()
    tof_sec = (arrive_date - depart_date) * 24 * 60 * 60
    transfer = Orbit.from_Lambert(r_planet_depart, r_planet_arrive, tof_sec, mu=mu_sun)
    
    r_planet_depart = depart_orbit.position_velocity()[0]
    r_depart = transfer.position_velocity()[0]
    assert util.norm(r_depart - r_planet_depart) < 1e-3*util.norm(r_planet_depart)
    
    r_planet_arrive = arrive_orbit.position_velocity()[0]
    r_arrive = transfer.position_velocity(tof_sec)[0]
    assert util.norm(r_arrive - r_planet_arrive) < 1e-3*util.norm(r_planet_arrive)
    
    return depart_orbit, arrive_orbit, transfer, tof_sec

def compute_rp(v_in, v_out, mu):
    assert v_in.shape == (3,)
    assert v_out.shape == (3,)
    try:
        assert abs(util.norm(v_in) - util.norm(v_out)) < 0.1, abs(util.norm(v_in) - util.norm(v_out)) # FIXME should the difference be this big?
    except AssertionError:
        # import pdb; pdb.set_trace()
        raise
    v_inf_sq = util.norm(v_in) * util.norm(v_out)
    turn_angle = np.arccos(np.dot(v_in, v_out) / v_inf_sq)
    r_p = mu / v_inf_sq * (1/np.cos((np.pi-turn_angle)/2) - 1)
    return r_p

def calc_pcp(depart_min, depart_max, depart_planet, depart_step,
             arrive_min, arrive_max, arrive_planet, arrive_step,
             max_tof=float('inf')):
    '''Does the computation for a porkchop plot, but doesn't plot anything
    
    Dates should be JDE'''
    # NOTE: can end up slightly overshooting the endpoint if (max-min) is not an even multiple of step size
    departure_dates = np.arange(depart_min, depart_max+depart_step, depart_step) # Use arange for integer increments
    arrival_dates = np.arange(arrive_min, arrive_max+arrive_step, arrive_step) # Use arange for integer increments
    
    # Calculate the v_inf values across the date ranges
    arrive_v_inf = np.zeros((len(arrival_dates), len(departure_dates)))
    arrive_v_inf_vec = {} # Use a dict instead of a very high-D array
    depart_v_inf = np.zeros((len(arrival_dates), len(departure_dates)))
    depart_v_inf_vec = {} # Use a dict instead of a very high-D array
    c3 = np.zeros((len(arrival_dates), len(departure_dates)))
    tof_days = np.zeros((len(arrival_dates), len(departure_dates)))
    for y, arrival_date in enumerate(arrival_dates):
        if y%10 == 0: print(f"{y}: ", end='')
        arrive_v_inf_vec[y] = {}
        depart_v_inf_vec[y] = {}
        for x, departure_date in enumerate(departure_dates):
            if y%10 == 0 and x%10 == 0: print(x, end=' ')
            if arrival_date - departure_date > max_tof:
                depart_v_inf_vec[y][x] = None
                depart_v_inf[y][x] = None
                arrive_v_inf_vec[y][x] = None
                arrive_v_inf[y][x] = None
                tof_days[y][x] = None
                continue
            
            try:
                depart_orbit, arrive_orbit, transfer, tof_sec = _calc(departure_date, arrival_date, depart_planet, arrive_planet)
            except:
                depart_v_inf_vec[y][x] = None
                depart_v_inf[y][x] = None
                arrive_v_inf_vec[y][x] = None
                arrive_v_inf[y][x] = None
                tof_days[y][x] = None
                continue
            
            # Departure
            r_planet_depart, v_planet_depart = depart_orbit.position_velocity()
            r_depart, v_depart = transfer.position_velocity()
            assert util.norm(r_depart - r_planet_depart) < 1e-3*util.norm(r_planet_depart)
            v_inf = v_depart - v_planet_depart
            depart_v_inf_vec[y][x] = v_inf
            v_inf = util.norm(v_inf)
            depart_v_inf[y][x] = v_inf
            assert depart_v_inf.dtype == np.float64
            c3[y][x] = v_inf**2
            
            # Arrival
            r_planet_arrive, v_planet_arrive = arrive_orbit.position_velocity()
            r_arrive, v_arrive = transfer.position_velocity(tof_sec)
            msg = util.norm(r_arrive - r_planet_arrive) / util.norm(r_planet_arrive)
            assert util.norm(r_arrive - r_planet_arrive) < 1e-3*util.norm(r_planet_arrive), msg
            v_inf = v_arrive - v_planet_arrive
            arrive_v_inf_vec[y][x] = v_inf
            v_inf = util.norm(v_inf)
            arrive_v_inf[y][x] = v_inf
            
            tof_days[y][x] = arrival_date - departure_date
            
        if y%10 == 0: print()
    
    # For plot: Departure dates on X, arrival dates on Y
    x = [departure_date - depart_min for departure_date in departure_dates]
    y = [arrival_date - arrive_min for arrival_date in arrival_dates]
    
    # Check for failure
    if np.isnan(tof_days).all():
        raise ValueError("No valid transfers found")
    
    return x, y, depart_v_inf, depart_v_inf_vec, c3, arrive_v_inf, arrive_v_inf_vec, tof_days

def plot_pcp(x, y,
             depart_date, depart_vals, depart_levels,
             arrive_date, arrive_vals, arrive_levels,
             tof_days, tof_levels,
             fig=None, ax=None):
        lines = []
        if fig is None and ax is None:
            fig, ax = plt.subplots(figsize=(11,8))
        elif fig is not None and ax is not None:
            pass # Use as given
        else: raise ValueError()
        ax.set_xlabel(f'Departure: Days past {depart_date}')
        ax.set_ylabel(f'Arrival: Days past {arrive_date}')
        
        if depart_vals is not None:
            CS = ax.contour(x, y, depart_vals, levels=depart_levels, colors=['red'], linewidths=.5)
            lines.append(CS.legend_elements()[0][0])
            ax.clabel(CS, CS.levels, inline=True, fontsize=10, inline_spacing=1)
        
        if arrive_vals is not None:
            CS = ax.contour(x, y, arrive_vals, levels=arrive_levels, colors=['blue'], linewidths=.5)
            lines.append(CS.legend_elements()[0][0])
            ax.clabel(CS, CS.levels, inline=True, fontsize=10, inline_spacing=1)
        
        CS = ax.contour(x, y, tof_days, levels=tof_levels, colors=['black'], linewidths=1)
        lines.append(CS.legend_elements()[0][0])
        ax.clabel(CS, CS.levels, inline=True, fontsize=10, inline_spacing=1)
        ax.grid()
        
        return fig, ax, lines
    
