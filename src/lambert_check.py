# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 23:44:45 2022

@author: Brian
"""

import numpy as np

import util
from toolbox import Orbit

def do_test(inputs, checks):
    tof_days = inputs.arrive - inputs.launch
    tof_sec = tof_days * 24 * 60 * 60
    o = Orbit.from_lambert(inputs.p_launch, inputs.p_arrive, tof=tof_sec, mu=util.mu['Sun'])
    
    def test(orbit):
        p_launch, v_launch = o.position_velocity()
        assert all(p_launch - checks.p_launch < 1e-6), (p_launch, checks.p_launch)
        assert all(v_launch - checks.v_launch < 1e-6), (v_launch, checks.v_launch)
    test(o)
    
    # Launch planet orbit
    fig, ax = Orbit.from_pvt(checks.p_launch, checks.v_planet_launch, mu=util.mu['Sun']).plot()
    # Destination planet orbit
    Orbit.from_pvt(checks.p_arrive, checks.v_planet_arrive, mu=util.mu['Sun']).plot(fig=fig)
    # Test Case orbit
    checks.orbit.plot(0, tof_sec, mode='time', fig=fig, plot_kwargs={'ls':'--'})
    o.plot(0, tof_sec, mode='time', fig=fig)
    
    while fig.legends:
        fig.legends[0].remove()
    fig.legend(['Central Body', 'Launch Planet', 'Arrival Planet', 'Correct Trajectory', 'Calculated Trajectory'])
    util.square_3d(ax)
    
#%% Earth to Venus - single-rev
def earth_venus_single():
    class inputs():
        launch = 2455450
        arrive = 2455610
        p_launch = np.array([147084764.9, -32521190, 467.1900914])
        p_arrive = np.array([-88002509.16, -62680223, 4220331.525])
    
    class checks():
        p_launch = inputs.p_launch
        p_arrive = inputs.p_arrive
        v_planet_launch = np.array([5.94623924, 28.9746412, -0.0007159])
        v_planet_arrive = np.array([20.0705936, -28.68983, -1.5512918])
        
        v_launch = np.array([4.651443497, 26.0824144, -1.39306043])
        v_arrive = np.array([16.79262045, -33.351675, 1.523021504])
        
        orbit = Orbit.from_pvt(p_launch, v_launch, mu=util.mu['Sun'], t0=0)
        
    do_test(inputs, checks)
earth_venus_single()
        
#%% Mars to Jupiter - single-rev
def mars_jupiter_single():
    class mj_inputs():
        launch = 2456300
        arrive = 2457500
        p_launch = np.array([170145121.3, -117637192.8, -6642044.272])
        p_arrive = np.array([-803451694.7, 121525767.1, 17465211.78])
        
    class mj_checks():
        p_launch = mj_inputs.p_launch
        p_arrive = mj_inputs.p_arrive
        v_planet_launch = np.array([14.70149986, 22.00292904, 0.100109562])
        v_planet_arrive = np.array([-2.110465959, -12.31199244, 0.098198408])
        
        v_launch = np.array([13.7407773, 28.83099312, 0.691285008])
        v_arrive = np.array([-0.883933069, -7.983627014, -0.240770598])
        
        orbit = Orbit.from_pvt(p_launch, v_launch, mu=util.mu['Sun'], t0=0)

    do_test(mj_inputs, mj_checks)
mars_jupiter_single()
    
#%% Earth to Venus - multi-rev type 3
def earth_venus_multi():
    class inputs():
        launch = 2460545
        arrive = 2460919
        p_launch = np.array([130423562.1, -76679031.85, 3624.816561])
        p_arrive = np.array([19195371.67, 106029328.4, 348953.802])
        
    class checks():
        p_launch = inputs.p_launch
        p_arrive = inputs.p_arrive
        v_planet_launch = np.array([14.61294123, 25.56747613, -0.001503446])
        v_planet_arrive = np.array([-34.57913611, 6.064190776, 2.078550651])
        
        v_launch = np.array([12.76771134, 22.79158874, 0.090338826])
        v_arrive = np.array([-37.30072389, -0.176853447, -0.066693083])
        
        orbit = Orbit.from_pvt(p_launch, v_launch, mu=util.mu['Sun'], t0=0)
    
    do_test(inputs, checks)
earth_venus_multi()
