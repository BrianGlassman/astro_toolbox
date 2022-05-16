# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 13:40:58 2021

@author: Brian
"""

import numpy as np
from scipy.optimize import root
from matplotlib import pyplot as plt
# 3D Plotting methodology from:
# https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html
from mpl_toolkits import mplot3d

from . import util

# TODO
# - Refactor to use the Universal Variable approach in Prussing/Conway
#   This should mean that the same approaches can be used regardless of elliptical vs. hyperbolic

class Orbit():
    @staticmethod
    def __calc_tau(e, a, f, t0, mu, do_checks=True):
        # Calculate time of periapsis passage
        # NOTE: e < 0 is caught above
        # NOTE: circular and/or un-inclined orbits are handled above so that f is
        #   always defined, so no special handling needed
        if e < 0:
            raise ValueError(f"Impossible e: {e}")
        elif 0 <= e < 1:
            # Elliptical or Circular
            # Try to find the first periapsis passage with tau > 0
            # NOTE: may be slightly negative due to floating point errors
            E = 2 * np.arctan2(np.sqrt((1-e)/(1+e)) * np.tan(f/2), 1) # Eccentric anomaly
            if E < 0:
                E = E + 2*np.pi
            tau = t0 - np.sqrt(a**3 / mu) * (E - e * np.sin(E))
            
            T = 2*np.pi*np.sqrt(a**3 / mu)
            while tau > T:
                tau -= T
                
            if tau < 0:
                if tau < -T*.01:
                    tau += T
                else:
                    pass # Approximately zero
                
            if do_checks: assert -.01*T <= tau < T
        elif e >= 1:
            # Hyperbolic or Parabolic
            # Tau is fixed. Tau > t0 if inbound, tau < t0 if outbound
            # Source: http://control.asu.edu/Classes/MAE462/462Lecture05.pdf
            #E = np.arccosh((np.cos(f) + e) / (1 + e*np.cos(f))) # Eccentric anomaly
            H = 2*np.arctanh(np.sqrt((e-1)/(e+1))*np.tan(f/2)) # Hyperbolic anomaly <-- tanh(H/2) = sqrt((e-1)/(e+1))*tan(f/2)
            if abs(H) > np.pi:
                raise RuntimeError(f"probably bad? H = {H}")
            
            M = e*np.sinh(H) - H # Mean Anomaly = sqrt(mu / (-a)**3) * (t0 - tau)
            tau = t0 - M*np.sqrt((-a)**3 / mu)
        else:
            raise RuntimeError("Fall-through shouldn't be possible")
        return tau
    
    @staticmethod
    def __normalize_f(f, e):
        if e < 1:
            # Force positive true anomaly for elliptical/circular orbits
            while f < 0:
                f += 2*np.pi
            while f > 2*np.pi:
                f -= 2*np.pi
            assert 0 <= f <= 2*np.pi, f"f = {f}"
        else:
            assert -np.pi < f <= np.pi, f
        return f
    
    def __init__(self, e=0, a=util.AU, i=0, LAN=0, AoP=0, f=0, t=0, mu='Sun', do_checks=True, time_mode='t0'):
        if isinstance(mu, str):
            mu = util.mu[mu]
        if i > np.pi and i <= 180:
            # Assume given in degrees
            i = util.rad(i)
        if do_checks: assert e >= 0, f"e = {e}"
        if do_checks:
            if not (0 <= i <= np.pi):
                # print(f"i = {i} out of range, normalizing")
                i = util.normalize_angle(i)
        if do_checks:
            if not (0 <= AoP <= 2*np.pi):
                # print(f"AoP = {AoP} out of range, normalizing")
                AoP = util.normalize_angle(AoP)
        if do_checks:
            f = self.__normalize_f(f, e)
        # assert 0 <= tau <= 2*np.pi # not sure I should force-check tau, and this is wrong anyway
        if do_checks: assert mu > 0
        
        time_mode = time_mode.lower()
        if time_mode == 'tau':
            # Set tau directly
            tau = t
        elif time_mode == 't0':
            # Set current time rather than setting periapsis time
            tau = self.__calc_tau(e, a, f, t, mu)
        
        self.e = e
        self.a = a
        self.i = i
        self.LAN = LAN
        self.AoP = AoP
        self.f = f
        self.tau = tau
        self.mu = mu
        
    @classmethod
    def from_orbit(cls, source, e=None, a=None, i=None, LAN=None, AoP=None, f=None, t=None,
                   do_checks=True, time_mode='tau'):
        """
        Duplicate an existing orbit, with overrides
        """
        if e is None: e = source.e
        if a is None: a = source.a
        if i is None: i = source.i
        if LAN is None: LAN = source.LAN
        if AoP is None: AoP = source.AoP
        if t is None:
            assert time_mode == 'tau'
            t = source.tau
        return cls(e=e, a=a, i=i, LAN=LAN, AoP=AoP, f=f, t=t,
                   mu=source.mu, do_checks=do_checks, time_mode=time_mode)
        
    @classmethod
    def from_pvt(cls, position, velocity, mu, t0=0, overrides=None, do_checks=True):
        """
        Compute the orbital elements for a 2 Body Problem orbit from position,
        velocity, and time.
        
        Equations from class notes: orbit_elements.pdf
        Order of operations from 9/9/2021 lecture
    
        Parameters
        ----------
        position : 3-vector of numerics
            The position vector in 3 space at the given time (units: km)
        velocity : 3-vector of numerics
            The velocity vector in 3 space at the given time (units: km)
        time : numeric
            The initial time at which position and velocity are given (units: seconds)
        mu : numeric
            The value to use for mu (units: km^3/s^2)
        t0 : numeric, optional
            The current time (defaults to 0)
        overrides : dict
            Orbital elements to override. Ex. {"e": 0} forces eccentricity to zero
        do_checks : boolean, optional
            If True (default) check for bad input/output
            If False, allow anything
    
        Returns
        -------
        an Orbit object that includes the given parameters
        """
        if isinstance(mu, str):
            mu = util.mu[mu]
        if overrides is None:
            overrides = {}
        
        # Calculate h to be used later
        h_vector = np.cross(position, velocity)
        h_hat = util.unit(h_vector)
        h = util.norm(h_vector)
        # Always uniquely and well defined
        
        # Can't use "p" because I'm using that for position magnitude
        orbit_param = h**2 / mu
        
        # Calculate i to return
        if 'i' in overrides:
            i = overrides['i']
        else:
            i = np.arccos(util.unit(h_vector)[2]) # h_hat dot z_hat = h_hat_z
        # Always uniquely and well defined
        if do_checks: assert 0 <= i <= np.pi, f"i = {i}"
        
        # Calculate LAN via the line of nodes, to return and used later
        if 'LAN' in overrides:
            raise NotImplementedError()
        # LAN = capital omega
        if i == 0 or i == np.pi:
            # Use x axis as fictitious line of nodes
            nodes = np.array([1, 0, 0])
        else:
            nodes = util.unit(np.cross([0,0,1], h_hat)) # Unit vector along line of nodes
        # Note use of arctan2 to give four-quadrant arctan
        LAN = np.arctan2(nodes[1], # y dot n = n_y
                         nodes[0]) # x dot n = n_x
        
        # Calculate e to return and use later
        if 'e' in overrides:
            e = overrides['e']
            if e != 0:
                raise NotImplementedError()
                # eccentricity vector is used to calculate AoP if e != 0
                # eccentricity unit vector is used to calculate f if e != 0
        else:
            eccentricity = 1/mu * np.cross(velocity, h_vector) - util.unit(position)
            e_hat = util.unit(eccentricity)
            e = util.norm(eccentricity)
            # Always uniquely and well defined
        if do_checks: assert e >= 0, f"Negative eccentricity: e = {e}"
            
        # Calculate the AoP via the in-plane line perpendicular to the line of nodes, to return
        # AoP = small omega
        n_perp = np.cross(h_hat, nodes)
        
        # Note use of arctan2 to give four-quadrant arctan
        # if e == 0 or n == 0:
        #     raise NotImplementedError("AoP when e == 0 or line of nodes is undefined?")
        # elif i == 0 or i == np.pi:
        #     raise NotImplementedError("AoP when LAN undefined")
        if e == 0:
            AoP = 0
        else:
            AoP = np.arctan2(np.dot(eccentricity, n_perp),
                             np.dot(eccentricity, nodes))
            # Undefined when eccentricity or line of nodes vectors are zero or undefined
        if AoP < 0:
            AoP += 2*np.pi
        if do_checks: assert 0 <= AoP <= 2*np.pi, f"AoP = {AoP}"
            
        if 'a' in overrides:
            a = overrides['a']
        else:
            if e == 1:
                a = float('inf')
            else:
                a = orbit_param / (1 - e**2)
                # Undefineed when e == 1
        
        # if e == 0 or i == 0 or i == np.pi:
        #     raise NotImplementedError("f is undefined")
        # else:
        if e == 0:
            # Circular. No periapsis, use line of nodes instead
            # NOTE: if un-inclined orbit, X axis is already used as line of nodes
            #   so no special handling required
            f = np.arctan2(np.dot(position, np.cross(h_hat, util.unit(nodes))),
                           np.dot(position, util.unit(nodes)))
        else:
            # Non-circular. Use periapsis (from eccentricity vector) as usual
            f = np.arctan2(np.dot(position, np.cross(h_hat, e_hat)),
                           np.dot(position, e_hat))
        f = cls.__normalize_f(f, e)
    
        # Time of periapsis passage
        if 'tau' in overrides:
            tau = overrides['tau']
        else:
            tau = cls.__calc_tau(e=e, a=a, f=f, t0=t0, mu=mu)
        
        return cls(e=e, a=a, i=i, LAN=LAN, AoP=AoP, f=f, t=tau,
                   time_mode='tau', mu=mu, do_checks=do_checks)

    @classmethod
    def from_Lambert(cls, position_1, position_2, tof, mu,
                     orbit_type=None, t0=0, overrides=None, do_checks=True, verbose=False):
        """Solves Lambert's problem to generate an orbit from two positions and
        the time of flight
        
        Methodology according to Lambert Handout from 6008, Kate Davis
    
        Parameters
        ----------
        position_1 : 3-vector of numerics
            The initial position vector in 3 space (units: km)
        position_2 :  3-vector of numerics
            The final position vector in 3 space (units: km)
        tof : numeric
            The time of flight between the two positions (units: seconds)
        mu : numeric
            The value to use for mu (units: km^3/s^2)
        orbit_type : None, 1, or 2, optional
            Allows forcing Type I or Type II orbits. Leaving as None will select
            automatically to ensure a prograde solution.
        t0: numeric, optional
            Offset the initial position by this amount, defaults to 0 (units: seconds)
        overrides : dict
            Orbital elements to override. Ex. {"e": 0} forces eccentricity to zero
        do_checks : boolean, optional
            If True (default) check for bad input/output
            If False, allow anything
    
        Returns
        -------
        an Orbit object that includes the given parameters
        """
        if isinstance(mu, str):
            mu = util.mu[mu]
        ### Special syntax to auto-calculate the positions
        if isinstance(position_1, tuple) and isinstance(position_2, tuple) and tof is None:
            _planet, _date1 = position_1
            _orbit = Orbit.from_meeus(_planet, _date1)
            position_1 = _orbit.position_velocity()[0]
            _planet, _date2 = position_2
            _orbit = Orbit.from_meeus(_planet, _date2)
            position_2 = _orbit.position_velocity()[0]
            tof = util.days_to_secs(_date2 - _date1)
            
        
        # Call the argument tof for convenience, but rename for ease of coding
        tof_target = tof ; del tof
        
        if orbit_type is None:
            def calc_TA(p_vector):
                # Assume orbits are approximately co-planar in the XY plane, ignore Z
                TA = np.arctan2(p_vector[1], p_vector[0])
                if TA < 0: TA += 2*np.pi
                if TA > 2*np.pi: TA -= 2*np.pi
                return TA
            dTA = calc_TA(position_2) - calc_TA(position_1)
            dTA = dTA % (2*np.pi)
            if dTA < np.pi:
                orbit_type = 1
            else:
                orbit_type = 2
            
        # DM (direction of motion)
        if orbit_type % 2 == 1:
            DM = +1
        else:
            DM = -1
            
        r1 = util.norm(position_1)
        r2 = util.norm(position_2)
            
        cosTA = np.dot(position_1, position_2) / (r1 * r2)
        A = DM*np.sqrt(r1*r2*(1+cosTA))
        
        if A == 0:
            raise ValueError("Trajectory can't be computed")
            
        c2 = 1/2
        c3 = 1/6
        psi = 0
        psi_high = 4*np.pi**2
        psi_low = -4*np.pi
        
        def check_min_tof():
            '''Check if the minimum value of psi is below the target TOF.
            If not, transfer is impossible'''
            psi = psi_low
            if psi > 1e-6:
                c2 = (1 - np.cos(np.sqrt(psi))) / psi
                c3 = (np.sqrt(psi) - np.sin(np.sqrt(psi))) / np.sqrt(psi**3)
            elif psi < -1e-6:
                # NOTE: hyperbolic trig
                c2 = (1 - np.cosh(np.sqrt(-psi))) / psi
                c3 = (np.sinh(np.sqrt(-psi)) - np.sqrt(-psi)) / np.sqrt((-psi)**3)
            else:
                c2 = 1/2
                c3 = 1/6
            y = r1 + r2 + A*(psi*c3 - 1) / np.sqrt(c2)
            if A > 0 and y < 0:
                # Readjust until y > 0
                while y < 0:
                    psi += 0.1
                    y = r1 + r2 + A*(psi*c3 - 1) / np.sqrt(c2)
            
            x = np.sqrt(y / c2)
            tof = (c3*x**3 + A*np.sqrt(y)) / np.sqrt(mu)
            if tof_target < tof:
                raise ValueError(f"Target TOF < TOF for psi_low ({tof:0.3e}). Transfer can't be computed")
            
        check_min_tof()
        
        tof_threshold = 1e-5
        count = 0
        while True and count < 1000:
            y = r1 + r2 + A*(psi*c3 - 1) / np.sqrt(c2)
            if A > 0 and y < 0:
                # Readjust until y > 0
                while y < 0:
                    psi += 0.1
                    y = r1 + r2 + A*(psi*c3 - 1) / np.sqrt(c2)
            
            x = np.sqrt(y / c2)
            tof = (c3*x**3 + A*np.sqrt(y)) / np.sqrt(mu)
            if verbose: print(f"psi ({psi:0.3f}): TOF ({tof:0.3f}) - target ({tof_target:0.3f}) = {tof - tof_target:0.3f}")
            
            if tof < tof_target:
                psi_low = psi
            elif tof == tof_target:
                break # Somehow hit it perfectly
            else:
                psi_high = psi
            
            psi = (psi_high + psi_low) / 2
            if psi > 1e-6:
                c2 = (1 - np.cos(np.sqrt(psi))) / psi
                c3 = (np.sqrt(psi) - np.sin(np.sqrt(psi))) / np.sqrt(psi**3)
            elif psi < -1e-6:
                # NOTE: hyperbolic trig
                c2 = (1 - np.cosh(np.sqrt(-psi))) / psi
                c3 = (np.sinh(np.sqrt(-psi)) - np.sqrt(-psi)) / np.sqrt((-psi)**3)
            else:
                c2 = 1/2
                c3 = 1/6
            
            if abs(tof - tof_target) < tof_threshold:
                break
            
            count += 1
        else:
            raise ValueError("Max count reached")
        
        f = 1 - y / r1
        g_dot = 1 - y / r2
        g = A*np.sqrt(y / mu)
        
        v1 = (position_2 - f*position_1) / g
        v2 = (g_dot*position_2 - position_1) / g
                
        obj = cls.from_pvt(position_1, v1,
                           mu=mu, t0=t0, overrides=overrides, do_checks=do_checks)
        obj.tof = tof
        return obj
    from_lambert = from_Lambert # Alias
    
    @classmethod
    def from_meeus(cls, planet, JDE):
        """Generates an orbit from Meeus coefficients
        
        Parameters
        ----------
        planet : string
            Which planet to generate an orbit for
        JDE : numeric
            The JDE time to evaluate at
    
        Returns
        -------
        an Orbit object that includes the given parameters
        """
        from . import meeus
        T = meeus.JDE_to_T(JDE)
        ele = meeus.get_elements(planet, T, angle_units='rad')
        ele['t'] = JDE
        obj = cls(**ele, mu='Sun', time_mode='t0')
        obj.planet = planet
        obj.color = util.colors[planet]
        return obj
    
    @classmethod
    def from_hohmann(cls, r_p, r_a, mu='Sun'):
        if isinstance(r_p, str):
            r_p = util.planetary_approx_orbit[r_p]
            r_p = r_p * util.AU
        if isinstance(r_a, str):
            r_a = util.planetary_approx_orbit[r_a]
            r_a = r_a * util.AU
        assert r_a >= r_p, "r_a must be larger than r_p"
        a = (r_p + r_a) / 2
        e = (r_a - r_p) / (r_a + r_p)
        return cls(e=e, a=a, mu=mu)

    @property
    def ele(self):
        '''Return the orbital elements as a namedtuple'''
        return util.Orbital_Elements(self.e, self.a, self.i, self.LAN, self.AoP, self.f, self.tau)
    
    @property
    def ele_gmat(self):
        '''Return the orbital elements as a namedtuple using GMAT's units'''
        return util.Orbital_Elements(self.e, self.a,
                                     util.deg(self.i), util.deg(self.LAN), util.deg(self.AoP), util.deg(self.f),
                                     self.tau)
    
    @property
    def T(self):
        '''Orbital Period (seconds)'''
        return 2*np.pi*np.sqrt(self.a**3 / self.mu)
    
    @property
    def r_p(self):
        '''Radius of periapsis'''
        # TODO add an argument vector=False which calculates the position of periapsis, not just radius
        # r_p = (h**2/mu) / (1+e) ... h**2/mu = a*(1-e**2)
        num = self.a * (1 - self.e**2)
        den = 1 + self.e
        return num / den
    
    @property
    def r_a(self):
        '''Radius of apoapsis'''
        # TODO add an argument vector=False which calculates the position of apoapsis, not just radius
        assert self.a > 0, "Apoapsis is undefined for parabolic/hyperbolic orbits"
        # FIXME I think r_a actually is defined for hyperbolic, just weirdly
        # r_a = (h**2/mu) / (1-e) ... h**2/mu = a*(1-e**2)
        num = self.a * (1 - self.e**2)
        den = 1 - self.e
        return num / den
    
    @property
    def E(self):
        '''Eccenctric Anomaly (radians)'''
        cosE = (self.e + np.cos(self.f)) / (1 + self.e*np.cos(self.f))
        sinE = (np.sqrt(1-self.e**2)*np.sin(self.f)) / (1 + self.e*np.cos(self.f))
        E = np.arctan2(sinE, cosE)
        return E
    
    def M(self, time):
        if self.e < 1:
            # Elliptical or circular
            M = np.sqrt(self.mu/self.a**3) * (time - self.tau)
        else:
            # Parabolic or hyperbolic
            M = np.sqrt(self.mu/(-self.a)**3) * (time - self.tau)
        return M
    
    def __str__(self):
        return str(self.ele).replace(util.Orbital_Elements.__name__, "Orbit")
    
    def __repr__(self):
        return str(self)
    
    def predict(self, time, angle_units='rad'):
        """
        Analytically solve Kepler's equation for the given orbit at the given time
        
        Equations from class notes: orbit_elements.pdf
        Order of operations from 9/9/2021 lecture
    
        Parameters
        ----------
        time : numeric
            Time to solve at (in seconds)
        angle_units : string, optional, either "rad" or "deg"
            What units to use for the returned angle values
            "rad" for radians (default)
            "deg" for degrees
    
        Raises
        ------
        ValueError
            Unknown angle units given
    
        Returns
        -------
        f : numeric
            True anomaly for the given orbit at the given time (in the chosen units)
        E/F : numeric
            Eccentric anomaly for the given orbit at the given time (in the chosen units)
            Note: is actuall hyperbolic eccentric anomaly (F) if orbit is hyperbolic
        """
        e = self.e # alias for readability
        # Sources for root finding:
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html#scipy.optimize.OptimizeResult
        # https://stackoverflow.com/questions/43047463/find-root-of-a-transcendental-equation-with-python
        if e < 1:
            # Equations primarily from 5052 Appendix A
            M = self.M(time)
            def fun(E):
                return M - (E - e * np.sin(E))
            ans = root(fun, M) # Find the root of the transcendental, using M as the starting point for E
            assert ans.status == 1, "Failed to solve"
            E = ans.x[0]
            
            # Equation A.5 from Appendix A
            num = np.cos(E) - e
            den = 1 - e*np.cos(E)
            cosf = num / den
            
            # Equation A.6 from Appendix A
            num = np.sqrt(1-e**2) * np.sin(E)
            den = 1 - e*np.cos(E)
            sinf = num / den
            
            # Run both equations through arctan2 to handle quadrants
            f = np.arctan2(sinf, cosf)
        else:
            # Equations primarily from https://en.wikipedia.org/wiki/Hyperbolic_trajectory
            # and https://space.stackexchange.com/a/27604
            M = self.M(time)
            def fun(F):
                return M - (e*np.sinh(F) - F)
            ans = root(fun, M) # Find the root of the transcendental, using M as the starting point for F
            assert ans.status == 1, "Failed to solve"
            F = ans.x[0]
            
            # tanh(F/2) = sqrt((e-1)/(e+2))tan(f/2)
            f = 2*np.arctan(np.tanh(F/2) * np.sqrt((e+1)/(e-1)))
            
            E = F
            
        if angle_units == 'rad':
            ans = f, E
        elif angle_units == 'deg':
            ans = util.deg(f), util.deg(E)
        else:
            raise ValueError(f"Unknown angle_units: {angle_units}")
        return ans
    
    kepler = predict # Alias
    
    def position_velocity(self, time=None):
        """
        Calculates the position and velocity vectors in intertial space for a given
        set of orbital elements
        
        Calls Orbit.predict to calculate f and E, then calculates position and
        velocity in Cartesian coordinates
    
        Parameters
        ----------
        time : numeric (optional)
            if given, return position and velocity at the given absolute time
    
        Returns
        -------
        position : 3-vector of numerics
            The position vector in 3 space at the given time (units: km)
        velocity : 3-vector of numerics
            The velocity vector in 3 space at the given time (units: km)
        """
        
        # Equations from class notes: Appendix_A
        ele = self.ele # rename for readability, and so it can be modified safely
        if time is not None:
            f, E = self.predict(time)
            ele = ele._replace(f=f)
        
        # Equation A.11
        r = ele.a * (1 - ele.e**2) / (1 + ele.e * np.cos(ele.f))
        
        # Equation A.16
        x = np.sin(ele.LAN)*np.sin(ele.i)
        y = -np.cos(ele.LAN)*np.sin(ele.i)
        z = np.cos(ele.i)
        H_hat = [x, y, z]
        # Equation A.17
        x = np.cos(ele.AoP)*np.cos(ele.LAN) - np.cos(ele.i)*np.sin(ele.AoP)*np.sin(ele.LAN)
        y = np.cos(ele.AoP)*np.sin(ele.LAN) + np.cos(ele.i)*np.sin(ele.AoP)*np.cos(ele.LAN)
        z = np.sin(ele.AoP)*np.sin(ele.i)
        e_hat = np.array([x, y, z])
        # Equation A.18
        e_perp = np.cross(H_hat, e_hat)
        
        # Equation A.20
        position = r*(np.cos(ele.f) * e_hat + np.sin(ele.f) * e_perp)
        
        # Equation A.21
        velocity = np.sqrt(self.mu / (ele.a*(1-ele.e**2))) * (-np.sin(ele.f)*e_hat + (ele.e + np.cos(ele.f))*e_perp)
        
        return position, velocity

    def plot(self, start=0, stop=1, steps=100, mode='period', threeD=True, fig=None, plot_point=False, fig_kwargs={}, plot_kwargs={}):
        # Avoid persistence between calls
        if fig_kwargs == {}: fig_kwargs = {}
        if plot_kwargs == {}: plot_kwargs = {}
        
        if stop == 'tof':
            stop = self.tof
            mode = 'time'
            assert start == 0
        
        mode = mode.lower()
        if mode not in ['time', 'times', 'period', 'periods']:
            raise ValueError("Mode must be time/times or period/periods")
            
        new_plot = fig is None
        
        # Normalize to using absolute time
        if mode.startswith('period'):
            T = self.T
            start = start * T
            stop = stop * T
            
        # Setup the 3D figure
        if new_plot:
            fig = plt.figure(**fig_kwargs)
            if threeD:
                ax = plt.axes(projection='3d')
            else:
                ax = plt.axes()
        else:
            assert len(fig.axes) == 1
            ax = fig.axes[0]
            threeD = ax.name == "3d" # https://stackoverflow.com/a/43563348/14501840
        
        # Draw a dot at (0,0,0) for convenience
        if new_plot:
            ax.plot([0], [0], 'k.', label='Central Body')
        else:
            pass # Assume it was already drawn
        
        ### Set values for different things if none are specified
        # Use the planet name as the label
        if hasattr(self, "planet"): plot_kwargs.setdefault('label', self.planet)
        # Use the default color
        if hasattr(self, "color"): plot_kwargs.setdefault('color', self.color)
        
            
        p = [self.position_velocity(time=t)[0] for t in np.linspace(start, stop, steps)]
        p = np.array(p)
        if threeD:
            ax.plot(p[:,0], p[:,1], p[:,2], **plot_kwargs)
        else:
            threshold = abs(p[:,0:2]).max() # Use max instead of min because of axis crossings
            if any(abs(p[:,2]) > threshold*.05):
                print(f"WARNING: Orbit Z value = {p[:,2].max()/threshold * 100:0.3f}% of X/Y")
            ax.plot(p[:,0], p[:,1], **plot_kwargs)
        
        if plot_point:
            self.plot_point(threeD=threeD, fig=fig, plot_kwargs=plot_kwargs)
    
        fig.legend(['Central Body', 'Orbit'])
        ax.set_xlabel("X (km)")
        ax.set_ylabel("Y (km)")
        if threeD:
            ax.set_zlabel("Z (km)")
        
        return fig, ax
    
    def plot_point(self, time=0, mode='period', threeD=True, fig=None, fig_kwargs={}, plot_kwargs={}):
        """Plot the planet's position at the given time
        Time can be either absolute or a multiple of orbit period"""
        # Avoid persistence between calls
        if fig_kwargs == {}: fig_kwargs = {}
        if plot_kwargs == {}: plot_kwargs = {}
        
        new_plot = fig is None
        
        # Normalize to using absolute time
        if mode.startswith('period'):
            T = self.T
            time = time * T
        
        # Setup the 3D figure
        if new_plot:
            fig = plt.figure(**fig_kwargs)
            if threeD:
                ax = plt.axes(projection='3d')
            else:
                ax = plt.axes()
        else:
            assert len(fig.axes) == 1
            ax = fig.axes[0]
            threeD = ax.name == "3d" # https://stackoverflow.com/a/43563348/14501840
        
        # Draw a dot at (0,0,0) for convenience
        if new_plot:
            ax.plot([0], [0], 'k.')
        else:
            pass # Assume it was already drawn
            
        if 'marker' not in plot_kwargs:
            plot_kwargs['marker'] = 'o'
            
        p = self.position_velocity(time=time)[0]
        if threeD:
            ax.plot(p[0], p[1], p[2], **plot_kwargs)
        else:
            threshold = abs(p[0:2]).max() # Use max instead of min because of axis crossings
            if abs(p[2]) > threshold*.05:
                print(f"WARNING: Point Z value = {p[2]/threshold * 100:0.3f}% of X,Y")
            ax.plot(p[0], p[1], **plot_kwargs)
    
        fig.legend(['Central body', 'Orbit'])
        ax.set_xlabel("X (km)")
        ax.set_ylabel("Y (km)")
        if threeD:
            ax.set_zlabel("Z (km)")
        
        return fig, ax
            

if __name__ == "__main__":
    def test():
        # Circular on X axis
        r = [1,0,0]
        v = [0,1,0]
        o = Orbit.from_pvt(position=r, velocity=v, mu=1)
        assert o.ele == util.Orbital_Elements(e=0, a=1, i=0, LAN=0, AoP=0,
                                              f=0, tau=0)
        
        # Circular on Y axis
        r = [0,1,0]
        v = [-1,0,0]
        o = Orbit.from_pvt(position=r, velocity=v, mu=1, t0=np.pi/2)
        assert o.ele == util.Orbital_Elements(e=0, a=1, i=0, LAN=0, AoP=0,
                                              f=np.pi/2, tau=0)
        
        # Elliptical
        r = [1,0,0]
        v = [0,np.sqrt(1.5),0]
        o = Orbit.from_pvt(position=r, velocity=v, mu=1)
        target = util.Orbital_Elements(e=0.5, a=2, i=0, LAN=0, AoP=0,
                                       f=0, tau=0)
        assert all(abs(oe - te) < 1e-6 for oe,te in zip(o.ele, target))
        
        # Parabolic
        # r = [1,0,0]
        # v = [0,np.sqrt(2)/2,0]
        # mu = 1
        # o = orbit.from_pvt(position=r, velocity=v, mu=mu)
        
        # Hyperbolic
    test()
    
    
    # o = orbit.from_pvt(position=[6e3, 6e4, 6e3], velocity=[-5,5,0], mu=4e5)