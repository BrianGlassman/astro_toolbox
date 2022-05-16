# -*- coding: utf-8 -*-
"""
Created on Sun May 15 14:03:30 2022

@author: Brian
"""

from functools import partial
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

if __name__ == '__main__':
    import os
    os.chdir('../..')

from astro_toolbox import Orbit, util
from astro_toolbox import meeus # For checking for valid planets

def gridconfigure(obj, rw=None, cw=None):
    """
    Convenience function to configure the grid for a TKinter object

    Parameters
    ----------
    obj : TKinter object
        The object to configure
    rw : list, optional
        List of row weights. If None, make one row of weight 1
    cw : list, optional
        List of column weights. If None, make one column of weight 1
    """
    if rw is None:
        rw = [1]
    if cw is None:
        cw = [1]
    
    for i, weight in enumerate(rw):
        obj.rowconfigure(i, weight=weight)
    for i, weight in enumerate(cw):
        obj.columnconfigure(i, weight=weight)
        
class Model:
    def __init__(self, *, master):
        self.planet_orbits = {}
        self.planet_lines = {}
        self.planet_dots = {}
        self.planet_positions = {}
        self.state = 'free' # For locking the model so only one window can edit at a time
        
        self.last_date = None # Initialize to None so that first update call runs
        self.current_date = tk.DoubleVar(value=2440422.5, name='Current Date')
        
        ### Create the orbit plot
        # https://matplotlib.org/3.1.0/gallery/user_interfaces/embedding_in_tk_sgskip.html
        fig = Figure(figsize=(5, 5), dpi=100)
        fig.add_subplot(111)
        canvas = FigureCanvasTkAgg(fig, master=master)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.fig = fig ; self.canvas = canvas
        
    def init_planet(self, planet):
        date = self.current_date.get()
        orbit = Orbit.from_meeus(planet, date)
        self.planet_orbits[planet] = orbit
        self.planet_positions[planet] = {date: orbit.position_velocity(util.days_to_secs(date))[0]}
        
        fig, ax = orbit.plot(fig=self.fig, plot_point=True)
        self.planet_lines[planet] = ax.lines[-2]
        self.planet_dots[planet] = ax.lines[-1]
        
        self._update_plot()
        
    def destroy_planet(self, planet):
        print("TODO destroy " + planet)
        
    def change_date(self, new_date):
        if new_date == self.current_date:
            return
        
        self._update_plot()
        
    def _update_planet(self, planet):
        date = self.current_date.get()
        # Get position for this date (generate it if needed)
        pos_dict = self.planet_positions[planet]
        position = pos_dict.get(date)
        if position is None:
            position = self.planet_orbits[planet].position_velocity(util.days_to_secs(date))[0]
            pos_dict[date] = position
        # Update the plot data
        self.planet_dots[planet].set_data(position[0], position[1])
        
    def _update_plot(self):
        ### Update planet positions
        for planet in self.planet_dots.keys():
            self._update_planet(planet)
        
        self.canvas.draw_idle()
        
class Window:
    def __init__(self, *args, model, title='', bg='yellow', **kwargs):
        self.tk = tk.Toplevel(bg=bg, *args, **kwargs)
        self.model = model
        
        self.tk.withdraw()
        self.tk.title(title)
        
        # Set action when X is clicked
        # (Doesn't seem to change what happens when Root is closed)
        self.tk.protocol('WM_DELETE_WINDOW', self.onDestroy)
            
    def onDestroy(self):
        self.tk.withdraw()
Placeholder = Window
        
class TargetsWindow(Window):
    def __init__(self, *args, model, title, bg='red', **kwargs):
        super().__init__(*args, model=model, title=title, bg=bg, **kwargs)
        
        ### Text box
        text = tk.Text(self.tk, height=8)
        text.insert('0.0', 'Earth')
        text.pack()
        self.text = text
        
        ### Button
        button = ttk.Button(self.tk)
        button.pack()
        self.button = button
        self.old_vals = []
        self.unlock()
        self.lock()
        
    def unlock(self):
        '''Unlock the window, lock the model'''
        assert self.model.state == 'free' # Make sure nothing else is being edited already
        self.model.state = 'editing targets' # Lock to this window
        
        ### Change state
        self.text['state'] = 'normal'
        self.text['bg'] = 'white'
        self.button['text'] = 'Save'
        self.button['command'] = self.lock
        
    def lock(self):
        '''Save and publish the updates, unlock the model'''
        assert self.model.state == 'editing targets'
        
        old_vals = self.old_vals
        new_vals = self.text.get('0.0', tk.END)
        new_vals = [x.strip() for x in new_vals.split('\n')]
        new_vals = [x for x in new_vals if x]
        
        ### Check for invalid planets
        # Return early (without changing state) if any invalid planets found
        if new_vals != old_vals:
            invalid = [val for val in new_vals if val not in meeus.meeus.keys()]
            if invalid:
                msg = 'Invalid planet'
                if len(invalid) > 1:
                    msg += 's'
                msg += ': '
                msg += ', '.join(invalid)
                tk.messagebox.showerror(title='Invalid Planet',
                                        message=msg)
                return
        
        ### Update with the new values
        if new_vals != old_vals:
            # Initialize any new planets
            for val in new_vals:
                if val not in old_vals:
                    self.model.init_planet(val)
            # Destroy any removed planets
            for val in old_vals:
                if val not in new_vals:
                    self.model.destroy_planet(val)
            
            self.old_vals = new_vals
        
        ### Change state
        self._disable()
    def _disable(self):
        '''Helper function to handle the state change itself'''
        self.text['state'] = 'disabled'
        self.text['bg'] = 'light grey'
        self.button['text'] = 'Edit'
        self.button['command'] = self.unlock
        
        self.model.state = 'free'
        
class DatesWindow(Window):
    def __init__(self, *args, model, title, bg='blue', **kwargs):
        super().__init__(*args, model=model, title=title, bg=bg, **kwargs)
        
        ### Text box
        text = tk.Text(self.tk, height=8)
        text.insert('0.0', 'Earth')
        text.pack()
        self.text = text
        
        ### Button
        button = ttk.Button(self.tk)
        button.pack()
        self.button = button
        self.old_vals = []
        self.lock()
        
    def unlock(self):
        '''Unlock the window, lock the model'''
        assert self.model.state == 'free' # Make sure nothing else is being edited already
        self.model.state = 'editing dates' # Lock to this window
        
        ### Change state
        # self.text['state'] = 'normal'
        # self.text['bg'] = 'white'
        self.button['text'] = 'Save'
        self.button['command'] = self.lock
        
    def lock(self):
        '''Save and publish the updates, unlock the model'''
        assert self.model.state == 'editing dates'
        
        ### Change state
        self._disable()
    def _disable(self):
        '''Helper function to handle the state change itself'''
        # self.text['state'] = 'disabled'
        # self.text['bg'] = 'light grey'
        self.button['text'] = 'Edit'
        self.button['command'] = self.unlock
        
        self.model.state = 'free'
        
class PlotsWindow:
    def __init__(self, *args, bg='black', **kwargs):
        self.frame = ttk.Frame(bg=bg, *args, **kwargs)

class DetailsWindow:
    def __init__(self, *args, bg='green', **kwargs):
        self.frame = ttk.Frame(bg=bg, *args, **kwargs)

# TODO
'''
Separate window for 3D plot
connection to Horizons
'''

class Main:
    def __init__(self, title='', width=600, height=600, *args, **kwargs):
        self.window = tk.Tk()
        gridconfigure(self.window, [3, 2], [2, 3])
        self.window.title(title)
        self.window.geometry(f'{width}x{height}')
        
        ### Model
        self.model = Model(master=self.window)
        
        ### Top menu bar
        self.menuBar = tk.Menu(master=self.window)
        self.windowMenu = tk.Menu(self.menuBar, tearoff=0)
        self.contents = {}
        for tag, cls in {'Targets': TargetsWindow,
                         'Dates': Placeholder, 'Plots': Placeholder, 'Details': Placeholder}.items():
            self.windowMenu.add_command(label=tag, command=partial(self.open_window, window=tag))
            self.init_window(cls, tag)
        self.menuBar.add_cascade(label="Windows", menu=self.windowMenu)
        self.window.configure(menu=self.menuBar)
        
        ### Time slider
        value = self.model.current_date.get()
        slider = ttk.Scale(master=self.window, variable=self.model.current_date, from_=value-180, to=value+180)
        slider['command'] = self.model.change_date
        slider.pack(anchor='s', fill='x')
        self.slider = slider
        
        
        
        ###
        self.window.mainloop()
        
    def init_window(self, cls, title):
        self.contents[title] = cls(title=title, model=self.model)
        
    def open_window(self, window):
        self.contents[window].tk.deiconify()

main = Main()

if __name__ == '__main__':
    os.chdir('astro_toolbox/Vizualizer')
