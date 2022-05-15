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
        self.planet_lines = {}
        self.planet_orbits = {}
        self.state = 'free' # For locking the model so only one window can edit at a time
        
        ### Create the orbit plot
        # https://matplotlib.org/3.1.0/gallery/user_interfaces/embedding_in_tk_sgskip.html
        fig = Figure(figsize=(5, 5), dpi=100)
        fig.add_subplot(111)
        canvas = FigureCanvasTkAgg(fig, master=master)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.fig = fig ; self.canvas = canvas
        
    def init_planet(self, planet, date=None):
        print("Init " + planet)
        if date is None:
            date = 2440422.5 # Moon Landing
            
        orbit = Orbit.from_meeus(planet, date)
        self.planet_orbits[planet] = orbit
        
        orbit.plot(fig=self.fig)
        
        self.update_plot()
        
    def destroy_planet(self, planet):
        print("TODO destroy " + planet)
        
    def update_plot(self):
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
        self.lock()
        
    def unlock(self):
        '''Unlock the list for editing'''
        assert self.model.state == 'free' # Make sure nothing else is being edited
        
        ### Change state
        self.text['state'] = 'normal'
        self.text['bg'] = 'white'
        self.button['text'] = 'Save'
        self.button['command'] = self.lock
        
    def lock(self):
        '''Save the list and publish the updates'''
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
        
class DatesWindow:
    def __init__(self, *args, bg='blue', **kwargs):
        self.frame = ttk.Frame(bg=bg, *args, **kwargs)
        
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
        
        
        
        ### TEMP
        self.open_window('Targets')
        
        
        
        self.window.mainloop()
        
    def init_window(self, cls, title):
        self.contents[title] = cls(title=title, model=self.model)
        
    def open_window(self, window):
        self.contents[window].tk.deiconify()

main = Main()

if __name__ == '__main__':
    os.chdir('astro_toolbox/Vizualizer')
