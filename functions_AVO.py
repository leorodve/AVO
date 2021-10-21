"""
Free to use

@author: Leonardo Rodriguez
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from obspy.io.segy.core import _read_segy
from obspy import *
from scipy.interpolate import CubicSpline

#Upload and read files defined by the user
def upload_segy():
    global segy_filename
    segy_filename = filedialog.askopenfilename(initialdir="/", title="Select A SEG-Y File", filetypes=(("SEG-Y files", "*.segy, *.sgy"),("all files", "*.*")))

def upload_velocity():
    global velocity_filename
    velocity_filename = filedialog.askopenfilename(initialdir="/", title="Select A Velocity File", filetypes=(("Velocity files", "*.txt"),("all files", "*.*")))

def read_vel_file():
    load_file = open(velocity_filename, "r")
    read_file = load_file.read()
    txt_splitted = read_file.splitlines()
    x = []
    vel_values_zero = np.zeros([ns,ng])
    for i, line in enumerate(txt_splitted):
        if line.startswith("HANDVEL"):
                x.append(txt_splitted[i+1].split())
    k=0
    for i, nv in enumerate(x):
        a = (len(nv) // 2) - 1
        b = 0
        c = 0
        if a > 0:
            while c <= a:
                j=(int(nv[b])//sr)-1
                vel_values_zero[j,k]=nv[b+1]
                c += 1
                b += 2
            k += 1
    return vel_values_zero

def create_vel_function():
    vel_values_file = read_vel_file()
    for i, a in enumerate(vel_values_file):
        for j, b in enumerate(a):
            if b == 0:
                vel_values_file[i,j] = np.nan
    vel_values_file = pd.DataFrame(vel_values_file)
    vel_values_file = vel_values_file.interpolate(method='linear', limit_direction='both', axis=0)
    vel_values_file = vel_values_file.to_numpy()
    vel_values_file = np.nan_to_num(vel_values_file)
    return vel_values_file

#Apply NMO correction to gather
def nmo_correction(cmp_gather, dt, offsets, velocities): 
    nmo = np.zeros_like(cmp_gather) 
    nsamples = cmp_gather.shape[0] 
    times = np.arange(0, nsamples*dt, dt) 
    for i, t0 in enumerate(times):
        for j, x in enumerate(offsets):
            t = reflection_time(t0, x, velocities[i])
            amplitude = sample_trace(cmp_gather[:, j], t, dt)
            if amplitude is not None: 
                nmo[i, j] = amplitude 
    return nmo

def reflection_time(t0, x, vnmo): 
    t = np.sqrt(t0**2 + x**2/vnmo**2) 
    return t

def sample_trace(trace, time, sr): 
    before = int(np.floor(time/sr)) 
    N = trace.size 
    a = before - 1
    b = before + 3
    samples = np.arange(a, b) 
    if any(samples < 0) or any(samples >= N): 
        amplitude = None 
    else: 
        times = sr*samples 
        amps = trace[samples] 
        interpolator = CubicSpline(times, amps) 
        amplitude = interpolator(time) 
    return amplitude

#Calculate AVO attributes
def rp_g(gather, vp, cmp_number):
    global nmo_gather, angle
    i = 0
    j = 0
    angle = np.zeros([1,fold])
    Intercept = np.zeros([ns,1])
    Gradient = np.zeros([ns,1])
    Offset = offset_calc(cmp_number)
    nmo_gather = nmo_correction(gather, sr/1000, Offset, vp)
    while i < ns:
        while j < fold:
            angle[0, j] = np.arcsin(Offset[j] / (vp[i] * ns * sr / 1000))
            j += 1
        x = np.power(np.sin(angle[0,:].reshape((-1,1))),2)
        model = LinearRegression().fit(x, nmo_gather[i,:])
        Intercept[i,0] = model.intercept_
        Gradient[i,0] = model.coef_
        j = 0
        i += 1
    return Intercept[:,0], Gradient[:,0], nmo_gather

def offset_calc(cmp_no):
    i = cmp_no * fold - fold
    offset = np.zeros(fold)
    j = 0
    while j < fold:
        offset[j] = st.traces[i].stats.segy.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group
        i+=1
        j+=1
    return offset

#Display Intercept and Gradient
def plot_data(plot1, plot2, vm, vm1, window):    
    figure = Figure(figsize=(10, 4), dpi=100)
    
    ax1 = figure.add_subplot(1, 2, 1)
    ax1.imshow(plot1, extent=[0,ng,(ns-1)*sr,0], cmap='brg', vmin=-vm, vmax=vm, aspect='auto')
#    ax1.imshow(plot1, extent=[0,ng,(ns-1)*sr,0], cmap='Greys', vmin=-vm, vmax=vm, aspect='auto')
    ax1.set_xlim(0,ng)
    ax1.set_ylim(ns*sr,0)
    ax1.set_title('Intercept (Rp)')
    ax1.set_xlabel('CMP Number')
    ax1.set_ylabel('Time [ms]')
    
    ax2 = figure.add_subplot(1, 2, 2)
    ax2.imshow(plot2, extent=[0,ng,(ns-1)*sr,0], cmap='brg', vmin=-vm1, vmax=vm1, aspect='auto')
    ax1.set_xlim(0,ng)
    ax1.set_ylim(ns*sr,0)
    ax2.set_title('Gradient (G)')
    ax2.set_xlabel('CMP Number')
    ax2.set_ylabel('Time [ms]')
    
    canvas = FigureCanvasTkAgg(figure, window)
    return canvas

#User Interface
def display_data_window():
    global ng, nt, sr, ns, fold, st
    for widget in root.winfo_children():
        widget.destroy()

    #   Read segy with obspy
    st = _read_segy(segy_filename)

    nt = len(st.traces)                      #Number of traces in SEGY file
    sr = int(st.traces[0].stats.segy.trace_header.sample_interval_in_ms_for_this_trace/1000) #Sample rate (ms)
    ns = len(st.traces[0].data)             #Number of samples per trace
    fold = fold_entry.get()                 #Fold given by the user 
    ng = int(nt / fold)
    Section = np.zeros((ns, fold))
    nmo_section = np.zeros((ns, fold))
    vel = create_vel_function()
    stack = np.zeros([ns,ng])
    Rp = np.zeros([ns,ng])
    G = np.zeros([ns,ng])
    i = 0
    while i <= nt-1:
        k = 0
        j = int(i // fold) + 1
        while j == int(i // fold) + 1:
            Section[:,k] += st.traces[i].data[:]
            i += 1
            k += 1
        Intercept (Rp) and Gradient (G)
        Rp[:,j-1], G[:,j-1], nmo_section = rp_g(Section, vel[:,j-1], j)
        for a, amplitude in enumerate(nmo_section):
            stack[a,j-1] = sum(amplitude)
    vm = np.percentile(Rp, 99)
    vm1 = np.percentile(G, 99)
    first_plot = plot_data(Rp, G, vm, vm1, root)
    first_plot.get_tk_widget().grid(row=0, column=0, columnspan=2)

    return root    

def welcome_window():
    #labels and grid
    ttk.Button(root, text="Please select a CMP ordered SEG-Y file", command=upload_segy).grid(column=0, row=1, columnspan=3)
    ttk.Button(root, text="Please select a Velocity file", command=upload_velocity).grid(column=0, row=2, columnspan=3)
    ttk.Label(root, text="Fold").grid(column=0, row=3)
    third_entry = ttk.Entry(root, textvariable=fold_entry)
    third_entry.grid(column=1, row=3, columnspan=2)
    ttk.Button(root, text="Next", command=display_data_window).grid(column=0, row=4, columnspan=3)
    
    return root

#window widget asking the user for file/parameters
root = Tk()
root.title("AVO Analysis Program")

#variables (default values)
fold_entry = IntVar()
fold_entry.set("10")

welcome_frame = welcome_window()
root.mainloop()