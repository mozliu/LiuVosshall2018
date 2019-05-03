# Revisiting magnotether data
# Originally used to detect saccades, located at saccades.py
# Molly Liu, 5/14/15

from neo.io import WinEdrIO
from copy import deepcopy
import numpy as np
from scipy.io import loadmat
#from colours import *
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colors
from scipy.stats import circmean, circstd
import scipy.signal as signal
from scipy.optimize import leastsq
#from plotting_help import keep_axes
import sys, os
import circular
import datetime

""" Tentative directory of .EDR data (as of 5/14/15):
Ch.0: orientation (range 0 - 4000mV corresponding to 0-360)
Ch.1: xpos of the arena, from MATLAB. Detects position of pattern on arena
Ch.2: ypos of the arena, from MATLAB. Pattern ID on random contours, stimulus size on loom
Ch.3: stim, from MATLAB. Coded for stimulus shape.
[Ch.4: tachometer --not present for experiments so far]
"""

def read_winEDR(edr_filename):
    """ Reads in edr_filename as a dictionary with keys indicating the channel """
    fh = WinEdrIO(filename=edr_filename)
    seg = fh.read_segment(lazy=False, cascade=True,)
    
    analog_signals_dict = {}
    for analog_signal in seg.analogsignals:
        analog_signals_dict[analog_signal.name.lower()] = analog_signal
    
    # Format of analog_signals_dict: {'ch.0': AnalogSignal(array([value],dtype=float32)*mV,[start,end],sampling rate: rate),...}
    return analog_signals_dict

def data2deg(orientation,xpos):
    """ Maps orientation and xpos to matching degrees. Modified 5/15/15 """
    xpos_deg = xpos*360./10000
    orientation_deg = (382 - orientation*(360./4000))%360
    return orientation_deg, xpos_deg
    
def simplify(analogsignal):
    return (analogsignal/analogsignal.units).simplified

def my_sin(x,freq,phase,amp,offset):
    return np.sin(x*freq + phase)*amp + offset

def find_nearest_i(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def xpos_params(xpos,t, possible_freq=np.array([0,0.1,0.2,0.5,1.0,2.0,4.0,6.0,8.0]), possible_deg=np.array([0,15,60,120])):
    """ Finds speed information """
    xpos = cts_trace(xpos)
    t = t-t[0]
    maxpos = max(xpos[20:-10])
    minpos = min(xpos[20:-10])
    guess_off = (maxpos+minpos)/2
    if maxpos-minpos/2:
        return 0.0,-1,0,guess_off
    else:
        # amp, offset, phase all trivial
        guess_amp = (maxpos-minpos)/2
        x0 = minpos+guess_amp
        guess_phase = -(np.sign(np.mean(np.diff(xpos[:20])))-1)*(np.pi/2)

        # frequency: done by counting peaks, then dividing time by it. have to multiply by 6 for some unknown reason
        peaks = np.where(xpos > maxpos-1)[0]
        unique_peak = peaks[np.where(np.diff(peaks) > 1)[0]]
        troughs = np.where(xpos < minpos+1)[0]
        unique_trough = troughs[np.where(np.diff(troughs) > 1)[0]]
        guess_peak = round(6*(len(unique_peak)+1)/t[-1])
        guess_trough = round(6*(len(unique_trough)+1)/t[-1])
        guess_freq = (guess_peak+guess_trough)/2
        
        # fitting ended up too finicky, so switched 6/12 to picking from bank of possible frequencies and degrees. max_speed goes away this way
    #     optimize_func = lambda x: np.sin(t*x[0] + x[1])*x[2] + x[3] - xpos
    #     p0 = [guess_freq, guess_phase, guess_amp, guess_offset]
    #     fit = leastsq(optimize_func,p0)[0]
    #     est_freq, est_phase, est_amp, _ = fit
    #     
    #     # graphs for debugging
    #     first_sine = my_sin(t,*p0)
    #     fit_sine = my_sin(t, *fit)
    #     plt.plot(t,xpos, c='g',lw=2,label='xpos')
    #     plt.plot(t,first_sine, c='c',label='est')
    #     plt.plot(t,fit_sine, c='b',label='fit')
    #     plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #     
    #     # find max speed of fitted curve
    #     fit_sine = my_sin(t,*fit)
    #     max_speed = max(np.diff(fit_sine))

        # adjustments
        adj_freq = possible_freq[find_nearest_i(possible_freq,guess_freq/6)]
        phase_sign = np.sign(guess_phase)
        deg_range = possible_deg[find_nearest_i(possible_deg,abs(guess_amp*2))]

        return adj_freq, phase_sign, deg_range, guess_off

def cts_trace(trace_orig):
    """ Converts wraparounds to a continuous trace for use """
    trace = deepcopy(trace_orig)
    jumps = np.diff(trace)
    d_jumps = np.where(abs(jumps) > 340)[0].tolist()
    d_jumps = [x for x in d_jumps if x!=0]
    segs = np.split(trace,d_jumps)
    
    status = 0
    for i in range(len(d_jumps)):
        if jumps[d_jumps[i]] > 0:
            status -= 1
        elif jumps[d_jumps[i]] < 0:
            status += 1
        segs[i+1] += status * 360
        try:
            segs[i+1][0] = segs[i+1][1]
            segs[i+1][-1] = segs[i+1][-2]
        except IndexError: # occurs when seg consists of one point
            segs[i+1][0] = segs[i][-1]
        
    trace_cts = np.concatenate(segs)
    #plt.scatter(trace_orig.times, trace_orig)
    #plt.scatter(trace_cts.times, trace_cts)
    
    return trace_cts
    

def get_velocity(trace):
    """ Part of the saccade detection function. Returns a continuous orientation trace and a filtered velocity for use in saccade detection (not useful for saccade characteristics) """
    
    # make trace continuous to avoid dealing with 2*pi jumps
    trace_cts = cts_trace(trace)
    
    # apply 8th-order Butterworth filter
    N = 8
    Wn = .25
    B, A = signal.butter(N,Wn,output='ba')
    trace_f = signal.filtfilt(B,A,trace_cts)
    
    # take derivative of orientation to get angular velocity
    d_trace = np.gradient(trace_f)
    
    #f,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)
    #f.set_size_inches(30,9)
    #ax1.plot(trace)
    #ax2.plot(trace_f)
    #ax2.set_ylim([-360,0])
    #ax3.plot(abs(d_trace))
    #ax3.set_ylim([0,5])
    #ax3.set_xlim([10000,30000])
    
    return d_trace, trace_cts, trace_f
    
def get_offset(trace):
    offset = trace['orientation'] - trace['xpos']
    offset[offset > 180] = offset[offset > 180] - 360
    offset[offset < -180] = offset[offset < -180] + 360
    return offset
    
def offset_mean(offset):
    radoff = np.radians(offset)
    avgrad = np.arctan2(sum(np.sin(radoff)),sum(np.cos(radoff)))
    degoff = np.degrees(avgrad)
    return degoff

def id_contour(id,ax=None):
    # convert to binary
    get_bin = lambda x, n: x >= 0 and str(bin(x))[2:].zfill(n) or "-" + str(bin(x))[3:].zfill(n)
    teststr = get_bin(id,8)
    cid = np.empty(8,dtype='int')
    for i in range(8):
        cid[i] = int(teststr[7-i])
    newcid = 1-(np.reshape(cid,(4,2),order='F'))
    if ax is None:
        ax = plt.gca()
    ax.imshow(newcid,cmap='Greys',interpolation='nearest')

class ShapeTrial():
    """ trial of a magnotether run """
    def __init__(self, t, shape, orientation, xpos):#, pretrial):
        self.t0 = t[0]
        self.t = t - t[0]
        self.shape = shape[0]
        self.xpos = xpos
        self.orientation = orientation
            
        
        freq, phase, deg_range, self.guess_off = xpos_params(xpos,t)
        self.freq = freq
        self.phase = phase
        self.deg_range = deg_range
        
        offset = orientation - self.guess_off
        offset[offset > 180] = offset[offset > 180] - 360
        offset[offset < -180] = offset[offset < -180] + 360
        self.offset = offset
    
    def sin_xpos(self):
        return my_sin(self.t,self.freq*6,-(np.sign(np.mean(np.diff(self.xpos[:20])))-1)*(np.pi/2),self.deg_range/2,0)
    
    def plot(self,ax=None):
        if ax is None:
            ax = plt.gca()
        ax.plot(self.t,self.orientation,c='b')
        ax.plot(self.t,self.xpos,c='g')
        ax.set_ylim(0,360)
        ax.set_title('t0 = '+str(self.t0)+'s, '+str(self.freq)+'Hz, '+str(self.deg_range))
        
class ContourTrial():
    def __init__(self, t, orientation, xpos, ypos):
        self.t0 = t[0]
        self.t = t - t[0]
        self.xpos = xpos
        self.orientation = orientation
        
        freq, phase, deg_range, guess_off = xpos_params(xpos,t)
        self.freq = freq
        self.phase = phase
        self.deg_range = deg_range
        
        offset = orientation - guess_off
        offset[offset > 180] = offset[offset > 180] - 360
        offset[offset < -180] = offset[offset < -180] + 360
        self.offset = offset
        
        ckey = np.arange(0,10000,39.21176470588235)
        self.id = int(np.round(np.mean([ find_nearest_i(ckey,yp) for yp in ypos ]),0))
        
    def plot(self,ax0,ax1):
        id_contour(self.id,ax0)
        otherside = self.xpos+30
        otherside[otherside > 360] = otherside[otherside > 360] - 360
        ax1.plot(self.t,self.xpos,c='g')
        ax1.plot(self.t,otherside,c='g')
        ax1.plot(self.t,self.orientation,c='b')
        ax1.set_ylim(0,360)
        ax1.set_title(str(self.id)+': t0 = '+str(self.t0)+'s, '+str(self.freq)+'Hz, '+str(self.deg_range))

class LoomTrial():
    """ loom stimulus """
    def __init__(self, t, orientation, xpos, ypos):
        self.t0 = t[0]
        self.t = t - t[0]
        self.orientation = orientation
        
        # determining xposition of the loom
        lcr = [120,240,360]
        id = int(np.round(np.mean([ find_nearest_i(lcr,x) for x in xpos ]),0))
        self.xpos = (id+1)*24/96.*360
        
        offset = orientation - xpos
        offset[offset > 180] = offset[offset > 180] - 360
        offset[offset < -180] = offset[offset < -180] + 360
        self.offset = offset
        
        # determining duration of the loom: (all speeds at 56)
        self.ypos = ypos
        loomt = np.where(np.diff(ypos) > 50)[0]
        self.loomt = [self.t[loomt[0]],self.t[loomt[-1]]]
    
    def plot(self,ax=None):
        if ax is None:
            ax = plt.gca()
        ax.axvspan(*self.loomt,color='0.5',alpha=0.5)
        ax.axhline(0,color='k',alpha=0.5)
        ax.set_ylim(-180,180)
        ax.set_title('t0 = '+str(self.t0))
        ax.plot(self.t, self.offset)
        #ax.set_xlim(0,15)

class GratingTrial():
    """ vertical grating """
    def __init__(self, t, orientation, xpos):
        self.t0 = t[0]
        self.t = t - t[0]
        self.orientation = orientation
        self.xpos = xpos
    
    def line_xpos(self):
    	slope, intercept = np.polyfit(self.t,self.xpos,1)
    	line = self.t*slope + intercept
    	return slope, line
        
    def plot(self,ax=None):
        if ax is None:
            ax = plt.gca()
        ax.plot(self.t,self.orientation,c='b')
        ax.plot(self.t,self.xpos,c='g')
        ax.set_ylim(0,360)
        ax.set_title('t0 = '+str(self.t0))

class Magno():
    """ Each MagnOrientation object corresponds to one experiment """
    # ----- Processing functions ----- #
    
    def __init__(self, edr_filename, name=None):
        """ Initializes Magno """
        self.edr_fname = edr_filename
        self.basename = ''.join(self.edr_fname.split('.')[:-1])
        self.name = name
    
    def open_edr(self,stim_key={'bar':0,'bar_l':1084,'mbar':2169,'mbar_l':3254,'spot':4339,'spot_l':5423,'contour':6508,'contour_l':7593,'blank':8678,'blank_l':9763,'grating':-3050,'loom':-6101,'loom_l':-9152},sampling_rate=100,stim_margin=50,manual=False,split=[],ythres=-1000,end=True):
        """ Reads in EDR data """
        self.edr_raw = read_winEDR(self.edr_fname)
        print("Reading in "+self.basename+"...")
                
        # subsamples data to 100 Hz (WinEDR records at 10000 Hz--this method should be robust to changes in recording framerate)
        self.raw_sampling_rate = float(self.edr_raw['orientation'].sampling_rate)
        self.sampling_rate = sampling_rate
        self.subsampling_step = int(self.raw_sampling_rate/sampling_rate)
        
        # subsample orientation, xpos, ypos, stim
        orientation = simplify(self.edr_raw['orientation'][::self.subsampling_step])
        xpos = simplify(self.edr_raw['xpos'][::self.subsampling_step])
        ypos = simplify(self.edr_raw['ypos'][::self.subsampling_step])
        stim = simplify(self.edr_raw['stim'][::self.subsampling_step])
        self.t = np.array(self.edr_raw['orientation'].times[::self.subsampling_step])
        
        # categorize stimulus based on stim_key
        shape = np.array(['']*len(orientation),dtype='a10')
        for key,voltage in stim_key.items():
            in_range = (stim > voltage-stim_margin) & (stim < voltage+stim_margin)
            shape[in_range] = key
        
        dt = np.dtype([('t','f4'),('orientation','f4'),('xpos','f4'),('ypos','f4'),('shape','a10')])
        data = np.zeros((len(orientation),),dtype=dt)
        data['t'] = self.t
        data['orientation'], data['xpos'] = data2deg(orientation,xpos)
        data['ypos'] = ypos
        data['shape'] = shape
        
        self.data = data
        self.shape, self.contour, self.grating, self.loom = self.by_trial(manual,split,ythres,end)
        
        # potential for improvement: once you get it going, save as trials (list of structured arrays) instead of long structured array. the trials are what you will eventually care about
    
    def append(self,otherMagno):
        self.data = np.append(self.data,otherMagno.data)
        self.shape.extend(otherMagno.shape)
        self.loom.extend(otherMagno.loom)
        self.grating.extend(otherMagno.grating)
    
    def by_trial(self,manual=False,split=[],ythres=-1000,end=True):
        """ Splits trials based on changes in ypos """
        
        if not manual:
            ypos0 = np.where(self.data['ypos'] < ythres)[0]
            split = ypos0[np.where(np.diff(ypos0) > 1)[0]]
            
        trials = np.split(self.data,split)
        #pretrial = split[1:] - 200
        
        # end points are used to mark trial transition, so trim
        trials = [ t[2:-10] for t in trials if len(t) > 100 ]
        # cut off last trial (presumably shut off)
        if end:
            trials = trials[:-1]
        
        shape = []
        contour = []
        grating = []
        loom = []
        
        for t in trials:
            # skip if shape is not uniform
            if not (t['shape'] == t['shape'][0]).all():
                continue
            elif b'contour' in t['shape'][0]:
                contour.append(ContourTrial(t['t'], t['orientation'], t['xpos'], t['ypos']))
            elif b'grating' in t['shape'][0]:
                grating.append(GratingTrial(t['t'], t['orientation'], t['xpos']))
            elif b'loom' in t['shape'][0]:
                loom.append(LoomTrial(t['t'], t['orientation'], t['xpos'], t['ypos']))
            else:
                shape.append(ShapeTrial(t['t'], t['shape'], t['orientation'], t['xpos']))
        
        return shape, contour, grating, loom

# ------ Manipulating lists of trial objects ------ #

# trials are stored as a list of trial objects    
def get_trials(magno,original=[],type='shape'):
    """ Create or extend a list of Trial objects.
        magno: list of Magno objects to extract trials from
        original: give a list to extend
        type: the kind of trial desired """
    if type is 'shape':
        original.extend( [m.shape for m in magno] )
    elif type is 'grating':
        original.extend( [m.grating for m in magno] )
    elif type is 'loom':
        original.extend( [m.loom for m in magno] )
    else:
        print('Invalid type!')
    return original

# for shape trials, retrieve by different properties
def find_unique(shapetrial,dec=[1,0,1]):
    freq = [ s.freq for s in shapetrial ]
    deg = [ s.deg_range for s in shapetrial ]
    speed = [ s.max_speed for s in shapetrial ]
    
    freq = np.unique([round(f,dec[0]) for f in freq])
    deg = np.unique([round(d,dec[1]) for d in deg])
    speed = np.unique([round(s,dec[2]) for s in speed])
    
    return freq, deg, speed

def shape_parse(shapetrial,shape='all',freq='all',deg='all',speed='all',prev=0):
    original = shapetrial
    if shape is not 'all':
        try:
            original = get_shape(original,shape=shape,prev=prev)
        except ValueError:
            print('Warning: Shape is not in the trials given!')
    if freq is not 'all':
        try:
            original = get_freq(original,freq=freq)
        except ValueError:
            print('Warning: Frequency is not in the trials given!')
    if deg is not 'all':
        try:
            original = get_deg_range(original,deg=deg)
        except ValueError:
            print('Warning: Degree range is not in the trials given!')
    if speed is not 'all':
        try:
            original = get_max_speed(shapetrial,speed=speed)
        except ValueError:
            print('Warning: Max speed is not in the trials given!')
            
    return original

# parse by shape
def get_shape(shapetrial,shape,prev=0):
    original = []
    for i in range(len(shapetrial)):
        if shape == 'notblank' and shapetrial[i].shape not in ['bg', 'blank','grating','loom']:
            if prev == 0:
                original.extend([shapetrial[i]])
            elif (prev > 0) & (i > 0):
                newtrial = ShapeTrial(shapetrial[i].t,shapetrial[i].shape,shapetrial[i].orientation,shapetrial[i].xpos)
                prevt = shapetrial[i-1].t[-prev*100:] - shapetrial[i-1].t[-1]
                prevoff = shapetrial[i-1].orientation[-prev*100:] - shapetrial[i].guess_off
                newtrial.t = np.concatenate([prevt,newtrial.t])
                newtrial.offset = np.concatenate([prevoff,newtrial.offset])
                original.extend([newtrial])
            elif prev < 0:
                newtrial = ShapeTrial(shapetrial[i].t,shapetrial[i].shape,shapetrial[i].orientation,shapetrial[i].xpos)
                newtrial.t = newtrial.t[-prev*100:]
                newtrial.offset = newtrial.offset[-prev*100:]
                original.extend([newtrial])
        if shapetrial[i].shape == shape:
            if prev == 0:
                original.extend([shapetrial[i]])
            elif (prev > 0) & (i > 0):
                newtrial = ShapeTrial(shapetrial[i].t,shapetrial[i].shape,shapetrial[i].orientation,shapetrial[i].xpos)
                prevt = shapetrial[i-1].t[-prev*100:] - shapetrial[i-1].t[-1]
                prevoff = shapetrial[i-1].orientation[-prev*100:] - shapetrial[i].guess_off
                newtrial.t = np.concatenate([prevt,newtrial.t])
                newtrial.offset = np.concatenate([prevoff,newtrial.offset])
                original.extend([newtrial])
            elif prev < 0:
                newtrial = ShapeTrial(shapetrial[i].t,shapetrial[i].shape,shapetrial[i].orientation,shapetrial[i].xpos)
                newtrial.t = newtrial.t[-prev*100:]
                newtrial.offset = newtrial.offset[-prev*100:]
                original.extend([newtrial])
    return original

# frequency 
def get_freq(shapetrial,freq,dec=1):
    original=[]
    for t in shapetrial:
        if round(t.freq,dec) == freq:
            original.extend([t])
    return original

# degree range
def get_deg_range(shapetrial,deg,dec=-1):
    original=[]
    for t in shapetrial:
        if round(t.deg_range,dec) == deg:
            original.extend([t])
    return original

# max speed
def get_max_speed(shapetrial,speed,dec=1):
    original=[]
    for t in shapetrial:
        if round(t.max_speed,dec) == speed:
            original.extend([t])
    return original
    
def contour_parse(contourtrial,id,freq='all',deg='all',speed='all'):
    original = []
    for t in contourtrial:
        if t.id == id:
            original.extend([t])
    if freq is not 'all':
        try:
            original = get_freq(original,freq=freq)
        except ValueError:
            print('Warning: Frequency is not in the trials given!')
    if deg is not 'all':
        try:
            original = get_deg_range(original,deg=deg)
        except ValueError:
            print('Warning: Degree range is not in the trials given!')
    if speed is not 'all':
        try:
            original = get_max_speed(shapetrial,speed=speed)
        except ValueError:
            print('Warning: Max speed is not in the trials given!')
            
    return original
    
# Plotting functions
# raw traces
def see_raw(fly):
    i = datetime.datetime.now()
    tstamp = i.strftime('%Y%m%d')
    
    shape = set([ s.shape for f in fly for s in f.shape ])
    shape = [ s for s in shape if s != '' ]
    for s in shape:
        strial = np.concatenate([ shape_parse(f.shape,s) for f in fly ])
        f,ax = plt.subplots(len(strial),sharex=True,sharey=True)
        f.set_size_inches(6,2.5*len(strial))
        for i in range(len(strial)):
            strial[i].plot(ax[i])
        savestr = tstamp+'_'+s+'.pdf'
        f.savefig(savestr,bbox_inches='tight')
        plt.close(f)
    
    ctrial = np.concatenate([ f.contour for f in fly ])
    f = plt.figure(1,[8,2.5*len(ctrial)])
    for i in range(len(ctrial)):
        ax0 = plt.subplot2grid((len(ctrial),4),(i,0))
        ax1 = plt.subplot2grid((len(ctrial),4),(i,1),colspan=3)
        ctrial[i].plot(ax0,ax1)
    savestr = tstamp+'_contour.pdf'
    f.savefig(savestr,bbox_inches='tight')
    plt.close(f)
    
    gtrial = np.concatenate([ f.grating for f in fly ])
    f,ax = plt.subplots(len(gtrial),sharex=True,sharey=True)
    f.set_size_inches(6,2.5*len(gtrial))
    for i in range(len(gtrial)):
        gtrial[i].plot(ax[i])
    savestr = tstamp+'_grating.pdf'
    f.savefig(savestr,bbox_inches='tight')
    plt.close(f)
        
    ltrial = np.concatenate([ f.loom for f in fly ])
    f,ax = plt.subplots(len(ltrial),sharex=True,sharey=True)
    f.set_size_inches(6,2.5*len(ltrial))
    for i in range(len(ltrial)):
        ltrial[i].plot(ax[i])
    savestr = tstamp+'_loom.pdf'
    f.savefig(savestr,bbox_inches='tight')
    plt.close(f)

def see_raw_v(fly):
    i = datetime.datetime.now()
    tstamp = i.strftime('%Y%m%d')
    
    strial = np.concatenate([ shape_parse(f.shape,shape='bar') for f in fly ])
    f,ax = plt.subplots(len(strial)*2)
    f.set_size_inches(6,5*len(strial))
    for i in range(len(strial)):
        ov, _, of = get_velocity(strial[i].orientation)
        xv, _, xf = get_velocity(strial[i].sin_xpos())
        ax[2*i].plot(strial[i].t,of,c='b')
        ax[2*i].plot(strial[i].t,xf,c='g')
        ax[2*i].set_title(str(i)+', t0 = '+str(strial[i].t0)+'s, '+str(strial[i].freq)+'Hz, '+str(strial[i].deg_range))
        ax[2*i+1].plot(strial[i].t,ov,c='b')
        ax[2*i+1].plot(strial[i].t,xv,c='g')
        ax[2*i+1].set_title(str(i)+' velocity')
    savestr = tstamp+'_bar_v.pdf'
    f.savefig(savestr,bbox_inches='tight')
    plt.close(f)
    
    gtrial = np.concatenate([ f.grating for f in fly ])
    f,ax = plt.subplots(len(gtrial)*2)
    f.set_size_inches(6,5*len(gtrial))
    for i in range(len(gtrial)):
        ov, _, of = get_velocity(gtrial[i].orientation)
        xv, xf = gtrial[i].line_xpos()
        ax[2*i].plot(gtrial[i].t,of,c='b')
        ax[2*i].plot(gtrial[i].t,xf,c='g')
        ax[2*i].set_title(str(i)+', t0 = '+str(gtrial[i].t0)+'s')
        ax[2*i+1].plot(gtrial[i].t,ov,c='b')
        ax[2*i+1].plot(gtrial[i].t,xv,c='g')
        ax[2*i+1].set_title(str(i)+' velocity')
    savestr = tstamp+'_grating_v.pdf'
    f.savefig(savestr,bbox_inches='tight')
    plt.close(f)

# hist2d of offset    
def hist2d_time(fly,shape='all',freq='all',plot=True,ax=None,prev=3,stimwidth=15.0,tend=10):
    img = np.empty([len(fly),prev+tend,24])
    for k in range(len(fly)):
        bars = shape_parse(fly[k].shape,shape,freq,prev=prev)
        t = [a.t for a in bars]
        offset = [np.radians(a.offset+stimwidth/2) for a in bars]
        allimg = np.empty([len(offset),prev+tend,24])
        for i in range(len(t)):
            offset[i][offset[i] > np.pi] -= 2*np.pi
            offset[i][offset[i] < -np.pi] += 2*np.pi
            a,_,_ = np.histogram2d(t[i],offset[i],bins=[tend+prev,24],range=[[-1*prev,tend],[-np.pi,np.pi]])
            for j in range(len(a)):
                a[j] = a[j] / sum(a[j])
            allimg[i,:,:] = a
        img[k,:,:] = np.nanmean(allimg,axis=0)
    imgtoplot = np.rot90(np.nanmean(img,axis=0))
    
    if plot:
        if ax is None:
            ax = plt.gca()
        h = ax.imshow(imgtoplot,cmap='inferno')
        ax.axvline(prev-0.5,color='w',linewidth=1,alpha=1)
        ax.axhline(24*(3/8)-0.5,color='w',ls='--',alpha=1,linewidth=1)
        ax.axhline(24*(5/8)-0.5,color='w',ls='--',alpha=1,linewidth=1)
        ax.set_xticks(np.arange(-0.5,prev+tend,3))
        ax.set_xticklabels(np.arange(-1*prev,tend,3))
        ax.set_yticks([-0.5,24*(3/8)-0.5,24*(5/8)-0.5,24-0.5])
        ax.set_yticklabels([180,45,-45,180])
        #plt.colorbar()

    return imgtoplot,h

def hist2d_listgrid(mozlist,shape=[b'bar',b'mbar',b'spot',b'blank'],prev=5,tend=15,save=False,fname='',vmax=None):
    n = len(shape)
    m = len(mozlist)
    if vmax is None:
        setv = True
        vmax = 0
    else:
        setv = False

    f,ax = plt.subplots(m,n,sharex=True,sharey=True)
    if m==1:
        ax = [ax]
    h = np.empty([m,n],dtype='O')
    f.set_size_inches(2*n,1.5*m)
    for s in shape:
        for i in range(len(mozlist)):
            imgtoplot,h[i][shape.index(s)] = hist2d_time(mozlist[i],shape=s,ax=ax[i][shape.index(s)],prev=prev,tend=tend)
            if setv & (np.max(imgtoplot) > vmax):
                vmax = np.max(imgtoplot)

        #ax[0,shape.index(s)].set_xlabel(s)
        
    #ticks = [-np.pi,-np.pi/2.,0,np.pi/2.,np.pi]
    #ticklabels = ['-180$^\circ$','-90$^\circ$','0$^\circ$','+90$^\circ$','+180$^\circ$']
    norm = colors.Normalize(vmin=0,vmax=vmax)
    for i in range(len(ax)):
        for j in range(len(ax[i])):
            h[i][j].set_norm(norm)
            if (j == 0) & (i == m-1):
                ax[i][j].spines['top'].set_visible(False)
                ax[i][j].spines['right'].set_visible(False)
                ax[i][j].spines['bottom'].set_visible(False)
                ax[i][j].spines['left'].set_visible(False)
                #ax[i][j].set_yticks(ticks)
                #ax[i][j].set_yticklabels(ticklabels)
                ax[i][j].set_ylabel('offset')
                ax[i][j].set_xlabel('time (s)')
            else:
                ax[i][j].set_axis_off()
    ax[0][0].set_xlim([-0.5,19.5])
    ax[0][0].set_ylim([-0.5,23.5])
    plt.colorbar(h[-1][-1])
    #plt.tight_layout()
    
    if save:
        plt.savefig(fname)

    return vmax

# residence in one zone over time
def latency(fly,alloff,shape='bar',freq='all',deg='all',thres=None,tstep=0.5,stdthres=1):
    #f,ax = plt.subplots(2,sharex=True)
    #print 'run'
    #determine threshold
    if thres is None:
        thres = circular.std(alloff)*stdthres
    tbins = np.arange(0,15+tstep,tstep)
    if 'bar' in shape:
        obins = [-np.pi,-thres,thres,np.pi]
    elif 'spot' in shape:
        obins = [-np.pi,-np.pi+thres,np.pi-thres,np.pi]
    elif 'blank' in shape:
        obins = [-np.pi,-thres,thres,np.pi]
    
    # computing the percentage of traces within the threshold for each
    perin = np.empty([len(fly),len(tbins)-1])
    for i in range(len(fly)):
        allbar = shape_parse(fly[i].shape,shape=shape,freq=freq,deg=deg)
        if len(allbar)==0:
            perin[i].fill(np.nan)
        else:
            t = [a.t for a in allbar]
            t = [item for sublist in t for item in sublist]
            offset = [np.radians(a.offset) for a in allbar]
            offset = [item for sublist in offset for item in sublist]
            hist2d = np.histogram2d(t,offset,bins=[tbins,obins])
            if 'bar' in shape:
                oin = [ hist2d[0][j][1] for j in range(len(hist2d[0])) ]
            elif 'spot' in shape:
                oin = [ hist2d[0][j][0] + hist2d[0][j][2] for j in range(len(hist2d[0])) ]
            elif 'blank' in shape:
                oin = [ hist2d[0][j][1] for j in range(len(hist2d[0])) ]
                
            perin[i] = oin/np.sum(hist2d[0],axis=1)
    perin_mean = np.nanmean(perin,axis=0)
    
    return perin_mean, hist2d[1][:-1], obins, perin

# residence in multiple zones over time
def hist_segment(fly,shape='all',freq='all',deg='all',thres=np.pi/3,tstep=0.1):
    tbins = np.arange(0,15+tstep,tstep)
    obins = [-np.pi,-np.pi+thres,-thres,thres,np.pi-thres,np.pi]
    
    fix = np.empty([len(fly),len(tbins)-1])
    antifix = np.empty([len(fly),len(tbins)-1])
    neither = np.empty([len(fly),len(tbins)-1])
    
    for i in range(len(fly)):
        alloff = shape_parse(fly[i].shape,shape=shape,freq=freq,deg=deg)
        if len(alloff)==0:
            fix[i].fill(np.nan)
            antifix[i].fill(np.nan)
            neither[i].fill(np.nan)
        else:
            t = [a.t for a in alloff]
            t = [item for sublist in t for item in sublist]
            offset = [np.radians(a.offset) for a in alloff]
            offset = [item for sublist in offset for item in sublist]
            hist2d = np.histogram2d(t,offset,bins=[tbins,obins])
            fix_i = [ hist2d[0][j][2] for j in range(len(hist2d[0])) ]
            antifix_i = [ hist2d[0][j][0] + hist2d[0][j][4] for j in range(len(hist2d[0])) ]
            neither_i = [ hist2d[0][j][1] + hist2d[0][j][3] for j in range(len(hist2d[0])) ]
            
            total = np.sum(hist2d[0],axis=1)
            fix[i] = fix_i/total
            antifix[i] = antifix_i/total
            neither[i] = neither_i/total
    
    fix_mean = np.nanmean(fix,axis=0)
    antifix_mean = np.nanmean(antifix,axis=0)
    neither_mean = np.nanmean(neither,axis=0)
    
    return fix_mean,neither_mean,antifix_mean,hist2d[1][:-1]

# an example plot    
def generic_latency_plot(fly):
    tstep = 0.1
    stdthres = [0.5,1,1.5]
    linestyles = ['-','--',':']
    f,ax = plt.subplots(3,sharex=True)
    f.set_size_inches([6,9])
    freq = np.arange(0.1,2.0,0.1)
    baroff = np.concatenate([hist2d_time(fly,'bar',freq=fr,plot=False)[0] for fr in freq])
    spotoff = np.concatenate([hist2d_time(fly,'spot',freq=fr,plot=False)[0] for fr in freq])
    tbar = np.concatenate([hist2d_time(fly,'bar',freq=fr,plot=False)[1] for fr in freq])
    tspot = np.concatenate([hist2d_time(fly,'spot',freq=fr,plot=False)[1] for fr in freq])
    ax[0].hist2d(tbar,baroff,bins=[150,90],range=[[0,15],[-np.pi,np.pi]])
    ax[1].hist2d(tspot,spotoff,bins=[150,90],range=[[0,15],[-np.pi,np.pi]])
    for i in range(len(stdthres)):
        perin_bar,t,obins,barper = magno.latency(fly,baroff,shape='bar',tstep=tstep,stdthres=stdthres[i])
        ax[2].plot(t,perin_bar,ls=linestyles[i],c='b',drawstyle='steps')
        ax[0].axhline(obins[1],ls=linestyles[i],c='w')
        ax[0].axhline(obins[2],ls=linestyles[i],c='w')
        ax[2].axhline(obins[2]/np.pi,c='0.5',ls=linestyles[i],label=str(stdthres[i])+r' $\sigma$')
    
        perin_spot,t,obins,spotper = magno.latency(fly,spotoff,shape='spot',tstep=tstep,stdthres=stdthres[i])
        ax[2].plot(t,perin_spot,ls=linestyles[i],c='r',drawstyle='steps')
        ax[1].axhline(obins[1],ls=linestyles[i],c='w')
        ax[1].axhline(obins[2],ls=linestyles[i],c='w')
    ax[0].set_yticklabels(['-180$^\circ$','-180$^\circ$','-120$^\circ$','-60$^\circ$','0$^\circ$','+60$^\circ$','+120$^\circ$','+180$^\circ$'])
    ax[1].set_yticklabels(['-180$^\circ$','-180$^\circ$','-120$^\circ$','-60$^\circ$','0$^\circ$','+60$^\circ$','+120$^\circ$','+180$^\circ$'])
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0.)
    plt.savefig('offset-time-latency-moving_20150607.pdf')

# average histogram by trial in one fly
def hist_avg(shapetrial,shape='all',freq='all',deg='all',speed='all',nbin=72,center=True,plot=True,ax=None,c='0.25'):
    """ histogram of the average degree by shape """
    st = shape_parse(shapetrial,shape,freq,deg,speed)
    if len(st) == 0:
        avghist = np.zeros(nbin+1)
    else:
        rad = [ np.radians(s.offset) for s in st ]
        allhist = np.empty([len(rad),nbin+1])
        for i in range(len(rad)):
            # shifts bins so that you graph hist in the middle
            if center:
                rad[i] += np.pi/nbin
                rad[i][rad[i]>np.pi] -= 2*np.pi
            hist = np.histogram(rad[i],nbin,(-np.pi,np.pi),density=True)
            allhist[i] = np.append(hist[0],[hist[0][0]])
        
        avghist = np.mean(allhist,axis=0)
        
        if plot:
            if ax is None:
                ax = plt.gca()
            ax.set_xticklabels(['-90$^\circ$','-45$^\circ$','0$^\circ$','+45$^\circ$','+90$^\circ$','+135$^\circ$','$\pm$180$^\circ$','-135$^\circ$'])
            ax.set_yticklabels([''])
            ax.spines['polar'].set_visible(False)
            ax.plot(hist[1]+np.pi/2,avghist,c=c)
        
    return avghist

# average histogram of trial over multiple flies
def hist_avgavg(fly,shape='all',freq='all',deg='all',nbin=72,center=True,plot=True,ax=None,c='k',save=False,fname=''):
    
    avghist = np.empty([len(fly),nbin+1])
    for i in range(len(fly)):
        avghist[i] = hist_avg(fly[i].shape,shape,freq,deg,nbin=nbin,center=center,plot=plot,ax=ax,c='0.5')
    avgavg = np.nanmean(avghist,axis=0)
    avgavg *= 1/(sum(avgavg*2*np.pi/nbin))
    stepsize = 2*np.pi/nbin
    steps = np.arange(-np.pi,np.pi+stepsize,stepsize)
    if plot:
        if ax is None:
            ax = plt.gca()
            
        ax.plot(steps+np.pi/2,avgavg,lw=2,c=c)
        ax.set_ylim(0,1)
    
    if save:
        plt.savefig(fname)
        
    return avgavg, avghist


# c = 'bgrcmykb'
# for i in range(3):
#     f = figure(i+1,[4,6])
#     ax = plt.subplot(111,polar=True)
#     #ax.set_title(shapes[i])
#     #avg = np.radians(magno.offset_mean(offset[i]))
#     for j in range(len(freq)):
#         try:
#             div = magno.get_freq(st[i],freq[j])
#             radoff = np.radians(np.concatenate([d.offset for d in div]))
#             hist = np.histogram(radoff,72,(-np.pi,np.pi),density=True)
#             ax.plot(hist[1][:-1]+np.pi/2,hist[0],color=c[j],label=str(freq[j]))
#         except ValueError:
#             continue
#     ax.set_xticklabels(['-90$^\circ$','-45$^\circ$','0$^\circ$','+45$^\circ$','+90$^\circ$','+135$^\circ$','$\pm$180$^\circ$','-135$^\circ$'])
#     ax.set_yticklabels([''])
#     ax.spines['polar'].set_visible(False)
#     plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#     savestr = 'offset_by_shape_and_freq-'+str(i)+'.pdf'
#     savefig(savestr)

# shapes = ['bar','mbar','spot','blank']
# st = [0]*len(shapes)
# for i in range(len(shapes)):
#     st[i] = magno.get_shape(fly.shape,shapes[i])
# speed_names = ['0hz','0.1 15','0.1 60','0.1 120','0.5 15','0.5 60','0.5 120',
#                '1 15','1 60','1 120','??']
# c = ['b','r','g','c','m','y','k','0.75','0.5','0.25','burlywood','chartreuse','skyblue']
# for i in range(len(shapes)):
#     f = figure(i+1,[6,4])
#     ax = plt.subplot(111)
#     ax.set_title(shapes[i])
#     for j in st[i]:
#         if round(j.freq,1) == 0:
#             n = 0
#         elif round(j.freq,1) == 0.1:
#             n = 1
#         elif round(j.freq,1) == 0.5:
#             n = 4
#         elif round(j.freq,1) == 1.0:
#             n = 7
#         else: n = 10
#         
#         if round(abs(j.deg_range),-1) == 60.0:
#             n += 1
#         elif round(abs(j.deg_range),-1) == 120:
#             n += 2
#         
#         ax.plot(j.t,j.offset,'.',c=c[n], label=speed_names[n])
#     
#     plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


# ------ Global plotting functions ------ #

def polarhist_shape(data_to_plot,shapes=['bar','mbar','spot','blank'],nbins=180):
    """ Polar histogram of orientations in data, separated by shape """
    
    allshapes = {}.fromkeys(shapes)
    for s in allshapes.iterkeys():
        allshapes[s] = sum([ f.by_shape(shape=[s])[s] for f in data_to_plot ])
    
    for i in range(len(shapes)):
        total = np.concatenate(allshapes[shapes[i]])
        toffset = offset(total)
        hist = np.histogram(toffset*np.pi/180,180,(-np.pi,np.pi),normed=True)
        ax = subplot(len(shapes),1,i+1,polar=True)
        ax.plot(hist[1][:-1],hist[0])
        ax.set_title(shapes[i])

def dotplot(ff,mc='k',pc=['0.25'],width=0.75,pos=None,ax=None,xlim=None,jitter=0.10,fixline=True):
    try:
        n,_ = ff.shape
    except ValueError:
        n = 1
    except AttributeError:
        n = len(ff)
    if ax is None:
        plt.figure(1,[n,2])
    else:
        plt.sca(ax)
    
    # default colors
    if len(mc) < n:
        mc = [mc]*n
    if len(pc) < n:
        pc = [pc]*n
    if pos is None:
        pos = range(n)
    if xlim is None:
        xlim = n-0.5
    
    # add guidelines
    #plt.axhline(0,c='k',ls='--')
    #plt.axhline(1,c='k',ls='--')
    if fixline:
        plt.axhline(0.25,c='0.5',ls=':')
    
    for i in range(n):
        # strip np.nan from dataset (messes up the graph)
        y = ff[i]
        y = y[~np.isnan(y)]
        m = np.median(y)

        # plot medians
        plt.plot([pos[i]-width/2,pos[i]+width/2],[m,m],c=mc[i],linewidth=3)
        
        # plot points
        x = np.random.normal(pos[i],jitter,len(y))
        plt.scatter(x,y,c=pc[i],alpha=0.5,edgecolor='None')
    
    # cosmetic
    plt.xlim(-0.5,xlim)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_yticks([0,.2,0.4,0.6,0.8,1])
    ax.set_xticks([])
    
    return ax
    
# def plot_mayer_v_hist(trace,title,save=False,name=''):
#     """ Mayer-style angular velocity histogram for MagnOrientation objects in trace """
#     figure(1,[5,5])
#     velocity = [0]*len(trace)
#     for i in range(len(trace)):
#         trimmed_flight = trace[i].trim_bouts(min_bout_len_s=600)
#         velocity[i] = trace[i].velocity[trimmed_flight==1]
#     all_v = np.concatenate(velocity)
#     n, bins, patches = hist(all_v*2, normed=1, bins = range(-50,51,1), histtype='stepfilled')
#     setp(patches, 'facecolor', '0.75', 'alpha', 0.75)
#     yscale('log', nonposy='clip')
#     title(title)
#     ylabel('Normalized frequency')
#     xlabel('Angular displacement size (deg/20ms)')
#     xlim(-50,50)
#     
#     if save:
#         savefig(name)
# 
# def plot_all_orientations(trace,title='',save=False,name=''):
#     """ Orientation histogram including all MagnOrientation objects in trace """
#     o = [0]*len(trace)
#     no = [0]*len(trace)
#     for i in range(len(o)):
#         trimmed = trace[i].trim_bouts(min_bout_len_s=60)
#         o[i] = trace[i].flight['orientation'][trimmed==1]
#         no[i] = trace[i].flight['orientation'][trimmed==0]
#     total_flight = (numpy.concatenate(o) * np.pi/180)%(2*np.pi)
#     #total_nonflight = (numpy.concatenate(no) * np.pi/180)%(2*np.pi)
#     
#     nbins = 360
#     histogram = np.histogram(total_flight,nbins,(0,2*np.pi))
#     bounds = max(histogram[0])
#     
#     figure(1,[5,5])
#     ax = plt.subplot(111,polar=True)
#     ax.set_ylim(0,bounds)
#     ax.fill(histogram[1][:-1],[(sum(histogram[0])/nbins) for k in range(histogram[0].size)],ec='none',fc='0.75',label="No bias") # if orientation was perfectly unbiased
#     ax.plot(histogram[1][:-1],histogram[0],label="Orientation")
#     #ax.plot([np.mean(total_flight),np.mean(total_flight)],[0,bounds],c='b')
#             
#     # graph "preferred orientation"
#     """for total_nonflight in no:
#         preferred_orientation = np.mean(total_nonflight)
#         ax.plot([preferred_orientation,preferred_orientation],[0,bounds],c='0')"""
#             
#     # graph possible points of fixation
#     tachometer = 1640 * np.pi/2000
#     leds1 = 2600 * np.pi/2000
#     leds2 = 3600 * np.pi/2000
#     leds3 = 600 * np.pi/2000
#     leds4 = 1600 * np.pi/2000
#     ax.plot([tachometer,tachometer],[0,bounds],linestyle='--',c='0',label="Tachometer")
#     ax.plot([leds1,leds1],[0,bounds],linestyle='--',c='r',label="LED")
#     ax.plot([leds3,leds3],[0,bounds],linestyle='--',c='r')
#     ax.plot([leds2,leds2],[0,bounds],linestyle='--',c='r')
#     ax.plot([leds4,leds4],[0,bounds],linestyle='--',c='r')
#     legend(bbox_to_anchor=(1.05, .9), loc=2, borderaxespad=0.)
#     title(title)
#     tight_layout()
#     
#     if save:
#         savefig(name)
# 
# def plot_flight_length(traces,labels=None,flight=1,times='all',stim='all',save=False):
#     """ Plots length of flight bout vs. time of bout initiation """
#     all_bouts = [0]*len(traces)
#     p = [0]*len(traces)
#     colors = 'bgrcmyk'
#     plt.figure(1,[8,6])
#     ax = plt.subplot(1,1,1)
#     for i in range(len(traces)):
#             all_bouts[i] = np.concatenate([ t.flight_bouts(flight,times,stim) for t in traces[i] ])
#             p[i] = ax.scatter(all_bouts[i]['start'],all_bouts[i]['duration'])
#             
#     if labels is not None:
#             ax.legend(p,labels)
#     ax.set_xlabel('Time of bout initiation (min)')
#     ax.set_ylabel('Length of flight bout (min)')
#     ax.set_title('Flight bouts by stimulus and genotype')
#             
#     if save:
#             plt.savefig(traces[0].basename+'_flight-length.pdf')
# 
#     return [ np.mean(ab['duration']) for ab in all_bouts ]
# 
# def get_all_saccades(trace):
#     """ Concatenates saccade arrays from a list of MagnOrientation """
#     turns = [ t.saccades for t in trace ]
#     all_s = np.concatenate(turns)
#     
#     # removes nan values
#     with_nan = []
#     for i in range(len(all_s)):
#         for col in all_s[i]:
#             if np.isnan(col):
#                 with_nan.append(i)
#     
#     all_s = np.delete(all_s,with_nan)
#     
#     return all_s
# 
# def amplitude_vs_duration(all_s,title='',bins=[30,100],range=[[0,400],[0,200]],save=False,name=''):
#     """ Histogram of saccade amplitudes given an array of saccades """
#     plt.figure(1,[5,5])
#     plt.hist2d(abs(all_s['duration'])*1000,abs(all_s['amplitude']),bins=bins,range=range)
#     #plt.yscale('log', nonposy='clip')
#     #plt.title(title)
#     plt.xlabel('duration (ms)')
#     plt.ylabel('amplitude (deg)')
#     if save:
#         plt.savefig(name)
#         
# def amplitude_vs_velocity(all_s,title='',bins=[30,100],range=[[0,400],[0,200]],save=False,name=''):
#     """ Histogram of saccade amplitudes given an array of saccades """
#     plt.figure(1,[5,5])
#     plt.hist2d(abs(all_s['max_v'])*100,abs(all_s['amplitude']),bins=bins,range=range)
#     #plt.yscale('log', nonposy='clip')
#     #plt.title(title)
#     plt.xlabel('|peak velocity| (deg/s)')
#     plt.ylabel('amplitude (deg)')
#     if save:
#         plt.savefig(name)
# 
# def hist_amplitude(all_s,title='',save=False,name=''):
#     """ Histogram of saccade amplitudes given an array of saccades """
#     plt.figure(1,[5,5])
#     n, bins, patches = plt.hist(all_s['amplitude'], normed=1, bins = arange(0,200,5), histtype='stepfilled')
#     plt.setp(patches, 'facecolor', '0.75', 'alpha', 0.75)
#     #plt.yscale('log', nonposy='clip')
#     #plt.title(title)
#     plt.xlabel('Amplitude (deg)')
#     plt.ylabel('Normalized frequency')
#     if save:
#         plt.savefig(name)
# 
# def hist_duration(all_s,title='',save=False,name=''):
#     """ Histogram of saccade durations given an array of saccades """
#     plt.figure(1,[5,5])
#     n, bins, patches = plt.hist(all_s['duration']*1000, normed=1, bins = arange(0,200,5), histtype='stepfilled')
#     plt.setp(patches, 'facecolor', '0.75', 'alpha', 0.75)
#     #plt.yscale('log', nonposy='clip')
#     #plt.title(title)
#     plt.xlabel('Duration (ms)')
#     plt.ylabel('Normalized frequency')
#     if save:
#         plt.savefig(name)
# 
# def hist_velocity(all_s,title='',save=False,name=''):
#     """ Histogram of saccade velocities given an array of saccades """
#     plt.figure(1,[5,5])
#     n, bins, patches = plt.hist(all_s['amplitude'], normed=1, bins = arange(0,200,5), histtype='stepfilled')
#     plt.setp(patches, 'facecolor', '0.75', 'alpha', 0.75)
#     #plt.yscale('log', nonposy='clip')
#     #plt.title(title)
#     plt.xlabel('Amplitude (deg)')
#     plt.ylabel('Normalized frequency')
#     plt.title(title)
#     
#     if save:
#         plt.savefig(name)
# 
# def get_isi(trace, min_bout_len_s_s=60):
#     """ Gets intersaccade intervals for a single trace """
#     bouts = trace.separate_bouts()
#     intervals = []
#     for i in range(len(bouts)):
#         if len(bouts[i]) > min_bout_len_s_s*100:
#             start_t, stop_t = bouts[i]['times'][0], bouts[i]['times'][-1]
#             f = np.where((trace.saccades['times'] > start_t) & (trace.saccades['times'] < stop_t))[0]
#             intervals.append([ abs(trace.saccades['times'][j+1] - trace.saccades['times'][j]) for j in range(f[0],f[-2]) ])
#     isi = np.concatenate(intervals)
#     return isi
# 
# def get_all_isi(traces,mbl=60):
#     """ Gets intersaccade intervals for all MagnOrientation in traces """
#     intervals = [0]*len(traces)
#     for i in range(len(traces)):
#         intervals[i] = get_isi(traces[i],min_bout_len_s=mbl)
#     all_isi = np.concatenate(intervals)
#     return all_isi
# 
# def hist_isi(all_isi,title='',save=False,name=''):
#     plt.figure(1,[5,5])
#     n, bins, patches = plt.hist(all_isi, normed=1, bins = arange(0,10,0.1), histtype='stepfilled')
#     plt.setp(patches, 'facecolor', '0.75', 'alpha', 0.75)
#     #plt.yscale('log', nonposy='clip')
#     #plt.xscale('log', nonposx='clip')
#     #plt.title(title)
#     plt.xlabel('intervals (s)')
#     if save:
#         plt.savefig(name)
# 
# def plot_fix_hist(trace,fix,save=False,name=''):
#     """ Orientation trace overlaid with position of a stripe, and a 2D histogram plotting the position of the stripe vs. the angle of the fly relative to the stripe """
#     plt.figure(1,[12,2])
#     ax = plt.subplot2grid((1,3),(0,0),colspan=2)
#     trace.plot(ax)
#     ax.plot(trace.flight['times'],trace.flight['stim'],'g')
#     ax1 = plt.subplot2grid((1,3),(0,2))
#     ax1.hist2d(fix['distance'],fix['stripe'],bins=[10,8],range=[[0,180],[0,360]])
#     ax1.set_ylim([0,360])
#     ax1.set_xlabel('angle from stripe (deg)')
#     ax1.axes.get_yaxis().set_ticklabels([])
#     #ax1.axes.get_xaxis().set_ticklabels([])
#     
#     if save:
#         plt.savefig(name)
#         
#     plt.show()
