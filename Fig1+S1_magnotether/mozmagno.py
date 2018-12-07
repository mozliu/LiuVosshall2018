# Mosquito magnotether scripts
# Adapted from magno.py, which was adapted from saccades.py
# Molly Liu, 2/25/16

from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colors
from scipy import stats
import scipy.signal as signal
from scipy.optimize import leastsq
import sys, os
import circular
import datetime
import pandas as pd
from traxio import readtrax
from glob import glob

# ------ basic utility functions ------ #

def cts_trace(trace_orig):
    """ Converts wraparounds to a continuous trace for use """
    trace = deepcopy(trace_orig)
    jumps = np.diff(trace)
    d_jumps = np.where(abs(jumps) > 200)[0].tolist()
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

def angles_per_s(df,smooth='100ms',add='s'):
    dfnew = pd.Series(cts_trace(df['orientation']), index = df.index)
    a = dfnew.diff().resample(smooth,how='sum')
    b = abs(a).resample(add,how='sum')
    return b

def get_velocity(trace_cts):
    # apply 8th-order Butterworth filter
    N = 8
    Wn = .2
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
    
    return d_trace, trace_f

def find_nearest_i(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
    
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

def fill_between_steps(ax, x, y1, y2=0, step_where='pre', **kwargs):
    ''' fill between a step plot and 

    Parameters
    ----------
    ax : Axes
       The axes to draw to

    x : array-like
        Array/vector of index values.

    y1 : array-like or float
        Array/vector of values to be filled under.
    y2 : array-Like or float, optional
        Array/vector or bottom values for filled area. Default is 0.

    step_where : {'pre', 'post', 'mid'}
        where the step happens, same meanings as for `step`

    **kwargs will be passed to the matplotlib fill_between() function.

    Returns
    -------
    ret : PolyCollection
       The added artist

    '''
    if step_where not in {'pre', 'post', 'mid'}:
        raise ValueError("where must be one of {{'pre', 'post', 'mid'}} "
                         "You passed in {wh}".format(wh=step_where))

    # make sure y values are up-converted to arrays 
    if np.isscalar(y1):
        y1 = np.ones_like(x) * y1

    if np.isscalar(y2):
        y2 = np.ones_like(x) * y2

    # temporary array for up-converting the values to step corners
    # 3 x 2N - 1 array 

    vertices = np.vstack((x, y1, y2))

    # this logic is lifted from lines.py
    # this should probably be centralized someplace
    if step_where == 'pre':
        steps = np.zeros((3, 2 * len(x) - 1), np.float)
        steps[0, 0::2], steps[0, 1::2] = vertices[0, :], vertices[0, :-1]
        steps[1:, 0::2], steps[1:, 1:-1:2] = vertices[1:, :], vertices[1:, 1:]

    elif step_where == 'post':
        steps = np.zeros((3, 2 * len(x) - 1), np.float)
        steps[0, ::2], steps[0, 1:-1:2] = vertices[0, :], vertices[0, 1:]
        steps[1:, 0::2], steps[1:, 1::2] = vertices[1:, :], vertices[1:, :-1]

    elif step_where == 'mid':
        steps = np.zeros((3, 2 * len(x)), np.float)
        steps[0, 1:-1:2] = 0.5 * (vertices[0, :-1] + vertices[0, 1:])
        steps[0, 2::2] = 0.5 * (vertices[0, :-1] + vertices[0, 1:])
        steps[0, 0] = vertices[0, 0]
        steps[0, -1] = vertices[0, -1]
        steps[1:, 0::2], steps[1:, 1::2] = vertices[1:, :], vertices[1:, :]
    else:
        raise RuntimeError("should never hit end of if-elif block for validated input")

    # un-pack
    xx, yy1, yy2 = steps

    # now to the plotting part:
    return ax.fill_between(xx, yy1, y2=yy2, **kwargs)

# ------ raw data into pd.Dataframe ------ #

def get_orientation(trxname):
    _,results = readtrax(trxname,True)
    orientation = pd.DataFrame({'orientation': (results['slope'][:,0])*180/np.pi},
                              index = pd.to_datetime(results['timestamp'][:,0]))
    return orientation

def get_stim(stimname):
    # xstim scaled and offset to match orientation:
    # moz facing trailing edge of dark bar
    stim = pd.read_pickle(stimname)
    stim['xstim'] = -72*stim['xstim']+340.07
    stim['xstim'][stim['xstim'] < 0] += 360
    return stim

def get_data(trxname,stimname):
    orientation = get_orientation(trxname)
    try:
        stim = get_stim(stimname)
        data = pd.concat([orientation,stim],axis=1)
    except IOError:
        data = orientation
    return data

# ------ Trial class (each mosquito is a list of Trial objects) ------ #

class Trial():
    """ trial of a magnotether run """
    def __init__(self, data, shape, name='', movingshapes=['bar','mbar','spot','blank','lbar','dbar','lspot','dspot','choice','gradient'],possible=np.array([0,0.1,0.5,1.0])):#, pretrial):
        self.name = name
        data = data.dropna()
        self.data = data
        self.data['orientation'][data['orientation'] > 360] -= 360
        self.data['orientation'][data['orientation'] < 0] += 360
        self.t = (data.index - data.index[0]).astype('timedelta64[10ms]')/100.
        self.shape = shape
        if shape in movingshapes:
            self.freq, self.guess_off = self.find_freq(possible)
            offset = data['orientation'] - self.guess_off
            offset[offset > 180] -= 360
            offset[offset < -180] += 360
            self.offset = offset
        
    def find_freq(self,possible=np.array([0,1.0])):
        xstim = cts_trace(self.data['xstim'])
        #plt.plot(xstim)
        maxpos = max(xstim[20:-10])
        minpos = min(xstim[20:-10])
        guess_off = (maxpos+minpos)/2
        if maxpos-minpos < 15:
            return 0.0, guess_off
        else:
            peaks = np.where(xstim[20:] > maxpos-10)[0]
            unique_peak = peaks[np.where(np.diff(peaks) > 1)[0]]
            troughs = np.where(xstim[20:] < minpos+10)[0]
            unique_trough = troughs[np.where(np.diff(troughs) > 1)[0]]
            guess_peak = round((len(unique_peak)+1)/self.t[-1],1)
            guess_trough = round((len(unique_trough)+1)/self.t[-1],1)
            guess_freq = (guess_peak+guess_trough)/2
            guess_freq = possible[find_nearest_i(possible,guess_freq)]
            return guess_freq, guess_off
    
    def save(self,prefix):
        savestr = prefix+self.shape
        self.data.to_pickle(savestr)

# ------ Loading and saving lists of Trials ------ #
    
def split_trials(fname,chsplit='ystim',method='diff',difft=0.5,off=173.1,tlength=1550.,manual=False,manualkey=[],shapes=['grating','bar','mbar','spot','blank']):
    # in update, ystim keys by shape.
    # 0: grating
    # 1: bar
    # 2: mbar
    # 3: spot
    # 4: blank
    
    data = get_data(fname[0],fname[1])
    #data = fname
    data = data.resample('10ms',how='mean').dropna(how='any')
    
    splitter = data[chsplit].fillna(method='pad')
    
    try:
        # splits by changes in chsplit
        if method == 'diff':
            thres = np.where(abs(np.diff(splitter)) > difft)[0] + 1
            thres = np.concatenate([thres[0:1],thres[1:][np.diff(thres) > 1]]) # removes consecutive indices to avoid arrays with length 1
        elif method == 'thres':
            thres = np.where(np.around(splitter,0)==off)[0]
            thres = thres[:-1][np.diff(thres) > tlength-50]
        data = np.split(data,thres)
        #else:
            #print 'Error: method has to be diff or thres'
    except KeyError:
        print('Error: chsplit key is not in data!')
    
    # above fails to catch consecutive presentations of shapes
    # so rough approach: since trials are same length, split trials that are too long into trial-length pieces (tlength)
    tooshort = []
    for i in range(len(data)):
        try:
            data[i] = np.array_split(data[i],int(np.around(len(data[i])/tlength,0)))
        # errors result from trials that are less than half-length
        except ValueError:
            tooshort.append(i)
            print(str(i)+' is too short')
    # delete trials shorter than half-length
    if len(tooshort) > 0:
        tooshort.reverse()
        for t in tooshort:
            del data[t]
        
    # above results in nested lists of data, so flatten list
    data = [item for sublist in data for item in sublist]
    # put in order of time
    data.sort(key=lambda d: d.index[0])
    
    # puts data into trial objects, annotated by shape
    if 'loom' in shapes:
        for i in range(len(data)):
            s = 'loom'
            if np.mean(data[i]['ystim']) < 1:
                s = 'grating'
            data[i] = Trial(data[i],s)          
    else:
        data = key_by_shape(data,manual,manualkey,shapes)
    
    return data

def key_by_shape(data,manual=False,manualkey=[],shapes=['grating','bar','mbar','spot','blank']):
    if manual:
        shapekey = manualkey
    else:
        shapekey = [ int(np.round(np.mean(d['ystim'])*len(shapes)/5.0)) for d in data ]
        
    trialdata = [0]*len(shapekey)
    for i in range(len(shapekey)):
        n = shapekey[i]
        if len(shapes) < 5:
            n -= 1
        trialdata[i] = Trial(data[i],shapes[n])
    return trialdata
    
def splitCO2(data,manual=False,manualsplit=-1):
    #data.sort(key=lambda d: d.data.index[0])
    if manual:
        split = manualsplit
    else:
        split = np.nonzero([ np.count_nonzero(np.around(t.data['CO2'],0)) for t in data ])[0][0:1]
    airco2 = np.split(data,split)
    return airco2

def save_trials(trialdata,prefix):
    for i in range(len(trialdata)):
        trialdata[i].save(prefix+'-'+str(i)+'-')
        
def load_trials(path,prefix,shapes=['bar','mbar','spot','blank','grating','loom','lbar','dbar','lspot','dspot','choice','gradient','30d','22.5d','15d','11.25d','bg']):
    loadstr = path+prefix+'*'
    files = glob(loadstr)
    trialdata = []
    for f in files:
        for s in shapes:
            if '-'+s in f:
                data = pd.read_pickle(f)
                trialdata.append(Trial(data,s,name=prefix))
    trialdata.sort(key=lambda d: d.data.index[0])
    return trialdata
    
def check_splits(data,explength=1550):
    toolong = []
    for i in np.where(np.array([ len(d) for d in data ]) > explength)[0]:
        toolong.append(i)
        data[i].plot()
    return toolong
    
    ## strategies for manually checking splits
    #np.where(np.around(testdata[toolong[6]]['xstim'],-1) == 150)
    #np.where(np.diff(testdata[toolong[8]]['xstim']) > 20)
    #plt.plot(range(len(testdata[toolong[7]]['xstim'])),testdata[toolong[7]]['xstim'])
    #plt.axvline(1524,c='r')
    #plt.axvline(3048,c='r')
    #plt.axvline(4572,c='r')
    
def resplit_trials(data,toolong,splitkey):
    for i in range(len(toolong)-1):
        if len(splitkey[i]) > 0:
            data[toolong[i]] = np.split(data[toolong[i]],splitkey[i])
            for j in range(1,len(data[toolong[i]])):
                data.append(data[toolong[i]][j])
            data[toolong[i]] = data[toolong[i]][0]
    data.sort(key=lambda d: d.index[0])
    return data
    
# ------ parsing lists of Trials ------ #

def shape_parse(shapetrial,shape='all',freq='all',prev=0):
    original = shapetrial
    if shape is not 'all':
        try:
            original = get_shape(original,shape=shape,prev=prev)
            prev = 0
        except ValueError:
            print('Warning: Shape is not in the trials given!')
    if freq is not 'all':
        try:
            original = get_freq(original,freq=freq,prev=prev)
        except ValueError:
            print('Warning: Frequency is not in the trials given!')
            
    return original

# parse by shape
def get_shape(shapetrial,shape,prev=0):
    original = []
    for i in range(len(shapetrial)):
        if shape == 'notblank' and shapetrial[i].shape not in ['bg', 'blank','grating','loom']:
            if prev == 0:
                original.extend([shapetrial[i]])
            elif (prev > 0) & (i > 0):
                newtrial = Trial(shapetrial[i].data,shapetrial[i].shape)
                prevt = shapetrial[i-1].t[-prev*100:] - shapetrial[i-1].t[-1]
                prevoff = shapetrial[i-1].data['orientation'][-prev*100:] - shapetrial[i].guess_off
                prevoff[prevoff > 180] -= 360
                prevoff[prevoff < -180] += 360
                newtrial.t = np.concatenate([prevt,newtrial.t])
                newtrial.offset = np.concatenate([prevoff,newtrial.offset])
                newtrial.data = pd.concat([shapetrial[i-1].data[-prev*100:],shapetrial[i].data])
                original.extend([newtrial])
            elif prev < 0:
                newtrial = Trial(shapetrial[i].data,shapetrial[i].shape)
                newtrial.data = newtrial.data[-prev*100:]
                newtrial.t = newtrial.t[-prev*100:]
                newtrial.offset = newtrial.offset[-prev*100:]
                original.extend([newtrial])
        if shapetrial[i].shape == shape:
            if prev == 0:
                original.extend([shapetrial[i]])
            elif (prev > 0) & (i > 0):
                newtrial = Trial(shapetrial[i].data,shapetrial[i].shape)
                prevt = shapetrial[i-1].t[-prev*100:] - shapetrial[i-1].t[-1]
                prevoff = shapetrial[i-1].data['orientation'][-prev*100:] - shapetrial[i].guess_off
                prevoff[prevoff > 180] -= 360
                prevoff[prevoff < -180] += 360
                newtrial.t = np.concatenate([prevt,newtrial.t])
                newtrial.offset = np.concatenate([prevoff,newtrial.offset])
                newtrial.data = pd.concat([shapetrial[i-1].data[-prev*100:],shapetrial[i].data])
                original.extend([newtrial])
            elif prev < 0:
                newtrial = Trial(shapetrial[i].data,shapetrial[i].shape)
                newtrial.t = newtrial.t[-prev*100:]
                newtrial.data = newtrial.data[-prev*100:]
                newtrial.offset = newtrial.offset[-prev*100:]
                original.extend([newtrial])
    return original

# frequency 
def get_freq(shapetrial,freq,dec=1,prev=0):
    original=[]
    for i in range(len(shapetrial)):
        tfreq = round(shapetrial[i].freq,dec)
        if freq == 'moving' and tfreq != 0.0:
            if prev == 0:
                original.extend([shapetrial[i]])
            elif (prev > 0) & (i > 0):
                newtrial = Trial(shapetrial[i].data,shapetrial[i].shape)
                prevt = shapetrial[i-1].t[-prev*100:] - shapetrial[i-1].t[-1]
                prevoff = shapetrial[i-1].data['orientation'][-prev*100:] - shapetrial[i].guess_off
                newtrial.t = np.concatenate([prevt,newtrial.t])
                newtrial.offset = np.concatenate([prevoff,newtrial.offset])
                original.extend([newtrial])
            elif prev < 0:
                newtrial = Trial(shapetrial[i].data,shapetrial[i].shape)
                newtrial.t = newtrial.t[-prev*100:]
                newtrial.offset = newtrial.offset[-prev*100:]
                original.extend([newtrial])
        if tfreq == freq:
            if prev == 0:
                original.extend([shapetrial[i]])
            elif (prev > 0) & (i > 0):
                newtrial = Trial(shapetrial[i].data,shapetrial[i].shape)
                prevt = shapetrial[i-1].t[-prev*100:] - shapetrial[i-1].t[-1]
                prevoff = shapetrial[i-1].data['orientation'][-prev*100:] - shapetrial[i].guess_off
                newtrial.t = np.concatenate([prevt,newtrial.t])
                newtrial.offset = np.concatenate([prevoff,newtrial.offset])
                original.extend([newtrial])
            elif prev < 0:
                newtrial = Trial(shapetrial[i].data,shapetrial[i].shape)
                newtrial.t = newtrial.t[-prev*100:]
                newtrial.offset = newtrial.offset[-prev*100:]
                original.extend([newtrial])
    return original

# ------ statistical tests ------ #

def get_means(moz,type='abs',shape='all',freq='all'):
    theta = []
    for m in moz:
        st = shape_parse(m,shape=shape,freq=freq,prev=0)
        if type == 'abs':
            moztheta = [ circular.mean(np.radians(t.data['orientation']),radian=True)[0] for t in st ]
        elif type == 'off':
            moztheta = [ circular.mean(np.radians(t.offset),radian=True)[0] for t in st ]
        else:
            print('Type must be abs or off!')
            return 0
        theta.append(circular.mean(moztheta,radian=True)[0])
    return np.array(theta)

def fix(m,shape,freq,fixrange=37.5,taft=3,tend=15,plothist=False,ax=None):
    st = shape_parse(m,shape=shape,freq=freq)
    fix = [ sum((s.offset > -1*fixrange) & (s.offset < fixrange) & (s.t > taft) & (s.t <= tend))/
           float(sum((s.t > taft) & (s.t <= tend))) for s in st ]
    if plothist:
        if ax is None:
            ax = plt.gca()
        ax.hist(fix)
    return fix
    
def kuiper(moz,type='abs',shape='all',freq='all',alpha=0):
    theta = get_means(moz,type,shape,freq)
    V,test = circular.kuiper(theta,alpha)
    return V,test
    
def rayleigh(moz,type='abs',shape='all',freq='all'):
    theta = get_means(moz,type,shape,freq)
    P = circular.rayleigh(theta)
    return P
    
# ------ plots ------ #

# raw traces
def all_raw(trials,fwidth=0,title='',save=False):
    f,ax = plt.subplots(2,sharex=True)
    if fwidth == 0:
        fwidth = len(trials)
    f.set_size_inches(fwidth,6)
    f.suptitle(title)
    trials.sort(key=lambda d: d.data.index[0])
    co = 'brykm'
    i = 0
    
    for t in trials:
        ax[0].scatter(t.data.index,t.data['xstim'],c='g',s=2,edgecolor='None')
        ax[0].scatter(t.data.index,t.data['orientation'],c=co[i%len(co)],s=2,edgecolor='None')
        ax[0].text(t.data.index[0],180,str(i))
        ax[1].scatter(t.data.index,t.data['ystim'],c='r',edgecolor='none',lw=3)
        i += 1
        
        if 'CO2' in t.data.columns.values:
            ax[1].plot(t.data.index,t.data['CO2'],c='m',lw=3)
    
    # add lines to indicate stimuli
    ystim = list(set(np.around([ np.mean(d.data['ystim']) for d in trials ],1)))
    ystim.sort()
    cy = 'kbgryk'
    for i in range(len(ystim)):
        ax[1].axhline(ystim[i],c=cy[i%len(cy)],ls='--')
    ax[0].set_xlim(trials[0].data.index[0],trials[-1].data.index[-1])
    ax[0].set_ylim(0,360)
    if save:
        plt.savefig(title+'.png')
        plt.close(f)

def plot_trial(trial,ax=None):
    if ax is None:
        ax = plt.gca()
    ax.axvline(0,color='0.5',linestyle='-')
    ax.plot(trial.t,trial.data['xstim'],c='b',s=2,edgecolor='None')
    ax.scatter(trial.t,trial.data['orientation'],c='k',s=2,edgecolor='None')
    ax.set_xlim(0,15)
    ax.set_ylim(0,360)
    ax.set_ylabel('orientation')
    ax.set_xlabel('time (s)')
    ax.set_yticklabels(['0$^\circ$','60$^\circ$','120$^\circ$','180$^\circ$','240$^\circ$','300$^\circ$','360$^\circ$'])
    ax.set_yticks([0,60,120,180,240,300,360])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
def plot_offset(trial,ax=None,c='k',alpha=1.0):
    if ax is None:
        ax = plt.gca()
    ax.axvline(0,color='0.5',linestyle='-')
    ax.axhline(-45,color='0.5',linestyle='--')
    ax.axhline(45,color='0.5',linestyle='--')
    ax.scatter(trial.t,trial.offset,c=c,alpha=alpha,s=2,edgecolor='None')
    ax.set_xlim(-3,10)
    ax.set_ylim(-180,180)
    ax.set_ylabel('offset')
    ax.set_yticklabels(['-180$^\circ$','-45$^\circ$','45$^\circ$','180$^\circ$'])
    ax.set_yticks([-180,-45,45,180])
    ax.set_xlabel('time (s)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

def checkcat(d,start=None,stop=None,freq=[0.0,0.1,0.5,1.0]):
    f,ax = plt.subplots(1,len(freq),sharex=True,sharey=True)
    f.set_size_inches(3*len(freq),2)
    for i in range(len(freq)):
        bar0 = shape_parse(d[start:stop],shape='notblank',freq=freq[i],prev=0)
        for b in bar0:
            plot_trial(b,ax=ax[i])
    ax[0].set_xlim(0,10)

def dfix(st,fixrange=60,taft=3,tend=10):
    fix = [ sum((s.offset > -1*fixrange) & (s.offset < fixrange) & (s.t > taft) & (s.t <= tend))/
               float(sum((s.t > taft) & (s.t <= tend))) for s in st ]
    blocks = [ s.block/10 for s in st ]
    pairs = range(0,blocks[-1]+1,2)
    d = np.empty(len(pairs))
    j = 0
    for i in pairs:
        try:
            d[j] = fix[blocks.index(i)] - fix[blocks.index(i+1)]
            j += 1
        except ValueError:
            d[j] = np.nan
            j += 1
    return d

def calcblocks(moz):
    for m in moz:
        t0 = m[0].data.index[0]
        tt = (np.array([ mo.data.index[0] for mo in m ]) - t0).astype('timedelta64[s]').astype('int')
        tt += 20 - ((m[0].data.index[-1] - m[0].data.index[0])+(m[1].data.index[-1] - m[1].data.index[0])).seconds # adjusting for shorter beginning periods
        blockn = tt/10
        for i in range(len(m)):
            m[i].block = blockn[i]

def calcfix(moz,shapes=['dbar','lbar','dspot','lspot','blank'],freq=[0.0,1.0],ntrial=12,fixrange=60,taft=3,tend=10):
    nrow = len(shapes)*len(freq)
    if nrow == 10:
        nrow = 9
    
    f = np.empty([len(moz),nrow,ntrial])
    f[:] = np.nan
    for m in range(len(moz)):
        for n in range(nrow):
            i = n/2
            j = n%2
            st = shape_parse(moz[m],shape=shapes[i],freq=freq[j],prev=0)
            ff = [ sum((s.offset > -1*fixrange) & (s.offset < fixrange) & (s.t > taft) & (s.t <= tend))/
           float(sum((s.t > taft) & (s.t <= tend))) for s in st ]
            blocks = [ s.block/10 for s in st ]
            f[m,n,blocks] = ff
    return f
    
def calcdfix(moz,shapes=['dbar','lbar','dspot','lspot','blank'],freq=[0.0,1.0],ntrial=6,save=False,fname=''):
    nrow = len(shapes)*len(freq)
    if nrow == 10: # adjusting for no blank
        nrow = 9
    
    d = np.empty([len(moz),nrow,ntrial])
    d[:] = np.nan
    for m in range(len(moz)):
        for n in range(nrow):
            i = n/2
            j = n%2
            st = shape_parse(moz[m],shape=shapes[i],freq=freq[j],prev=0)
            df = dfix(st)
            d[m,n,:len(df)] = df
    if save:
        np.save(fname,d)
    m = np.empty([nrow,len(moz)])
    for i in range(nrow):
        m[i] = np.nanmean(d[:,i,1:],axis=1)
    return d,m

# violinplot -- assumes that you already have fixation computed
def violinplot(ff,fc=None,ec=None,pos=None,ax=None,xlim=None,mc='k',pc='0.25',jitter=0.02,fixline=True):
    try:
    	n,_ = ff.shape
    except ValueError:
    	n = 1
    if ax is None:
    	plt.figure(1,[n,2])
    else:
    	plt.sca(ax)
    
    # default colors
    if fc is None:
    	fc = ['0.5']*n
    if ec is None:
    	ec = ['0.5']*n
    if pos is None:
    	pos = range(n)
    if xlim is None:
    	xlim = n-0.5
    
    # add guidelines
    #plt.axhline(0,c='k',ls='--')
    #plt.axhline(1,c='k',ls='--')
    if fixline:
    	plt.axhline(0.33,c='0.5',ls=':')
    
    for i in range(n):
        # strip np.nan from dataset (messes up the graph)
        y = ff[i]
        y = y[~np.isnan(y)]
        
        # plot violins
        a = plt.violinplot(y,positions=[pos[i]],widths=0.75,showextrema=False,showmedians=True,showmeans=False)
        for bod in a['bodies']:
            bod.set(facecolor=fc[i],edgecolor=ec[i],linewidth=2,alpha=0.75)
        a['cmedians'].set(color=mc,linewidth=1)
        
        # plot points
        x = np.random.normal(pos[i],jitter,len(y))
        plt.scatter(x,y,c=pc,alpha=0.25,edgecolor='None')
    
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

# violinplot -- assumes that you already have fixation computed
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

def allfol(moz,vint=300,thres=0.3,maxtrials=10):
	mfol = np.empty([len(moz),maxtrials])
	mv = np.empty([len(moz),maxtrials])
	mfol[:] = np.nan
	mv[:] = np.nan
	for i in range(len(moz)):
		bg = shape_parse(moz[i],shape='bg')
		v = np.empty(len(bg))
		fol = np.empty(len(bg))
		v[:] = np.nan
		fol[:] = np.nan
		for j in range(len(bg)):
			fol[j],v[j] = following(bg[j],vint,thres)
		mfol[i,:len(fol)] = fol
		mv[i,:len(v)] = v
	return mfol,mv

def following(trial,vint,thres,plot=False,textx=7,ax=None):
	xstim = cts_trace(trial.data['xstim'])
	ori = cts_trace(trial.data['orientation'])
	t = trial.t
	a,_,_,_,_ = stats.linregress(t,xstim) # slope of xstim
	
	v = np.empty(len(t))
	v[:vint/2] = [ (ori[i+vint/2]-ori[0])/(t[i+vint/2]-t[0]) for i in range(vint/2) ]
	v[vint/2:len(t)-vint/2] = [ (ori[i+vint/2]-ori[i-vint/2])/(t[i+vint/2]-t[i-vint/2]) for i in range(vint/2,len(t)-vint/2) ]
	v[len(t)-vint/2:] = [ (ori[-1]-ori[i-vint/2])/(t[-1]-t[i-vint/2]) for i in range(len(t)-vint/2,len(t)) ]
	
	if a >= 0:
		fol = (v > (1-thres)*a) & (v < (1+thres)*a)
	else:
		fol = (v < (1-thres)*a) & (v > (1+thres)*a)
	
	result = sum(fol)/float(len(t))
	
	if plot:
		if ax is None:
			ax = plt.gca()
		ax.scatter(t,trial.data['xstim'],c='g',s=2,edgecolor='none',label='stimulus')
		ax.scatter(t,trial.data['orientation'],c='b',s=2,edgecolor='none',label='mosquito')
		ax.scatter(t[fol],trial.data['orientation'][fol],c='r',s=15,alpha=0.1,edgecolor='none',label='following')
		
		ax.set_xlim(min(t),max(t))
		ax.set_xlabel('time (s)')
		ax.set_ylabel('orientation')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.text(textx,10,str(round(result,2)),color='r')
		ax.set_ylim(0,360)
		ax.set_yticklabels(['0$^\circ$','60$^\circ$','120$^\circ$','180$^\circ$','240$^\circ$','300$^\circ$','360$^\circ$'])
		ax.set_yticks([0,60,120,180,240,300,360])
		#plt.legend()
	
	return sum(fol)/float(len(t)),a

def co2change_raw(m,shapes=['dbar','lbar','dspot','lspot','blank'],freq=[0.0,1.0],save=False,fname=''):
    try:
        col = (max([mo.block for mo in m])+1)/2
    except AttributeError:
        t0 = m[0].data.index[0]
        tt = (np.array([ mo.data.index[0] for mo in m ]) - t0).astype('timedelta64[s]').astype('int')
        blockn = (tt+6)/100
        for i in range(len(m)):
            m[i].block = blockn[i]
        col = (max([mo.block for mo in m])+1)/2

    d = np.empty([9,col])
    d[:] = np.nan
    f,ax = plt.subplots(len(shapes)*len(freq),col,sharex=True,sharey=True)
    f.set_size_inches(1.5*col,9)
    for n in range(9):
        # order of the rows is dbar 0.0, dbar 1.0, lbar 0.0 ...
        i = n/2
        j = n%2
        st = shape_parse(m,shape=shapes[i],freq=freq[j],prev=0)
        for s in st:
            if s.block%2 == 0:
                c = 'm'
            elif s.block%2 == 1:
                c = 'k'
            plot_offset(s,ax=ax[n,s.block/2],c=c,alpha=0.75)
            ax[n,s.block/2].axhline(60,color='0.5',ls='--')
            ax[n,s.block/2].axhline(-60,color='0.5',ls='--')
        df = dfix(st)
        d[n,:len(df)] = df
            
        for k in range(len(d[n])):
            ax[n,k].text(7,-180,str(round(d[n][k],2)),color='m')
        
    ax[0,0].set_xlim(0,10)
    for i in range(len(ax)):
        for j in range(len(ax[i])):
            if (j == 0) & (i == 8):
                ax[i][j].spines['top'].set_visible(False)
                ax[i][j].spines['right'].set_visible(False)
                ax[i][j].spines['bottom'].set_visible(False)
                ax[i][j].spines['left'].set_visible(False)
            else:
                ax[i][j].set_axis_off()
    
    if save:
        plt.savefig(fname)
    
    return d

# plot circular means of trials over time for an indiv mosquito
def circmean_time(m,shape='all',freq='all',data='off',taft=2,plot=True,ax=None,c='k'):
    trials = shape_parse(m,shape,freq,prev=-1*taft)
    if data == 'off':
        circmean = [ abs(circular.mean(st.offset,radian=False)[0]) for st in trials ]
    elif data == 'abs':
        circmean = [ abs(circular.mean(st.data['orientation'],radian=False)[0]) for st in trials ]
    else:
        'Data type must be off (for offset) or abs (for absolute position)!'
        pass
    t = [ (st.data.index[0] - m[0].data.index[0])/pd.Timedelta(minutes=1) for st in trials ]
    
    if plot:
        if ax is None:
            ax = plt.gca()
        ax.plot(t,circmean,linestyle=':',marker='o',c=c,mfc=c,mec='none')
        ax.set_ylim(0,180)
        ax.set_xlim(0,26)
    return t,circmean

def circmean_multi(moz,shape='all',freq='all',taft=2,ax=None):
    if ax is None:
        ax = plt.gca()
    colormap = plt.cm.viridis
    c = [ colormap(i) for i in np.linspace(0,0.9,len(moz)) ]
    i = 0
    for m in moz:
        _,_ = circmean_time(m,shape=shape,freq=freq,taft=taft,plot=True,ax=ax,c=c[i])
        i += 1

# hist2d of offset    
def hist2d_time(fly,shape='all',freq='all',plot=True,ax=None,prev=3,stimwidth=15.0,tend=10):
    img = np.empty([len(fly),prev+tend,24])
    for k in range(len(fly)):
        bars = shape_parse(fly[k],shape,freq,prev=prev)
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

# grid of hist2d split by frequency and shape
def hist2d_grid(moz,shape=['bar','mbar','spot','blank','notblank'],freq=[0.0,0.1,0.5,1.0,'moving','all'],prev=5,save=False,fname=''):
    n = len(shape)
    m = len(freq)
    f,ax = plt.subplots(m,n,sharex=True,sharey=True)
    f.set_size_inches(2*n,1.5*m)
    for s in shape:
        for fr in freq:
            _ = hist2d_time(moz,shape=s,freq=fr,ax=ax[shape.index(s)],prev=prev)
            if shape.index(s) is n-1:
                ax[0].set_ylabel(fr)
        #ax[shape.index(s)].set_xlabel(s)
    plt.tight_layout()
    if save:
        plt.savefig(fname)
    return ax
        
def hist2d_listgrid(mozlist,shape=['bar','mbar','spot','blank'],prev=3,save=False,fname='',vmax=None,tend=10):
    n = len(shape)
    m = len(mozlist)
    if vmax is None:
        setv = True
        vmax = 0
    else:
        setv = False

    f,ax = plt.subplots(m,n,sharex=True,sharey=True)
    h = np.empty([m,n],dtype='O')
    f.set_size_inches(2*n,1.5*m)
    for s in shape:
        for i in range(len(mozlist)):
            imgtoplot,h[i][shape.index(s)] = hist2d_time(mozlist[i],shape=s,ax=ax[i,shape.index(s)],prev=prev,tend=tend)
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
    ax[0][0].set_xlim([-0.5,12.5])
    ax[0][0].set_ylim([-0.5,23.5])
    plt.colorbar(h[-1][-1])
    #plt.tight_layout()
    
    if save:
        plt.savefig(fname)

    return vmax

# residence in one zone over time
def latency(fly,shape='bar',freq='all',thres=15./72.*np.pi,tstep=0.5,plot=False,ax=None,ci=0.95,c='b',ec='b',label='',ms=6,prev=5,stimwidth=15.0):
    #f,ax = plt.subplots(2,sharex=True)
    #print 'run'
    #determine threshold
    #if thres is None:
    #    thres = circular.std(alloff)*stdthres
    tbins = np.arange(-1*prev,15+tstep,tstep)
    obins = [-np.pi,-thres,thres,np.pi]
    
    # computing the percentage of traces within the threshold for each
    perin = np.empty([len(fly),len(tbins)-1])
    for i in range(len(fly)):
        allbar = shape_parse(fly[i],shape=shape,freq=freq,prev=prev)
        if len(allbar)==0:
            perin[i].fill(np.nan)
        else:
            # residency of all trials within a fly (not finding the mean of anything first)
            t = [a.t for a in allbar]
            t = [item for sublist in t for item in sublist]
            offset = [np.radians(a.offset+stimwidth/2) for a in allbar]
            offset = [item for sublist in offset for item in sublist]
            offset[offset > np.pi] -= np.pi
            hist2d = np.histogram2d(t,offset,bins=[tbins,obins])
            # finds fixation in the center (between -thres and thres)
            oin = [ hist2d[0][j][1] for j in range(len(hist2d[0])) ]
            # finding percentage (normalizing so sum of all is one)
            perin[i] = oin/np.sum(hist2d[0],axis=1)
    # averages across all flies
    perin_mean = np.nanmean(perin,axis=0)
    perin_ci = stats.t.interval(ci,len(fly)-1)[1] * np.nanstd(perin,axis=0)/np.sqrt(len(fly))
    
    if plot:
        if ax is None:
            ax = plt.gca()
        
        ax.errorbar(tbins[:-1],perin_mean,yerr=perin_ci,fmt='o',color=c,ecolor=c,mec=ec,ms=ms,label=label)
        ax.axhline(thres/np.pi,c='0.5')
        ax.axvline(0,color='0.5')

# antifixation
def antilatency(fly,shape='bar',freq='all',thres=15./72.*np.pi,tstep=0.5,plot=False,ax=None,ci=0.95,c='b',ec='b',label='',ms=6,prev=5):
    #f,ax = plt.subplots(2,sharex=True)
    #print 'run'
    #determine threshold
    #if thres is None:
    #    thres = circular.std(alloff)*stdthres
    tbins = np.arange(-1*prev,15+tstep,tstep)
    obins = [-np.pi,-np.pi+thres,np.pi-thres,np.pi]
    
    # computing the percentage of traces within the threshold for each
    perin = np.empty([len(fly),len(tbins)-1])
    for i in range(len(fly)):
        allbar = shape_parse(fly[i],shape=shape,freq=freq,prev=prev)
        if len(allbar)==0:
            perin[i].fill(np.nan)
        else:
            # residency of all trials within a fly (not finding the mean of anything first)
            t = [a.t for a in allbar]
            t = [item for sublist in t for item in sublist]
            offset = [np.radians(a.offset) for a in allbar]
            offset = [item for sublist in offset for item in sublist]
            hist2d = np.histogram2d(t,offset,bins=[tbins,obins])
            # finds fixation in the center (between -thres and thres)
            oin = [ hist2d[0][j][0] + hist2d[0][j][2] for j in range(len(hist2d[0])) ]
            # finding percentage (normalizing so sum of all is one)
            perin[i] = oin/np.sum(hist2d[0],axis=1)
    # averages across all flies
    perin_mean = np.nanmean(perin,axis=0)
    perin_ci = stats.t.interval(ci,len(fly)-1)[1] * np.nanstd(perin,axis=0)/np.sqrt(len(fly))
    
    if plot:
        if ax is None:
            ax = plt.gca()
        ax.errorbar(tbins[:-1],perin_mean,yerr=perin_ci,fmt='o',color=c,ecolor=c,mec=ec,ms=ms,label=label)
        ax.axhline(thres/np.pi,c='0.5')
        ax.axvline(0,color='0.5')

# combined fixation and antifixation latency graph
#def latanti(fly,shape='bar',freq='all',thres=15./72.*np.pi,tstep=1.0,plot=False,ax=None,ci=0.95,c=['b','0.5'],ec='b',label='',ms=6,prev=5):
    #tbins = np.arange(-1*prev,15+tstep,tstep)
    #obins = [-np.pi,-np.pi+thres,-thres,thres,

# traditional heat map and latency graph
def hist2d_latency(mozlist,color=['0.5','k'],shape=['bar','mbar','spot','blank'],freq='all',c='bgrk',tstep=1,save=False,fname=''):
    ticks = [-np.pi,-np.pi/2.,0,np.pi/2.,np.pi]
    ticklabels = ['-180$^\circ$','-90$^\circ$','0$^\circ$','+90$^\circ$','+180$^\circ$']
    
    nrows = len(mozlist)+2
    f,ax = plt.subplots(nrows,len(shape),sharex=True)
    f.set_size_inches(3*len(shape),2*nrows)
    
    for i in range(len(shape)):
        for j in range(len(mozlist)):
            _ = hist2d_time(mozlist[j],shape=shape[i],freq=freq,plot=True,ax=ax[j,i])
            latency(mozlist[j],shape=shape[i],freq=freq,plot=True,ax=ax[-2,i],tstep=tstep,c=color[j],ec=c[i])
            antilatency(mozlist[j],shape=shape[i],freq=freq,plot=True,ax=ax[-1,i],tstep=tstep,c=color[j],ec=c[i])
        ax[-2,i].set_ylim(0,1)
        ax[-1,i].set_ylim(0,1)
    ax[j,0].set_ylabel('offset')
    ax[-2,0].set_ylabel('% fix')
    ax[-1,0].set_ylabel('% antifix')
    ax[-1,0].set_xlabel('time (s) since appearance')
    
    for i in range(len(ax)):
        for j in range(len(ax[i])):
            if j == 0:
                ax[i][j].spines['top'].set_visible(False)
                ax[i][j].spines['right'].set_visible(False)
                ax[i][j].spines['bottom'].set_visible(False)
                ax[i][j].spines['left'].set_visible(False)
            else:
                ax[i][j].set_axis_off()
            if i < len(mozlist):
                ax[i][j].set_yticks(ticks)
                ax[i][j].set_yticklabels(ticklabels)    
    plt.tight_layout()
    if save:
        plt.savefig(fname)                  

# average histogram by trial in one fly
def polhist_indiv(m,shape='all',freq='all',prev=15,taft=2,nbin=24,plot=False,ax=None,c='k',save=False,fname='',stimwidth=15.0):
    """ histogram of the average degree by shape """
    st = shape_parse(m,shape=shape,freq=freq,prev=prev)
    allbef = np.empty([len(st),nbin+1])
    allaft = np.empty([len(st),nbin+1])
    
    for i in range(len(st)):
        # separate into before and after and shift to center the histogram
        split = [find_nearest_i(st[i].t,0),find_nearest_i(st[i].t,2)]
        bef,_,aft = np.split(st[i].offset,split)
        bef = np.radians(bef)
        aft = np.radians(aft)
        bef += np.pi/nbin
        bef[bef>np.pi] -= 2*np.pi
        aft += np.pi/nbin
        aft[aft>np.pi] -= 2*np.pi
    
        histbef = np.histogram(bef,nbin,(-np.pi,np.pi))
        histaft = np.histogram(aft,nbin,(-np.pi,np.pi))
    
        allbef[i] = np.append(histbef[0],[histbef[0][0]])/float(np.sum(histbef[0]))
        allaft[i] = np.append(histaft[0],[histaft[0][0]])/float(np.sum(histaft[0]))
    
    avgbef = np.mean(allbef,axis=0)
    avgaft = np.mean(allaft,axis=0)
    
    if plot:
        stepsize = 2*np.pi/nbin
        steps = np.arange(-np.pi,np.pi+stepsize,stepsize)
        steps += np.pi/2+np.radians(stimwidth/2)
        
        if ax is None:
            ax = plt.gca()
        ax.spines['polar'].set_visible(False)
        ax.set_xticks([np.pi/2,3*np.pi/2])
        ax.set_xticklabels(['0$^\circ$','$\pm$180$^\circ$'])
        
        ax.plot(steps,avgaft,c=c,lw=1,drawstyle='steps-post')
        
        if save:
            plt.savefig(fname)
        
    return avgbef,avgaft

# average histogram of trial over multiple flies showing individual trials
def polhist_indivall(moz,shape='all',freq='all',taft=2,nbin=24,plot=True,ax=None,col='k',save=False,fname='',stimwidth=15.0):
    colormap = plt.cm.gist_rainbow
    c = [ colormap(i) for i in np.linspace(0,0.9,len(moz)) ]
    all = np.empty([len(moz),nbin+1])
    for i in range(len(moz)):
        _,all[i] = polhist_indiv(moz[i],shape,freq,prev=0,taft=taft,nbin=nbin,plot=True,ax=ax,c=c[i])
    avg = np.nanmean(all,axis=0)
    
    if plot:
        if ax is None:
            ax = plt.gca()
        ax.set_xticks([np.pi/2,3*np.pi/2])
        ax.set_xticklabels(['0$^\circ$','$\pm$180$^\circ$'])
        ax.spines['polar'].set_visible(False)
        stepsize = 2*np.pi/nbin
        steps = np.arange(-np.pi,np.pi+stepsize,stepsize)
        steps += np.pi/2+np.radians(stimwidth/2)
        ax.plot(steps,avg,lw=2,c=col,drawstyle='steps-post')
        
        if save:
            plt.savefig(fname)

# average histogram of trial over multiple flies comparing trials to previous trial
def polhist(moz,shape='all',freq='all',prev=13,taft=2,nbin=24,plot=False,rmax=0,ci=True,ax=None,c='k',save=False,fname=''):
    
    allbef = np.empty([len(moz),nbin+1])
    allaft = np.empty([len(moz),nbin+1])
    for i in range(len(moz)):
        allbef[i],allaft[i] = polhist_indiv(moz[i],shape,freq,prev,taft,nbin)
    avgbef = np.nanmean(allbef,axis=0)
    avgaft = np.nanmean(allaft,axis=0)
    
    stepsize = 2*np.pi/nbin
    steps = np.arange(-np.pi,np.pi+stepsize,stepsize)
    if plot:
        if ax is None:
            ax = plt.gca()
        ax.set_xticks([np.pi/2,3*np.pi/2])
        ax.set_xticklabels(['0$^\circ$','$\pm$180$^\circ$'])
        ax.spines['polar'].set_visible(False)
        
        lw = 1
        if ci:
            cibef = stats.t.interval(0.95,len(moz)-1)[1]*np.nanstd(allbef,axis=0)/np.sqrt(len(moz))
            ciaft = stats.t.interval(0.95,len(moz)-1)[1]*np.nanstd(allaft,axis=0)/np.sqrt(len(moz))
            fill_between_steps(ax,steps+np.radians(7.5)+np.pi/2,avgbef-cibef,avgbef+cibef,step_where='post',alpha=0.5,edgecolor='None',facecolor='0.5')
            fill_between_steps(ax,steps+np.radians(7.5)+np.pi/2,avgaft-ciaft,avgaft+ciaft,step_where='post',alpha=0.5,edgecolor='None',facecolor=c)
            lw = 2
            if rmax == 0:
                rmax = max(np.amax(avgaft+ciaft),np.amax(avgbef+cibef))
        else:
            if rmax == 0:
                rmax = max(np.amax(avgaft),np.amax(avgbef))
        
        ax.plot([np.radians(-7.5)+np.pi/2,np.radians(-7.5)+np.pi/2],[0,rmax],ls='--',c='k')
        ax.plot([np.radians(7.5)+np.pi/2,np.radians(7.5)+np.pi/2],[0,rmax],ls='--',c='k')
        ax.plot([np.radians(-37.5)+np.pi/2,np.radians(-37.5)+np.pi/2],[0,rmax],ls='--',c='0.5')
        ax.plot([np.radians(37.5)+np.pi/2,np.radians(37.5)+np.pi/2],[0,rmax],ls='--',c='0.5')
        ax.plot(steps+np.radians(7.5)+np.pi/2,avgbef,c='0.5',lw=lw,drawstyle='steps-post') # shift such that 0 is facing the very center of the shape
        ax.plot(steps+np.radians(7.5)+np.pi/2,avgaft,c=c,lw=lw,drawstyle='steps-post')
        ax.set_ylim(0,rmax)
        ax.set_yticks([rmax])
        ax.set_yticklabels([str(rmax)[:4]])
        
        if save:
            plt.savefig(fname)
        
    return avgbef, avgaft

def polhist_comp(mozlist,shape='all',freq='all',color='km',ec='None',taft=2,nbin=24,stimwidth=15.0,rmax=0.5,ci=0.95,ax=None,save=False,fname=''):
    if ax is None:
        ax = plt.gca()
    
    stepsize = 2*np.pi/nbin
    steps = np.arange(-np.pi,np.pi+stepsize,stepsize)
    steps += np.pi/2+np.radians(stimwidth/2)
    
    # plot average and ci for each moz in mozlist
    c = 0
    for moz in mozlist:
        all = np.empty([len(moz),nbin+1])
        for i in range(len(moz)):
            _,all[i] = polhist_indiv(moz[i],shape,freq,0,taft,nbin)
        avg = np.nanmean(all,axis=0)
        cip = stats.t.interval(ci,len(moz)-1)[1]*np.nanstd(all,axis=0)/np.sqrt(len(moz))
        fill_between_steps(ax,steps,avg-cip,avg+cip,step_where='post',alpha=0.5,edgecolor=ec,facecolor=color[c])
        ax.plot(steps,avg,lw=2,c=color[c],drawstyle='steps-post')
        c += 1
    
    # making the plot pretty
    ax.spines['polar'].set_visible(False)
    ax.set_xticks([np.pi/2,3*np.pi/2])
    ax.set_xticklabels(['0$^\circ$','$\pm$180$^\circ$'])
    # lines indicating stimulus width and movement
    ax.plot([np.radians(-stimwidth/2)+np.pi/2,np.radians(-stimwidth/2)+np.pi/2],[0,rmax],ls='--',c='k')
    ax.plot([np.radians(stimwidth/2)+np.pi/2,np.radians(stimwidth/2)+np.pi/2],[0,rmax],ls='--',c='k')
    ax.plot([np.radians(-stimwidth/2-30)+np.pi/2,np.radians(-stimwidth/2-30)+np.pi/2],[0,rmax],ls='--',c='0.5')
    ax.plot([np.radians(stimwidth/2+30)+np.pi/2,np.radians(stimwidth/2+30)+np.pi/2],[0,rmax],ls='--',c='0.5')
    ax.set_ylim(0,rmax)
    ax.set_ylim(0,rmax)
    ax.set_yticks([rmax])
    ax.set_yticklabels([str(rmax)[:4]])
    
    if save:
        plt.savefig(fname)  

def polhist_compgrid(mozlist,shape=['bar','mbar','spot','blank'],freq=['all',0.0,0.1,0.5,1.0],c='km',ec='bgrk',taft=2,nbin=24,stimwidth=15.0,ci=0.95,rmax=0.15,rticks=[0.04,0.08,0.12],save=False,fname=''):
    f,ax = plt.subplots(len(freq),len(shape),subplot_kw=dict(projection='polar'))
    f.set_size_inches(4*len(shape),4*len(freq))
    for i in range(len(freq)):
        for j in range(len(shape)):
            polhist_comp(mozlist,shape=shape[j],freq=freq[i],color=c,ec=ec[j],taft=taft,nbin=nbin,stimwidth=stimwidth,rmax=rmax,ci=ci,ax=ax[i,j])
            ax[i,j].set_yticks(rticks)
    
    if save:
        plt.savefig(fname)

# average histogram of absolute position by trial in one mosquito
def abshist_indiv(m,shape='all',freq='all',nbin=24,plot=False,ax=None,pltpos=0.22,c='0.5'):
    """ histogram of the average degree by shape """
    st = shape_parse(m,shape=shape,freq=freq,prev=0)
    all = np.empty([len(st),nbin+1])
    xpos = []
    
    for i in range(len(st)):
        # separate into before and after and shift to center the histogram
        pos = np.radians(st[i].data['orientation'])
        pos += np.pi/nbin
        pos[pos>np.pi] -= 2*np.pi
    
        hist = np.histogram(pos,nbin,(-np.pi,np.pi))
    
        all[i] = np.append(hist[0],[hist[0][0]])/float(np.sum(hist[0]))
        try:
            xpos.append(st[i].guess_off)
        except AttributeError:
            pass
    
    avg = np.nanmean(all,axis=0)
    
    if plot:
        if ax==None:
            ax = plt.gca()
        stepsize = 2*np.pi/nbin
        steps = np.arange(-np.pi,np.pi+stepsize,stepsize)
        steps += np.pi/2
        ax.plot(steps,avg,c=c,drawstyle='steps-post',lw=1)
        ax.scatter(xpos,[pltpos]*len(xpos),c=c,alpha=0.15,edgecolors='none')
    
    return avg, xpos
    
def abshist(moz,shape='all',freq='all',color='k',ec='None',nbin=24,ax=None,save=False,fname=''):
    if ax is None:
        ax = plt.gca()
    
    stepsize = 2*np.pi/nbin
    steps = np.arange(-np.pi,np.pi+stepsize,stepsize)
    steps += np.pi/2
    
    # plot average and ci for each moz in mozlist
    all = np.empty([len(moz),nbin+1])
    #xpos = np.empty(0)
    
    colormap = plt.cm.gist_rainbow
    c = [ colormap(i) for i in np.linspace(0,0.9,len(moz)) ]
    for i in range(len(moz)):
        all[i],x = abshist_indiv(moz[i],shape,freq,nbin,plot=True,ax=ax,c=c[i],pltpos=(0.2+0.005*i))
        #xpos = np.concatenate([xpos,np.radians(x)])
        
    #xhist = np.histogram(xpos,nbin,(0,2*np.pi))[0]
    #xhist = np.append(xhist,xhist[0])/(float(np.sum(xhist))*2)
    #ax.plot(steps,xhist,lw=2,c='g',drawstyle='steps-post')
    
    avg = np.nanmean(all,axis=0)
    ax.plot(steps,avg,lw=2,c=color,drawstyle='steps-post')
    
    # making the plot pretty
    ax.set_xticklabels(['270$^\circ$','315$^\circ$','0$^\circ$',
                    '45$^\circ$','90$^\circ$','135$^\circ$',
                    '180$^\circ$','225$^\circ$'])
    ax.set_ylim(0,0.2+0.005*i)
    ax.spines['polar'].set_visible(False)
    #ax.set_yticks([])
    
    if save:
        plt.savefig(fname)
