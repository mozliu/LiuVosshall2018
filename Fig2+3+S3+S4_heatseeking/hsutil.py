import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy
import scipy.signal as signal
import scipy.stats as stats
from scipy.interpolate import spline
import scipy.io as sio
from statsmodels.robust.scale import mad

def getmat(path,fn):
    #dat = sio.loadmat(path+fn+'\\'+fn+'_topy.mat')
    dat = sio.loadmat(path+fn+'_topy.mat')
    Rtot = dat['Rtot'][0]
    try:
        Rtemp = np.transpose(dat['Rtemp'])[0]
        return Rtot, Rtemp
    except KeyError:
        return Rtot    

def itergetmat(path,flist,length=9180):
    n = len(flist)
    Rtot = np.empty([n,length])
    Rtemp = np.empty([n,length])
    for i in range(n):
        a = getmat(path,flist[i])
        if len(a) == 2:
            Rtot[i] = a[0]
            Rtemp[i] = a[1]
        else:
            Rtot[i] = a
    return Rtot,Rtemp

def plot_trace(t,tot,c='0.5',ax=None,label=None,indiv=True,med=False):
    if ax is None:
        ax = plt.gca()
    avg = np.nanmedian(tot,axis=0)
    if indiv:
        for trace in tot:
            ax.plot(t,trace,c=c,alpha=0.8,lw=0.5)
    if med:
        medabs = mad(tot,c=1,axis=0)
        ax.fill_between(t,avg-medabs,avg+medabs,alpha=0.4,facecolor=c,edgecolor=None)
    ax.plot(t,avg,alpha=1,c=c,label=label,lw=0.75)
    return ax

def plot_stacked_traces(Rtemp,dRtot,cRtot,c=['orange','0.5'],indiv=True,length=9180,start=540,sv=0,svstr=''):
    t = np.array(range(length-start)) / 60.

    f,ax = plt.subplots(3,sharex=True)
    f.set_size_inches(7,3)

    nd = len(dRtot)
    nc = len(cRtot)

    if indiv:
        sem = True
    else:
        sem = False

    ax[0] = plot_trace(t,Rtemp[:,start:],c='r',ax=ax[0],indiv=indiv,sem=sem)
    ax[0].set_xlim(min(t),max(t))
    ax[0].set_xticks(range(0,144,12))
    ax[0].set_ylabel('Peltier T ($^\circ$C)')

    ax[1] = plot_trace(t,dRtot[:,start:]*100,c=c[0],ax=ax[1],indiv=indiv,sem=sem)
    ax[1].set_ylim(0,50)
    ax[1].set_yticks(np.arange(0,60,10))
    ax[1].set_ylabel('% Pelt')
    #ax[1].text(136,0.52,'n='+str(nd),color='orange')

    ax[2] = plot_trace(t,cRtot[:,start:]*100,c=c[1],ax=ax[2],indiv=indiv,sem=sem)
    ax[2].set_ylim(0,50)
    ax[2].set_yticks(np.arange(0,60,20))
    ax[2].set_yticklabels(['']*4)
    ax[2].set_xlabel('time (min)')
    #ax[2].text(136,0.52,'n='+str(nc),color='0.5')

    for a in ax:
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.spines['left'].set_visible(False)
        a.spines['bottom'].set_visible(False)
        a.xaxis.grid(True)
        a.yaxis.grid(False)

    if sv:
        plt.tight_layout()
        plt.savefig(svstr)

def avgheat(cRtot,sampling=np.arange(630,9180,720),interval=90):

    n = len(sampling)
    cRpts = np.empty([len(cRtot),n])
    cRpts[:] = np.nan
    for i in range(len(cRtot)):
        cRpts[i] = [ np.mean(cRtot[i,j:j+interval]) for j in sampling ]
    cRpts = np.transpose(cRpts)
    return cRpts

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
        y = y[~np.isnan(y) & ~np.isinf(y)]
        
        # plot violins
        a = plt.violinplot(y,positions=[pos[i]],widths=0.75,showextrema=False,showmedians=True,showmeans=False)
        for bod in a['bodies']:
            bod.set(facecolor=fc[i],edgecolor=ec[i],linewidth=1,alpha=0.75)
        a['cmedians'].set(color=mc,linewidth=3)
        
        # plot points
        x = np.random.normal(pos[i],jitter,len(y))
        plt.scatter(x,y,c=pc,alpha=0.6,s=15,edgecolor='None')
    
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

def dotplot(ff,mc='k',pc='0.25',width=0.75,pos=None,ax=None,xlim=None,jitter=0.10,fixline=True):
    try:
        n,_ = ff.shape
    except ValueError:
        n = 1
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
        plt.axhline(5.42,c='0.5',ls=':')
    
    for i in range(n):
        # strip np.nan from dataset (messes up the graph)
        y = ff[i]
        y = y[~np.isnan(y)]
        m = np.median(y)

        # plot medians
        plt.plot([pos[i]-width/2,pos[i]+width/2],[m,m],c=mc[i],linewidth=2)
        
        # plot points
        x = np.random.normal(pos[i],jitter,len(y))
        plt.scatter(x,y,c=pc[i],s=20,edgecolor='None')
    
    # cosmetic
    plt.xlim(-0.5,xlim)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    #ax.set_yticks([0,.2,0.4,0.6,0.8,1])
    ax.set_xticks([])
    
    return ax

def dotplot_match(fa,fb,ax=None,peltLines=1,mc='rb',pc='rb'):
    if ax is None:
        f,ax = plt.subplots(1)
        f.set_size_inches(7,2)

    try:
        n_stim,n_trials = fa.shape
    except AttributeError:
        n_stim = len(fa)

    # draw a line between the points
    if peltLines:
        pos = np.array(list(range(0,n_stim*3,3))+list(range(1,n_stim*3,3)))
        sort_ind = np.argsort(pos)
        for i in range(n_trials):
            #ax.plot(pos,np.concatenate([posL[:,i],posR[:,i]])[sort_ind],c='0.5',alpha=0.2)
            for j in range(n_stim):
                ax.plot([pos[j],pos[j]+1],[fa[j,i],fb[j,i]],c='0.5',alpha=1,lw=0.5)

    ax = dotplot(fa,mc=mc[0],pos=range(0,n_stim*3,3),pc=pc[0],fixline=False,ax=ax)
    ax = dotplot(fb,mc=mc[1],pos=range(1,n_stim*3,3),pc=pc[1],fixline=False,ax=ax)
    ax.set_xlim(-1.5,n_stim*3)
    #ax.set_ylim(0,20)
    #ax.set_yticks(np.arange(0,21,5))
    ax.set_ylabel('mean %')
    plt.tight_layout()

    return ax

def friedman(fda):
    print(stats.friedmanchisquare(*(d for d in fda))[1])
    ntests = 0
    for i in range(fda.shape[0]):
        printstr = ''
        for j in range(i+1,fda.shape[0]):
            printstr = printstr + str(round(stats.ranksums(fda[i],fda[j])[1],4)) + '\t'
            ntests += 1
        print(printstr)
    print(0.05/ntests)
    
def multiMWU(fa,fb,alpha=0.05):
    n,_ = fa.shape
    ntests = 0
    for i in range(n):
        print(stats.mannwhitneyu(fa[i],fb[i]))
        ntests += 1
    print('p_adj = '+str(alpha/ntests))

def multiWilcoxon(fa,fb,alpha=0.05):
    n,_ = fa.shape
    ntests = 0
    for i in range(n):
        print(stats.ranksums(fa[i],fb[i]))#,zero_method='pratt'))
        ntests += 1
    print('p_adj = '+str(alpha/ntests))