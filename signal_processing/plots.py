import numpy as np
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot2D_shot(data,pclip=95,cmap='gray',aspect='auto',size=(12,6),out=False,name='data',labels=['Data','Traces','Time, sec.']):
    """
    Plot 2D shot gather using plt.imshow
    
    data      seismic data   [traces,times]
    pclip     PercentileRank [vmin,vmax]
    labels[title,xlabel,ylabel]
    """
    plt.rcParams.update({'font.size': 14})
    pclip = np.percentile(data,pclip)
    fig = plt.figure(figsize=size)
    plt.imshow(data.T,cmap=cmap,aspect=aspect,vmin=-pclip,vmax=pclip)
    plt.title(labels[0])
    plt.xlabel(labels[1])
    plt.ylabel(labels[2])
    plt.colorbar()
    if out:
        fig.savefig(name + '.png',dpi=300,transparent=True)

##2 Picks
def plot_picks(picks,x,xlabel='Offsets [km]',ylabel='Time [sec]',
        out=False,name='QC_picks',colors=['r','b--','b'],
        figsize=(10,10),
        labels=['observed','initial','updated'],ylim=(0,20),xlim=(0,0)):
    """
    plots first arrivals
    picks[observed,initials..,final picks]
    """

    plt.rc('lines', linewidth=2)
    fig,ax = plt.subplots(1,figsize=figsize)
    
    for i in range(len(picks)):
        ax.plot(x,picks[i],colors[i],label=labels[i])
    
    ax.set(ylim=ylim)
    if xlim[1] > 0:
        ax.set(xlim=xlim)
    ax.invert_yaxis()
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlabel(xlabel,fontsize=14)
    ax.set_ylabel(ylabel,fontsize=14)
    ax.legend(fontsize=14,frameon=False)
    if out:
        fig.savefig(name + '.png',dpi=300,transparent=True)
    plt.tight_layout()
##3
def plot_model(data,
               figsize = (14,9),
               labels  = ['Horizontal distance, km','Depth, km','title'],
               extent  =[0,500,20,0],
               out=None,cmap='jet',v=[1450,8000],
               levels=np.array([1500,2000,3000,4000,5000,6000,7000,8000]),
               aspect='auto'):
    """
    plot models
    data[time,traces]

    """
    np.set_printoptions(precision=1)
    
    fig,ax = plt.subplots(figsize=figsize)
    if levels is not None:
        levels = levels
        CS = ax.contour(data.T, levels, colors='k', origin='upper', extent=extent,alpha=0.2)
        ax.clabel(CS, inline=True, fontsize=10,fmt='%1.0f')
    
    im = ax.imshow(data.T,extent=extent,aspect=aspect,cmap=cmap,vmin=v[0],vmax=v[1])#'turbo')#'jet')#'coolwarm')

    cbar = plt.colorbar(im,ax=ax, orientation='horizontal',pad=0.1)
    cbar.set_label('Velocity m/s',fontsize=14)
    
    im.axes.set_xlabel(labels[0],fontsize=14)
    im.axes.set_ylabel(labels[1],fontsize=14)
    if out is not None:
        fig.savefig(out + '.png',dpi=500,transparent=True)
    plt.show()
###
def ModelVelPlot(model,velcoor=[],sources=None,velticks=[1000,3000,5000],
                vellimit=[700,6500],v = [800,6000],cmap='jet',
                xticks=np.arange(0,200,25),
                out = None,title='Velocit model'):

    fig, axs = plt.subplots(ncols = 2, nrows = 1,sharey=True,
                        gridspec_kw = {'width_ratios': [3, 1]},
                       figsize=(10,4))
    plt.tight_layout()
    #  plot model
    ax = axs[0]
    im = ax.imshow(model,aspect='auto',vmin=v[0],vmax=v[1],cmap=cmap)
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right", size="15%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.plot(velcoor[0],velcoor[1],'r')
    # sources
    if sources is not None:
        sx = sources[0]
        sz = sources[0]
        ax.plot(sx,sz,'r.',markersize=8,label='sources')
        ax.legend(loc='lower left',fontsize=14)
    im.axes.set_xlabel('Horizontal distance [m]');im.axes.set_ylabel('Depth [m]')
    ax.set_xticks(xticks)
    # velocity model plotting aside
    axs[1].plot(velcoor[-1],velcoor[1],'r')
    axs[1].grid(axis='x')
    axs[1].set_xticks(velticks)
    axs[1].set_xlim(vellimit)
    axs[1].set_xlabel('Velocity [m/s]')
    plt.subplots_adjust(wspace=0)
    if out is not None:
        fig.savefig(out + '.jpg',dpi=600,bbox_inches='tight',format='jpeg')
###
def PickDiffPlot(data=[],ylim1=[0,0.08],ylim2=[-0.001,0.016],
                yticks1 = [0.00,0.02,0.04,0.06],yticks2= [-0.001,0.005,0.010,0.015],
                label=['true','initial'],out=None,
                axislabels=['Time [sec]','Horizontal distance [m]','Time [sec]']):

    fig, axs = plt.subplots(2, 1,figsize=(8,6),
                            gridspec_kw={'height_ratios': [3, 1]},
                           sharex=True)

    axs[0].plot(offsets,p0,'r-',label=label[0])
    axs[0].plot(offsets,p1,'b-',label=label[1])
    axs[0].set_ylim(ylim1)
    axs[0].set_yticks(yticks1)
    axs[0].legend()
    axs[0].axes.invert_yaxis()
    axs[0].set_ylabel(axislabels[0])

    axs[1].plot(offsets,p1-p0,'.r')
    #axs[1].set_yticks(np.round(np.linspace(min(p1-p0),max(p1-p0),4),3))
    axs[1].set_yticks(yticks2)
    axs[1].set_ylim(ylim2)
    axs[1].set_xlabel(axislabels[1])
    axs[1].set_ylabel(axislabels[2])

    axs[1].grid(axis='both')
    axs[0].grid(axis='both')

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()

    if out is not None:
        fig.savefig(out + '.jpg',dpi=600,bbox_inches='tight',format='jpeg')

