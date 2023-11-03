#================================PlotFuncs.py==================================#
# Created by Ciaran O'Hare 2020

# Description:
# This file has many functions which are used throughout the project, but are
# all focused around the bullshit that goes into making the plots

#==============================================================================#

from numpy import *
from numpy.random import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib import colors
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
#from scipy.stats import norm
import matplotlib.patheffects as pe

pltdir_png = '../plots/'
pltdir = pltdir_png+'pdfs/'


def UpperFrequencyAxis(ax,N_Hz=1,tickdir='out',xtick_rotation=0,labelsize=25,xlabel=r"$\nu_a$ [Hz]",lfs=40,tick_pad=8,tfs=25,xlabel_pad=10):
    m_min,m_max = ax.get_xlim()
    ax2 = ax.twiny()
    ax2.set_xlim([m_min*241.8*1e12/N_Hz,m_max*241.8*1e12/N_Hz])
    ax2.set_xlabel(xlabel,fontsize=lfs,labelpad=xlabel_pad)
    ax2.set_xscale('log')
    plt.xticks(rotation=xtick_rotation)
    ax2.tick_params(labelsize=tfs)
    ax2.tick_params(which='major',direction=tickdir,width=2.5,length=13,pad=tick_pad)
    ax2.tick_params(which='minor',direction=tickdir,width=1,length=10)
    locmaj = mpl.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=50)
    locmin = mpl.ticker.LogLocator(base=10.0, subs=arange(2, 10)*.1,numticks=100)
    ax2.xaxis.set_major_locator(locmaj)
    ax2.xaxis.set_minor_locator(locmin)
    ax2.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    plt.sca(ax)




#==============================================================================#


class DarkPhoton():
    def FigSetup(xlab=r'Dark photon mass [eV]',ylab='Kinetic mixing',\
             chi_min = 1.0e-18,chi_max = 1.0e0,\
             m_min = 3e-18,m_max = 7e5,\
             lw=2.5,lfs=40,tfs=25,tickdir='out',\
             Grid=False,Shape='Rectangular',mathpazo=True,\
             TopAndRightTicks=False,FrequencyAxis=False,FrequencyLabels=False,UnitAxis=True,f_rescale=1,\
            tick_rotation = 20,width=20,height=10,upper_tickdir='out', reduced = False):
        
        if reduced:
            chi_min = 1.0e-13
            chi_max = 1.0e-7
            m_min = 1e-3
            m_max = 1e0

        plt.rcParams['axes.linewidth'] = lw
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=tfs)

        if mathpazo:
            plt.rcParams.update({"text.usetex": True,"font.family": "serif","font.serif": ["Palatino"],})


        if Shape=='Wide':
            fig = plt.figure(figsize=(16.5,5))
        elif Shape=='Rectangular':
            fig = plt.figure(figsize=(16.5,11))
        elif Shape=='Custom':
            fig = plt.figure(figsize=(width,height))

        ax = fig.add_subplot(111)

        ax.set_xlabel(xlab,fontsize=lfs)
        ax.set_ylabel(ylab,fontsize=lfs)

        ax.tick_params(which='major',direction=tickdir,width=2.5,length=13,right=TopAndRightTicks,top=TopAndRightTicks,pad=7)
        ax.tick_params(which='minor',direction=tickdir,width=1,length=10,right=TopAndRightTicks,top=TopAndRightTicks)


        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim([m_min,m_max])
        ax.set_ylim([chi_min,chi_max])

        locmaj = mpl.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=50)
        locmin = mpl.ticker.LogLocator(base=10.0, subs=arange(2, 10)*.1,numticks=100)
        ax.xaxis.set_major_locator(locmaj)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

        locmaj = mpl.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
        locmin = mpl.ticker.LogLocator(base=10.0, subs=arange(2, 10)*.1,numticks=100)
        ax.yaxis.set_major_locator(locmaj)
        ax.yaxis.set_minor_locator(locmin)
        ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

        if Shape=='Rectangular':
            plt.xticks(rotation=tick_rotation)

        if Grid:
            ax.grid(zorder=0)

        if FrequencyAxis:
            ax2 = ax.twiny()



            ax2.set_xscale('log')
            ax2.tick_params(which='major',direction=upper_tickdir,width=2.5,length=13,pad=7)
            ax2.tick_params(which='minor',direction=upper_tickdir,width=1,length=10)
            locmaj = mpl.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=50)
            locmin = mpl.ticker.LogLocator(base=10.0, subs=arange(2, 10)*.1,numticks=100)
            ax2.xaxis.set_major_locator(locmaj)
            ax2.xaxis.set_minor_locator(locmin)
            ax2.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

            if FrequencyLabels:
                ax2.set_xticks([1e-3,1e0,1e3,1e6,1e9,1e12,1*241.8*1e12,1000*241.8*1e12])
                ax2.set_xticklabels(['mHz','Hz','kHz','MHz','GHz','THz','eV','keV'])
            ax2.set_xlim([m_min*241.8*1e12/f_rescale,m_max*241.8*1e12/f_rescale])

            plt.sca(ax)
        return fig,ax


    def Haloscopes(ax,fs=17,projection=True,text_on=True,col='darkred'):
        y2 = ax.get_ylim()[1]
        zo = 0.3

        # ADMX
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/ADMX.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1,lw=3)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/ADMX2018.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/ADMX2019_1.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/ADMX2019_2.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/ADMX2021.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/ADMX_Sidecar.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)

        # HAYSTAC
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/HAYSTAC.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/HAYSTAC_2020.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/HAYSTAC_2022.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)

        # CAPP
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/CAPP-1.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/CAPP-2.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/CAPP-3.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/CAPP-4.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/CAPP-5.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/CAPP-6.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)
        dat = loadtxt("limit_data/DarkPhoton/Rescaled/CAST-CAPP.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)

        dat = loadtxt("limit_data/DarkPhoton/Rescaled/TASEH.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.1)

        dat = loadtxt("limit_data/DarkPhoton/Rescaled/ORGAN-1a.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,facecolor=col,zorder=0.21)

        dat = loadtxt("limit_data/DarkPhoton/Rescaled/QUAX.txt")
        dat[0,1] = 1e0
        plt.plot(dat[:,0],dat[:,1],zorder=0.2,color=col,lw=2)

        dat = loadtxt("limit_data/DarkPhoton/Rescaled/QUAX2.txt")
        dat[0,1] = 1e0
        plt.plot(dat[:,0],dat[:,1],zorder=0.2,color=col,lw=2)

        dat = loadtxt("limit_data/DarkPhoton/Rescaled/QUAX3.txt")
        dat[0,1] = 1e0
        plt.plot(dat[:,0],dat[:,1],zorder=0.2,color=col,lw=2)

        if text_on:
            plt.text(1.4e-6,0.5e-14,r'{\bf ADMX}',fontsize=fs,color=col,rotation=90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.text(0.8e-5,0.1e-13,r'{\bf CAPP}',fontsize=fs-2,color=col,rotation=90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.text(0.19e-4,3e-15,r'{\bf HAYSTAC}',fontsize=fs-5,color=col,rotation=90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.text(0.47e-4,3e-12,r'{\bf QUAX}',fontsize=fs-8,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center',clip_on=True)

        return


    def StellarBounds(ax,fs=19,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        # Stellar physics constraints

        # Globular clusters
        HB_col = [0.01, 0.75, 0.24]
        HB = loadtxt("limit_data/DarkPhoton/RG.txt")
        plt.fill_between(HB[:,0],HB[:,1],y2=y2,edgecolor=None,facecolor=HB_col,zorder=0.9)
        plt.plot(HB[:,0],HB[:,1],color='k',alpha=1,zorder=0.9,lw=lw)

        # Globular clusters
        HB_col = 'DarkGreen'
        HB = loadtxt("limit_data/DarkPhoton/HB.txt")
        plt.fill_between(HB[:,0],HB[:,1],y2=y2,edgecolor=None,facecolor=HB_col,zorder=0.95)
        plt.plot(HB[:,0],HB[:,1],color='k',alpha=1,zorder=0.95,lw=lw)

        # Solar bound
        Solar_col = 'ForestGreen'
        Solar = loadtxt("limit_data/DarkPhoton/Solar.txt")
        plt.fill_between(Solar[:,0],Solar[:,1],y2=y2,edgecolor=None,facecolor=Solar_col,zorder=1.02,label="Solar")
        plt.plot(Solar[:,0],Solar[:,1],color='k',alpha=1,zorder=1.02,lw=lw)

        Solar = loadtxt("limit_data/DarkPhoton/Solar-Global.txt")
        plt.fill_between(Solar[:,0],Solar[:,1]/Solar[:,0],y2=y2,edgecolor=None,facecolor=Solar_col,zorder=1.021)
        plt.plot(Solar[:,0],Solar[:,1]/Solar[:,0],color='k',alpha=1,zorder=1.021,lw=lw)

        if text_on:
            plt.text(0.8e2,1.5e-14,r'{\bf Solar}',fontsize=fs,color='w',rotation=-41,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(1e3,0.7e-14,r'{\bf HB}',fontsize=fs,color='w',rotation=-38,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(0.7e4,0.68e-14,r'{\bf RG}',fontsize=fs,color='w',rotation=-37,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
        return


    def Xenon(ax,col='crimson',fs=23,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Xenon1T.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)

        plt.fill_between(1e3*dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)
        plt.plot(1e3*dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.5,lw=lw)

        dat = loadtxt("limit_data/DarkPhoton/Xenon1T_S1S2.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)

        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.5,lw=lw)


        dat = loadtxt("limit_data/DarkPhoton/XENON1T_SE.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.5,lw=lw)


        dat = loadtxt("limit_data/DarkPhoton/XENON1T_Solar_SE.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.302,label="XENON1T")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.302,lw=lw)


        dat = loadtxt("limit_data/DarkPhoton/XENONnT.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.5,lw=lw)

        if text_on:
            plt.text(8e2,2.5e-17,r'{\bf XENON}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            #plt.text(0.65e-3,2.4e-11,r'{\bf XENON1T}',color='w',rotation=-41,fontsize=15,path_effects=line_background(1,'k'),clip_on=True)


        return





    def DAMIC(ax,col='salmon',fs=21,text_on=True,lw=1.5):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/DAMIC.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)

        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]
        dat[-1,1] = y2[-1]
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.001)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.001,lw=lw)

        if text_on:
            plt.text(4e-1,1.3e-14,r'{\bf DAMIC}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.plot([4e0,1e1],[3e-14,6e-14],'-',lw=2.5,color=col,path_effects=line_background(3.5,'k'))
        return

    def MuDHI(ax,col='#400927',fs=15,text_on=True,lw=1.5):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/MuDHI.txt")

        y2 = interp(dat[:,0],m1,y1)
        dat[dat[:,1]>y2,1] = y2[dat[:,1]>y2]
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=10)
        plt.plot(dat[::2,0],dat[::2,1],color='k',alpha=1,zorder=10,lw=lw)

        if text_on:
            plt.text(0.18e-2,1e-11,r'{\bf MuDHI}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.plot([1.2e-2,1.5e0],[1e-11,4e-11],'-',lw=2.5,color=col,path_effects=line_background(3.5,'k'))
        return


    def FUNK(ax,col='red',fs=21,text_on=True,lw=1.5):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/FUNK.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(2/3/0.27)

        y2 = interp(dat[:,0],m1,y1)
        #dat[0,1] = y2[0]
        dat[dat[:,1]>y2,1] = y2[dat[:,1]>y2]
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.3)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.3,lw=lw)

        if text_on:
            plt.text(1.9e-1,0.8e-13,r'{\bf FUNK}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.plot([7e-1,3e0],[2.5e-13,1e-12],'-',lw=2.5,color=col,path_effects=line_background(3.5,'k'))
        return

    def SENSEI(ax,col='firebrick',fs=21,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SENSEI.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)

        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1,lw=lw)

        if text_on:
            plt.text(1.7e0,1e-15,r'{\bf SENSEI}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.plot([7e0,1e1],[3e-15,9e-15],'-',lw=2.5,color=col,path_effects=line_background(3.5,'k'))
        return

    def SuperCDMS(ax,col=[0.4,0,0],fs=18,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SuperCDMS.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.6)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.6,lw=lw)

        if text_on:
            plt.text(0.5e1,1.5e-16,r'{\bf SuperCDMS}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.plot([5e1,0.8e2],[3e-16,9e-16],'-',lw=2.5,color=col,path_effects=line_background(3.5,'k'))
        return

    def Nanowire(ax,col='pink',fs=22,text_on=True,lw=1.5):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/WSi_Nanowire.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]/1.1
        dat[-1,1] = y2[-1]/1.1
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.3)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.3,lw=lw)

        if text_on:
            plt.text(5e-4,1e-10,r'{\bf WSi Nanowire}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.plot([9e-3,3e-3],[3e-10,9e-10],'-',lw=2.5,color=col)
        return


    def SQMS(ax,col='#02734b',fs=17,text_on=True,lw=0.5,ms=10):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SQMS.txt")
        dat[:,1] = dat[:,1]*sqrt(1/3/0.019)
        plt.plot(dat[:,0],dat[:,1],lw=lw,color=col,alpha=1,zorder=0.0)
        if text_on:
            plt.text(5.7e-6,0.65e-14,r'{\bf SQMS}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
        return

    def LAMPOST(ax,col='#471710',fs=15,text_on=True,lw=1.5):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/LAMPOST.txt")
        dat[:,1] = dat[:,1]*sqrt(0.4/0.45)*sqrt(2/3/0.27)

        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]/1.1
        dat[-1,1] = y2[-1]/1.1
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0)
        plt.plot(dat[:,0],dat[:,1],color=col,alpha=1,zorder=0,lw=lw)

        if text_on:
            plt.text(0.3e-1,4.5e-13,r'{\bf LAMPOST}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.plot([3e-1,0.6e0],[6e-13,1e-12],'-',lw=1.5,color=col,path_effects=line_background(2,'k'))
        return

    def Tokyo(ax,col='darkred',fs=15,text_on=True,lw=1.5):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/Tokyo-Dish.txt")
        dat[:,1] = dat[:,1]*sqrt(2/3/0.6)
        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]
        dat[-1,1] = y2[-1]
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.4)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.4,lw=lw)


        dat = loadtxt("limit_data/DarkPhoton/Tokyo-Knirck.txt")
        dat[:,1] = dat[:,1]*sqrt(1/3/0.175)
        plt.fill_between(dat[:,0],dat[:,1],y2=1e0,edgecolor='k',facecolor=col,zorder=1.09)

        dat = loadtxt("limit_data/DarkPhoton/Tokyo-Tomita.txt")
        plt.plot([dat[1,0],dat[1,0]],[dat[1,1],1e0],'-',color=col,lw=3,zorder=0.2)
        if text_on:
            #plt.text(2e-4,1e-10,r'{\bf Tokyo-3}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.45e-3,3e-8,r'{\bf Tokyo-2}',fontsize=fs-2,color='k',rotation=90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.text(0.03e-1,2e-12,r'{\bf Tokyo-1}',fontsize=fs+4,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.plot([0.3e-1,3e0],[2e-12,8e-12],'-',lw=2.5,color=col,path_effects=line_background(3.5,'k'))
        return

    def FAST(ax,col='tomato',fs=10,text_on=True,lw=1.5,edge_on=False,zorder=0.11):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/FAST.txt")
        dat[:,1] = dat[:,1]*sqrt(2/3/0.6)
        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]/1.1
        dat[-1,1] = y2[-1]/1.1
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=zorder)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=zorder,lw=lw)

        if text_on:
            plt.text(7e-6,4e-12,r'{\bf FAST}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
        return

    def LOFAR(ax,col='red',fs=10,text_on=True,lw=1.5,edge_on=False,zorder=0.11):
        # Solar corona bound
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/LOFAR.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]/1.1
        dat[-1,1] = y2[-1]/1.1
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=zorder)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=zorder,lw=lw)

        if text_on:
            plt.text(1.95e-7,3e-14,r'{\bf LOFAR (Sun)}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
        return

    def Jupiter(ax,col='Green',fs=15,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Jupiter.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=2)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=2,lw=lw)
        if text_on:
            plt.text(0.1e-14,4.5e-1,r'{\bf Jupiter}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1,'k'),clip_on=True)
        return

    def Earth(ax,col='DarkGreen',fs=17,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Earth.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.9)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.9,lw=lw)
        if text_on:
            plt.text(0.4e-13,2e-1,r'{\bf Earth}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1,'k'),clip_on=True)
        return


    def Crab(ax,col=[0.1,0.4,0.1],fs=17,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Crab.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=2)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=2,lw=lw)

    #     dat = loadtxt("limit_data/DarkPhoton/Crab_2.txt")
    #     plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.9,lw=lw)
    #     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.9)
        if text_on:
            plt.text(0.5e-6,3e-1,r'{\bf Crab}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1,'k'),clip_on=True)
            plt.text(0.8e-6,0.9e-1,r'{\bf nebula}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1,'k'),clip_on=True)

        return


    def QUALIPHIDE(ax,col='r',fs=9,text_on=True,edge_on=False,lw=0.8,zorder=0):
        # data file is for randomly polarised case and 0.3 GeV/cm^3
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/QUALIPHIDE.txt")
        dat[:,1] = dat[:,1]*sqrt(1/3/0.13)*sqrt(0.3/0.45)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='k',facecolor=col,zorder=zorder,lw=0)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],'k-',lw=lw,zorder=zorder)
        if text_on:
            plt.text(3.5e-5,0.13e-12,r'{\bf QUALIPHIDE}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
        return
        
    def SHUKET(ax,col='maroon',fs=13,text_on=False,edge_on=False,lw=0.8):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SHUKET.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(1/3/0.038)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.2)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],'k-',lw=lw,zorder=0.2)
        if text_on:
            plt.text(3.5e-5,0.13e-12,r'{\bf SHUKET}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
        return

    def DarkEfield(ax,col='darkred',fs=17,text_on=True,edge_on=False,lw=0.8):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/DarkEfield.txt")
        dat[:,1] = dat[:,1]*sqrt(1.64/5) # convert from 5 sigma CL to 95%
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(1/3/0.129)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.2)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],'k-',lw=lw,zorder=0.2)
        if text_on:
            plt.text(0.8e-7/1.2,0.2e-12,r'{\bf Dark}',fontsize=fs,color=col,rotation=90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.text(2e-7/1.2,0.2e-12,r'{\bf E-field}',fontsize=fs,color=col,rotation=90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
        return

    def ORPHEUS(ax,col='darkred',fs=10,text_on=True,edge_on=False,lw=0.8):
        # data file is for randomly polarised case
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/ORPHEUS.txt")
        dat[:,1] = dat[:,1]*sqrt(1/3/0.01944939)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='k',facecolor=col,zorder=0.1,lw=0)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],'k-',lw=lw,zorder=0.2)
        if text_on:
            plt.text(6.5e-5,0.5e-13,r'{\bf ORPHEUS}',color=col,rotation=-90,fontsize=fs,clip_on=True)
        return

    def WISPDMX(ax,col='crimson',fs=12,text_on=True,edge_on=False,lw=0.8):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/WISPDMX.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(1/3/0.23)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.201)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.202,lw=lw)

        if text_on:
            plt.text(9e-7,4.1e-12/1.2,r'{\bf WISP}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.text(9e-7,1.8e-12/1.2,r'{\bf DMX}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)

        return

    def DOSUE(ax,col='red',fs=9,text_on=True,edge_on=False,lw=0.8):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/DOSUE-RR.txt")
        dat[:,1] = dat[:,1]*sqrt(2/3/0.29377804)*sqrt(0.39/0.45)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.201)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.202,lw=lw)

        if text_on:
            plt.text(90e-6,0.26e-10,r'{\bf DOSUE-RR}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
        return


    def SQuAD(ax,col=[0.7,0,0],fs=12,text_on=True,lw=0.5,point_on=False,ms=10):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SQuAD.txt")
        dat[:,1] = dat[:,1]*sqrt(0.4/0.45)*sqrt(1/3/0.019)
        plt.plot([dat[0,0],dat[0,0]],[y2,dat[0,1]],lw=lw,color=col,alpha=1,zorder=0.2)
        if point_on:
            plt.plot(dat[0,0],dat[0,1],'o',mfc=col,mec='k',mew=lw+1,zorder=0.2,markersize=ms)
        if text_on:
            plt.text(36e-6,0.25e-14,r'{\bf SQuAD}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
        return


    def DMPathfinder(ax,col='pink',fs=13,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/DM-Pathfinder.txt")
        dat[:,1] = dat[:,1]*sqrt(1/0.075)
        plt.plot([dat[0,0],dat[0,0]],[y2,dat[0,1]],lw=2,color=col,alpha=1,zorder=0.6)
        if text_on:
            plt.text(2.1e-9,0.5e-8/1.9,r'{\bf DM}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)
            plt.text(2.1e-9,0.2e-8/1.9,r'{\bf Pathfinder}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center',clip_on=True)

        return


    def QuantumCyclotron(ax,col='orangered',fs=13,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/QuantumCyclotron.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
        plt.plot([dat[0,0],dat[0,0]],[y2,dat[0,1]],lw=2,color=col,alpha=1,zorder=0.6,path_effects=line_background(2.5,'k'))
        if text_on:
            plt.text(0.95e-3,1e-10,r'{\bf QC}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center',clip_on=True)
        return

    def DarkMatter(ax,Witte_col='royalblue',Caputo_col='dodgerblue',Arias_col='navy',fs=20,projection=True,text_on=True):
        y2 = ax.get_ylim()[1]
        zo = 0.3
        pek=[pe.Stroke(linewidth=7, foreground='k'), pe.Normal()]

        # Combined limits
        dat = loadtxt("limit_data/DarkPhoton/DM_combined.txt")
        plt.plot(dat[:,0],dat[:,1],'-',color='w',alpha=1,zorder=zo+0.1,lw=2.5,path_effects=pek)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor='lightgray',zorder=zo,alpha=1.0)
        plt.plot([1e-16,dat[0,0]],[dat[0,1],dat[0,1]],'--',color='w',alpha=1,zorder=zo+0.1,lw=2.5,path_effects=pek)
        plt.fill_between([1e-16,dat[0,0]],[dat[0,1],dat[0,1]],y2=y2,edgecolor=None,facecolor='lightgray',zorder=zo+0.1,alpha=1.0)
        plt.plot(dat[40:,0],dat[40:,1],'--',color='w',alpha=1,lw=2.5,zorder=1000,solid_capstyle='round')

        # Individual limits
        dat2 = loadtxt("limit_data/DarkPhoton/Cosmology_Witte_inhomogeneous.txt")
        dat4 = loadtxt("limit_data/DarkPhoton/Cosmology_Caputo_HeII.txt",delimiter=',')
        dat5 = loadtxt("limit_data/DarkPhoton/Cosmology_Arias.txt")

        plt.fill_between(dat2[:,0],dat2[:,1],y2=y2,edgecolor='k',facecolor=Witte_col,zorder=0.305,alpha=0.8)
        plt.fill_between(dat4[:,0],dat4[:,1],y2=y2,edgecolor='k',facecolor=Caputo_col,zorder=0.305,alpha=0.8)
        plt.fill_between(dat5[:,0],dat5[:,1],y2=y2,edgecolor='k',facecolor=Arias_col,zorder=0.306,alpha=1)

        if text_on:
            plt.gcf().text(0.295,0.42-0.04,r'{\bf DPDM} HeII',fontsize=15,color='w',ha='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.gcf().text(0.295,0.4-0.04,r'Reionisation',fontsize=15,color='w',ha='center',clip_on=True)
            plt.gcf().text(0.295,0.38-0.04,r'(Caputo et al.)',fontsize=13,color='w',ha='center',clip_on=True)

            plt.gcf().text(0.365,0.37,r'{\bf DPDM}',fontsize=17,color='w',ha='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.gcf().text(0.365,0.35,r'(Witte et al.)',fontsize=13,color='w',ha='center',clip_on=True)

            plt.gcf().text(0.485,0.43,r'{\bf DPDM}',rotation=21.5,fontsize=18,color='w',va='center',ha='center',path_effects=line_background(1.5,'k'),rotation_mode='anchor',clip_on=True)
            plt.gcf().text(0.49,0.41,r'(Arias et al.)',rotation=21.5,fontsize=16,color='w',va='center',ha='center',path_effects=line_background(1,'k'),rotation_mode='anchor',clip_on=True)
    
        return

    def COBEFIRAS(ax,col=[0.1,0.2,0.5],text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat3 = loadtxt("limit_data/DarkPhoton/COBEFIRAS.txt",delimiter=',')
        plt.fill_between(dat3[:,0],dat3[:,1],y2=y2,edgecolor='k',facecolor=col,zorder=0.5,alpha=1)
        plt.plot(dat3[:,0],dat3[:,1],'k-',lw=lw,zorder=0.5)
        if text_on:
            plt.gcf().text(0.29,0.70,r'{\bf COBE/FIRAS}',fontsize=22,color='w',ha='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.gcf().text(0.29,0.67,r'$\gamma \rightarrow X$',fontsize=22,color='w',ha='center',path_effects=line_background(1,'k'),clip_on=True)
        return


    def LSW(ax,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SPring-8.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.45, 0.05, 0.1],zorder=1.1001)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.1001,lw=lw)

        dat = loadtxt("limit_data/DarkPhoton/ALPS.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.55, 0.0, 0.16],zorder=1.091)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.091,lw=lw)

        dat = loadtxt("limit_data/DarkPhoton/LSW_UWA.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.6, 0.0, 0.2],zorder=1.09)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.09,lw=lw)

        dat = loadtxt("limit_data/DarkPhoton/LSW_ADMX.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.65, 0.1, 0.24],zorder=1.089)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.089,lw=lw)

    #     dat = loadtxt("limit_data/DarkPhoton/LSW_CERN.txt")
    #     plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.089,lw=2)
    #     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.65, 0.15, 0.2],zorder=1.089)

        dat = loadtxt("limit_data/DarkPhoton/CROWS.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.7, 0.2, 0.2],zorder=1.08)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.08,lw=lw)

        dat = loadtxt("limit_data/DarkPhoton/DarkSRF.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.5, 0.2, 0.2],zorder=1.06)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.06,lw=lw)

   

        if text_on:
            plt.text(0.4e-6,0.15e-3,r'{\bf LSW-ADMX}',fontsize=17,color='w',rotation=-58,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(1e-5,5e-5,r'{\bf LSW-UWA}',fontsize=14,color='w',rotation=-56,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(0.55e0,0.9e-4,r'{\bf LSW-SPring-8}',fontsize=13,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(1.2e-4,0.9e-5,r'{\bf ALPS}',fontsize=25,color='w',rotation=-56,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(0.75e-7,9.9e-5,r'{\bf CROWS}',fontsize=24,color='w',rotation=-56,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(8.2e-7,0.4e-8,r'{\bf DarkSRF}',fontsize=17,color='w',rotation=-42,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)

        return

    def Coulomb(ax,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Cavendish.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.7,0,0],zorder=1.07)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.07,lw=lw)

        dat = loadtxt("limit_data/DarkPhoton/PlimptonLawton.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor='crimson',zorder=1.071)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.071,lw=lw)

        dat = loadtxt("limit_data/DarkPhoton/Spectroscopy.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.4, 0.0, 0.13],zorder=1.11)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.11,lw=lw)

        dat = loadtxt("limit_data/DarkPhoton/AFM.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.4, 0.2, 0.2],zorder=1.5)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.5,lw=lw)
        if text_on:
            plt.text(2.5e-10,0.35e-1,r'{\bf Plimpton-Lawton}',fontsize=15,color='w',rotation=-38,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(3e1,3e-1,r'{\bf AFM}',fontsize=20,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(0.5e-8,4e-6,r'{\bf Cavendish-Coulomb}',fontsize=23,color='w',rotation=-38,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(0.2e2,1e-3,r'{\bf Spectroscopy}',fontsize=23,color='w',rotation=-34,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)

        return

    def NeutronStarCooling(ax,col='#004d00',fs=18,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/NeutronStarCooling.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.1001)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.1001,lw=lw)
        if text_on:
            plt.text(0.9e4,0.4e-6,r'{\bf Neutron stars}',fontsize=fs,color='w',rotation=-45,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1,'k'),clip_on=True)
        return

    def CAST(ax,col='maroon',fs=19,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/CAST.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.1)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.1,lw=lw)
        if text_on:
            plt.text(0.95e-3,6e-6,r'{\bf CAST}',fontsize=fs,color='w',rotation=-59,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
        return

    def HINODE(ax,col='#700606',fs=16,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/HINODE.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.1001)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.1001,lw=lw)
        if text_on:
            plt.text(5e-3,0.3e-5,r'{\bf HINODE}',fontsize=fs,color='w',rotation=-59,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
        return

    def SHIPS(ax,col='indianred',fs=20,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SHIPS.txt")
        dat[:,1] = dat[:,1]/dat[:,0]
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.09)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.09,lw=lw)

        if text_on:
            plt.text(0.6e-1,0.08e-8,r'{\bf SHIPS}',fontsize=fs,color='w',rotation=-32,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
        return

    def TEXONO(ax,col=[0.5, 0.0, 0.13],fs=15,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/TEXONO.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.101)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.101,lw=lw)
        if text_on:
            plt.text(0.25e2,0.1e-4,r'{\bf TEXONO}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
        return

    def IGM(ax,col='seagreen',fs=18,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/IGM.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.49)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.49,lw=lw)

        if text_on:
            plt.text(4e-12,0.03e-7,r'{\bf IGM}',fontsize=fs,color='w',rotation=-39,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.gcf().text(0.233,0.565,r'{\bf DPDM heating}',color='w',fontsize=23,path_effects=line_background(1.5,'k'))

        return

    def LeoT(ax,col='mediumseagreen',fs=18,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/LeoT.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.3061)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.3062,lw=lw)

        if text_on:
            plt.text(7e-13,0.2e-9,r'{\bf Leo T}',fontsize=fs,color='w',rotation=-39,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
        return

    def GasClouds(ax,col='#00cc66',fs=18,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/GasClouds.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.306)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.307,lw=lw)

        if text_on:
            plt.text(0.86e-13,1e-10,r'{\bf Gas clouds}',fontsize=fs,color='w',rotation=-38.5,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
        return

    def SuperMAG(ax,col='#b5403e',fs=18,text_on=True,lw=1.5):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SuperMAG.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1,lw=lw)

        if text_on:
            plt.text(1.5e-17,1e-1/1.4,r'{\bf Super}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)
            plt.text(1.5e-17,0.2e-1/1.4,r'{\bf MAG}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center',path_effects=line_background(1.5,'k'),clip_on=True)

        return

##############################################################################################################################################################
############################# IAXO ###########################################################################################################################
##############################################################################################################################################################

    def IAXO(ax,col='magenta',fs=30,text_on=True,lw=2.5,pureL=False):
        y2 = ax.get_ylim()[1]
        
        suffix = "-10keV"
        #suffixGas = "-tPlasmon-newerE-gas"
        suffixGas = "AtlasGas-100eV-full"
        #"newstats-5yr-2-babyIAXO-tPlasmonGas.dat"
        #newstats-10keVbaselineIAXO-tPlasmonGas

        if pureL:
            #col = 'yellow' stats-10keVagainbaselineIAXO-tPlasmonGas.dat
            datGas = loadtxt("../data/limits/babyIAXO{}-pureL.dat".format(suffix))
            #datGas = loadtxt("../data/limits/baselineIAXO-tPlasmonGas.dat".format(suffix))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=301,lw=lw)
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.3, alpha=1.)

            datGas = loadtxt("../data/limits/baselineIAXO{}-pureL.dat".format(suffix))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=301,lw=lw, ls='-')
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.3, alpha=1.)
        
            datGas = loadtxt("../data/limits/upgradedIAXO{}-pureL.dat".format(suffix))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=301,lw=lw, ls='-')
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.3, alpha=0.6)


            datGas = loadtxt("../data/limits/stats-70eVagainbaselineIAXO-tPlasmonGas.dat".format(suffix))
            plt.plot(datGas[:,0],datGas[:,1],color='magenta', ls='--',alpha=1,zorder=0.301,lw=lw)

            if text_on:
                plt.text(1e-1,5e-11,r'{\bf babyIAXO}',fontsize=25,color='black',rotation=-33,rotation_mode='anchor',ha='center',va='center', zorder=105.5)


        else:
    #		babyIAXO
            datGas = loadtxt("../data/limits/stats-babyIAXO-{}.dat".format(suffixGas))
            #datGas = loadtxt("../data/limits/{}babyIAXO-tPlasmonGas.dat".format(suffixGas))
            plt.plot(datGas[:,0],datGas[:,1],color='cyan',alpha=1,zorder=100.301,lw=lw,ls='-',label="babyIAXO")
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor='cyan',zorder=0.3, alpha=1.)

    #		baselineIAXO
            datGas = loadtxt("../data/limits/stats-baselineIAXO-{}.dat".format(suffixGas))
            #datGas = loadtxt("../data/limits/{}baselineIAXO-tPlasmonGas.dat".format(suffixGas))
            plt.plot(datGas[:,0],datGas[:,1],color='orange',alpha=1,zorder=100.301,lw=lw,ls='--',label="IAXO")
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor='cyan',zorder=0.3, alpha=1.)
            
    #		upgradedIAXO
            datGas = loadtxt("../data/limits/stats-upgradedIAXO-{}.dat".format(suffixGas))
            #datGas = loadtxt("../data/limits/{}upgradedIAXO-tPlasmonGas.dat".format(suffixGas))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=100.301,lw=lw, ls='-.',label="IAXO+")
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor='cyan',zorder=0.3, alpha=1.)
            plt.vlines(datGas[-1,0], datGas[-1,1], 1, zorder=0.301,lw=lw, color='cyan')


            if text_on:
                plt.text(1e-1,5e-11,r'{\bf IAXO}',fontsize=fs,color='white',rotation=-32,rotation_mode='anchor',ha='center',va='center', zorder=105.5)

        return


##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################


    def SHIPS(ax,col='indianred',fs=20,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SHIPS.txt")
        dat[:,1] = dat[:,1]/dat[:,0]
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.09,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.09)
        if text_on:
            plt.text(0.6e-1*(1-0.05),0.08e-8*(1+0.1),r'{\bf SHIPS}',fontsize=fs,color='k',rotation=-32,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.6e-1,0.08e-8,r'{\bf SHIPS}',fontsize=fs,color='w',rotation=-32,rotation_mode='anchor',ha='center',va='center')
        return

    def TEXONO(ax,col=[0.5, 0.0, 0.13],fs=15,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/TEXONO.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.101,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.101)
        if text_on:
            plt.text(0.25e2*(1-0.01),0.1e-4*(1+0.08),r'{\bf TEXONO}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.25e2,0.1e-4,r'{\bf TEXONO}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
        return

    def IGM(ax,col='seagreen',fs=18,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/IGM.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.49)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.6,zorder=0.49,lw=2)

        if text_on:
            plt.text(4e-12*(1-0.05),0.03e-7*(1+0.07),r'{\bf IGM}',fontsize=fs,color='k',rotation=-39,rotation_mode='anchor',ha='center',va='center')
            plt.text(4e-12,0.03e-7,r'{\bf IGM}',fontsize=fs,color='w',rotation=-39,rotation_mode='anchor',ha='center',va='center')
            plt.gcf().text(0.233*(1-0.005),0.565*(1+0.003),r'{\bf DPDM heating}',color='k',fontsize=23)
            plt.gcf().text(0.233,0.565,r'{\bf DPDM heating}',color='w',fontsize=23)

        return

    def LeoT(ax,col='mediumseagreen',fs=18,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/LeoT.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.3061)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.68,zorder=0.3062,lw=2)

        if text_on:
            plt.text(7e-13*(1-0.05),0.2e-9*(1+0.07),r'{\bf Leo T}',fontsize=fs,color='k',rotation=-39,rotation_mode='anchor',ha='center',va='center')
            plt.text(7e-13,0.2e-9,r'{\bf Leo T}',fontsize=fs,color='w',rotation=-39,rotation_mode='anchor',ha='center',va='center')
        return

    def GasClouds(ax,col='#00cc66',fs=18,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/GasClouds.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.306)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.6,zorder=0.307,lw=2)

        if text_on:
            plt.text(0.86e-13*(1-0.07),1e-10*(1+0.07),r'{\bf Gas clouds}',fontsize=fs,color='k',rotation=-38,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.86e-13,1e-10,r'{\bf Gas clouds}',fontsize=fs,color='w',rotation=-38,rotation_mode='anchor',ha='center',va='center')
        return

    def SuperMAG(ax,col='#b5403e',fs=18,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SuperMAG.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1)

        if text_on:
            plt.text(1.5e-17*(1-0.05),1e-1*(1+0.05)/1.4,r'{\bf Super}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(1.5e-17,1e-1/1.4,r'{\bf Super}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(1.5e-17*(1-0.05),0.2e-1*(1+0.05)/1.4,r'{\bf MAG}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(1.5e-17,0.2e-1/1.4,r'{\bf MAG}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')

        return


#==============================================================================#
def MySaveFig(fig,pltname,pngsave=True):
    fig.savefig(pltdir+pltname+'.pdf',bbox_inches='tight')
    if pngsave:
        fig.set_facecolor('w') # <- not sure what matplotlib fucked up in the new version but it seems impossible to set png files to be not transparent now
        fig.savefig(pltdir_png+pltname+'.png',bbox_inches='tight',transparent=False)

def line_background(lw,col):
    return [pe.Stroke(linewidth=lw, foreground=col), pe.Normal()]
#==============================================================================#
