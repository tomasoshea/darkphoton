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


def FigSetup(xlab=r'$m_a$ [eV]',ylab='',\
                 g_min = 1.0e-19,g_max = 1.0e-6,\
                 m_min = 1.0e-12,m_max = 1.0e7,\
                 lw=2.5,lfs=45,tfs=25,tickdir='out',\
                 Grid=False,Shape='Rectangular',\
                 mathpazo=False,TopAndRightTicks=False,\
                xtick_rotation=20.0,tick_pad=8,x_labelpad=10,y_labelpad=10,\
             FrequencyAxis=False,N_Hz=1,upper_xlabel=r"$\nu_a$ [Hz]",**freq_kwargs):



    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)

    if mathpazo:
            plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Palatino"],
            })

    if Shape=='Wide':
        fig = plt.figure(figsize=(16.5,5))
    elif Shape=='Rectangular':
        fig = plt.figure(figsize=(16.5,11))
    elif Shape=='Square':
        fig = plt.figure(figsize=(14.2,14))

    ax = fig.add_subplot(111)

    ax.set_xlabel(xlab,fontsize=lfs,labelpad=x_labelpad)
    ax.set_ylabel(ylab,fontsize=lfs,labelpad=y_labelpad)

    ax.tick_params(which='major',direction=tickdir,width=2.5,length=13,right=TopAndRightTicks,top=TopAndRightTicks,pad=tick_pad)
    ax.tick_params(which='minor',direction=tickdir,width=1,length=10,right=TopAndRightTicks,top=TopAndRightTicks)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([m_min,m_max])
    ax.set_ylim([g_min,g_max])

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

    plt.xticks(rotation=xtick_rotation)

    if Grid:
        ax.grid(zorder=0)

    if FrequencyAxis:
        UpperFrequencyAxis(ax,N_Hz=N_Hz,tickdir='out',\
                           xtick_rotation=xtick_rotation,\
                           xlabel=upper_xlabel,\
                           lfs=lfs/1.3,tfs=tfs,tick_pad=tick_pad-2,**freq_kwargs)

    return fig,ax

#==============================================================================#


class DarkPhoton():
    def FigSetup(xlab=r'Dark photon mass [eV]',ylab='Kinetic mixing',\
             chi_min = 1.0e-18,chi_max = 1.0e0,\
             m_min = 3e-18,m_max = 1e5,\
             lw=2.5,lfs=40,tfs=25,tickdir='out',\
             Grid=False,Shape='Rectangular',mathpazo=True,\
             TopAndRightTicks=False,FrequencyAxis=True,FrequencyLabels=True,UnitAxis=True,f_rescale=1,\
            tick_rotation = 20,width=20,height=10,upper_tickdir='out'):

        """             chi_min = 1.0e-18,chi_max = 1.0e0,\
             m_min = 3e-18,m_max = 1e5,\ """
        
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

    def FigSetup_reduced(xlab=r'Dark photon mass [eV]',ylab='Kinetic mixing',\
             chi_min = 1.0e-13,chi_max = 1.0e-9,\
             m_min = 1e-2,m_max = 0.5e1,\
             lw=2.5,lfs=40,tfs=25,tickdir='out',\
             Grid=False,Shape='Rectangular',mathpazo=True,\
             TopAndRightTicks=False,FrequencyAxis=True,FrequencyLabels=True,UnitAxis=True,f_rescale=1,\
            tick_rotation = 20,width=20,height=10,upper_tickdir='out'):

        """             chi_min = 1.0e-18,chi_max = 1.0e0,\
             m_min = 3e-18,m_max = 1e5,\ """
        
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
        

    def Haloscopes(ax,col=[0.75, 0.2, 0.2],fs=17,projection=True,text_on=True):
        y2 = ax.get_ylim()[1]
        zo = 0.3

        HAYSTAC_col = 'indianred'
        CAPP_col = 'crimson'
        QUAX_col = 'r'
        ADMX_col = 'firebrick'

        # ADMX
        costh = sqrt(0.019)
        B = 7.6
        dat = loadtxt("limit_data/AxionPhoton/ADMX.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=ADMX_col,zorder=0.1,lw=3)

        B = 6.8
        dat = loadtxt("limit_data/AxionPhoton/ADMX2018.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=ADMX_col,zorder=0.1)

        B = 7.6
        dat = loadtxt("limit_data/AxionPhoton/ADMX2019_1.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=ADMX_col,zorder=0.1)

        B = 7.6
        dat = loadtxt("limit_data/AxionPhoton/ADMX2019_2.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=ADMX_col,zorder=0.1)

        B = 7.6
        dat = loadtxt("limit_data/AxionPhoton/ADMX2021.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=ADMX_col,zorder=0.1)

    #     B = 3.11
    #     dat = loadtxt("limit_data/AxionPhoton/ADMX_Sidecar_AC.txt")
    #     dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    #     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=ADMX_col,facecolor=ADMX_col,zorder=0.1)

    #     B = 5.0
    #     dat = loadtxt("limit_data/AxionPhoton/ADMX_SLIC.txt")
    #     dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    #     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=ADMX_col,facecolor=ADMX_col,zorder=100)



        B = 9
        dat = loadtxt("limit_data/AxionPhoton/HAYSTAC_highres.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=HAYSTAC_col,zorder=0.1)
        dat = loadtxt("limit_data/AxionPhoton/HAYSTAC_2020_highres.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=HAYSTAC_col,zorder=0.1)


        # CAPP
        B = 7.3
        dat = loadtxt("limit_data/AxionPhoton/CAPP-1.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=CAPP_col,zorder=0.1)

        B = 7.8
        dat = loadtxt("limit_data/AxionPhoton/CAPP-2.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=CAPP_col,zorder=0.1)

        B = 7.9
        dat = loadtxt("limit_data/AxionPhoton/CAPP-3.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=CAPP_col,zorder=0.1)

        # CAPP-3 [KSVZ]
        dat_min = dat[argmin(dat[:,1]),:]
        dat_min[1] = dat_min[1]*costh/sqrt(0.2)
        plt.plot([dat_min[0],dat_min[0]],[1e-10,dat_min[1]],'-',color=CAPP_col,lw=1.5,zorder=0.1)


        B = 8.1
        costh = sqrt(0.03)
        dat = loadtxt("limit_data/AxionPhoton/QUAX.txt")
        dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*0.023*dat[:,0]))
        plt.fill_between([dat[0,0],dat[0,0]],[y2,dat[0,1]],y2=y2,color=QUAX_col,zorder=0.1)


        if text_on:
            plt.text(1.4e-6,0.5e-14,r'{\bf ADMX}',fontsize=fs,color=ADMX_col,rotation=90,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.8e-5,0.1e-13,r'{\bf CAPP}',fontsize=fs-2,color=CAPP_col,rotation=90,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.19e-4,3e-15,r'{\bf HAYSTAC}',fontsize=fs-5,color=HAYSTAC_col,rotation=90,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.47e-4,3e-12,r'{\bf QUAX}',fontsize=fs-8,color=QUAX_col,rotation=-90,rotation_mode='anchor',ha='center',va='center')

        return


    def StellarBounds(ax,fs=19,text_on=True):
        y2 = ax.get_ylim()[1]
        # Stellar physics constraints

        # Globular clusters
        HB_col = [0.01, 0.75, 0.24]
        HB = loadtxt("limit_data/DarkPhoton/RG.txt")
        plt.plot(HB[:,0],HB[:,1],color='k',alpha=0.5,zorder=0.9,lw=2)
        plt.fill_between(HB[:,0],HB[:,1],y2=y2,edgecolor=None,facecolor=HB_col,zorder=0.9)

        # Globular clusters
        HB_col = 'DarkGreen'
        HB = loadtxt("limit_data/DarkPhoton/HB.txt")
        plt.plot(HB[:,0],HB[:,1],color='k',alpha=0.5,zorder=0.95,lw=2)
        plt.fill_between(HB[:,0],HB[:,1],y2=y2,edgecolor=None,facecolor=HB_col,zorder=0.95)

        # Solar bound
        Solar_col = 'ForestGreen'
        Solar = loadtxt("limit_data/DarkPhoton/Solar.txt")
        plt.plot(Solar[:,0],Solar[:,1],color='k',alpha=0.5,zorder=1.02,lw=2)
        plt.fill_between(Solar[:,0],Solar[:,1],y2=y2,edgecolor=None,facecolor=Solar_col,zorder=1.02)

        Solar = loadtxt("limit_data/DarkPhoton/Solar-Global.txt")
        plt.plot(Solar[:,0],Solar[:,1]/Solar[:,0],color='k',alpha=0.5,zorder=1.021,lw=2)
        plt.fill_between(Solar[:,0],Solar[:,1]/Solar[:,0],y2=y2,edgecolor=None,facecolor=Solar_col,zorder=1.021)

        if text_on:
            plt.text(0.8e2*(1-0.01),1.5e-14*(1+0.05),r'{\bf Solar}',fontsize=fs,color='k',rotation=-41,rotation_mode='anchor',ha='center',va='center')
            plt.text(1e3*(1-0.01),0.7e-14*(1+0.05),r'{\bf HB}',fontsize=fs,color='k',rotation=-38,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.8e4*(1-0.01),0.7e-14*(1+0.05),r'{\bf RG}',fontsize=fs,color='k',rotation=-37,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.8e2,1.5e-14,r'{\bf Solar}',fontsize=fs,color='w',rotation=-41,rotation_mode='anchor',ha='center',va='center')
            plt.text(1e3,0.7e-14,r'{\bf HB}',fontsize=fs,color='w',rotation=-38,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.8e4,0.7e-14,r'{\bf RG}',fontsize=fs,color='w',rotation=-37,rotation_mode='anchor',ha='center',va='center')
        return


    def Xenon(ax,col='crimson',fs=23,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Xenon1T.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)

        plt.plot(1e3*dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.5,lw=2)
        plt.fill_between(1e3*dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)

        dat = loadtxt("limit_data/DarkPhoton/Xenon1T_S1S2.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)

        plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.5,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)


        dat = loadtxt("limit_data/DarkPhoton/XENON1T_SE.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.5,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)


        dat = loadtxt("limit_data/DarkPhoton/XENON1T_Solar_S2.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.5,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)


        dat = loadtxt("limit_data/DarkPhoton/XENONnT.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.5,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)

        if text_on:
            plt.text(8e2,2.5e-17,r'{\bf XENON}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')

        return


    def DAMIC(ax,col='salmon',fs=21,text_on=True):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/DAMIC.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)

        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]
        dat[-1,1] = y2[-1]
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.001,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.001)
        if text_on:
            plt.text(6e-1,1.3e-14,r'{\bf DAMIC}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.plot([5e0,1e1],[3e-14,6e-14],'-',lw=2.5,color=col)
        return

    def MuDHI(ax,col='#a82844',fs=15,text_on=True):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/MuDHI.txt")

        y2 = interp(dat[:,0],m1,y1)
        dat[dat[:,1]>y2,1] = y2[dat[:,1]>y2]
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=10)
        plt.plot(dat[::2,0],dat[::2,1],color='k',alpha=1,zorder=10,lw=1.2)

        if text_on:
            plt.text(0.42e-2,1.8e-11,r'{\bf MuDHI}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.plot([2.8e-2,1.5e0],[2e-11,4e-11],'-',lw=2,color=col,path_effects=line_background(1,'k'))
        return


    def FUNK(ax,col='red',fs=21,text_on=True):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/FUNK.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(2/3/0.27)

        y2 = interp(dat[:,0],m1,y1)
        #dat[0,1] = y2[0]
        dat[dat[:,1]>y2,1] = y2[dat[:,1]>y2]
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.3,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.3)
        if text_on:
            plt.text(2.6e-1,1e-13,r'{\bf FUNK}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.plot([9e-1,3e0],[3e-13,1e-12],'-',lw=2.5,color=col)
        return

    def SENSEI(ax,col='firebrick',fs=21,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SENSEI.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)

        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1)
        if text_on:
            plt.text(1.7e0,1e-15,r'{\bf SENSEI}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.plot([7e0,1e1],[3e-15,9e-15],'-',lw=2.5,color=col)
        return

    def SuperCDMS(ax,col=[0.4,0,0],fs=18,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SuperCDMS.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.6,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.6)

        if text_on:
            plt.text(0.5e1,1.5e-16,r'{\bf SuperCDMS}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.plot([5e1,0.8e2],[3e-16,9e-16],'-',lw=2.5,color=col)
        return

    def Nanowire(ax,col='pink',fs=22,text_on=True):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/WSi_Nanowire.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]/1.1
        dat[-1,1] = y2[-1]/1.1
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.3,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.3)
        if text_on:
            plt.text(5e-4,1e-10,r'{\bf WSi Nanowire}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.plot([9e-3,3e-3],[3e-10,9e-10],'-',lw=2.5,color=col)
        return


    def SQMS(ax,col='#02734b',fs=17,text_on=True,lw=0.5,ms=10):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SQMS.txt")
        dat[:,1] = dat[:,1]*sqrt(1/3/0.019)
        plt.plot(dat[:,0],dat[:,1],lw=lw,color=col,alpha=1,zorder=0.0)
        if text_on:
            plt.text(5.7e-6,0.65e-14,r'{\bf SQMS}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
        return

    def LAMPOST(ax,col='red',fs=15,text_on=True):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/LAMPOST.txt")
        dat[:,1] = dat[:,1]*sqrt(0.4/0.45)*sqrt(2/3/0.27)

        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]/1.1
        dat[-1,1] = y2[-1]/1.1
        plt.plot(dat[:,0],dat[:,1],color=col,alpha=1,zorder=0,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0)
        if text_on:
            plt.text(0.3e-1,5e-13,r'{\bf LAMPOST}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.plot([3e-1,0.6e0],[6e-13,1e-12],'-',lw=1.5,color=col)
        return

    def Tokyo(ax,col='darkred',fs=15,text_on=True):
        m1,y1 = loadtxt("limit_data/DarkPhoton/DM_combined.txt",unpack=True)
        dat = loadtxt("limit_data/DarkPhoton/Tokyo-Dish.txt")
        dat[:,1] = dat[:,1]*sqrt(2/3/0.6)
        y2 = interp(dat[:,0],m1,y1)
        dat[0,1] = y2[0]
        dat[-1,1] = y2[-1]
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.4,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.4)


        dat = loadtxt("limit_data/DarkPhoton/Tokyo-Knirck.txt")
        dat[:,1] = dat[:,1]*sqrt(1/3/0.175)
        plt.fill_between(dat[:,0],dat[:,1],y2=1e0,edgecolor='k',facecolor=col,zorder=1.09)

        dat = loadtxt("limit_data/DarkPhoton/Tokyo-Tomita.txt")
        plt.plot([dat[1,0],dat[1,0]],[dat[1,1],1e0],'-',color=col,lw=3,zorder=0.2)
        if text_on:
            plt.text(2.3e-4,2.5e-10,r'{\bf Tokyo-3}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.45e-3,3e-8,r'{\bf Tokyo-2}',fontsize=fs-2,color='k',rotation=90,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.2e-1,4e-12,r'{\bf Tokyo-1}',fontsize=fs+4,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.plot([2.05e-1,4e0],[5e-12,8e-12],'-',lw=2.5,color=col)
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
            plt.text(7e-6,4e-12,r'{\bf FAST}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
        return

    def Jupiter(ax,col='Green',fs=15,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Jupiter.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=2,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=2)
        if text_on:
            plt.text(0.1e-14*(1-0.02),4.5e-1*(1+0.07),r'{\bf Jupiter}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.1e-14,4.5e-1,r'{\bf Jupiter}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
        return

    def Earth(ax,col='DarkGreen',fs=17,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Earth.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.9,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.9)
        if text_on:
            plt.text(0.4e-13*(1-0.01),2e-1*(1+0.05),r'{\bf Earth}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.4e-13,2e-1,r'{\bf Earth}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
        return


    def Crab(ax,col=[0.1,0.4,0.1],fs=17,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Crab.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.09999,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.09999)

    #     dat = loadtxt("limit_data/DarkPhoton/Crab_2.txt")
    #     plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.9,lw=2)
    #     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.9)
        if text_on:
            plt.text(0.5e-6*(1-0.02),3e-1*(1+0.07),r'{\bf Crab}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.5e-6,3e-1,r'{\bf Crab}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')

            plt.text(0.8e-6*(1-0.02),0.9e-1*(1+0.07),r'{\bf nebula}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.8e-6,0.9e-1,r'{\bf nebula}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')

        return


    def SHUKET(ax,col='maroon',fs=13,text_on=True,edge_on=False,lw=0.8):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SHUKET.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(1/3/0.038)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.2)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],'k-',lw=lw,zorder=0.2)
        if text_on:
            plt.text(3.5e-5,0.13e-12,r'{\bf SHUKET}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
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
            plt.text(0.8e-7/1.2,0.2e-12,r'{\bf Dark}',fontsize=fs,color=col,rotation=90,rotation_mode='anchor',ha='center',va='center')
            plt.text(2e-7/1.2,0.2e-12,r'{\bf E-field}',fontsize=fs,color=col,rotation=90,rotation_mode='anchor',ha='center',va='center')
        return

    def ORPHEUS(ax,col='darkred',fs=10,text_on=True,edge_on=False,lw=0.8):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/ORPHEUS.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='k',facecolor=col,zorder=0.1,lw=0)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],'k-',lw=lw,zorder=0.2)
        if text_on:
            plt.text(6.5e-5,1e-13,r'{\bf ORPHEUS}',color=col,rotation=-90,fontsize=fs)
        return

    def WISPDMX(ax,col='crimson',fs=12,text_on=True,edge_on=False,lw=0.8):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/WISPDMX.txt")
        dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(1/3/0.23)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.201)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.202,lw=lw)

        if text_on:
            plt.text(9e-7,4.1e-12/1.2,r'{\bf WISP}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(9e-7,1.8e-12/1.2,r'{\bf DMX}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')

        return

    def DOSUE(ax,col='red',fs=9,text_on=True,edge_on=False,lw=0.8):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/DOSUE-RR.txt")
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.201)
        if edge_on:
            plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.202,lw=lw)

        if text_on:
            plt.text(90e-6,0.26e-10,r'{\bf DOSUE-RR}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
        return


    def SQuAD(ax,col=[0.7,0,0],fs=12,text_on=True,lw=0.5,point_on=False,ms=10):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SQuAD.txt")
        dat[:,1] = dat[:,1]*sqrt(0.4/0.45)*sqrt(1/3/0.019)
        plt.plot([dat[0,0],dat[0,0]],[y2,dat[0,1]],lw=lw,color=col,alpha=1,zorder=0.2)
        if point_on:
            plt.plot(dat[0,0],dat[0,1],'o',mfc=col,mec='k',mew=lw+1,zorder=0.2,markersize=ms)
        if text_on:
            plt.text(36e-6,0.25e-14,r'{\bf SQuAD}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
        return


    def DMPathfinder(ax,col='pink',fs=13,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/DM-Pathfinder.txt")
        dat[:,1] = dat[:,1]*sqrt(1/0.075)
        plt.plot([dat[0,0],dat[0,0]],[y2,dat[0,1]],lw=2,color=col,alpha=1,zorder=0.6)
        if text_on:
            plt.text(2.1e-9,0.5e-8/1.9,r'{\bf DM}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(2.1e-9,0.2e-8/1.9,r'{\bf Pathfinder}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')

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
            plt.gcf().text(0.295,0.42-0.04,r'{\bf DPDM} HeII',fontsize=15,color='w',ha='center')
            plt.gcf().text(0.295,0.4-0.04,r'Reionisation',fontsize=15,color='w',ha='center')
            plt.gcf().text(0.295,0.38-0.04,r'(Caputo et al.)',fontsize=13,color='w',ha='center')

            plt.gcf().text(0.365,0.37,r'{\bf DPDM}',fontsize=17,color='w',ha='center')
            plt.gcf().text(0.365,0.35,r'(Witte et al.)',fontsize=13,color='w',ha='center')

            plt.gcf().text(0.49,0.48,r'{\bf DPDM}',fontsize=18,color='w',ha='center')
            plt.gcf().text(0.49,0.46,r'(Arias et al.)',fontsize=16,color='w',ha='center')

        return

    def COBEFIRAS(ax,col=[0.1,0.2,0.5],text_on=True):
        y2 = ax.get_ylim()[1]
        dat3 = loadtxt("limit_data/DarkPhoton/COBEFIRAS.txt",delimiter=',')
        plt.fill_between(dat3[:,0],dat3[:,1],y2=y2,edgecolor='k',facecolor=col,zorder=0.5,alpha=1)
        if text_on:
            plt.gcf().text(0.29,0.70,r'{\bf COBE/FIRAS}',fontsize=22,color='w',ha='center')
            plt.gcf().text(0.29,0.67,r'$\gamma \rightarrow X$',fontsize=22,color='w',ha='center')
        return


    def LSW(ax,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/SPring-8.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.1001,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.45, 0.05, 0.1],zorder=1.1001)

        dat = loadtxt("limit_data/DarkPhoton/ALPS.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.091,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.55, 0.0, 0.16],zorder=1.091)

        dat = loadtxt("limit_data/DarkPhoton/LSW_UWA.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.09,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.6, 0.0, 0.2],zorder=1.09)

        dat = loadtxt("limit_data/DarkPhoton/LSW_ADMX.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.089,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.65, 0.1, 0.24],zorder=1.089)

    #     dat = loadtxt("limit_data/DarkPhoton/LSW_CERN.txt")
    #     plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.089,lw=2)
    #     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.65, 0.15, 0.2],zorder=1.089)

        dat = loadtxt("limit_data/DarkPhoton/CROWS.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.08,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.7, 0.2, 0.2],zorder=1.08)

        if text_on:
            plt.text(0.4e-6,0.15e-3,r'{\bf LSW-ADMX}',fontsize=17,color='w',rotation=-58,rotation_mode='anchor',ha='center',va='center')
            plt.text(1e-5,5e-5,r'{\bf LSW-UWA}',fontsize=14,color='w',rotation=-56,rotation_mode='anchor',ha='center',va='center')

            plt.text(0.55e0*(1-0.02),0.9e-4*(1+0.08),r'{\bf LSW-SPring-8}',fontsize=13,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.55e0,0.9e-4,r'{\bf LSW-SPring-8}',fontsize=13,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')


            plt.text(1.2e-4*(1-0.02),0.9e-5*(1+0.08),r'{\bf ALPS}',fontsize=25,color='k',rotation=-56,rotation_mode='anchor',ha='center',va='center')
            plt.text(1.2e-4,0.9e-5,r'{\bf ALPS}',fontsize=25,color='w',rotation=-56,rotation_mode='anchor',ha='center',va='center')

            plt.text(0.75e-7*(1-0.01),9.9e-5*(1+0.05),r'{\bf CROWS}',fontsize=24,color='k',rotation=-56,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.75e-7,9.9e-5,r'{\bf CROWS}',fontsize=24,color='w',rotation=-56,rotation_mode='anchor',ha='center',va='center')
        return

    def Coulomb(ax,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/Cavendish.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.07,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.7,0,0],zorder=1.07)

        dat = loadtxt("limit_data/DarkPhoton/PlimptonLawton.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.071,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor='crimson',zorder=1.071)

        dat = loadtxt("limit_data/DarkPhoton/Spectroscopy.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.11,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.4, 0.0, 0.13],zorder=1.11)

        dat = loadtxt("limit_data/DarkPhoton/AFM.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.5,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.4, 0.2, 0.2],zorder=1.5)
        if text_on:
            plt.text(2.5e-10*(1-0.02),0.35e-1*(1+0.08),r'{\bf Plimpton-Lawton}',fontsize=15,color='k',rotation=-38,rotation_mode='anchor',ha='center',va='center')
            plt.text(2.5e-10,0.35e-1,r'{\bf Plimpton-Lawton}',fontsize=15,color='w',rotation=-38,rotation_mode='anchor',ha='center',va='center')

            plt.text(3e1*(1-0.02),3e-1*(1+0.08),r'{\bf AFM}',fontsize=20,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
            plt.text(3e1,3e-1,r'{\bf AFM}',fontsize=20,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')

            plt.text(0.5e-8*(1-0.02),4e-6*(1+0.08),r'{\bf Cavendish-Coulomb}',fontsize=23,color='k',rotation=-38,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.5e-8,4e-6,r'{\bf Cavendish-Coulomb}',fontsize=23,color='w',rotation=-38,rotation_mode='anchor',ha='center',va='center')

            plt.text(0.2e2*(1-0.01),1e-3*(1+0.08),r'{\bf Spectroscopy}',fontsize=23,color='k',rotation=-34,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.2e2,1e-3,r'{\bf Spectroscopy}',fontsize=23,color='w',rotation=-34,rotation_mode='anchor',ha='center',va='center')

        return

    def NeutronStarCooling(ax,col='#004d00',fs=18,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/NeutronStarCooling.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.1001,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.1001)

        if text_on:
            plt.text(0.9e4*(1-0.03),0.4e-6*(1+0.05),r'{\bf Neutron stars}',fontsize=fs,color='k',rotation=-43,rotation_mode='anchor',ha='center',va='center')
            plt.text(0.9e4,0.4e-6,r'{\bf Neutron stars}',fontsize=fs,color='w',rotation=-43,rotation_mode='anchor',ha='center',va='center')
        return

    def CAST(ax,col='maroon',fs=27,text_on=True):
        y2 = ax.get_ylim()[1]
        dat = loadtxt("limit_data/DarkPhoton/CAST.txt")
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.1,lw=2)
        plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.1)
        if text_on:
            plt.text(4e-3*(1-0.01),0.8e-6*(1+0.08),r'{\bf CAST}',fontsize=fs,color='k',rotation=-59,rotation_mode='anchor',ha='center',va='center')
            plt.text(4e-3,0.8e-6,r'{\bf CAST}',fontsize=fs,color='w',rotation=-59,rotation_mode='anchor',ha='center',va='center')
        return

##############################################################################################################################################################
############################# IAXO ###########################################################################################################################
##############################################################################################################################################################

    def IAXO(ax,col='magenta',fs=30,text_on=True,lw=1.5,pureL=False):
        y2 = ax.get_ylim()[1]
        
        suffix = "-newE"
        suffixGas = "-newE"

        if pureL:
            col = 'yellow'
            datGas = loadtxt("../data/limits/babyIAXO{}-pureL.dat".format(suffix))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=0.301,lw=lw)
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.3, alpha=1.)

            datGas = loadtxt("../data/limits/baselineIAXO{}-pureL.dat".format(suffix))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=0.301,lw=lw, ls='-')
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.3, alpha=1.)#
            
            datGas = loadtxt("../data/limits/upgradedIAXO{}-pureL.dat".format(suffix))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=0.301,lw=lw, ls='-')
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.3, alpha=1.)

            if text_on:
                plt.text(1e-1,5e-11,r'{\bf Upgraded IAXO}',fontsize=20,color=col,rotation=-32,rotation_mode='anchor',ha='center',va='center', zorder=105.5)


        else:
    #		babyIAXO
            datGas = loadtxt("../data/limits/babyIAXO-tPlasmon{}.dat".format(suffixGas))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=100.301,lw=lw)
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor='cyan',zorder=0.3, alpha=1.)

    #		baselineIAXO
            datGas = loadtxt("../data/limits/baselineIAXO-tPlasmon{}.dat".format(suffixGas))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=100.301,lw=lw, ls='-')
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor='cyan',zorder=0.3, alpha=1.)
            
    #		upgradedIAXO
            datGas = loadtxt("../data/limits/upgradedIAXO-tPlasmon{}.dat".format(suffixGas))
            plt.plot(datGas[:,0],datGas[:,1],color='black',alpha=1,zorder=100.301,lw=lw, ls='-')
            plt.fill_between(datGas[:,0],datGas[:,1],y2=y2,edgecolor=None,facecolor='cyan',zorder=0.3, alpha=1.)
            plt.vlines(datGas[-1,0], datGas[-1,1], 1, zorder=100.301,lw=lw, color='black')

            if text_on:
                plt.text(1e-1,5e-11,r'{\bf IAXO}',fontsize=fs,color='cyan',rotation=-32,rotation_mode='anchor',ha='center',va='center', zorder=105.5)

        return

#####################################################################

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