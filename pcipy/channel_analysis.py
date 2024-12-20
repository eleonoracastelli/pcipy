#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:40:50 2024

@author: E Castelli, J Baker

These are various functions for producing diagnositic analyses of PCI channels
"""
from . import plotting
import matplotlib.pyplot as plt
import numpy as np

def stationarity_plots(filter_obj, ydata, title='',ndens=1000,split_orders=False,nchan=None,zero_mean=False,detrend=False):
    '''
    Make time-domain plots to explore PCI channel stationarity

    Arguments:
       filter_obj   ndarray   An ndarray containing the PCI filter instance or a
                    or LIST   list of those for side-by-side comparison.
       ydata        ndarray   ydata to be filtered
    Optional:
       title        str       Plot title
       ndens        int       target number of samples to show in each plot curve
                              or 0 for all samples
       
    '''
    filters=filter_obj
    if not isinstance(filters,list):
        filters=[filters]
    if not isinstance(title,list):
        title=[title]*len(filters)

    if nchan is None:
        nchan=filters[0].maxcompts
    
    plotting.plotconfig(lbsize=20, lgsize=16, fsize=18, 
                    ticklabelsize=20, style='publication',
                    fontfamily = 'STIXGeneral')
    fig1, axs = plt.subplots(ncols=len(filters),squeeze=False,sharey=True)
    axs=axs[0]
    for i,xf in enumerate(filters):
        print("i=",i)
        x=np.array(xf.apply_for_channels(ydata, n_channels=nchan,zero_mean=zero_mean,detrend=detrend))
        for ich in range(nchan):
            #xch=xf.channels[-(ich+1)]
            xch=x[-(ich+1)]
            
            if ndens==0:
                ev=1
            else:
                ev=len(xch)//ndens+1
            axs[i].plot(xch[::ev],label='c'+str(ich),zorder=nchan-ich)
            axs[i].set_title(title[i])
            #check var:
            chmean=np.mean(xch)
            print('mean:',chmean)
            print('variance: ',sum((xch-chmean)**2)/(len(xch)-1))
            axs[i].legend()
    plt.show()

    if split_orders:
        order=max([xf.order for xf in filters])
        plotting.plotconfig(lbsize=20, lgsize=16, fsize=18, 
                            ticklabelsize=20, style='publication',
                            fontfamily = 'STIXGeneral')
        nf=len(filters)
        fig1, axs = plt.subplots(nrows=order+2,ncols=len(filters),squeeze=False,figsize=[6.4*nf, 4.8*((order+2))],sharex=True,sharey=True)
        for ax in axs.flat:
            ax.label_outer()
            
        for i,xf in enumerate(filters):
            x=np.array(xf.apply_for_channels_split_orders(ydata, n_channels=nchan,zero_mean=zero_mean,detrend=detrend))
            xsum=np.sum(x,axis=0)
            ax=axs[0,i]
            for ich in range(x.shape[1]):
                ax.plot(xsum[-(ich+1)][::ev],label='c'+str(ich)+' net',zorder=nchan-ich)
            for iord in range(xf.order+1):
                ax=axs[iord+1,i]
                for ich in range(x.shape[1]):
                    ax.plot(x[iord,-(ich+1)][::ev],label='c'+str(ich)+',$T^'+str(iord)+'$ part',zorder=nchan-ich)
                    #t=title[i]
                    #if len(t)>0:t=t+": Partial at time-order "+str(iord)
                    #ax.set_title(t)
                    ax.legend()
            
        plt.show()



def temporal_variance_corr_plots(filter_obj, ydata, title='',ndens=1000,split_orders=False,nchan=None,zero_mean=False,detrend=False):
    filters=filter_obj
    if not isinstance(filters,list):
        filters=[filters]
    if not isinstance(title,list):
        title=[title]*len(filters)
    nplot=nchan
    nf=len(filters)
    
    if nchan is None:
        nchan=filters[0].maxcompts

    plotting.plotconfig(lbsize=20, lgsize=16, fsize=18, 
                    ticklabelsize=20, style='publication',
                    fontfamily = 'STIXGeneral')
    fig1, axs1 = plt.subplots(ncols=len(filters),nrows=1,squeeze=False,sharey=True,figsize=[6.4*nf, 4.8*((1))])
    fig1, axs2 = plt.subplots(ncols=len(filters),nrows=nplot,squeeze=False,sharey=True,figsize=[6.4*nf, 4.8*((nplot))])

    for i,xf in enumerate(filters):
        print("i=",i)
        x=np.array(xf.apply_for_channels(ydata, n_channels=nchan,zero_mean=zero_mean,detrend=detrend))
        #we compute channel covariance in batches of size nbatch
        print('plotting',nplot,'of',nchan,'channels')
        #trim the data to a multiple of nbatch
        nbatch=1000
        ndata=len(x.T)
        ntrim=(ndata//nbatch)*nbatch
        binned=x[:,:ntrim].reshape(nchan,-1,nbatch)
        if True:
            print("subtracting bin means!")
            binned=(binned.T-binned.mean(axis=2).T).T
        #print(binned.shape)
        print(ntrim//nbatch,'batches of size',nbatch)
        #print(binned[0,0:2])
        if ndens==0:
            ev=1
        else:
            ev=ntrim//nbatch//ndens+1

        #title='Channel stdevs for '+plots_label

        ax=axs1[0,i]
        stdevs=np.sqrt(np.mean(binned**2,axis=2))
        for j in range(nplot):
            jj=nchan-1-j
            print('stdevs shape',stdevs.shape,'ev',ev)
            print(stdevs[jj,::ev].shape)
            ax.semilogy(stdevs[jj,::ev], marker='.',linewidth=1, label=r'$\mathrm{ch}$[' + str(j+1) + ']')
        ax.legend(loc='upper left', ncol=2)
        ax.grid(linewidth=1, which='both', color='gray', linestyle='dotted')
        ax.set_xlabel("sample")
        ax.set_ylabel(r"stdev")
        ax.set_title(title)
        #plt.show()
        for k in range(nplot-1):                    
            kk=nchan-1-k                
            #title='Channel cross correlation vs channel '+str(k+1)+' for '+plots_label
            ax=axs2[k,i]
            for j in range(nplot-k):
                jj=kk-j
                covar=np.mean(binned[jj]*binned[kk],axis=1)
                #print(jj,kk,covar[1],stdevs[jj,1],stdevs[kk,1],stdevs[jj,1]*stdevs[kk,1])
                #print((binned[jj]*binned[kk])[:,1])
                #print(covar.shape)
                ax.plot((covar/stdevs[jj]/stdevs[kk])[::ev], marker='.',linewidth=1, label=r'$\mathrm{ch}$[' + str(nchan-kk)+','+str(nchan-jj)+ ']')
            ax.legend(loc='upper left', ncol=2)
            ax.grid(linewidth=1, which='both', color='gray', linestyle='dotted')
            ax.set_xlabel("sample")
            ax.set_ylabel(r"stdev")
            #axes.set_xlim([1e-3, 1])
            # axes.set_ylim([1e-21, 1e-16])
            ax.set_title(title)

