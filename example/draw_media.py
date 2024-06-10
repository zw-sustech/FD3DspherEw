"""
                    imgfmt=fignm.split('.')[1]
draw_media: Plot and save the media used to FD3DTopoEw forward

Author:
  Yuanhang Huo  (yhhuo@mail.ustc.edu.cn)
  Wei Zhang     (zhangwei@sustech.edu.cn)
Affiliation:
  University of Science and Technology of China
  Southern University of Science and Technology
Date:
  2019-09-06

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from fdpytool import *

def draw_media(runpath,varnm,subs,subc,subt,setaxis='scaled',flag_title=1,flag_emlast=1,\
               flag_km=1,flag_clb=1,clbtype='jet',clbrange=[None,None],figsize=[5,4],\
               figdpi=300,figshow=1,figsave=0,fignm='fd3dmedia.png'):

    """
    draw_media: Plot and save the media used to forward modeling
    Input:
        runpath: path of the run example
        varnm: variable name to plot and save. Option: 'Vp','Vs','rho','mu','lambda'
        subs: starting grid index to the media
        subc: counts of grid to the media
        subt: grid stride to the media
        setaxis: set option to the axis. Default: 'scaled'. \
                 Option: 'on','off','equal','scaled','tight','auto','normal','image','square'
        flag_title: title flag, display title in the figure or not. Default: 'True'
        flag_emlast: emlast flag. Default 'True'
        flag_km: axis units selection, if 'True', it's 'km'; if not, 'm'. Default: 'True'
        flag_clb: display colorbar in the figure or not. Default: 'True'
        clbtype: colormap type. Option: 'jet','viridis','plasma','magma','Greys','spring', \
                 'summer','coolwarm','bwr','seismic','hsv','rainbow'. Default: 'jet' \
                 See full to visit: https://matplotlib.org/tutorials/colors/colormaps.html
        clbrange: colorbar range. Default: None 
        figsize: figure size to show and save. Default: [5,4]
        figdpi: figure resolution to show and save. Default: 300
        figshow: show figure or not. Warning: 'show' can block the run. Default: 'True'
        figsave: save figure or not. Default: 'False'
        fignm: figure name including the format suffix when 'figsave' flag is True. \
               Default: 'fd3dmedia.png'
    Output:
        Figure windows and saved figures
    Usage:
        >>> from draw_media import draw_media
        >>> draw_media('../fz','Vp',[1,1,1],[-1,1,-1],[1,1,1])
        >>> draw_media('../example','rho',[1,1,1],[-1,-1,1],[1,1,1],setaxis='square',figsize=[3,3],figdpi=300,figshow=0,figsave=1,fignm='mediarho.eps')
        >>> draw_media('../run','Vs',[10,10,10],[1,-1,-1],[2,2,2],figdpi=600,figshow=1)
    """
    
    flag_jetwr=0
    
    fnm_conf,dir_coord,dir_metric,dir_media,dir_source,dir_station,dir_out=get_sim_path(runpath)
    
    snapinfo=locate_snap(fnm_conf,0,'start',subs,'count',subc,'stride',subt)
    
    x,y,z=gather_coord(snapinfo,'coorddir',dir_coord)
    nx=x.shape[0]; ny=x.shape[1]; nz=z.shape[2]
    
    str_unit='m'
    if flag_km:
        x=x/1e3; y=y/1e3; z=z/1e3
        str_unit='km'
    
    if varnm == 'Vp':
        rho=gather_media(snapinfo,'rho','mediadir',dir_media)
        mu=gather_media(snapinfo,'mu','mediadir',dir_media)
        lamd=gather_media(snapinfo,'lambda','mediadir',dir_media)
        v=((lamd+2*mu)/rho)**0.5
        v=v/1e3
    elif varnm == 'Vs':
        rho=gather_media(snapinfo,'rho','mediadir',dir_media)
        mu=gather_media(snapinfo,'mu','mediadir',dir_media)
        v=(mu/rho)**0.5
        v=v/1e3
    elif varnm == 'rho':
        v=gather_media(snapinfo,varnm,'mediadir',dir_media)
        v=v/1e3
    else:
        v=gather_media(snapinfo,varnm,'mediadir',dir_media)
    
    
    # Pcolor plot
    plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
    if nx == 1:
        if flag_emlast:
            Y=np.squeeze(y).transpose(1,0)[::-1,:]
            Y=np.append(Y,Y[-1,:].reshape(1,-1),axis=0)
            Y=np.append(Y,Y[:,-1].reshape(-1,1),axis=1)
            Z=np.squeeze(z).transpose(1,0)[::-1,:]
            Z=np.append(Z,Z[-1,:].reshape(1,-1),axis=0)
            Z=np.append(Z,Z[:,-1].reshape(-1,1),axis=1)
            V=np.squeeze(v).transpose(1,0)[::-1,:]
            V=np.append(V,V[-1,:].reshape(1,-1),axis=0)
            V=np.append(V,V[:,-1].reshape(-1,1),axis=1)
    
            plt.pcolor(Y,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
    
        else:
            Y=np.squeeze(y).transpose(1,0)
            Y=np.append(Y,Y[-1,:].reshape(1,-1),axis=0)
            Y=np.append(Y,Y[:,-1].reshape(-1,1),axis=1)
            Z=np.squeeze(z).transpose(1,0)
            Z=np.append(Z,Z[-1,:].reshape(1,-1),axis=0)
            Z=np.append(Z,Z[:,-1].reshape(-1,1),axis=1)
            V=np.squeeze(v).transpose(1,0)
            V=np.append(V,V[-1,:].reshape(1,-1),axis=0)
            V=np.append(V,V[:,-1].reshape(-1,1),axis=1)
    
            plt.pcolor(Y,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
    
        plt.xlabel('Y ' + '(' + str_unit + ')')
        plt.ylabel('Z ' + '(' + str_unit + ')')
        plt.axis(setaxis)
        plt.xlim(np.min(Y),np.max(Y))
        plt.ylim(np.min(Z),np.max(Z))
    
    elif ny == 1:
        if flag_emlast:
            X=np.squeeze(x).transpose(1,0)[::-1,:]
            X=np.append(X,X[-1,:].reshape(1,-1),axis=0)
            X=np.append(X,X[:,-1].reshape(-1,1),axis=1)
            Z=np.squeeze(z).transpose(1,0)[::-1,:]
            Z=np.append(Z,Z[-1,:].reshape(1,-1),axis=0)
            Z=np.append(Z,Z[:,-1].reshape(-1,1),axis=1)
            V=np.squeeze(v).transpose(1,0)[::-1,:]
            V=np.append(V,V[-1,:].reshape(1,-1),axis=0)
            V=np.append(V,V[:,-1].reshape(-1,1),axis=1)
    
            plt.pcolor(X,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
    
        else:
            X=np.squeeze(x).transpose(1,0)
            X=np.append(X,X[-1,:].reshape(1,-1),axis=0)
            X=np.append(X,X[:,-1].reshape(-1,1),axis=1)
            Z=np.squeeze(z).transpose(1,0)
            Z=np.append(Z,Z[-1,:].reshape(1,-1),axis=0)
            Z=np.append(Z,Z[:,-1].reshape(-1,1),axis=1)
            V=np.squeeze(v).transpose(1,0)
            V=np.append(V,V[-1,:].reshape(1,-1),axis=0)
            V=np.append(V,V[:,-1].reshape(-1,1),axis=1)
    
            plt.pcolor(X,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
    
        plt.xlabel('X ' + '(' + str_unit + ')')
        plt.ylabel('Z ' + '(' + str_unit + ')')
        plt.axis(setaxis)
        plt.xlim(np.min(X),np.max(X))
        plt.ylim(np.min(Z),np.max(Z))
        
    
    else:
        if flag_emlast:
            X=np.squeeze(x).transpose(1,0)[::-1,:]
            X=np.append(X,X[-1,:].reshape(1,-1),axis=0)
            X=np.append(X,X[:,-1].reshape(-1,1),axis=1)
            Y=np.squeeze(y).transpose(1,0)[::-1,:]
            Y=np.append(Y,Y[-1,:].reshape(1,-1),axis=0)
            Y=np.append(Y,Y[:,-1].reshape(-1,1),axis=1)
            V=np.squeeze(v).transpose(1,0)[::-1,:]
            V=np.append(V,V[-1,:].reshape(1,-1),axis=0)
            V=np.append(V,V[:,-1].reshape(-1,1),axis=1)
    
            plt.pcolor(X,Y,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
    
        else:
            X=np.squeeze(x).transpose(1,0)
            X=np.append(X,X[-1,:].reshape(1,-1),axis=0)
            X=np.append(X,X[:,-1].reshape(-1,1),axis=1)
            Y=np.squeeze(y).transpose(1,0)
            Y=np.append(Y,Y[-1,:].reshape(1,-1),axis=0)
            Y=np.append(Y,Y[:,-1].reshape(-1,1),axis=1)
            V=np.squeeze(v).transpose(1,0)
            V=np.append(V,V[-1,:].reshape(1,-1),axis=0)
            V=np.append(V,V[:,-1].reshape(-1,1),axis=1)
    
            plt.pcolor(X,Y,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
    
        plt.xlabel('X ' + '(' + str_unit + ')')
        plt.ylabel('Y ' + '(' + str_unit + ')')
        plt.axis(setaxis)
        plt.xlim(np.min(X),np.max(X))
        plt.ylim(np.min(Y),np.max(Y))
    
    
    if flag_clb:
        clb=plt.colorbar()
        if varnm == 'Vp' or varnm == 'Vs':
            clb.ax.set_ylabel('$\mathrm{\mathsf{(km/s)}}$')
        if varnm == 'rho':
            clb.ax.set_ylabel('$\mathrm{\mathsf{(g/cm^3)}}$')
    
    if flag_title:
        plt.title(varnm.title())
    
    if figsave:
        plt.savefig(fignm)
    if figshow:
        plt.show()
    
