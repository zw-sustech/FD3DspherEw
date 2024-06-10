"""
draw_snap: Plot and save the snapshot generated from FD3DTopoEw

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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imageio
import subprocess
import sys

#import utm
from subprocess import call
from sys import stdout,stderr
import math
import os
#from geopy.distance import geodesic
sys.path.append(".")
from fdpytool import *

def draw_snap(runpath,varnm,snapid,subs,subc,subt,tinds,tinde,tstep,setaxis='scaled',\
              flag_emlast=1,flag_km=1,flag_clb=1,clbtype='seismic',clbrange=[None,None],\
              figsize=[5,4],figdpi=300,figshow=0,imgsave=0,gifsave=0,fignm='fd3dsnap.png'):

    """
    draw_snap: Plot and save the snapshot generated from FD3DTopoEw
    Input:
        runpath: path of the run example
        varnm: variable name to plot and save.\
               Option: 'Vx','Vy','Vz','Txx','Txy','Txz','Tyy','Tyz','Tzz'
        snapid: snap index, see 'SeisFD3D.conf'
        subs: starting grid index to snapshot
        subc: counts of grid to snapshot
        subt: grid stride to snapshot
        tinds: starting time layer to plot and save
        tinde: ending time layer to plot and save
        tstep: stride of time layer to plot and save
        setaxis: set option to the axis. Default: 'scaled'.\
                 Option: 'on','off','equal','scaled','tight','auto','normal','image','square'
        flag_emlast: emlast flag, Default: 'True'
        flag_km: axis units selection, if 'True', it's 'km'; if not, 'm'. Default: 'True'
        flag_clb: display colorbar in the figure or not. Default: 'True'
        clbtype: colrmap tyep. Option: 'jet','viridis','plasma','magma','Greys','spring',\
                 'summer','coolwarm','bwr','seismic','hsv','rainbow'. Default: 'seismic' \
                 See full to visit: https://matplotlib.org/turorials/colors/colormaps.html
        clbrange: colorbar range. Default: None
        figsize: figure size to show and save. Default: [5,4]
        figdpi: figure resolution to show and save. Default: 300
        figshow: show figure or not. Warning: 'show' can block the run. Default: 'True'
        imgsave: save figures to every time layer or not. Default: 'False'
        gifsave: save snapshot GIF or not. Default: 'False'
        fignm: figure name including the format suffix when 'imgsave' or 'gifsave' \
               flag is True. Default: 'fd3dsnap.png'
    Output:
        Figure windows, saved figures or GIF animation
    Usage:
        >>> from draw_snap import draw_snap
        >>> draw_snap('../fz','Vz',[1,1,1],[-1,1,-1],[1,1,1],1,400,10,figshow=0,imgsave=1,gifsave=1,fignm='snapVz.jpg')
        >>> draw_snap('../example','Vx',[1,1,1],[-1,-1,1],[1,1,1],100,600,20,setaxis='square',flag_km=0,figsize=[6,6],figdpi=600)
        >>> draw_snap('../run','Tzz',[1,1,1],[1,-1,-1],[2,2,2],1,1000,100)
    """
    
    
    fnm_conf,dir_coord,dir_metric,dir_media,dir_source,dir_station,dir_out=get_sim_path(runpath)
    
    snapinfo=locate_snap(fnm_conf,snapid,'start',subs,'count',subc,'stride',subt)
    #print(snapinfo)
    x,y,z=gather_coord(snapinfo,'coorddir',dir_coord)
    
#    nx=x.shape[0]; ny=x.shape[1]; nz=z.shape[2]
    nx=x.shape[0]; ny=y.shape[0]; nz=z.shape[0]
    print(nx,ny,nz)
    
    str_unit='m'
#    if flag_km:
#        x=x/1e3; y=y/1e3; z=z/1e3
#        str_unit='km'
    if flag_km:
        x=x*180/math.pi; y=y*180/math.pi; z=z/1e3
        str_unit='km'
#    plt.ion()
    plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
    for nlayer in range(tinds,tinde+tstep,tstep):
        print('nlayer',nlayer)
        v,t=gather_snap(snapinfo,snapid,nlayer,varnm,'outdir',dir_out)
#        print("Drawing " + str(nlayer) + "th time layer (t = " + "{:.4f}".format(t) + " s)")
        print('v.shape',v.shape)
        figtitle="Snapshot of " + varnm + " at " "{:9.4f}".format(t) + " s"
        
        if nx == 1:
            if flag_emlast:
#                Y=np.squeeze(y).transpose(1,0)[::-1,:]
#                Y=np.append(Y,Y[-1,:].reshape(1,-1),axis=0)
#                Y=np.append(Y,Y[:,-1].reshape(-1,1),axis=1)
#                Z=np.squeeze(z).transpose(1,0)[::-1,:]
#                Z=np.append(Z,Z[-1,:].reshape(1,-1),axis=0)
#                Z=np.append(Z,Z[:,-1].reshape(-1,1),axis=1)
                #Y=y[::-1]
                Y=y
                Y=np.append(Y,Y[-1]) # lon
                #print(Y)
                Z=z[::-1]
                Z=np.append(Z,Z[-1]) # depth
                #print(Z)
                
                V=np.squeeze(v).transpose(1,0)[::-1,:]
                V=np.append(V,V[-1,:].reshape(1,-1),axis=0)
                V=np.append(V,V[:,-1].reshape(-1,1),axis=1)


                plt.clf()
                plt.cla()
                plt.pcolor(Y,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
                plt.yticks(fontsize=3)
                plt.xlabel('Y ' + '(' + str_unit + ')')
                plt.ylabel('Z ' + '(' + str_unit + ')')
                #plt.axis(setaxis)
                plt.title(figtitle)
                plt.colorbar()
                if figshow:
                    plt.pause(0.1)
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
    
                plt.clf()
                plt.cla()
                plt.pcolor(Y,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
        
                plt.xlabel('Y ' + '(' + str_unit + ')')
                plt.ylabel('Z ' + '(' + str_unit + ')')
                plt.axis(setaxis)
                plt.title(figtitle)
                plt.colorbar()
                if figshow:
                    plt.pause(0.1)
        
        elif ny == 1:
            if flag_emlast:
#                X=np.squeeze(x).transpose(1,0)[::-1,:]
#                X=np.append(X,X[-1,:].reshape(1,-1),axis=0)
#                X=np.append(X,X[:,-1].reshape(-1,1),axis=1)
#                Z=np.squeeze(z).transpose(1,0)[::-1,:]
#                Z=np.append(Z,Z[-1,:].reshape(1,-1),axis=0)
#                Z=np.append(Z,Z[:,-1].reshape(-1,1),axis=1)
                X=x[::-1]
                X=np.append(X,X[-1])
                
                Z=z[::-1]
                Z=np.append(Z,Z[-1])
                
                V=np.squeeze(v).transpose(1,0)[::-1,:]
                V=np.append(V,V[-1,:].reshape(1,-1),axis=0)
                V=np.append(V,V[:,-1].reshape(-1,1),axis=1)
    
                plt.clf()
                plt.cla()
                plt.pcolor(X,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
        
                plt.xlabel('X ' + '(' + str_unit + ')')
                plt.ylabel('Z ' + '(' + str_unit + ')')
                plt.axis(setaxis)
                plt.title(figtitle)
                plt.colorbar()
                if figshow:
                    plt.pause(0.1)
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
    
                plt.clf()
                plt.cla()
                plt.pcolor(X,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
        
                plt.xlabel('X ' + '(' + str_unit + ')')
                plt.ylabel('Z ' + '(' + str_unit + ')')
                plt.axis(setaxis)
                plt.title(figtitle)
                plt.colorbar()
                if figshow:
                    plt.pause(0.1)
    
        
        else:
            if flag_emlast:
                #print(x)
                #print(np.squeeze(x))
                #print(len(x))
#                X=np.squeeze(x).transpose(1,0)[::-1,:]
#                X=np.append(X,X[-1,:].reshape(1,-1),axis=0)
#                X=np.append(X,X[:,-1].reshape(-1,1),axis=1)
#                Y=np.squeeze(y).transpose(1,0)[::-1,:]
#                Y=np.append(Y,Y[-1,:].reshape(1,-1),axis=0)
#                Y=np.append(Y,Y[:,-1].reshape(-1,1),axis=1)
                X=x[::-1]
                X=np.append(X,X[-1])
          
                Y=y[::-1]
                Y=np.append(Y,Y[-1])
                
                V=np.squeeze(v).transpose(1,0)[::-1,:]
                V=np.append(V,V[-1,:].reshape(1,-1),axis=0)
                V=np.append(V,V[:,-1].reshape(-1,1),axis=1)
  

                plt.clf()
                plt.cla()
                plt.pcolor(X,Y,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
                plt.xticks(fontsize=2)
                plt.xlabel('X ' + '(' + str_unit + ')')
                plt.ylabel('Y ' + '(' + str_unit + ')')
                plt.axis(setaxis)
                plt.title(figtitle)
                plt.colorbar()
                if figshow:
                    plt.pause(0.1)
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
    
                plt.clf()
                plt.cla()
                plt.pcolor(X,Y,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
        
                plt.xlabel('X ' + '(' + str_unit + ')')
                plt.ylabel('Y ' + '(' + str_unit + ')')
                plt.axis(setaxis)
                plt.title(figtitle)
                plt.colorbar()
                if figshow:
                    plt.pause(0.1)
        
        if imgsave or gifsave:
            imgnm=fignm[:-(len(fignm.split('.')[-1])+1)]
            imgfmt=fignm.split('.')[-1]
            imgfullnm='{}_nlayer_{}.{}'.format(imgnm,nlayer,imgfmt)
            plt.savefig(imgfullnm)

    if gifsave:
        frames=[]
        for nlayer in range(tinds,tinde+tstep,tstep):
            imgfullnm='{}_nlayer_{}.{}'.format(imgnm,nlayer,imgfmt)
            frames.append(imageio.imread(imgfullnm))
        imageio.mimsave('{}.gif'.format(imgnm),frames,'GIF',duration=0.5)
            
    if imgsave == 0 and gifsave == 1:
        subprocess.call('rm {}_nlayer_*.{}'.format(imgnm,imgfmt),shell=True)

    
#    plt.ioff()
#    plt.show()
    if figshow:
        plt.show()
#        plt.show(block=False)
#    plt.pause(6)
#    plt.close()
    #print(X)
    #print(Z)
    #print(x[0])
    

