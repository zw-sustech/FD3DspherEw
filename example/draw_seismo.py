"""
draw_seismo: Plot and save the seismograms or the stress tensor for single receiver or multi-receivers

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
from scipy import integrate
from scipy import signal
import sys
sys.path.append(".")
from fdpytool import *


def draw_seismo(runpath,lineid,srec,erec,comp,figsize=[8,4],figdpi=128,figshow=1,figsave=0,fignm='fd3dseismo.png',filtflag=0,Fbp=[1,10]):

    """
    draw_seismo: Plot and save the seismograms or the stress tensor for single receiver or multi-receivers
    Input:
        runpath: path of the run example
        lineid: line index. Set '0' to isolated multi receivers 
        srec: starting receiver index to plot 
        erec: ending receiver index to plot
        comp: variables (velocity or stress) to plot. Option: 'Vx','Vy','Vz','Txx','Tyy','Tzz','Txy','Txz','Tyz','Ux','Uy','Uz','Ax','Ay','Az'
        figsize: figure size to show and save. Default: [8,4]
        figdpi: figure resolution to show and save. Default: 128
        figshow: show figure or not, 'show' can block the run. Default: 'True'
        figsave: save figure or not. Default: 'False'
        fignm: figure name including the format suffix when 'figsave' flag is true. Default: 'fd3dseismo.png'
        filtflag: butterworth bandpass filter or not. Default: 'False'
        Fbp: pass-band, only valid when 'filtflag' flag is true. Default: [1,10]
    Output:
        Figure windows and saved figures
    Usage:
        >>> from draw_seismo import draw_seismo
        >>> draw_seismo('../example',0,1,1,['Vx','Vy','Vz'],figdpi=300,filtflag=1,Fbp=[1,5])
        >>> draw_seismo('../fz',3,1,10,['Uz'],figsize=[10,4],figdpi=300,figshow=0,figsave=1,fignm='seismovz.eps')
        >>> draw_seismo('../run',0,3,8,['Tzz'])
        >>> draw_seismo('../run',2,1,6,['Vx','Vy'],filtflag=1,Fbp=[10,50])
    """

    fnm_conf,dir_coord,dir_metric,dir_media,dir_source,dir_station,dir_out=get_sim_path(runpath)

    comp0=[]
    for compnm in comp:
        if compnm[0] == 'U':
            comp0.append('U0')
        elif compnm[0] == 'V':
            comp0.append('V0')
        elif compnm[0] == 'A':
            comp0.append('A0')
        else:
            comp0.append('T0')

    ryid=[]; stid=[]
    for n in range(srec,erec+1):
        ryid.append(lineid)
        stid.append(n)

    ncomp=len(comp)
    nrecv=len(stid)

    coord_grid=np.empty((nrecv,3))
    for n in range(nrecv):
        seismoinfo=locate_seismo(fnm_conf,ryid[n],stid[n],dir_station)
        coord_grid[n,:]=retrieve_station(seismoinfo,'grid','stationdir',dir_station)
    
        VVx  =retrieve_seismo(seismoinfo,'Vx','outdir',dir_out)
        VVy  =retrieve_seismo(seismoinfo,'Vy','outdir',dir_out)
        VVz  =retrieve_seismo(seismoinfo,'Vz','outdir',dir_out)
    
        TTxx =retrieve_seismo(seismoinfo,'Txx','outdir',dir_out)
        TTxy =retrieve_seismo(seismoinfo,'Txy','outdir',dir_out)
        TTxz =retrieve_seismo(seismoinfo,'Txz','outdir',dir_out)
        TTyy =retrieve_seismo(seismoinfo,'Tyy','outdir',dir_out)
        TTyz =retrieve_seismo(seismoinfo,'Tyz','outdir',dir_out)
        TTzz =retrieve_seismo(seismoinfo,'Tzz','outdir',dir_out)
    
        T    =retrieve_seismo(seismoinfo,'time','outdir',dir_out)
        stept=T[-1]-T[-2]
    
        UUx=integrate.cumtrapz(VVx,dx=stept)
        UUy=integrate.cumtrapz(VVy,dx=stept)
        UUz=integrate.cumtrapz(VVz,dx=stept)
    
        AAx=np.gradient(VVx,stept)
        AAy=np.gradient(VVy,stept)
        AAz=np.gradient(VVz,stept)
    
        # Butterworth filter or not
        if filtflag:
            Noder=2
            Fnyq=1/stept/2
            fb,fa=signal.butter(Noder,np.array(Fbp)/Fnyq,btype='bandpass')
            VVx=signal.filtfilt(fb,fa,VVx)
            VVy=signal.filtfilt(fb,fa,VVy)
            VVz=signal.filtfilt(fb,fa,VVz)
            UUx=signal.filtfilt(fb,fa,UUx)
            UUy=signal.filtfilt(fb,fa,UUy)
            UUz=signal.filtfilt(fb,fa,UUz)
            AAx=signal.filtfilt(fb,fa,AAx)
            AAy=signal.filtfilt(fb,fa,AAy)
            AAz=signal.filtfilt(fb,fa,AAz)
    
        if n == 0:
            Vx=VVx.reshape(1,-1);   Vy=VVy.reshape(1,-1);   Vz=VVz.reshape(1,-1)
            Ux=UUx.reshape(1,-1);   Uy=UUy.reshape(1,-1);   Uz=UUz.reshape(1,-1)
            Ax=AAx.reshape(1,-1);   Ay=AAy.reshape(1,-1);   Az=AAz.reshape(1,-1)
            Txx=TTxx.reshape(1,-1); Txy=TTxy.reshape(1,-1); Txz=TTxz.reshape(1,-1)
            Tyy=TTyy.reshape(1,-1); Tyz=TTyz.reshape(1,-1); Tzz=TTzz.reshape(1,-1)
        else:
            Vx=np.append(Vx,VVx.reshape(1,-1),axis=0)
            Vy=np.append(Vy,VVy.reshape(1,-1),axis=0)
            Vz=np.append(Vz,VVz.reshape(1,-1),axis=0)
            Ux=np.append(Ux,UUx.reshape(1,-1),axis=0)
            Uy=np.append(Uy,UUy.reshape(1,-1),axis=0)
            Uz=np.append(Uz,UUz.reshape(1,-1),axis=0)
            Ax=np.append(Ax,AAx.reshape(1,-1),axis=0)
            Ay=np.append(Ay,AAy.reshape(1,-1),axis=0)
            Az=np.append(Az,AAz.reshape(1,-1),axis=0)
            Txx=np.append(Txx,TTxx.reshape(1,-1),axis=0)
            Txy=np.append(Txy,TTxy.reshape(1,-1),axis=0)
            Txz=np.append(Txz,TTxz.reshape(1,-1),axis=0)
            Tyy=np.append(Tyy,TTyy.reshape(1,-1),axis=0)
            Tyz=np.append(Tyz,TTyz.reshape(1,-1),axis=0)
            Tzz=np.append(Tzz,TTzz.reshape(1,-1),axis=0)
    
    
    # Maximum for normalization show
    Vx0=np.abs(Vx).max();   Vy0=np.abs(Vy).max();   Vz0=np.abs(Vz).max()
    Ux0=np.abs(Ux).max();   Uy0=np.abs(Uy).max();   Uz0=np.abs(Uz).max()
    Ax0=np.abs(Ax).max();   Ay0=np.abs(Ay).max();   Az0=np.abs(Az).max()
    Txx0=np.abs(Txx).max(); Txy0=np.abs(Txy).max(); Txz0=np.abs(Txz).max()
    Tyy0=np.abs(Tyy).max(); Tyz0=np.abs(Tyz).max(); Tzz0=np.abs(Tzz).max()
    
    V0=max(Vx0,Vy0,Vz0);  U0=max(Ux0,Uy0,Uz0);  A0=max(Ax0,Ay0,Az0)
    T0=max(Txx0,Txy0,Txz0,Tyy0,Tyz0,Tzz0)
    
    
    # Plot seismograms
    # Normalization show for multi-receiver
    if nrecv > 1:
        for icomp in range(ncomp):
            plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
            for irec in range(nrecv):
                Wshow=eval(comp[icomp] + '[' + str(irec) + ',:' + ']')/eval(comp0[icomp])

                Wshow1 = eval(comp[icomp] + '[' + str(irec-srec) + ',:' + ']')
                np.savetxt(dir_out+'/'+str(irec)+'.txt',Wshow1)
                with open ('Vz_rec'+str(irec)+'.txt','w') as f:
                    f.write(str(len(Wshow1))+'\n')
                    for i in range(len(Wshow1)):
                        f.write("{:.2f}  {:.5e} \n".format(stept*(i+1),Wshow1[i]))



                Tshow=T[0:Wshow.size]
                plt.plot(Tshow,Wshow+irec+1,'k-',linewidth=1.2)
                plt.title(comp[icomp])
                plt.xlabel('Time (s)')
                plt.ylabel('Receiver Index')
                plt.xlim(Tshow[0],Tshow[-1])
            #locs,labs=plt.yticks()
            plt.yticks(range(1,nrecv+1),range(1,nrecv+1))
            if figsave:
                if ncomp > 1:
                    imgnm=fignm[:-(len(fignm.split('.')[-1])+1)]
                    imgfmt=fignm.split('.')[-1]
                    imgfullnm='{}_{}.{}'.format(imgnm,comp[icomp],imgfmt)
                    plt.savefig(imgfullnm)
                else:
                    plt.savefig(fignm)
        if figshow:
            plt.show()

    
    # NO Normalization show for single receiver
    else:
        for icomp in range(ncomp):
            for irec in range(nrecv):
                Wshow=eval(comp[icomp] + '[' + str(irec) + ',:' + ']')

                Wshow1 = eval(comp[icomp] + '[' + str(irec-srec) + ',:' + ']')
                np.savetxt(dir_out+'/'+str(irec)+'.txt',Wshow1)
                with open ('Vz_rec'+str(irec)+'.txt','w') as f:
                    f.write(str(len(Wshow1))+'\n')
                    for i in range(len(Wshow1)):
                        f.write("{:.2f}  {:.5e} \n".format(stept*(i+1),Wshow1[i]))




                Tshow=T[0:Wshow.size]
                plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
                plt.plot(Tshow,Wshow,'k-',linewidth=1.2)
                plt.title(comp[icomp])
                plt.xlabel('Time (s)')
                plt.ylabel('Amplitude')
                plt.xlim(Tshow[0],Tshow[-1])
            if figsave:
                if ncomp > 1:
                    imgnm=fignm[:-(len(fignm.split('.')[-1])+1)]
                    imgfmt=fignm.split('.')[-1]
                    imgfullnm='{}_{}.{}'.format(imgnm,comp[icomp],imgfmt)
                    plt.savefig(imgfullnm)
                else:
                    plt.savefig(fignm)
        if figshow:
            plt.show()
    
