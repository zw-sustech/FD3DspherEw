"""
fd3dtopoew.plot.save.py: Python scripts to plot and save the figures for seismo, media and snapshot

Author:
  Yuanhang Huo  (yhhuo@mail.ustc.edu.cn)
  Wei Zhang     (zhangwei@sustech.edu.cn)
Affiliation:
  University of Science and Technology of China
  Southern University of Science and Technology
Date:
  2019-09-04
"""

import os
import time
import datetime
import subprocess
import matplotlib
matplotlib.use('Agg')
import sys
sys.path.append('.')    # path of python subfunctions
from draw_seismo import draw_seismo
from draw_media import draw_media
from draw_snap import draw_snap

# **************************  Set path and directory ****************************** #
PROJDIR='./'
FIGDIR='./figs'
if not os.path.exists(FIGDIR):
    os.makedirs(FIGDIR)

# Starting time show
print('Starting Time: ' + str(datetime.datetime.now()) + '\n')

# cancel warning of max figures
matplotlib.rcParams.update({'figure.max_open_warning': 0})


'''
# ********************************************************************************* #
#                                QC: Plot the model
# ********************************************************************************* #
print('Start to plot and save the model images ... ')
ts=time.perf_counter()

# ----------------------------- Parameters Box Start ------------------------------ #
# variables name to plot
varplot=['Vp','Vs','rho']          # Option: 'Vp','Vs','rho','lambda','mu'. One or more
# for all slices init
subsall=[]
subcall=[]
subtall=[]
# for the first media slice (XOY)
subs=[1,1,1]            # starting grid index to the media
subc=[-1,-1,1]          # counts of grid to media
subt=[1,1,1]            # grid stride to the media
subsall.append(subs)
subcall.append(subc)
subtall.append(subt)
# for the second media slice (XOZ)
subs=[1,1,1]
subc=[-1,1,-1]
subt=[1,1,1]
subsall.append(subs)
subcall.append(subc)
subtall.append(subt)
## for the third media slice (YOZ)
subs=[1,1,1]
subc=[1,-1,-1]
subt=[1,1,1]
subsall.append(subs)
subcall.append(subc)
subtall.append(subt)
# figure and axis
axisopt='scaled'        # axis scaling
flagtitle=1             # if show figure title, '1' to 'True' and '0' to 'False'
flagemlast=1            # if emlast, '1' to 'True' and '0' to 'False'
flagkm=1                # axis units selection, '1' to 'km' and '0' to 'm'
flagclb=1               # if show figure colorbar, '1' to 'True' and '0' to 'False'
clbtype='jet'           # colormap type
clbrange=[None,None]    # colorbar range
figsize=[6,4]           # figure size to show and save 
figdpi=300              # figure resolution to show and save 
figshow=0               # if show figure window, '1' to 'True' and '0' to 'False'. \
                        # Warning: 'True' can block the script run
figsave=1               # save figure or not, '1' to 'True' and '0' to 'False'
fignm='fd3dmedia.png'   # figure name including format suffix when 'figsave' flag is 'True'
# ----------------------------- Parameters Box End -------------------------------- #

# plot
for isub in range(len(subsall)):
    for varnm in varplot:
        if len(subsall) == 1 and len(varplot) == 1:
            draw_media(PROJDIR,varnm,subsall[0],subcall[0],subtall[0],setaxis=axisopt,\
                    flag_title=flagtitle,flag_emlast=flagemlast,flag_km=flagkm,flag_clb=flagclb,\
                    clbtype=clbtype,clbrange=clbrange,figsize=figsize,figdpi=figdpi,figshow=0,\
                    figsave=figsave,fignm=fignm)
            if figsave:
                subprocess.call('mv {} {}'.format(fignm,FIGDIR),shell=True)
        else:
            imgnm=fignm[:-(len(fignm.split('.')[-1])+1)]
            imgfmt=fignm.split('.')[-1]
            imgfullnm='{}_slice{}_{}.{}'.format(imgnm,isub+1,varnm,imgfmt)
            draw_media(PROJDIR,varnm,subsall[isub],subcall[isub],subtall[isub],setaxis=axisopt,\
                    flag_title=flagtitle,flag_emlast=flagemlast,flag_km=flagkm,flag_clb=flagclb,\
                    clbtype=clbtype,clbrange=clbrange,figsize=figsize,figdpi=figdpi,figshow=0,\
                    figsave=figsave,fignm=imgfullnm)
            if figsave:
                subprocess.call('mv {} {}'.format(imgfullnm,FIGDIR),shell=True)

# running time show
te=time.perf_counter()
runtime=te-ts
runtm,runts=divmod(runtime,60)
runth,runtm=divmod(runtm,60)
print('   Time cost: {:0>2d}:{:0>2d}:{:0>2d}'.format(int(runth),int(runtm),int(runts)))
# ********************************************************************************* #

'''

# ********************************************************************************* #
#                                QC: Plot the seismo
# ********************************************************************************* #
print('Start to plot and save the seismo images ... ')
ts=time.perf_counter()

# ----------------------------- Parameters Box Start ------------------------------ #
lineid=0                # line index, '0' to isolated multi receivers, '1','2'... to #line
srec=1                  # starting receiver index to plot
erec=6                 # ending receiver index to plot
#comp=['Vx','Vy','Vz']
comp=['Vy']                        # variable names to plot. Option: 'Vx','Vy','Vz','Ux','Uy','Uz',\
                        # 'Ax','Ay','Az','Txx','Txy','Txz','Tyy','Tyz','Tzz'
figsize=[10,5]           # figure size to show and save
figdpi=800              # figure resolution to show and save
figshow=0               # show figure window or not, '1' to 'True' and '0' to 'False'. \
                        # Warning: 'True' can block the script run
figsave=1               # save figure or not, '1' to 'True' and '0' to 'False'
fignm='fd3dseismo1.png'  # figure name including format suffix when 'figsave' flag is 'True'
filtflag=0              # butterworth bandpass filter or not, '1' to 'True' and '0' to 'False'
Fbp=[10,50]              # pass-band, only valid when 'filtflag' flag is 'True'
# ----------------------------- Parameters Box End -------------------------------- #

# plot
draw_seismo(PROJDIR,lineid,srec,erec,comp,figsize=figsize,figdpi=figdpi,figshow=figshow,figsave=figsave,fignm=fignm,filtflag=filtflag,Fbp=Fbp)

# move to FIGDIR
if figsave:
    subprocess.call('mv {}*.{} {}'.format(fignm.split('.')[0],fignm.split('.')[1],FIGDIR),shell=True)

# running time show
te=time.perf_counter()
runtime=te-ts
runtm,runts=divmod(runtime,60)
runth,runtm=divmod(runtm,60)
print('   Time cost: {:0>2d}:{:0>2d}:{:0>2d}'.format(int(runth),int(runtm),int(runts)))
# ********************************************************************************* #





# ********************************************************************************* #
#                                QC: Plot the snapshot
# ********************************************************************************* #
print('Start to plot and save the snapshot images ... ')
ts=time.perf_counter()

# ----------------------------- Parameters Box Start ------------------------------ #
# variable names to plot
varplot=['Vy']    
                        # Option: 'Vx','Vy','Vz','Txx','Txy','Txz','Tyy','Tyz','Tzz'
# for all slices init
snapidall=[]
subsall=[]
subcall=[]
subtall=[]
# for the first snapshot slice
snapid=1          # snap slice index, see 'SeisFD3D.conf'
subs=[1,1,1]            # staring grid index to snapshot
subc=[-1,-1,-1]          # counts of grid to snapshot
subt=[1,1,1]            # grid stride to snapshot
snapidall.append(snapid)
subsall.append(subs)
subcall.append(subc)
subtall.append(subt)
# for the second snapshot slice
#snapid=2
#subs=[1,1,1]
#subc=[-1,-1,1]
#subt=[1,1,1]
#snapidall.append(snapid)
#subsall.append(subs)
#subcall.append(subc)
#subtall.append(subt)
# time layer, axis and figure
tinds=1                 # staring time layer to plot and save
tinde=10000           # ending time layer to plot and save 这个步数 乘以时间步长 为想要的时刻
tstep=100            # stride of time layer to plot and save
axisopt='scaled'        # axis scaling. Option: 'on','off','equal','scaled','tight','auto', \
                        # 'normal','image','square'
flagemlast=1            # if emlast, '1' to 'True' and '0' to 'False'
flagkm=1                # axis units selection, '1' to 'km' and '0' to 'm'
flagclb=1               # if show figure colorbar, '1' to 'True' and '0' to 'False'
clbtype='seismic'       # colormap type
clbrange=[-1e-11,1e-11]    # colorbar range
figsize=[5,4]           # figure size to show and save
figdpi=300              # figure resolution to show and save
figshow=0               # if show figure window, '1' to 'True' and '0' to 'False'
imgsave=1               # if save figures to every time layer, '1' to 'True' and '0' to 'False'
gifsave=0               # if save snapshot GIF, '1' to 'True' and '0' to 'False'
fignm='fd3dsnap_snap0.png'    # figure name including format suffix when 'imgsave' or 'gifsave' is 'True'
# ----------------------------- Parameters Box End -------------------------------- #

# plot
for isub in range(len(snapidall)):
    for varnm in varplot:
        if len(snapidall) == 1 and len(varplot) == 1:
            draw_snap(PROJDIR,varnm,snapidall[0],subsall[0],subcall[0],subtall[0],\
                    tinds,tinde,tstep,setaxis=axisopt,flag_emlast=flagemlast,flag_km=flagkm,\
                    flag_clb=flagclb,clbtype=clbtype,clbrange=clbrange,figsize=figsize,\
                    figdpi=figdpi,figshow=figshow,imgsave=imgsave,gifsave=gifsave,fignm=fignm)
            if imgsave or gifsave:
                subprocess.call('mv {}*.* {}'.format(fignm.split('.')[0],FIGDIR),shell=True)
        else:
            imgnm=fignm[:-(len(fignm.split('.')[-1])+1)]
            imgfmt=fignm.split('.')[-1]
            imgfullnm='{}_snap{}_{}.{}'.format(imgnm,isub+1,varnm,imgfmt)
            draw_snap(PROJDIR,varnm,snapidall[isub],subsall[isub],subcall[isub],subtall[isub],\
                    tinds,tinde,tstep,setaxis=axisopt,flag_emlast=flagemlast,flag_km=flagkm,\
                    flag_clb=flagclb,clbtype=clbtype,clbrange=clbrange,figsize=figsize,\
                    figdpi=figdpi,figshow=figshow,imgsave=imgsave,gifsave=gifsave,fignm=imgfullnm)
            if imgsave or gifsave:
                subprocess.call('mv {}*.* {}'.format(imgfullnm.split('.')[0],FIGDIR),shell=True)

# running time show
te=time.perf_counter()
runtime=te-ts
runtm,runts=divmod(runtime,60)
runth,runtm=divmod(runtm,60)
print('   Time cost: {:0>2d}:{:0>2d}:{:0>2d}'.format(int(runth),int(runtm),int(runts)))
# ********************************************************************************* #


print('\nEnding Time:   ' + str(datetime.datetime.now()))
