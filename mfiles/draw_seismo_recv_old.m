% retreive seismogram on individual receivers.
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
id=0; n=4;
flag_acc = 0; flag_vel = 0; flag_dis = 0;
flag_1canvas = 1;
flag_twin=0;
twin=[34 40];
%twin=[34.6 36.7];
%twin=[38.5 41.5];
%twin=[34.6 36.7];
flag_overlap = 1;
spec_line='r';

fnm_conf='../fy/SeisFD3D.conf';
pnm_info='../fy/input/';
pnm_out='../fy/output/';

%fnm_conf='../fy-a100-surf/SeisFD3D.conf';
%pnm_info='../fy-a100-surf/input/';
%pnm_out='../fy-a100-surf/output/';

% -------------------- load data --------------------------
[Vx,Vy,Vz,t,x,y,z] = ...
     retrieve_seismo(id,n,'confname',fnm_conf,'outdir',pnm_out,'infodir',pnm_info);
stept=t(2)-t(1);

if flag_twin
   indx=find(t>=twin(1) & t<=twin(2));
end

if flag_1canvas==1
   nplot=flag_acc+flag_vel+flag_dis;
else
   nplot=1;
end

ncol=0;

% -------------------- plot figures ------------------------
% -- plot Ax[yz] --
if flag_acc==1

Ax=gradient(Vx,stept);
Ay=gradient(Vy,stept);
Az=gradient(Vz,stept);

if flag_1canvas==1
   ncol=ncol+1;
else
   ncol=1;
end

if flag_overlap~=1
   figure
end

nrow=1;
subplot(3,nplot,ncol+(nrow-1)*nplot)
hold on;
plot(t,Ax,spec_line);
title('Ax');
if flag_twin
   hold on
   plot(t(indx),Ax(indx),'r');
end

nrow=2;
subplot(3,nplot,ncol+(nrow-1)*nplot)
hold on;
plot(t,Ay,spec_line);
title('Ay');
if flag_twin
   hold on
   plot(t(indx),Ay(indx),'r');
end

nrow=3;
subplot(3,nplot,ncol+(nrow-1)*nplot)
hold on;
plot(t,Az,spec_line);
title('Az');
if flag_twin
   hold on
   plot(t(indx),Az(indx),'r');
end

end

% -- plot Vx[yz] --
if flag_vel==1

if flag_1canvas==1
   ncol=ncol+1;
else
   ncol=1;
end

if flag_overlap~=1 | flag_1canvas~=1
   figure
end

nrow=1;
subplot(3,nplot,ncol+(nrow-1)*nplot)
hold on;
plot(t,Vx,spec_line);
title('Vx');
if flag_twin
   hold on
   plot(t(indx),Vx(indx),'r');
end

nrow=2;
subplot(3,nplot,ncol+(nrow-1)*nplot)
hold on;
plot(t,Vy,spec_line);
title('Vy');
if flag_twin
   hold on
   plot(t(indx),Vy(indx),'r');
end

nrow=3;
subplot(3,nplot,ncol+(nrow-1)*nplot)
hold on;
plot(t,Vz,spec_line);
title('Vz');
if flag_twin
   hold on
   plot(t(indx),Vz(indx),'r');
end

end

% -- plot Ux[yz] --
if flag_dis==1

Ux=cumtrapz(Vx)*stept;
Uy=cumtrapz(Vy)*stept;
Uz=cumtrapz(Vz)*stept;

if flag_1canvas==1
   ncol=ncol+1;
else
   ncol=1;
end

if flag_overlap~=1 | flag_1canvas~=1
   figure
end

nrow=1;
subplot(3,nplot,ncol+(nrow-1)*nplot)
hold on;
plot(t,Ux,spec_line);
title('Ux');
if flag_twin
   hold on
   plot(t(indx),Ux(indx),'r');
end

nrow=2;
subplot(3,nplot,ncol+(nrow-1)*nplot)
hold on;
plot(t,Uy,spec_line);
title('Uy');
if flag_twin
   hold on
   plot(t(indx),Uy(indx),'r');
end

nrow=3;
subplot(3,nplot,ncol+(nrow-1)*nplot)
hold on;
plot(t,Uz,spec_line);
title('Uz');
if flag_twin
   hold on
   plot(t(indx),Uz(indx),'r');
end

end

if 1
Ux=cumtrapz(Vx)*stept;
Uy=cumtrapz(Vy)*stept;
Uz=cumtrapz(Vz)*stept;

subplot(3,1,1)
hold on;
plot(t,Ux/max(max(abs(Ux))),spec_line);
title('Ux');

subplot(3,1,2)
hold on;
plot(t,Uy/max(max(abs(Uy))),spec_line);
title('Uy');

subplot(3,1,3)
hold on;
plot(t,Uz/max(max(abs(Uz))),spec_line);
title('Uz');

end
