% retreive seismogram on individual receivers.
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
id=0; n=4;
flag_acc = 0;
flag_vel = 1;
flag_dis = 0;
flag_normalize = 0; % 0:unnorm; 1:norm together; 2:norm seperate
flag_overlap = 1; spec_line='r';
flag_1canvas = 1;

flag_twin=0;
twin=[34 40];
%twin=[34.6 36.7];
%twin=[38.5 41.5];
%twin=[34.6 36.7];

%fnm_conf='../SeisFD3D.conf'; pnm_info='../input/'; pnm_out='../output/';

%fnm_conf='../fy/SeisFD3D.conf'; pnm_info='../fy/input/'; pnm_out='../fy/output/';
fnm_conf='../fx/SeisFD3D.conf'; pnm_info='../fx/input/'; pnm_out='../fx/output/';
%fnm_conf='../../d30-drinc/fy/SeisFD3D.conf'; pnm_info='../../d30-drinc/fy/input/'; pnm_out='../../d30-drinc/fy/output/';
%fnm_conf='../../d30-drinc/fx/SeisFD3D.conf'; pnm_info='../../d30-drinc/fx/input/'; pnm_out='../../d30-drinc/fx/output/';


% -------------------- load data --------------------------
[Vx,Vy,Vz,t,x,y,z] = ...
     retrieve_seismo(id,n,'confname',fnm_conf,'outdir',pnm_out,'infodir',pnm_info);
stept=t(2)-t(1);

if flag_twin; indx=find(t>=twin(1) & t<=twin(2)); end

if flag_1canvas==1; nplot=flag_acc+flag_vel+flag_dis; else nplot=1; end

ncol=0;

% -------------------- plot figures ------------------------
% -- plot Ax[yz] --
if flag_acc==1
   Ax=gradient(Vx,stept); Ay=gradient(Vy,stept); Az=gradient(Vz,stept);
   if flag_normalize==0
      Ax0=1;Ay0=1;Az0=1;
   elseif flag_normalize==1
      Ax0=max(max(abs([Ax Ay Az]))); Ay0=Ax0; Az0=Ax0;
   else
      Ax0=max(max(abs(Ax))); Ay0=max(max(abs(Ay))); Az0=max(max(abs(Az)));
   end
   
   if flag_1canvas==1; ncol=ncol+1; else ncol=1; end
   if flag_overlap~=1; figure; end
   
   nrow=1;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Ax/Ax0,spec_line); title('Ax');
   if flag_twin; hold on; plot(t(indx),Ax(indx)/Ax0,'r'); end
   
   nrow=2;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Ay/Ay0,spec_line); title('Ay');
   if flag_twin; hold on; plot(t(indx),Ay(indx)/Ay0,'r'); end
   
   nrow=3;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Az/Az0,spec_line); title('Az');
   if flag_twin; hold on; plot(t(indx),Az(indx)/Az0,'r'); end
end

% -- plot Vx[yz] --
if flag_vel==1
   if flag_normalize==0
      Vx0=1;Vy0=1;Vz0=1;
   elseif flag_normalize==1
      Vx0=max(max(abs([Vx Vy Vz]))); Vy0=Vx0; Vz0=Vx0;
   else
      Vx0=max(max(abs(Vx))); Vy0=max(max(abs(Vy))); Vz0=max(max(abs(Vz)));
   end
   
   if flag_1canvas==1; ncol=ncol+1; else ncol=1; end
   if flag_overlap~=1 | flag_1canvas~=1; figure; end
   
   nrow=1;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Vx/Vx0,spec_line); title('Vx');
   if flag_twin; hold on; plot(t(indx),Vx(indx)/Vx0,'r'); end
   
   nrow=2;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Vy/Vy0,spec_line); title('Vy');
   if flag_twin; hold on; plot(t(indx),Vy(indx)/Vy0,'r'); end
   
   nrow=3;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Vz/Vz0,spec_line); title('Vz');
   if flag_twin; hold on; plot(t(indx),Vz(indx)/Vz0,'r'); end
end

% -- plot Ux[yz] --
if flag_dis==1
   Ux=cumtrapz(Vx)*stept; Uy=cumtrapz(Vy)*stept; Uz=cumtrapz(Vz)*stept;
   if flag_normalize==0
      Ux0=1;Uy0=1;Uz0=1;
   elseif flag_normalize==1
      Ux0=max(max(abs([Ux Uy Uz]))); Uy0=Ux0; Uz0=Ux0;
   else
      Ux0=max(max(abs(Ux))); Uy0=max(max(abs(Uy))); Uz0=max(max(abs(Uz)));
   end
   
   if flag_1canvas==1; ncol=ncol+1; else ncol=1; end
   if flag_overlap~=1 | flag_1canvas~=1; figure; end
   
   nrow=1;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Ux/Ux0,spec_line); title('Ux');
   if flag_twin; hold on; plot(t(indx),Ux(indx)/Ux0,'r'); end
   
   nrow=2;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Uy/Uy0,spec_line); title('Uy');
   if flag_twin; hold on; plot(t(indx),Uy(indx)/Uy0,'r'); end
   
   nrow=3;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Uz/Uz0,spec_line); title('Uz');
   if flag_twin; hold on; plot(t(indx),Uz(indx)/Uz0,'r'); end
end

