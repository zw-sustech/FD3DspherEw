% retreive seismogram on individual receivers.
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
flag_acc = 0; flag_vel = 1; flag_dis = 0;
flag_normalize = 2; % 0:unnorm; 1:norm together; 2:norm seperate
flag_1canvas = 1;

if 1
OUTPUT_ROOT='../';
%OUTPUT_ROOT='/net/fs01/data/wzhang/Iran.FFTomo/station.ASAO/fx/';
%OUTPUT_ROOT='/net/fs01/data/wzhang/Iran.FFTomo/2005.01.10_18.47/';spec_line='g';
%OUTPUT_ROOT='/net/fs01/data/wzhang/Iran.FFTomo/2005.01.10_18.47.vlow/';spec_line='b';
%OUTPUT_ROOT='/net/fs01/data/wzhang/Iran.FFTomo/2005.01.10_18.47.charac/';spec_line='g';
%OUTPUT_ROOT='/net/fs01/data/wzhang/Iran.FFTomo/2005.01.10_18.47.vlow.high/';spec_line='r';
%OUTPUT_ROOT='/net/fs01/data/wzhang/Iran.FFTomo/2005.01.10_18.47.vzero/';spec_line='r';
%OUTPUT_ROOT='/net/fs01/data/wzhang/Iran.FFTomo/2005.01.10_18.47.vhoc/';spec_line='k';
%OUTPUT_ROOT='/net/fs01/data/wzhang/IRAN/2005.01.10_18.47_all_vext_i1100/';
%OUTPUT_ROOT='/net/fs01/data/wzhang/IRAN/2005.01.10_18.47/';

OUTPUT_ROOT='/net/fs01/data/wzhang/Iran.test/';spec_line='r';

fnm_conf=[OUTPUT_ROOT 'SeisFD3D.conf'];
pnm_out=[OUTPUT_ROOT 'output/'];
%pnm_out=[OUTPUT_ROOT 'output.ppval/'];
pnm_info=[OUTPUT_ROOT 'input/'];

id=0; n=1;
flag_overlap=1;
%spec_line='r';
enlarger=24;
flag_twin=0; twin=[34 40];
flag_filter =1 ; N=4; Fc= 0.1;
end

% -------------------- load data --------------------------
[Vxb,Vyb,Vzb,t,x,y,z] = ...
     retrieve_seismo(id,n,'confname',fnm_conf,'outdir',pnm_out,'infodir',pnm_info);
stept=t(2)-t(1);

if flag_filter==1
   Feff=1/stept/2;
   Wn=Fc/Feff;
   [b,a]=butter(N,Wn);
end

if flag_twin; indx=find(t>=twin(1) & t<=twin(2)); end

if flag_1canvas==1; nplot=flag_acc+flag_vel+flag_dis; else nplot=1; end

ncol=0;

% -------------------- plot figures ------------------------
% -- plot Ax[yz] --
if flag_acc==1
   Axb=gradient(Vxb,stept); Ayb=gradient(Vyb,stept); Azb=gradient(Vzb,stept);
   if flag_filter
      Ax=filter(b,a,Axb);Ay=filter(b,a,Ayb);Az=filter(b,a,Azb);
   else
      Ax=Axb;Ay=Ayb;Az=Azb;
   end

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

   if flag_filter
      %Vx=filter(b,a,Vxb);Vy=filter(b,a,Vyb);Vz=filter(b,a,Vzb);
      Vx=filtfilt(b,a,Vxb);Vy=filtfilt(b,a,Vyb);Vz=filtfilt(b,a,Vzb);
   else
      Vx=Vxb;Vy=Vyb;Vz=Vzb;
   end

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
   if flag_twin; hold on; plot(t(indx),Vy(indx)/Vy0,'r'); end
   
   nrow=2;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Vy/Vy0,spec_line); title('Vy');
   if flag_twin; hold on; plot(t(indx),Vx(indx)/Vx0,'r'); end
   
   nrow=3;
   subplot(3,nplot,ncol+(nrow-1)*nplot)
   hold on; plot(t,Vz/Vz0,spec_line); title('Vz');
   if flag_twin; hold on; plot(t(indx),Vz(indx)/Vz0,'r'); end
end

% -- plot Ux[yz] --
if flag_dis==1
   Uxb=cumtrapz(Vxb)*stept; Uyb=cumtrapz(Vyb)*stept; Uzb=cumtrapz(Vzb)*stept;

   if flag_filter
      Ux=filter(b,a,Uxb);Uy=filter(b,a,Uyb);Uz=filter(b,a,Uzb);
   else
      Ux=Uxb;Uy=Uyb;Uz=Uzb;
   end

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

% -------------------- plot figures ------------------------
if 0
%figure
Ux=cumtrapz(Vx)*stept; Uy=cumtrapz(Vy)*stept; Uz=cumtrapz(Vz)*stept;
%U0=max(max(abs([Ux,Uy,Uz]))); Ux0=U0;Uy0=U0;Uz0=U0;
%Ux0=max(max(abs(Ux]))); Uy0=max(max(abs(Uy]))); Uz0=max(max(abs(Uz])));
%Ux0=1;Uy0=1;Uz0=1;
Ux0=1/20;Uy0=1/20;Uz0=1/20;
subplot(3,1,1)
hold on
plot(t,Ux/Ux0,spec_line)
subplot(3,1,2)
hold on
plot(t,Uy/Uy0,spec_line);
subplot(3,1,3)
hold on
plot(t,Uz/Uz0,spec_line);
end

if 0
% ------------------- save sac -----------------------------
Vsac=sacb(t,Vx); sacw('syn.s001.BHN.SAC',Vsac);
Vsac=sacb(t,Vy); sacw('syn.s001.BHE.SAC',Vsac);
Vsac=sacb(t,Vz); sacw('syn.s001.BHZ.SAC',Vsac);
end
