% draw seismogram on inline receivers or recv points
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
flag_acc = 1; flag_vel = 1; flag_dis = 1;

enlarger=6;

%id=1; n1=7;n2=12;dn=1; % hill result
id=1; n1=10;n2=190;dn=10;
%id=1; n1=30;n2=80;dn=10;

fnm_conf='../SeisFD3D.conf';
pnm_info='../input/';
pnm_out='../output/';

% -------------------- load data --------------------------
npt=0;
for n=n1:dn:n2
npt=npt+1;
[Vx(:,npt),Vy(:,npt),Vz(:,npt),t,x(npt),y(npt),z(npt)] ...
    =retrieve_seismo(id,n,'confname',fnm_conf,'outdir',pnm_out,'infodir',pnm_info);
end

Lx=max(max(x))-min(min(x));
Ly=max(max(y))-min(min(y));
Lz=max(max(z))-min(min(z));
if Lz==max(max(Lx,Ly),Lz)
   L=z;
elseif Ly==max(max(Lx,Ly),Lz)
   L=y;
else
   L=x;
end
stept=t(2)-t(1);

% -------------------- plot figures ------------------------
% -- plot Ax[yz] --
if flag_acc==1

for n=1:npt
    Ax(:,n)=gradient(Vx(:,n),stept);
    Ay(:,n)=gradient(Vy(:,n),stept);
    Az(:,n)=gradient(Vz(:,n),stept);
end

a0=max(max(abs([Ax,Ay,Az])));
a0=a0*enlarger;

figure
subplot(1,3,1)
for n=1:npt
    plot(t,Ax(:,n)/a0*(L(npt)-L(1))+L(n));
    hold on
end
title('Ax');
subplot(1,3,2)
for n=1:npt
    plot(t,Ay(:,n)/a0*(L(npt)-L(1))+L(n));
    hold on
end
title('Ay');
subplot(1,3,3)
for n=1:npt
    plot(t,Az(:,n)/a0*(L(npt)-L(1))+L(n));
    hold on
end
title('Az');
end

% -- plot Vx[yz] --
if flag_vel==1

v0=max(max(abs([Vx,Vy,Vz])));
v0=v0*enlarger;

figure
subplot(1,3,1)
for n=1:npt
    plot(t,Vx(:,n)/v0*(L(npt)-L(1))+L(n));
    hold on
end
title('Vx');
subplot(1,3,2)
for n=1:npt
    plot(t,Vy(:,n)/v0*(L(npt)-L(1))+L(n));
    hold on
end
title('Vy');
subplot(1,3,3)
for n=1:npt
    plot(t,Vz(:,n)/v0*(L(npt)-L(1))+L(n));
    hold on
end
title('Vz');
end

% -- plot Ux[yz] --
if flag_dis==1
for n=1:npt
    Ux(:,n)=cumtrapz(Vx(:,n))*stept;
    Uy(:,n)=cumtrapz(Vy(:,n))*stept;
    Uz(:,n)=cumtrapz(Vz(:,n))*stept;
end
u0=max(max(abs([Ux,Uy,Uz])));
u0=u0*enlarger;

figure
subplot(1,3,1)
for n=1:npt
    plot(t,Ux(:,n)/u0*(L(npt)-L(1))+L(n));
    hold on
end
title('Ux');
subplot(1,3,2)
for n=1:npt
    plot(t,Uy(:,n)/u0*(L(npt)-L(1))+L(n));
    hold on
end
title('Uy');
subplot(1,3,3)
for n=1:npt
    plot(t,Uz(:,n)/u0*(L(npt)-L(1))+L(n));
    hold on
end
title('Uz');
end
