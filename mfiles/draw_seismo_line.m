% draw seismogram on inline receivers or recv points
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
flag_acc = 0; flag_vel = 1; flag_dis = 0;

if 1
OUTPUT_ROOT='../';
fnm_conf=[OUTPUT_ROOT 'SeisFD3D.conf'];
pnm_info=[OUTPUT_ROOT 'input/'];pnm_out=[OUTPUT_ROOT 'output/'];
id=1; n1=10;n2=190;dn=10;
spec_line='r'; enlarger=24;
end

% -------------------- load data --------------------------
npt=0;
for n=n1:dn:n2
npt=npt+1;
[Vx(:,npt),Vy(:,npt),Vz(:,npt),t,x(npt),y(npt),z(npt)] ...
    =retrieve_seismo(id,n,'confname',fnm_conf,'outdir',pnm_out,'infodir',pnm_info);
end

theta=x/pi*180; phi=y/pi*180; r=z;
Lx=max(max(theta))-min(min(theta));
Ly=max(max(phi))-min(min(phi));
Lz=max(max(r))-min(min(r));
if Lz==max(max(Lx,Ly),Lz)
   L=r;
elseif Ly==max(max(Lx,Ly),Lz)
   L=phi;
else
   L=theta;
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

A0=max(max(abs([Ax,Ay,Az]))); Ax0=A0;Ay0=A0;Az0=A0;
%Ax0=max(max(abs(Ax]))); Ay0=max(max(abs(Ay]))); Az0=max(max(abs(Az])));
Ax0=Ax0*enlarger; Ay0=Ay0*enlarger; Az0=Az0*enlarger;
%Ax0=1.0;Ay0=1.0;Az0=1.0;

figure
subplot(1,3,1)
for n=1:npt
    hold on
    plot(t,Ax(:,n)/Ax0*(L(npt)-L(1))+L(n),spec_line);
end
title('Ax');
subplot(1,3,2)
for n=1:npt
    hold on
    plot(t,Ay(:,n)/Ax0*(L(npt)-L(1))+L(n),spec_line);
end
title('Ay');
subplot(1,3,3)
for n=1:npt
    hold on
    plot(t,Az(:,n)/Ax0*(L(npt)-L(1))+L(n),spec_line);
end
title('Az');
end

% -- plot Vx[yz] --
if flag_vel==1

V0=max(max(abs([Vx,Vy,Vz]))); Vx0=V0;Vy0=V0;Vz0=V0;
%Vx0=max(max(abs(Vx]))); Vy0=max(max(abs(Vy]))); Vz0=max(max(abs(Vz])));
Vx0=Vx0*enlarger; Vy0=Vy0*enlarger; Vz0=Vz0*enlarger;
%Vx0=1.0;Vy0=1.0;Vz0=1.0;

figure
subplot(1,3,1)
for n=1:npt
    hold on
    plot(t,Vx(:,n)/Vx0*(L(npt)-L(1))+L(n),spec_line);
end
title('Vx');
subplot(1,3,2)
for n=1:npt
    hold on
    plot(t,Vy(:,n)/Vy0*(L(npt)-L(1))+L(n),spec_line);
end
title('Vy');
subplot(1,3,3)
for n=1:npt
    hold on
    plot(t,Vz(:,n)/Vz0*(L(npt)-L(1))+L(n),spec_line);
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
U0=max(max(abs([Ux,Uy,Uz]))); Ux0=U0;Uy0=U0;Uz0=U0;
%Ux0=max(max(abs(Ux]))); Uy0=max(max(abs(Uy]))); Uz0=max(max(abs(Uz])));
Ux0=Ux0*enlarger; Uy0=Uy0*enlarger; Uz0=Uz0*enlarger;
%Ux0=1;Uy0=1;Uz0=1;

figure
subplot(1,3,1)
for n=1:npt
    hold on
    plot(t,Ux(:,n)/Ux0*(L(npt)-L(1))+L(n),spec_line);
end
title('Ux');
subplot(1,3,2)
for n=1:npt
    hold on
    plot(t,Uy(:,n)/Uy0*(L(npt)-L(1))+L(n),spec_line);
end
title('Uy');
subplot(1,3,3)
for n=1:npt
    hold on
    plot(t,Uz(:,n)/Uz0*(L(npt)-L(1))+L(n),spec_line);
end
title('Uz');
end

% -------------------- plot figures ------------------------
if 0
for n=1:npt
    Ux(:,n)=cumtrapz(Vx(:,n))*stept;
    Uy(:,n)=cumtrapz(Vy(:,n))*stept;
    Uz(:,n)=cumtrapz(Vz(:,n))*stept;
end
%figure
%U0=max(max(abs([Ux,Uy,Uz]))); Ux0=U0;Uy0=U0;Uz0=U0;
%Ux0=max(max(abs(Ux]))); Uy0=max(max(abs(Uy]))); Uz0=max(max(abs(Uz])));
%Ux0=1;Uy0=1;Uz0=1;
Ux0=1/20;Uy0=1/20;Uz0=1/20;
n=1;
subplot(3,1,1)
hold on
plot(t,Ux(:,n)/Ux0,spec_line)
subplot(3,1,2)
hold on
plot(t,Uy(:,n)/Uy0,spec_line);
subplot(3,1,3)
hold on
plot(t,Uz(:,n)/Uz0,spec_line);
end
