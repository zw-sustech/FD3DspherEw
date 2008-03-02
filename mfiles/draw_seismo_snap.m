% draw seismogram on a speical receiver by using node snap nc files.
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
flag_acc = 0; flag_vel = 0; flag_dis = 1;

varnm='Vy';
id=1;
%subs=[1,1,1,1];subc=[-1,-1,-1,-1];subt=[1,1,1,5];
%subs=[1,80,1,1];subc=[-1,1,1,-1];subt=[1,1,1,1];

subs=[346,221,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];

enlarger=6;

fnm_conf='../SeisFD3D.conf';
pnm_metric='../input/';
pnm_out='../output/';

% -------------------- load data --------------------------
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
[X,Y,Z]=gather_coord(snapinfo,'metricdir',pnm_metric);
x=squeeze(X);y=squeeze(Y);z=squeeze(Z);

[V,t,tindx,varnm]=retrieve_seismo_snap(snapinfo,id,'Vx','outdir',pnm_out);
Vx=squeeze(permute(V,[4 3 2 1]));
[V,t,tindx,varnm]=retrieve_seismo_snap(snapinfo,id,'Vy','outdir',pnm_out);
Vy=squeeze(permute(V,[4 3 2 1]));
[V,t,tindx,varnm]=retrieve_seismo_snap(snapinfo,id,'Vz','outdir',pnm_out);
Vz=squeeze(permute(V,[4 3 2 1]));

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
npt=length(x);

% -------------------- plot figures ------------------------
% -- plot Ax[yz] --
if flag_acc==1

for n=1:npt
    Ax(:,n)=gradient(Vx(:,n),stept);
    Ay(:,n)=gradient(Vy(:,n),stept);
    Az(:,n)=gradient(Vz(:,n),stept);
end

a0=max(max(max(abs([Ax,Ay,Az]))));
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
subplot(3,1,1)
for n=1:npt
    %plot(t,Ux(:,n)/u0*(L(npt)-L(1))+L(n));
    plot(t,Ux(:,n)/max(max(max(abs(Ux(:,n))))));
    hold on
end
title('Ux');
subplot(3,1,2)
for n=1:npt
    plot(t,Uy(:,n)/max(max(max(abs(Uy(:,n))))));
    hold on
end
title('Uy');
subplot(3,1,3)
for n=1:npt
    plot(t,Uz(:,n)/max(max(max(abs(Uz(:,n))))));
    hold on
end
title('Uz');
end
