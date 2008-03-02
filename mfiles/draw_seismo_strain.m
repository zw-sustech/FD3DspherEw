% $Date$
% $Revision$
% $LastChangedBy$

%clear all

% ----------------------- parameter -----------------------
flag_acc = 1; flag_vel = 0; flag_dis = 0;

id=1;
%subs=[1,1,1,1];subc=[-1,-1,-1,-1];subt=[1,1,1,5];
%subs=[1,80,1,1];subc=[-1,1,1,-1];subt=[1,1,1,1];
%subs=[173,1,185,1];subc=[1,-1,1,-1];subt=[1,5,1,1];
%subs=[173,1,185,1];subc=[1,-1,1,1000];subt=[1,5,1,1];
%subs=[60,1,185,1];subc=[7,-1,1,1500];subt=[25,5,1,1];

%subs=[173,61,185,1];subc=[1,8,1,1500];subt=[1,25,1,1];
%subs=[355,221,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];

%subs=[346,471,10,1];subc=[1,1,1,3500];subt=[1,1,1,1]; % 5:5
%subs=[346,521,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 3:7
%subs=[346,571,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 1:9

%enlarger=12;
enlarger=24;
%M=[ 1 0 0; 0 1 0; 0 0 1 ];
%M=src_mech(0,90,0); %001
%M=src_mech(0,90,90); %002
%M=src_mech(0,90,180); %003
%M=src_mech(0,90,270); %004
%M=src_mech(0,90,45); %005
M=src_mech(45,45,45); %010
%M=src_mech(45,45,90); %011

pnm_metric='../input/';
pnm_media ='../input/';

%fnm_conf='../SeisFD3D.conf';
%pnm_metric='../input/';
%pnm_media ='../input/';
%%pnm_outx='../fx/output/';
%pnm_outx='../fx/output.old/';
%%pnm_outy='../fx/output.old/';
%%pnm_outz='../fx/output.old/';
%pnm_outy=pnm_outx;pnm_outz=pnm_outx;

%x0=-15+30*0.04; s1=(10-5-x0)/0.04+1
%subs=[346,346,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 5
%subs=[346,471,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 5
subs=[346,521,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 3
%subs=[346,571,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 1

if 0
fnm_conf='../fy-a100-surf/SeisFD3D.conf';
pnm_outx='../fy-a100-surf/output/';
pnm_outy='../fy-a100-surf/output/';
pnm_outz='../fy-a100-surf/output/';
end

if 1
fnm_conf='../fx-timg/SeisFD3D.conf';
pnm_outx='../fx-timg/output/';
pnm_outy='../fy/output/';
pnm_outz='../fz-timg/output/';
end

if 0
pnm_metric='../../r4km-newmedia/prem.model/';
pnm_media='../../r4km-newmedia/prem.model/';
subs=[126,126,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 3
%subs=[126,176,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 1
fnm_conf='../../r4km-newmedia/fx/SeisFD3D.conf';
pnm_outx='../../r4km-newmedia/fx/output/';
pnm_outy='../../r4km-newmedia/fy/output/';
pnm_outz='../../r4km-newmedia/fz/output/';
end

if 0
pnm_metric='../../r1h2-newmedia/prem.model/';
pnm_media ='../../r1h2-newmedia/prem.model/';
subs=[251,251,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 3
%subs=[251,351,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 1
fnm_conf='../../r1h2-newmedia/fy/SeisFD3D.conf';
pnm_outx='../../r1h2-newmedia/fy/output/';
pnm_outy='../../r1h2-newmedia/fy/output/';
pnm_outz='../../r1h2-newmedia/fy/output/';
end

if 1
pnm_metric='../../r0.5h2-newmedia/prem.model/';
pnm_media ='../../r0.5h2-newmedia/prem.model/';
subs=[251,251,20,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 3
%subs=[251,351,20,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 1
fnm_conf='../../r0.5h2-newmedia/fy/SeisFD3D.conf';
pnm_outx='../../r0.5h2-newmedia/fy/output/';
pnm_outy='../../r0.5h2-newmedia/fy/output/';
pnm_outz='../../r0.5h2-newmedia/fy/output/';
end

if 1
pnm_metric='../prem.model/';
pnm_media ='../prem.model/';
%id=2;subs=[125,43,200,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 20
subs=[1,126,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 20
%subs=[1,251,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 15
%subs=[1,376,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 10
%subs=[1,501,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 5
fnm_conf='../fy/SeisFD3D.conf';
pnm_outx='../fx/output/';
pnm_outy='../fy/output/';
pnm_outz='../fx/output/';
end

%subs=[346,221,10,1];subc=[1,1,1,3500];subt=[1,1,1,1]; % 15:-5
%subs=[346,471,10,1];subc=[1,1,1,3500];subt=[1,1,1,1]; % 5:5
%subs=[346,521,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 3:7
%subs=[346,571,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 1:9
%fnm_conf='../../../r1km/SeisFD3D.conf';
%pnm_metric='../../../r1km/input/';
%pnm_media ='../../../r1km/input/';
%pnm_outx='../../../r1km/output/';
%pnm_outy=pnm_outx;pnm_outz=pnm_outx;

%subs=[346,221,10,1];subc=[1,1,1,3500];subt=[1,1,1,1]; % 15:-5
%%subs=[346,571,10,1];subc=[1,1,1,3500];subt=[1,1,1,1];  % 1:9
%fnm_conf='../../../r1km-newmedia/SeisFD3D.conf';
%pnm_metric='../../../r1km-newmedia/input/';
%pnm_media ='../../../r1km-newmedia/input/';
%pnm_outx='../../../r1km-newmedia/output/';
%pnm_outy=pnm_outx;pnm_outz=pnm_outx;
% -------------------- load data --------------------------
if 1
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
[X,Y,Z]=gather_coord(snapinfo,'metricdir',pnm_metric);

[Ux0,tx]=retrieve_seismo_represent(snapinfo,id,M,'outdir',pnm_outx,'mediadir',pnm_media);
[Uy0,ty]=retrieve_seismo_represent(snapinfo,id,M,'outdir',pnm_outy,'mediadir',pnm_media);
[Uz0,tz]=retrieve_seismo_represent(snapinfo,id,M,'outdir',pnm_outz,'mediadir',pnm_media);

Ux0=permute(Ux0,[4 1 2 3]);
Uy0=permute(Uy0,[4 1 2 3]);
Uz0=permute(Uz0,[4 1 2 3]);

else
load Uxyz0.mat
flag_vel=0;
end

nt=min(min([length(tx), length(ty), length(tz)]));
Ux0(nt+1:end,:,:,:)=[];
Uy0(nt+1:end,:,:,:)=[];
Uz0(nt+1:end,:,:,:)=[];
t=tx(1:nt);

nx=size(Ux0,2);ny=size(Ux0,3);nz=size(Ux0,4);nt=size(Ux0,1);

% -------------------- convolv ------------------------
dt=t(2)-t(1);
%stffc=0.05; stft0=30;
stffc=0.025; stft0=20;
S=stf_ricker(t,stffc,stft0);

if 0
for k=1:nz
for j=1:ny
for i=1:nx
    Ux(:,i,j,k)=convolv_trapz(S,Ux0(:,i,j,k),dt,1,nt);
    Uy(:,i,j,k)=convolv_trapz(S,Uy0(:,i,j,k),dt,1,nt);
    Uz(:,i,j,k)=convolv_trapz(S,Uz0(:,i,j,k),dt,1,nt);
end
end
end
else
Ux=Ux0;Uy=Uy0;Uz=Uz0;
end

%stffc=5.575; stft0=18;
%S1=stf_gauss(t,stffc,stft0);

% -------------------- plot ------------------------
Ux=squeeze(Ux);Uy=squeeze(Uy);Uz=squeeze(Uz);
x=squeeze(X);y=squeeze(Y);z=squeeze(Z);
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
npt=length(theta);

%L=L+10;
%Uy=-Uy;

% -------------------- plot figures ------------------------
% -- plot Vx[yz] --
if flag_vel==1

for n=1:npt
    Vx(:,n)=gradient(Ux(:,n),stept);
    Vy(:,n)=gradient(Uy(:,n),stept);
    Vz(:,n)=gradient(Uz(:,n),stept);
end

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

u0=max(max(abs([Ux,Uy,Uz])));
u0=u0*enlarger;
u1=max(max(abs([Ux])));
u2=max(max(abs([Uy])));
u3=max(max(abs([Uz])));
u1=u1*enlarger; u2=u2*enlarger; u3=u3*enlarger;

figure
subplot(1,3,1)
for n=1:npt
    u1=max(max(abs(Ux(:,n))))*enlarger;
    plot(t,Ux(:,n)/u0*(L(npt)-L(1))+L(n));
    hold on
end
title('Ux');
subplot(1,3,2)
for n=1:npt
    u2=max(max(abs(Uy(:,n))))*enlarger;
    plot(t,Uy(:,n)/u0*(L(npt)-L(1))+L(n));
    hold on
end
title('Uy');
subplot(1,3,3)
for n=1:npt
    u3=max(max(abs(Uz(:,n))))*enlarger;
    plot(t,Uz(:,n)/u0*(L(npt)-L(1))+L(n));
    hold on
end
title('Uz');
end

% x0=-15+30*0.04=-13.8
% x=(([0.1 1 2 3 5 10 15 20]-10)-x0)/0.08
%   48.7500   60.0000   72.5000   85.0000  110.0000  172.5000  235.0000  297.5000

if 0
%figure
n=1;
subplot(3,1,1)
hold on
plot(t,Ux(:,n)/max(max(abs(Ux(:,n)))),'g')
subplot(3,1,2)
hold on
plot(t,Uy(:,n)/max(max(abs(Uy(:,n)))),'g');
subplot(3,1,3)
hold on
plot(t,Uz(:,n)/max(max(abs(Uz(:,n)))),'g');
end

if 0
%figure
u0=max(max(abs([Ux,Uy,Uz])));
u0=1;
u0=1/20;
n=1;
subplot(3,1,1)
hold on
plot(t,Ux(:,n)/u0,'g')
subplot(3,1,2)
hold on
plot(t,Uy(:,n)/u0,'g');
subplot(3,1,3)
hold on
plot(t,Uz(:,n)/u0,'g');
end


if 1
U5=Uy;
%figure
u0=max(max(abs([Ux,Uy,Uz])));
u0=1;
u0=1/20;
n=1;
subplot(4,1,4)
hold on
plot(t,Uy(:,n)/u0,'g');
end

