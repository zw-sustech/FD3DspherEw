% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
flag_vel = 0; flag_dis = 0;

%M=[ 1 0 0; 0 1 0; 0 0 1 ];
%M=src_mech(0,90,0); %001
%M=src_mech(0,90,90); %002
%M=src_mech(0,90,180); %003
%M=src_mech(0,90,270); %004
%M=src_mech(0,90,45); %005
%M=src_mech(45,45,45); %010
M=src_mech(45,45,90); %011

if 1
MODEL_ROOT='../'; pnm_metric=[MODEL_ROOT 'prem.model/']; pnm_media=[MODEL_ROOT 'prem.model/'];

OUTPUT_ROOT='../';
fnm_conf=[OUTPUT_ROOT 'fy-abs/SeisFD3D.conf'];
pnm_outx=[OUTPUT_ROOT 'fx-abs/output/'];pnm_outy=[OUTPUT_ROOT 'fy-abs/output/'];pnm_outz=[OUTPUT_ROOT 'fz-abs/output/'];
%id=2;subs=[125,43,200,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 20
id=1;
%subs=[1,126,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 20
subs=[1,251,10,1],subc=[1,1,1,10000];subt=[1,1,1,1];  % 15
%subs=[1,376,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 10
%subs=[1,501,10,1];subc=[1,1,1,10000];subt=[1,1,1,1];  % 5
spec_line='r'; enlarger=24;
end

% -------------------- load data --------------------------
if 1
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
[X,Y,Z]=gather_coord(snapinfo,'metricdir',pnm_metric);

[Ux0,tx]=retrieve_seismo_sgt(snapinfo,id,M,'outdir',pnm_outx,'mediadir',pnm_media);
[Uy0,ty]=retrieve_seismo_sgt(snapinfo,id,M,'outdir',pnm_outy,'mediadir',pnm_media);
[Uz0,tz]=retrieve_seismo_sgt(snapinfo,id,M,'outdir',pnm_outz,'mediadir',pnm_media);

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
if 1
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

