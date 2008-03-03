clear all;

src='0010'; % 45/45/45
%src='0011'; % 45/45/90
%src='0012'; % 8km,exp
%src='0013'; % 8km,0/90/0
%src='0014'; % 8km,0/90/45
%src='0015'; % 8km,0/45/45
%src='0016'; % 8km,45/45/45
%src='0017'; % 8km,45/45/90

rcv='03'; % 5deg
rcv='04'; % 10deg

spec_line='k';

fnm_out=['ELA/S' src '.R' rcv '.ASC'];

% ------------------ load data ------------------------
fid=fopen(fnm_out);
srcnm=fgetl(fid);
srcloc=fscanf(fid,'%f',[3]); str=fgetl(fid);
revnm=fgetl(fid);
revloc=fscanf(fid,'%f',[3]); str=fgetl(fid);
str=fgetl(fid);
str=fgetl(fid);
U=fscanf(fid,'%f',[4,inf]);
U=U';

U(1,:)=[];
t=U(:,1);
Uz0=U(:,2);
Uy0=U(:,3);
Ux0=U(:,4);

nt=size(U,1);
dt=t(2)-t(1);

% ------------------ convolv stf ------------------------
S0=stf_gauss(t,5.575,18);D0=stf_gaussderiv(t,5.575,18);
S1=stf_ricker(t,0.05,30);D1=stf_rickerderiv(t,0.05,30);

Ux=convolv_trapz(D1,Ux0,dt,1,nt);
Uy=convolv_trapz(D1,Uy0,dt,1,nt);
Uz=convolv_trapz(D1,Uz0,dt,1,nt);
%Ux=convolv_trapz(S1,Ux ,dt,1,nt);
%Uy=convolv_trapz(S1,Uy ,dt,1,nt);
%Uz=convolv_trapz(S1,Uz ,dt,1,nt);

%S1=stf_ricker(t,0.025,20);D1=stf_rickerderiv(t,0.05,30);
%Ux=convolv_trapz(S1,Ux ,dt,1,nt);
%Uy=convolv_trapz(S1,Uy ,dt,1,nt);
%Uz=convolv_trapz(S1,Uz ,dt,1,nt);

%Ux=Ux0;Uy=Uy0;Uz=Uz0;

% ------------------ plot figure ------------------------
if 0
subplot(3,1,1)
hold on
plot(t,-Ux/max(max(abs(Ux))),'r')
subplot(3,1,2)
hold on
plot(t,Uy/max(max(abs(Uy))),'r');
subplot(3,1,3)
hold on
plot(t,Uz/max(max(abs(Uz))),'r');
end

if 1
%figure
%U0=max(max(abs([Ux,Uy,Uz]))); Ux0=U0;Uy0=U0;Uz0=U0;
%Ux0=max(max(abs(Ux]))); Uy0=max(max(abs(Uy]))); Uz0=max(max(abs(Uz])));
Ux0=1;Uy0=1;Uz0=1;
subplot(3,1,1)
hold on; plot(t,-Ux/Ux0,spec_line)
subplot(3,1,2)
hold on; plot(t,Uy/Uy0,spec_line);
subplot(3,1,3)
hold on; plot(t,Uz/Uz0,spec_line);

subplot(3,1,1)
title(fnm_out)
end

