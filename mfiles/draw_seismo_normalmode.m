clear all;

%fnm_out='S0001.R01.ASC';

%fnm_out='S0000.R03.ASC';
%fnm_out='S0002.R03.ASC';
%fnm_out='S0003.R03.ASC';
%fnm_out='S0004.R03.ASC';
%fnm_out='S0005.R01.ASC';
%fnm_out='S0010.R01.ASC';
%fnm_out='S0011.R01.ASC';

%fnm_out='S0000.R01.ASC';
%fnm_out='S0000.R03.ASC';
%fnm_out='S0001.R03.ASC';
%fnm_out='S0002.R03.ASC';
%fnm_out='S0003.R03.ASC';
%fnm_out='S0004.R03.ASC';
%fnm_out='S0006.R03.ASC';
%fnm_out='S0010.R03.ASC';
%fnm_out='S0011.R03.ASC';

%fnm_out='S0010.R04.ASC'; %5
%fnm_out='S0010.R05.ASC'; %10
fnm_out='S0010.R06.ASC'; %15
%fnm_out='S0010.R07.ASC'; %20

%fnm_out='S0001.R04.ASC';
%fnm_out='S0002.R04.ASC';
%fnm_out='S0003.R04.ASC';
%fnm_out='S0004.R04.ASC';
%fnm_out='S0005.R04.ASC';
%fnm_out='S0006.R04.ASC';
%fnm_out='S0010.R04.ASC';
%fnm_out='S0011.R04.ASC';

%fnm_out='S0011.R06.ASC';

%fnm_out='S0001.R04.ASC';
%fnm_out='S0001.R06.ASC';
%fnm_out='S0000.R03.ASC';
%fnm_out='S0000.R04.ASC';
%fnm_out='S0000.R06.ASC';
%fnm_out='S0012.R04.ASC'; %8km,exp,5km
%fnm_out='S0012.R03.ASC'; %8km,exp,5km
%fnm_out='S0013.R01.ASC'; %8km,0/90/0
%fnm_out='S0013.R03.ASC';
%fnm_out='S0013.R04.ASC';
%fnm_out='S0013.R06.ASC';
%fnm_out='S0014.R01.ASC'; %8km,0/90/45
%fnm_out='S0014.R03.ASC';
%fnm_out='S0014.R04.ASC';
%fnm_out='S0015.R01.ASC'; %8km,0/45/45
%fnm_out='S0015.R03.ASC';
%fnm_out='S0015.R04.ASC';
%fnm_out='S0016.R01.ASC'; %8km,45/45/45
%fnm_out='S0016.R03.ASC';
%fnm_out='S0016.R04.ASC';
%fnm_out='S0017.R01.ASC'; %8km,45/45/90
%fnm_out='S0017.R03.ASC';
%fnm_out='S0017.R04.ASC';

fnm_out=['ELA/' fnm_out];

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

if 0
figure;
subplot(3,1,1)
plot(t,Ux)
subplot(3,1,2)
plot(t,Uy);
subplot(3,1,3)
plot(t,Uz);
end

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
u0=max(max(abs([Ux,Uy,Uz])));
u0=1;
%u0=1/8;
subplot(3,1,1)
hold on
plot(t,-Ux/u0,'r')
subplot(3,1,2)
hold on
plot(t,Uy/u0,'r');
subplot(3,1,3)
hold on
plot(t,Uz/u0,'r');
end

if 0
%figure
u0=max(max(abs([Ux,Uy,Uz])));
u0=1;
subplot(4,1,4)
hold on
plot(t,Uy/u0,'r');
end

