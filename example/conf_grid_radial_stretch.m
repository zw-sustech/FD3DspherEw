function check_stretch(fct)
% b=1.5;
% b >=1.0, more close to 1.0, more stretch

%gindx_vmap=[120 150];
%steph= 0.125;

%indx_layer = [150 120 100];
%z_layer = [0 -1.5 -3.355];
%b_layer=[10.0 1.30];

%nk    = 230
%npts = [ 70   30   30   30   30    39 ]
%thic = [ 70e3 48e3 60e3 90e3 120e3 250e3]
%bval = [ 10   1.45  10   1.49  3   1.5]

%npts = [ 40   30   15   30    30    39    15    30]
%thic = [ 40e3 48e3 30e3 90e3  130e3 230e3 110e3 330e3]
%bval = [ 10   1.45  10   1.49 2     1.9   2.5     1.6]

%npts = [ 24     86      10    30    39    10    30]
%thic = [ 24.4e3 195.6e3 35e3  130e3 230e3 75e3 330e3]
%bval = [ 10     1.2     3     2.1   1.9   2.5     1.6]

% for homogeneous model
%npts = [ 25     50    154]
%thic = [ 25e3   130e3  616e3]
%bval = [ 50     1.15   10 ]

% for valid.incdr1km
% npts 230
npts = [ 10   8    15   30    30    107   30]
thic = [ 10e3 10e3 28e3 80e3  105e3 500e3 170e3]
bval = [ 10   1.75 1.8  1.9   2.2   1.9   2.8]

sum(npts)

nlayer=length(npts);
for n=1:nlayer
    b=bval(n);

    nx=npts(n)+1;
    n2=sum(npts(1:n))+1;
    n1=n2-npts(n);

    x=0:nx-1; x=x/(nx-1);
    ymax=thic(n);

    y= ymax*( (b+1)-(b-1)*( (b+1)/(b-1) ).^(1-x) ) ...
         ./( ( (b+1)/(b-1) ).^(1-x) +1 );
    x=x*nx;
    z(n1:n2)=sum(thic(1:n))-thic(n)+y(1:nx);
end

%figure;set(gcf,'renderer','zbuffer')

subplot(1,3,1)
plot(z/1e3,'*');
set(gca,'ydir','reverse');
title('DepDisplay');
print('Dep.png','-dpng');

subplot(1,3,2)
dz=diff(z);
plot(dz,(z(1:end-1)+z(2:end))/2/1e3,'*');
set(gca,'ydir','reverse');
grid on;
title('DepDiff_AvgDep');
print('DepDiff_AvgDep.png','-dpng');


subplot(1,3,3)
%dz=gradient(z);
plot(dz,1:length(dz),'*');
set(gca,'ydir','reverse');
title('DepDiff_Indx');
print('DepDiff_Indx.png','-dpng');

r=6371e3 - fliplr(z);
fid=fopen('SeisGrid.radial.dat','w');
fprintf(fid,'%f\n',r)
fclose(fid)

