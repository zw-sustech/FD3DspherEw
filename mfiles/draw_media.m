% draw_model
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

path([pwd '/fun-spool'],path);
path([pwd '/saclab'],path);

% ----------------------- parameter -----------------------
flag_surf=0;
flag_pcolor=1;
flag_overlap = 0;

flag_km=1;
flag_subplot = 0;
flag_print = 0;

id=0;
%fnm_conf='../SeisFD3D.conf';
%pnm_metric='../input/'; pnm_media ='../input/';
%pnm_metric='../../prem.model/'; pnm_media ='../../prem.model/';
fnm_conf='/net/fs01/data/wzhang/Iran.FFTomo/2005.01.10_18.47/SeisFD3D.conf';
pnm_metric='/net/fs01/data/wzhang/Iran.FFTomo/model.iran.volume/';
pnm_media='/net/fs01/data/wzhang/Iran.FFTomo/model.iran.volume/';
%pnm_metric='/net/fs01/data/wzhang/IRAN/test_media_verpoly/input/';
%pnm_media='/net/fs01/data/wzhang/IRAN/test_media_verpoly/input/';

%vlon=39.00:0.02:65.98; nlon=1350
%vlat=22.00:0.02:43.98; nlat=1100

%lon=38:0.03:67.37E, 980
%lat=20.53:0.03:44.5N, 800

fnm_conf='/net/fs01/data/wzhang/Iran.FFTomo/test/2005.05.01.18.58.colat/SeisFD3D.conf';
pnm_metric='/net/fs01/data/wzhang/Iran.FFTomo/test/2005.05.01.18.58.colat/input/';
pnm_media='/net/fs01/data/wzhang/Iran.FFTomo/test/2005.05.01.18.58.colat/input/';

% [x,y,z] order
%subs=[ 1  1  1  ]; subc=[-1  1 -1  ]; subt=[ 2  2  2  ];
%subs=[ 50 50 21  ]; subc=[ 200 150 181]; subt=[ 2  2  1 ];
%subs=[ 991 1  1  ]; subc=[ 1   -1 -1]; subt=[ 1 1 1 ];
%subs=[ 351 1  1  ]; subc=[ 1   -1 -1]; subt=[ 1 1 1 ];
%subs=[ 1 1101  1  ]; subc=[ -1  1 -1]; subt=[ 1 1 1 ];
%subs=[ 1 100  1  ]; subc=[ -1 1 -1]; subt=[ 1 1 1 ];
subs=[ 1 1  250  ]; subc=[ -1 -1 1]; subt=[ 1   1  1 ];

%subs=[ 851 1  1  ]; subc=[ 1   -1 -1]; subt=[ 1 1 1 ];  % 39

%--------------------------------------------------------------------
% validation: 3D
%flag_overlap = 1;
%subs=[ ceil((29.5-22)/0.02)+1 1  1  ]; subc=[ 1   -1 -1]; subt=[ 1 1 1 ];  % 29.5
%subs=[ ceil((37.5-22)/0.02)+1 1  1  ]; subc=[ 1   -1 -1]; subt=[ 1 1 1 ];  % 37.5
%subs=[ 1 ceil((50-39)/0.02)+1 1]; subc=[-1 1 -1]; subt=[ 1 1 1 ];  % 50
%subs=[ 1 ceil((57-39)/0.02)+1 1]; subc=[-1 1 -1]; subt=[ 1 1 1 ];  % 57
%xlabel('longitude'); ylabel('latitude'); zlabel('radius');
%shading interp; axis tight;
%print -dpng iran.3D.Vs.png

%--------------------------------------------------------------------
% validation: sliced
%subs=[ 1 1  250  ]; subc=[ -1 -1 1]; subt=[ 1   1  1 ]; %surface
%subs=[ 1 1  246  ]; subc=[ -1 -1 1]; subt=[ 1   1  1 ]; %4km
%subs=[ 1 1  244  ]; subc=[ -1 -1 1]; subt=[ 1   1  1 ]; %6km
%subs=[ 1 1  242  ]; subc=[ -1 -1 1]; subt=[ 1   1  1 ]; %8km
%subs=[ round((34-22)/0.02)+1 1  1  ]; subc=[ 1   -1 -1]; subt=[ 1 1 1 ];  % 34N: 601 ,449
%subs=[ 1 round((57.02-39)/0.02)+1 1]; subc=[-1 1 -1]; subt=[ 1 1 1 ];  % 57.02E: 902, 635
%xlabel('longitude'); ylabel('latitude');
%shading interp; axis tight;
%print -dpng iran.3D.Vs.png

%--------------------------------------------------------------------
% validation: well 34N, 57.02
%subs=[ 1 round((57.02-39)/0.02)+1 1]; subc=[-1 1 -1]; subt=[ 1 1 1 ];  % 57.02E: 902, 635
%plot(squeeze(v(601,1,:)),6371-squeeze(Z(601,1,:)))
%set(gca,'ydir','reverse');

%--------------------------------------------------------------------
var_list=[]; scl_caxis_list=[];
%var_list{end+1}='rho';    %scl_caxis_list{end+1}=[0,5000];
%var_list{end+1}='mu';     %scl_caxis_list{end+1}=[0,1e10];
%var_list{end+1}='lambda'; %scl_caxis_list{end+1}=[0,1e10];
var_list{end+1}='Vp';     %scl_caxis_list{end+1}=[0,1e10];
%var_list{end+1}='Vs';     %scl_caxis_list{end+1}=[0,1e10];

% -------------------- load data --------------------------

[snapinfo]=locate_snap(fnm_conf,0,'start',subs,'count',subc,'stride',subt);

%-- get coord data
[X,Y,Z]=gather_coord(snapinfo,'coorddir',pnm_metric);
nx=size(X,1);ny=size(X,2);nz=size(X,3);
%X=permute(X,[2,1,3]); Y=permute(Y,[2,1,3]); Z=permute(Z,[2,1,3]);
%[x,y,z]=sph2cart(X,Y,Z);
[x,y,z]=sph2cart(Y,X-pi/2,Z);

str_unit='m';
if flag_km
   x=x/1e3;y=y/1e3;z=z/1e3;Z=Z/1e3; str_unit='km';
end

% ------------------------ plot model ----------------------------------
nvar=length(var_list);

for n=1:nvar

if (strcmp(var_list{n},'Vp') | strcmp(var_list{n},'Vs')) & (~exist('rho'))
   [rho,varnm]=gather_media(snapinfo,'rho','mediadir',pnm_media);
end
if (strcmp(var_list{n},'Vp') | strcmp(var_list{n},'Vs')) & (~exist('mu'))
   [mu,varnm]=gather_media(snapinfo,'mu','mediadir',pnm_media);
end
if strcmp(var_list{n},'Vp') & (~exist('lambda'))
   [lambda,varnm]=gather_media(snapinfo,'lambda','mediadir',pnm_media);
end

if strcmp(var_list{n},'Vp')
   v=( (lambda+2*mu)./rho ).^0.5;
elseif strcmp(var_list{n},'Vs')
   v=( mu./rho ).^0.5;
else
   [v,varnm]=gather_media(snapinfo,var_list{n},'mediadir',pnm_media);
   %v=permute(V,[2,1,3]);
   if strcmp(var_list{n},'rho')
      rho=v;
   elseif strcmp(var_list{n},'mu')
      mu=v;
   elseif strcmp(var_list{n},'lambda')
      lambda=v;
   end
end

if flag_overlap==1
   hold on
elseif flag_subplot==0 | n==1
   hid=figure;set(hid,'renderer','zbuffer');
end
if flag_subplot==1
   subplot(nvar,1,n)
elseif flag_subplot==2
   subplot(1,nvar,n)
end

if flag_surf==1
   %sid=surf(squeeze(x),squeeze(y),squeeze(z),squeeze(v));
   sid=surf(squeeze(Y)/pi*180,squeeze(X)/pi*180-90,squeeze(Z),squeeze(v));
   %View = [-37.5 70];
else
   if nx==1
      %hid=pcolor(squeeze(y),squeeze(z),squeeze(v));
      hid=pcolor(squeeze(Y/pi*180),squeeze(Z),squeeze(v));
   elseif ny==1
      %[x,y,z]=sph2cart(zeros(size(Y)),X-pi/2,Z);
      %hid=pcolor(squeeze(x),squeeze(z),squeeze(v));
      %hid=pcolor(squeeze(X/pi*180-90),squeeze(Z),squeeze(v));
      hid=pcolor(squeeze(X/pi*180),squeeze(Z),squeeze(v));
   else
      %hid=pcolor(squeeze(x),squeeze(y),squeeze(v));
      %hid=pcolor(squeeze(Y/pi*180),squeeze(X/pi*180-90),squeeze(v));
      hid=pcolor(squeeze(Y/pi*180),squeeze(X/pi*180),squeeze(v));
   end
   %view(-90,90)
end % flag_surf

%rotate(hid,[0 0 1],90)
%axis image
%shading interp;
shading flat;
if exist('scl_caxis_list') & ~ isempty(scl_caxis_list)
    caxis(scl_caxis_list{n});
end
colorbar('vert')
title(var_list{n})
%camlight

if flag_subplot==0 & flag_print==1
   print(gcf,'-dpng',[var_list{n} '_ak135.png']);
end

end % nvar

% -------------------- save figures ------------------------
if flag_subplot>0 & flag_print==1
   print -dpng model_ak135.png
end

% plot(squeeze(X(:,1,248))/pi*180,squeeze(v(:,1,248)),'k')
% plot(squeeze(X(:,1,248))/pi*180-90,squeeze(v(:,1,248)),'k')
% plot(squeeze(X(:,1,247))/pi*180-90,squeeze(v(:,1,247)),'k')
% plot(squeeze(X(:,1,247))/pi*180-90,squeeze(v(:,1,247)),'k')
% plot(squeeze(X(:,1,246))/pi*180-90,squeeze(v(:,1,246)),'k')
% draw_media
% plot(squeeze(X(:,1,248))/pi*180-90,squeeze(v(:,1,248)),'k')
% hold on
% plot(squeeze(X(:,1,247))/pi*180-90,squeeze(v(:,1,247)),'b')
% plot(squeeze(X(:,1,246))/pi*180-90,squeeze(v(:,1,246)),'r')

