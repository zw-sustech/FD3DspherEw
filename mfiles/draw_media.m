% draw_model
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
flag_surf=0;
flag_pcolor=1;
flag_overlap = 0;

flag_km=1;
flag_subplot = 2;
flag_print = 0;

id=0;
fnm_conf='../SeisFD3D.conf';
pnm_metric='../input/';
pnm_media ='../input/';

% [x,y,z] order
%subs=[ 1  1  1  ]; subc=[-1  1 -1  ]; subt=[ 2  2  2  ];
%subs=[ 50 50 21  ]; subc=[ 200 150 181]; subt=[ 2  2  1 ];
subs=[ 126 1  1  ]; subc=[ 1   -1 -1]; subt=[ 1   1  1 ];

var_list=[]; scl_caxis_list=[];
var_list{end+1}='rho';    %scl_caxis_list{end+1}=[0,5000];
var_list{end+1}='mu';     %scl_caxis_list{end+1}=[0,1e10];
var_list{end+1}='lambda'; %scl_caxis_list{end+1}=[0,1e10];
var_list{end+1}='Vp';     %scl_caxis_list{end+1}=[0,1e10];
var_list{end+1}='Vs';     %scl_caxis_list{end+1}=[0,1e10];

% -------------------- load data --------------------------

[snapinfo]=locate_snap(fnm_conf,0,'start',subs,'count',subc,'stride',subt);

%-- get coord data
[X,Y,Z]=gather_coord(snapinfo,'metricdir',pnm_metric);
nx=size(X,1);ny=size(X,2);nz=size(X,3);
X=permute(X,[2,1,3]); Y=permute(Y,[2,1,3]); Z=permute(Z,[2,1,3]);
[x,y,z]=sph2cart(X,Y,Z);

str_unit='m';
if flag_km
   x=x/1e3;y=y/1e3;z=z/1e3; str_unit='km';
end

% ------------------------ plot model ----------------------------------
nvar=length(var_list);

for n=1:nvar

if strcmp(var_list{n},'Vp')
   v=( (lambda+2*mu)./rho ).^0.5;
elseif strcmp(var_list{n},'Vp')
   v=( mu./rho ).^0.5;
else
   [V,varnm]=gather_media(snapinfo,var_list{n},'mediadir',pnm_media);
   v=permute(V,[2,1,3]);
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
   sid=surf(squeeze(x),squeeze(y),squeeze(z),squeeze(v));
   %View = [-37.5 70];
else
   if nx==1
      hid=pcolor(squeeze(y),squeeze(z),squeeze(v));
   elseif ny==1
      hid=pcolor(squeeze(x),squeeze(z),squeeze(v));
   else
      hid=pcolor(squeeze(x),squeeze(y),squeeze(v));
   end
   view(-90,90)
end % flag_surf

%rotate(hid,[0 0 1],90)
axis image
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

