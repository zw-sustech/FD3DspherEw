% draw_model
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
flag_km=0;
flag_slice=0;
flag_surf=0; flag_pcolor=1;
flag_subplot = 1;
flag_print = 0;
flag_jetwr = 0;

id=0;
fnm_conf='../SeisFD3D.conf';
pnm_metric='../input/';
pnm_media ='../input/';

%subs=[ 1  1  1  ]; % [x,y,z] order
%subc=[-1  1 -1  ];
%subt=[ 2  2  2  ];
subs=[ 50 50 21  ]; % [x,y,z] order
subc=[ 200 150 181];
subt=[ 2  2  1 ];
sx=[40 80 150]*1e3; sy=[60 100]*1e3; sz=[-50]*1e3;
subs=[ 1 50 1  ]; % [x,y,z] order
subc=[ -1 1 -1];
subt=[ 1  1  1 ];

var_list={'rho'};
%var_list={'rho','mu','lambda'};
%scl_caxis_list={ ...
%      [0,5000], ...
%      [0,1e10], ...
%      [0,1e10]  ...
%      };

% -------------------- load data --------------------------
[snapinfo]=locate_snap(fnm_conf,0,'start',subs,'count',subc,'stride',subt);

%-- get coord data
[X,Y,Z]=gather_coord(snapinfo,'metricdir',pnm_metric);
nx=size(X,1);ny=size(X,2);nz=size(X,3);
X=permute(X,[2,1,3]);
Y=permute(Y,[2,1,3]);
Z=permute(Z,[2,1,3]);
[x,y,z]=sph2cart(Y,pi/2-X,Z);

if nx==1
   indx_pcolor=2;
elseif ny==1
   indx_pcolor=1;
else
   indx_pcolor=3;
end

if flag_km
   x=x/1e3;y=y/1e3;z=z/1e3;
   str_unit='km';
   if exist('sx')
      sx=sx/1e3;sy=sy/1e3;sz=sz/1e3;
   end
end

% ------------------------ plot model ----------------------------------
nvar=length(var_list);

for n=1:nvar

[V,varnm]=gather_media(snapinfo,var_list{n},'mediadir',pnm_media);
v=permute(V,[2,1,3]);

if flag_subplot==0 | n==1
   hid=figure;set(hid,'renderer','zbuffer');
end
if flag_subplot==1
   subplot(nvar,1,n)
elseif flag_subplot==2
   subplot(1,nvar,n)
end

if flag_slice==1

sid=slice(x,y,z,v,sx,sy,sz);

elseif flag_surf==1

sid=surf(squeeze(x),squeeze(y),squeeze(z),squeeze(v))

else

if indx_pcolor==1
   pcolor(squeeze(x),squeeze(z),squeeze(v));
elseif indx_pcolor==2
   pcolor(squeeze(y),squeeze(z),squeeze(v));
else
   pcolor(squeeze(x),squeeze(y),squeeze(v));
end
view(-90,90)

end % flag_slice

axis image
%shading interp;
shading flat;
if flag_jetwr==1
   colormap('jetwr');
end
if exist('scl_caxis_list')
    caxis(scl_caxis_list{n});
end
colorbar('vert')
title(var_list{n})
View = [-37.5 70];
%camlight

if flag_subplot==0 & flag_print==1
   print(gcf,'-dpng',[var_list{n} '_ak135.png']);
end

end % nvar

% -------------------- save figures ------------------------
if flag_subplot>0 & flag_print==1
   print -dpng model_ak135.png
end

