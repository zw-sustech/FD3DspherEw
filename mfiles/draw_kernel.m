% draw_snap
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
flag_km=1;
flag_slice=1; sx=[80]*1e3; sy=[80]*1e3; sz=[0]*1e3;
flag_surf=0; flag_pcolor=0;
flag_subplot = 1;
flag_print = 0;
flag_jetwr = 0;

id=1;
subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
fnm_conf='../SeisFD3D.conf';
pnm_metric='../input/';
pnm_out='../tomo/';
str_unit='m';

var_list={'Kapz','Kaqz','Kbpz','Kbqz'};
scl_caxis_list={ ...
      [-1e-11,1e-11], ...
      [-2e-11,2e-11], ...
      [-1e-11,1e-11], ...
      [-2e-11,2e-11]  ...
      %[-1e-16,1e-16], ...  %Kap
      %[-2e-16,2e-16], ...  %Kaq
      %[-1e-16,1e-16], ...  %Kbp
      %[-2e-16,2e-16]  ...  %Kbq
      };

% -------------------- load data --------------------------
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
%-- get coord data
[X,Y,Z]=gather_coord(snapinfo,'metricdir',pnm_metric);
nx=size(X,1);ny=size(X,2);nz=size(X,3);
x=permute(X,[2,1,3]); y=permute(Y,[2,1,3]); z=permute(Z,[2,1,3]);

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

% ----------------------- plot kernel -----------------------------------
nvar=length(var_list);

for n=1:nvar

[V,varnm]=gather_dist(snapinfo,id,var_list{n},'outdir',pnm_out);
v=permute(V,[2,1,3]);

if flag_subplot==0 | n==1
   hid=figure;set(hid,'renderer','zbuffer');
end
if flag_subplot==1
   subplot(2,2,n)
elseif flag_subplot==2
   subplot(4,1,n)
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
%view(-90,90)

end % flag_slice

axis image
shading flat;
if flag_jetwr==1
   colormap('jetwr');
end
if exist('scl_caxis_list')
    caxis(scl_caxis_list{n});
end
if flag_slice==1
   view(36.5,-26)
   %View = [-37.5 70];
   camlight
end
colorbar('vert')

title(var_list{n})

if flag_subplot==0 & flag_print==1
   print(gcf,'-dpng',[var_list{n} '_ak135_ricker05.png']);
end

end % nvar

% -------------------- save figures ------------------------
if flag_subplot>0 & flag_print==1
   print -dpng kall_ak135_ricker05.png
end
