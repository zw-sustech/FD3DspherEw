% draw_snap
%
% $Date$
% $Revision$
% $LastChangedBy$

clear all

% ----------------------- parameter -----------------------
flag_km=1;
flag_slice=0; sx=[150]*1e3; sy=[80]*1e3; sz=[0]*1e3;
flag_surf=0; flag_pcolor=1;
flag_print = 0;
flag_avi = 0;

%varnm='Vy'; taut=1.0;
varnm='Tyy'; taut=1.0;
id = 2;
%subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
%subs=[346,1,10,1];subc=[1,-1,1,3500];subt=[1,1,1,1];
%subs=[346,1,1,1];subc=[1,-1,-1,3500];subt=[1,1,1,1];
subs=[126,1,1];subc=[1,-1,-1];subt=[1,1,1];
%n1=1000; n2=5000; dn=50;
n1=1000; n2=1500; dn=50;
%n1=1; n2=5000; dn=1;
%n1=1200; n2=n1; dn=1;
n1=858; n2=n1; dn=1;
%scl_caxis=[-1e-4 1e-4];
scl_caxis=[-1e4 1e4];
%scl_daspect=[10 10 1];
%scl_daspect=[1 1 1];

%fnm_conf='../SeisFD3D.conf';
%pnm_metric='../input/';
%pnm_out='../output/';
fnm_conf='../fy/SeisFD3D.conf';
%pnm_out='../fx/output/';
pnm_out='../fy/output/';
pnm_metric='../prem.model/';
str_unit='m';

% --------------- if using gather_snap_show ---------------
if 0
gather_snap_show('id',id,'n',[n1 n2 dn],vnm, ...
   'light','surf', ...
   'daspect',[1 1 1],'caxis',[-0.1 0.1],'outdir',pnm_out)
return
end

% -------------------- load data --------------------------
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
%-- get coord data
[X,Y,Z]=gather_coord(snapinfo,'metricdir',pnm_metric);
nx=size(X,1);ny=size(X,2);nz=size(X,3);
X=permute(X,[2,1,3]); Y=permute(Y,[2,1,3]); Z=permute(Z,[2,1,3]);
%[x,y,z]=sph2cart(Y,pi/2-X,Z);
[x,y,z]=sph2cart(X,Y,Z);

if nx==1
   indx_pcolor=1;
elseif ny==1
   indx_pcolor=2;
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

% ----------------------- plot figure -----------------------------------
% -- create new window --
if 1
scrsz = get(0,'ScreenSize');
%hid=figure('Position',[50 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
hid=figure('Position',[50 200 1024 468])
%hid=figure;
set(hid,'BackingStore','on');
set(hid,'renderer','zbuffer');
%set(hid,'menubar','none');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperUnits','points')
set(gcf,'PaperPosition',[0 0 1024 768])
end

if flag_avi
   aviid = avifile(['snap_' num2str(id,'%3.3i') '_' varnm '.avi']);
end

% -- time loop

for nlayer=n1:dn:n2

[V,t,ntime,varnm]=gather_snap(snapinfo,id,nlayer,varnm,'outdir',pnm_out);
v=permute(V,[2,1,3]);

disp([ '  draw ' num2str(nlayer) 'th layer(t=' num2str(t) ')']);

if flag_slice==1

sid=slice(x,y,z,v,sx,sy,sz);

elseif flag_surf==1

sid=surf(squeeze(x),squeeze(y),squeeze(z),squeeze(v));

else

if indx_pcolor==1
   hid=pcolor(squeeze(y),squeeze(z),squeeze(v));
elseif indx_pcolor==2
   hid=pcolor(squeeze(x),squeeze(z),squeeze(v));
else
   hid=pcolor(squeeze(x),squeeze(y),squeeze(v));
end
view(-90,90)

end % flag_slice

%rotate(hid,[0 0 1],90)
axis image
%caxis([-1e-2,1e-2])
%shading interp;
shading flat;
if exist('scl_caxis')
   caxis(scl_caxis);
end
if exist('scl_daspect')
   daspect(scl_daspect);
end
if flag_slice==1
   view(36.5,-26)
   camlight
end
if flag_surf==1
   View = [-37.5 70];
   camlight
end
colorbar('vert')

titlestr=['Snapshot of ' varnm ' at ' ...
          '{\fontsize{16}{\bf ' ...
          num2str(double(t),'%07.3f') ...
          '}}s'];
title(titlestr)

drawnow
pause(taut);

if flag_print==1
   fnm_out=[varnm '_ndim',num2str(nlayer,'%5.5i')];
   set(gca,'FontName','FixedWidth');
   print(gcf,'-dpng',[fnm_out '.png']);
end

if flag_avi==1
   F = getframe(gca);
   aviid = addframe(aviid,F);
end

end

if flag_avi==1
   aviid = close(aviid);
end

