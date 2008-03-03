function gather_snap_show(varargin)
% draw seismic wave snapshot
% Usage:
%   gather_snap_show('id',1,'n',[1 10 1],'delay',0.5,'Vz','topointerp');
%
% $Date$
% $Revision$
% $LastChangedBy$

disp(' ')

%-- parameters --
id = 1;
taut=0.0;
flag_nodisplay=0;
flag_colorbar=1;
flag_light=0;
flag_savefig=0;
flag_avi=0;
canv_size=[0 0 1024 768];

fnm_conf='../SeisFD3D.conf';

%-- flags --
n=1;
while n<=nargin

switch varargin{n}
case 'id'
    id=varargin{n+1};n=n+1;
case 'n'
    nbound=varargin{n+1}; n=n+1;
case 'confname'
    fnm_conf=varargin{n+1}; n=n+1;
case 'delay'
    taut=varargin{n+1}; n=n+1;
case 'noMenu'
    flag_noMenu = 1; disp('set flag_noMenu ...');
case 'slice'
    flag_slice = 1; disp('set flag_slice ...');
case 'surf'
    flag_surf = 1; disp('set flag_surf ...');
case 'flat'
    flag_flat = 1; disp('set flag_flat ...');
case 'interp'
    flag_interp = 1; disp('set flag_interp ...');
case 'caxis'
    scl_caxis=varargin{n+1}; n=n+1;
case 'daspect'
    scl_daspect=varargin{n+1}; n=n+1;
case 'light'
    flag_light=1;
case 'nocolorbar'
    flag_colorbar=0;
case 'nodisplay'
    flag_nodisplay=1;
case 'jpg'
    flag_savefig=1;
    fig_drv='-djpeg';
    fig_sub='jpg';
case 'png'
    flag_savefig=1;
    fig_drv='-dpng';
    fig_sub='png';
case 'eps'
    flag_savefig=1;
    fig_drv='-depsc2';
    fig_sub='eps';
case 'avi'
    flag_avi=1;
case 'keyboard'
    flag_keyboard = 1; disp('set flag_keyboard ...');
end
n=n+1;

end

% -- set prameters --
dn=1;
if length(nbound)==3
   dn=nbound(3);
end

[snapinfo]=locate_snap(fnm_conf,id,varargin{:});

% -- create new window --
if ~ flag_nodisplay
hid=figure;
set(hid,'BackingStore','on');
set(hid,'renderer','zbuffer');

if exist('flag_noMenu') & flag_noMenu==1
   set(hid,'menubar','none');
end
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperUnits','points')
set(gcf,'PaperPosition',canv_size)
end

if flag_avi
   aviid = avifile(['snap_' num2str(id,'%3.3i') '_' varnm '.avi']);
end

%-- get coord data
[x,y,z]=gather_coord(snapinfo,varargin{:});
nx=size(x,1);ny=size(x,2);nz=size(x,3);
x=permute(x,[2,1,3]); y=permute(y,[2,1,3]); z=permute(z,[2,1,3]);
[x,y,z]=sph2cart(x,y,z);

% -- time loop
for nlayer=nbound(1):dn:nbound(2)
    % window
    if flag_nodisplay
       close;
       hid=figure('pos',canv_size);
       %set(hid,'BackingStore','on');
       set(hid,'renderer','zbuffer');
       set(hid, 'PaperPositionMode', 'manual');
       set(hid,'PaperUnits','points')
       set(hid,'PaperPosition',canv_size)
    end

    [v,t,varnm]=gather_snap(snapinfo,id,nlayer,varargin{:});
    v=permute(v,[2,1,3]);

% -- plot --
    if exist('flag_slice') & flag_slice==1
       sid=slice(x,y,z,v,sx,sy,sz);
    elseif exist('flag_surf') & flag_surf==1
       sid=surf(squeeze(x),squeeze(y),squeeze(z),squeeze(v));
    else
       if nx==1
          sid=pcolor(squeeze(y),squeeze(z),squeeze(v));
       elseif ny==1
          sid=pcolor(squeeze(x),squeeze(z),squeeze(v));
       else
          sid=pcolor(squeeze(x),squeeze(y),squeeze(v));
       end
    end

    % set propeties
    if exist('flag_flat') & flag_flat==1
       shading flat;
    end
    if exist('flag_interp') & flag_interp==1
       shading interp;
    end
    if exist('scl_daspect')
       daspect(scl_daspect);
    end
    if exist('scl_caxis')
       caxis(scl_caxis);
    end
    if flag_colorbar
       cid=colorbar('vert');
       %cid=colorbar('horizontal','position',[0.13 0.627 0.1 0.03])
       %set(cid,'xticklabel',{' ','m/s','1'})
    end
    if flag_light==1
       axis tight
       %camlight
       %lightangle(-45,45,'infinite') 
       lightangle(-45,45) 
       %view(0,90)
    end

    titlestr=['Snapshot of ' varnm ' at ' ...
              '{\fontsize{16}{\bf ' ...
              num2str(double(t),'%07.3f') ...
              '}}s'];
    title(titlestr)

    drawnow
    pause(taut);

if flag_savefig
if flag_nodisplay
   fnm_out=[varnm '_ndim',num2str(nlayer,'%5.5i')];
   xt1=get(gca,'xtick'); yt1=get(gca,'ytick'); zt1=get(gca,'ztick');
   set(gca,'ztick',[],'zminortick','off');
   set(gca,'xtick',[]); set(gca,'ytick',[]);
   set(cid,'visible','off');
   set(gca,'FontName','FixedWidth');
   img=imcapture(gcf,'all');
   imwrite(img,[fnm_out '_skel.' fig_sub],fig_sub);

   set(gca,'xtick',xt1,'xgrid','off','xminortick','off','xminorgrid','off');
   set(gca,'ytick',yt1,'ygrid','off','yminortick','off','yminorgrid','off');
   set(sid,'visible','off');
   set(cid,'visible','on');
   print(gcf,fig_drv,[fnm_out '_mask.' fig_sub]);
else
   fnm_out=[varnm '_ndim',num2str(nlayer,'%5.5i')];
   set(gca,'FontName','FixedWidth');
   print(gcf,fig_drv,[fnm_out '.' fig_sub]);
end
end

if flag_avi
   F = getframe(gca);
   aviid = addframe(aviid,F);
end

end % loop_ntime

if flag_avi
   aviid = close(aviid);
end

disp(' ')

if exist('flag_keyboard')
   disp('press return to exit the gather_snap_show function');
   keyboard
end

