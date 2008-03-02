function [Vx,Vy,Vz,t,x,y,z]=retrieve_seismo(id,indx,varargin)
% retreive seismogram on individual receivers.
% id=0 means recv type, otherwise line type
%
% $Date$
% $Revision$
% $LastChangedBy$

%pnm_out='../wave.src1/output/';
pnm_out='../output/';
pnm_info='../input/';
fnm_conf='../SeisFD3D.conf';

%-- flags --
n=1;
while n<=nargin-2

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case 'outdir'
       pnm_out=varargin{n+1}; n=n+1;
   case 'infodir'
       pnm_info=varargin{n+1}; n=n+1;
   case 'confname'
       fnm_conf=varargin{n+1}; n=n+1;
   end
end
n=n+1;

end

% check
if ~ exist(fnm_conf,'file')
   error([mfilename ': file ' fnm_conf ' does not exist']);
end
if ~ exist(pnm_out,'dir')
   error([mfilename ': directory ' pnm_out ' does not exist']);
end
if ~ exist(pnm_info,'dir')
   error([mfilename ': directory ' pnm_info ' does not exist']);
end

[n_i,n_j,n_k,npt]=locate_seismo(fnm_conf,id,indx,pnm_info);

fnm_seismo=[pnm_out 'seismo_mpi' ...
            num2str(n_i,'%2.2i')...
            num2str(n_j,'%2.2i')...
            num2str(n_k,'%2.2i') ...
            '.nc'];
fnm_seismoinfo=[pnm_info 'seismoinfo_mpi' ...
            num2str(n_i,'%2.2i')...
            num2str(n_j,'%2.2i')...
            num2str(n_k,'%2.2i') ...
            '.nc'];

x=nc_varget(fnm_seismoinfo,'coord',[npt-1,0],[1 1]);
y=nc_varget(fnm_seismoinfo,'coord',[npt-1,1],[1 1]);
z=nc_varget(fnm_seismoinfo,'coord',[npt-1,2],[1 1]);

t=nc_varget(fnm_seismo,'time');

Vx=nc_varget(fnm_seismo,'Vx',[0 npt-1],[-1,1]);
Vy=nc_varget(fnm_seismo,'Vy',[0 npt-1],[-1,1]);
Vz=nc_varget(fnm_seismo,'Vz',[0 npt-1],[-1,1]);
