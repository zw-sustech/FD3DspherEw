function [x,y,z]=gather_coord(snapinfo,varargin)
% gather coordinate
%
% $Date$
% $Revision$
% $LastChangedBy$

pnm_metric='../input/';

%-- flags --
n=1;
while n<=nargin-1

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case 'metricdir'
       pnm_metric=varargin{n+1}; n=n+1;
   end
end

n=n+1;

end

% check
if ~ exist(pnm_metric,'dir')
   error([mfilename ': directory ' pnm_metric ' does not exist']);
end

nthd=length(snapinfo);
%-- get coord data
for n=1:nthd
    n_i=snapinfo(n).thisid(1);n_j=snapinfo(n).thisid(2);n_k=snapinfo(n).thisid(3);
    i1=snapinfo(n).indxs(1);j1=snapinfo(n).indxs(2);k1=snapinfo(n).indxs(3);
    i2=snapinfo(n).indxe(1);j2=snapinfo(n).indxe(2);k2=snapinfo(n).indxe(3);
    subs=snapinfo(n).wsubs;subc=snapinfo(n).wsubc;subt=snapinfo(n).wsubt;
    fnm_metric=get_fnm_metric(pnm_metric,n_i,n_j,n_k);
    if ~ exist(fnm_metric,'file')
       error([mfilename ': file ' fnm_metric 'does not exist']);
    end

    %subs=fliplr(subs);subc=fliplr(subc);subt=fliplr(subt);
    x1(            i1:i2)=nc_varget(fnm_metric,'x',subs(1)-1,subc(1),subt(1));
    y1(      j1:j2      )=nc_varget(fnm_metric,'y',subs(2)-1,subc(2),subt(2));
    z1(k1:k2            )=nc_varget(fnm_metric,'z',subs(3)-1,subc(3),subt(3));
end

nx=length(x1); ny=length(y1); nz=length(z1);
for k=1:nz
for j=1:ny
    x(k,j,1:nx)=x1;
end
end
for k=1:nz
for i=1:nx
    y(k,1:ny,i)=y1;
end
end
for j=1:ny
for i=1:nx
    z(1:nz,j,i)=z1;
end
end

x=permute(x,[3 2 1]);
y=permute(y,[3 2 1]);
z=permute(z,[3 2 1]);

