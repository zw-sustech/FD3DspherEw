function [V,varargout]=gather_media(snapinfo,varargin)
% gather modelinate
%
% $Date$
% $Revision$
% $LastChangedBy$

pnm_nc='../input/';
varnm='rho';

%-- flags --
nargs=nargin-1;

n=1;
while n<=nargs

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case {'rho','lambda','mu','Qs', ...
         'Px','Py','Pz', ...
         'Uxy','Uxz','Uyz' }
       varnm=varargin{n};
   case {'C11','C12','C13','C14','C15','C16', ...
               'C22','C23','C24','C25','C26', ...
                     'C33','C34','C35','C36', ...
                           'C44','C45','C46', ...
                                 'C55','C56', ...
                                       'C66'}
       varnm=varargin{n};
   case 'mediadir'
       pnm_nc=varargin{n+1}; n=n+1;
   end
end

n=n+1;

end

% check
if ~ exist(pnm_nc,'dir')
   error([mfilename ': directory ' pnm_nc ' does not exist']);
end

nthd=length(snapinfo);
%-- get model data
for n=1:nthd
    n_i=snapinfo(n).thisid(1);n_j=snapinfo(n).thisid(2);n_k=snapinfo(n).thisid(3);
    i1=snapinfo(n).indxs(1);j1=snapinfo(n).indxs(2);k1=snapinfo(n).indxs(3);
    i2=snapinfo(n).indxe(1);j2=snapinfo(n).indxe(2);k2=snapinfo(n).indxe(3);
    subs=snapinfo(n).wsubs;subc=snapinfo(n).wsubc;subt=snapinfo(n).wsubt;
    fnm_nc=get_fnm_media(pnm_nc,n_i,n_j,n_k);
    if ~ exist(fnm_nc,'file')
       error([mfilename ': file ' fnm_nc 'does not exist']);
    end

    V(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_nc,varnm, ...
          fliplr(subs)-1,fliplr(subc),fliplr(subt));
end

V=permute(V,[3,2,1]);

% pack varargout
if nargout>=2, varargout(1)={varnm}; end
