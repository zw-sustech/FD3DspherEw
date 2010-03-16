function set_mfiles_path
%
% set_mfiles_path: Set the path of mfiles function.
%
% Usage: set_mfiles_path
%
% Input:
%   none
%
% Output:
%   none
%
% Example:
%   >> set_mfiles_path()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date$
% $Revision$
% $LastChangedBy$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MFILE_ROOT='../../../mfiles'

path([MFILE_ROOT '/saclab'],path);
path([MFILE_ROOT '/fun-spool'],path);
%path([MFILE_ROOT '/fileexchange'],path);
%path([MFILE_ROOT '/others'],path);
