function [n_i,n_j,n_k,npt]=locate_seismo(fnm_conf,id,indx,pnm_info)
% locate snap index in mpi threads
%
% $Date$
% $Revision$
% $LastChangedBy$

% check
if ~ exist(fnm_conf,'file')
   error([mfilename ': file ' fnm_conf ' does not exist']);
end

LenFD=3;
str_cm='#';

fid=fopen(fnm_conf);
conf=textscan(fid,'%s','delimiter','\n','whitespace','');
fclose(fid);
nline=size(conf{1});

for n=1:nline

str=conf{1}{n};
if isempty(str)
   continue
end
npit=findstr(str,str_cm);
if ~ isempty(npit)
   str(npit(1):end)=[];
end
str=regexprep(str,{'=','\|'},' ');
[tag,s]=strtok(str);
if isempty(tag)
   continue
end

switch tag
case 'dims'
    dims=sscanf(s,'%f',3)';
case 'ni'
    ni=sscanf(s,'%f',1);
case 'nj'
    nj=sscanf(s,'%f',1);
case 'nk'
    nk=sscanf(s,'%f',1);
case 'nt'
    nt=sscanf(s,'%f',1);
case 'pnm_grid'
    pnm_grid=sscanf(s,'%s',1);
case 'pnm_media'
    pnm_media=sscanf(s,'%s',1);
case 'pnm_src'
    pnm_src=sscanf(s,'%s',1);
case 'pnm_output'
    pnm_output=sscanf(s,'%s',1);
case 'number_of_recv'
    num_recv=sscanf(s,'%f',1);
case 'number_of_inline'
    num_line=sscanf(s,'%f',1);
case 'number_of_snap'
    num_snap=sscanf(s,'%f',1);
end  % select

end

for n_k=0:dims(3)-1
for n_j=0:dims(2)-1
for n_i=0:dims(1)-1
fnm_seismoinfo=[pnm_info 'seismoinfo_mpi' ...
            num2str(n_i,'%2.2i')...
            num2str(n_j,'%2.2i')...
            num2str(n_k,'%2.2i') ...
            '.nc'];
recvid=nc_varget(fnm_seismoinfo,'id');
npt=find(recvid(:,1)==indx & recvid(:,2)==id);
if npt>0
   return
end
end
end
end

error('this receiver doesn''t exist')

