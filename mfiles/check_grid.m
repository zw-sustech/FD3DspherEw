clear all;

fnm_conf_grid='../SeisGrid.conf';

% check
if ~ exist(fnm_conf_grid,'file')
   error([mfilename ': file ' fnm_conf_grid ' does not exist']);
end

% read data
str_cm='#';
fid=fopen(fnm_conf_grid);
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
case 'vmap_nlayer'
    nlayer=sscanf(s,'%f',1)';
case 'vmap_ijmax'
    ijmax=sscanf(s,'%f',2)';
case '<vmap'
    m=n
    for j=1:ijmax(2)
    for i=1:ijmax(1)
        m=m+1;
        str=conf{1}{m};
        vgrid=sscanf(str,'%f',2+nlayer);
        x(i)=vgrid(1); y(j)=vgrid(2);
        zgrid(j,i,1:nlayer)=vgrid(3:nlayer+2);
    end
    end
    break
end  % select

end

% plot
hid=figure;
set(hid,'BackingStore','on');
set(hid,'renderer','zbuffer');
%set(hid,'menubar','none');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperUnits','points')
set(gcf,'PaperPosition',[0 0 1024 768])

for n=1:nlayer-1
    surf(x,y,zgrid(:,:,n));
    hold on
end

