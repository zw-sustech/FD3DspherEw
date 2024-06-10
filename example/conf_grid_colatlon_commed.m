%close all
clear all

%==== number of node per degree, number of cell will be one less that this value ====
nxnod = 35;
nynod = 35;

%==== number of cell or spacing per degree ====
nxcel=nxnod-1;
nycel=nynod-1;

%==== elementary coordinate for node per degree ====
xelem = linspace(0,1,nxnod);
yelem = linspace(0,1,nynod);

%====  coordinate of the first node in each colatitude and longitude segments===
XFirstNode  = [ 84:94 ];
YFirstNode  = [ 0:35 ];

%==== calculate coordinate of all nodes ====
x     = [];
y     = [];
for n = 1:length(XFirstNode)
    x((n-1)*nxcel+1:n*nxcel+1)=xelem(1:nxnod)+XFirstNode(n); 
    %-- avoid accumulated error if calculating from the first node
end
for n = 1:length(YFirstNode)
    y((n-1)*nycel+1:n*nycel+1)=yelem(1:nynod)+YFirstNode(n); 
end

%==== output ====
nx=length(x);
ny=length(y);

fmt_out='%f\n';
% x
fid=fopen('SeisGrid.colat.dat','w');
fprintf(fid,fmt_out,x);
fclose(fid);
% y
fid=fopen('SeisGrid.lon.dat','w');
fprintf(fid,fmt_out,y);
fclose(fid);

