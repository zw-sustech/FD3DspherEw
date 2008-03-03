clear all

nk=230;

r0=6371;

dz=1.0;
dh=0.05;

r(nk-30:nk)=[-30:0]*dz+r0;
for k=nk-31:-1:nk-130
    r(k)=r(k+1)-dz-dh*((nk-31)-k);
end
for k=nk-131:-1:1
    r(k)=r(k+1)-(r(nk-129)-r(nk-130))-dh*((nk-131)-k);
end
