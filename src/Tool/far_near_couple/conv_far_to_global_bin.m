% draw_snap_pcolor: Draw wavefield snapshot by using pcolor.

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Initial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date$
% $Revision$
% $LastChangedBy$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

set_mfiles_path
[fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
  dir_station,dir_out]=get_simul_path('root','../../../run/homo')

%-------------------------------------------------------------------------------
% parameter
%-------------------------------------------------------------------------------
    SIZ_INT=4;
    SIZ_FLOAT=4;

    n1=1; n2=6000; dn=1; NTWIN=floor( (n2-n1+1)/dn );
    %n1=1; n2=10; dn=1; NTWIN=floor( (n2-n1+1)/dn );

    %id = 2; subs=[1,1,1];subc=[6,106,104];subt=[1,1,1];
    %    fnm_out='wave_far_x1.bin'
    %    indx0=[27 27 27];
    %    indx_small=[indx0(1) indx0(1)+subc(1)-1 ...
    %                indx0(2) indx0(2)+subc(2)-1 ...
    %                indx0(3) indx0(3)+subc(3)-1 ...
    %                NTWIN];

    id = 3; subs=[1,1,1];subc=[6,106,104];subt=[1,1,1];
        fnm_out='wave_far_x2.bin'
        indx0=[127 27 27];
        indx_small=[indx0(1) indx0(1)+subc(1)-1 ...
                    indx0(2) indx0(2)+subc(2)-1 ...
                    indx0(3) indx0(3)+subc(3)-1 ...
                    NTWIN];

    %id = 4; subs=[1,1,1];subc=[106,6,104];subt=[1,1,1];
    %    fnm_out='wave_far_y1.bin'
    %    indx0=[27 27 27];
    %    indx_small=[indx0(1) indx0(1)+subc(1)-1 ...
    %                indx0(2) indx0(2)+subc(2)-1 ...
    %                indx0(3) indx0(3)+subc(3)-1 ...
    %                NTWIN];

    %id = 5; subs=[1,1,1];subc=[106,6,104];subt=[1,1,1];
    %    fnm_out='wave_far_y2.bin'
    %    indx0=[27 127 27];
    %    indx_small=[indx0(1) indx0(1)+subc(1)-1 ...
    %                indx0(2) indx0(2)+subc(2)-1 ...
    %                indx0(3) indx0(3)+subc(3)-1 ...
    %                NTWIN];

    %id = 6; subs=[1,1,1];subc=[106,106,6];subt=[1,1,1];
    %    fnm_out='wave_far_z1.bin'
    %    indx0=[27 27 27];
    %    indx_small=[indx0(1) indx0(1)+subc(1)-1 ...
    %                indx0(2) indx0(2)+subc(2)-1 ...
    %                indx0(3) indx0(3)+subc(3)-1 ...
    %                NTWIN];

%-- locate
    [snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

%-- create output file
    fid=fopen(fnm_out,'w','ieee-le');
%-- file header
    siz_pad= SIZ_INT * 7;
    fwrite(fid,siz_pad,'int32');    %- padding
    fwrite(fid,indx_small,'int32'); %- header info
    fwrite(fid,siz_pad,'int32');    %- padding

    W=zeros(subc(1),subc(2),subc(3),9);
%-- time loop
    for nlayer=n1:dn:n2
    
        [ Vx,t]=gather_snap(snapinfo,id,nlayer, 'Vx','outdir',dir_out);
        [ Vy,t]=gather_snap(snapinfo,id,nlayer, 'Vy','outdir',dir_out);
        [ Vz,t]=gather_snap(snapinfo,id,nlayer, 'Vz','outdir',dir_out);
        [Txx,t]=gather_snap(snapinfo,id,nlayer,'Txx','outdir',dir_out);
        [Tyy,t]=gather_snap(snapinfo,id,nlayer,'Tyy','outdir',dir_out);
        [Tzz,t]=gather_snap(snapinfo,id,nlayer,'Tzz','outdir',dir_out);
        [Tyz,t]=gather_snap(snapinfo,id,nlayer,'Tyz','outdir',dir_out);
        [Txz,t]=gather_snap(snapinfo,id,nlayer,'Txz','outdir',dir_out);
        [Txy,t]=gather_snap(snapinfo,id,nlayer,'Txy','outdir',dir_out);

        %W=[Vx; Vy; Vz; Txx; Tyy; Tzz; Tyz; Txz; Tyz];
        W(:,:,:,1)= Vx;
        W(:,:,:,2)= Vy;
        W(:,:,:,3)= Vz;
        W(:,:,:,4)=Txx;
        W(:,:,:,5)=Tyy;
        W(:,:,:,6)=Tzz;
        W(:,:,:,7)=Tyz;
        W(:,:,:,8)=Txz;
        W(:,:,:,9)=Txy;
        
        disp([ '  draw ' num2str(nlayer) 'th layer(t=' num2str(t) ')']);

    %-- write to output file 
        siz_pad= SIZ_FLOAT * prod(subc) *9;
        fwrite(fid,siz_pad,'int32');    %- padding
        fwrite(fid,W,'float32');        %- data
        fwrite(fid,siz_pad,'int32');    %- padding
    end
    
%-- close file
    fclose(fid);
