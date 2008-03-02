#!/bin/bash
 
# $Date$
# $Revision$
# $LastChangedBy$

vnm="Vz"
figfmt="png"
pnm_out="../output/"

n1=100
n2=4400
dn=100

matlab -nojvm -nodisplay <<EOF
gather_snap_show('id',1, ...
   'n',[${n1} ${n2} ${dn}], ...
   'nodisplay', ...
   'light', ...
   'flat', ... %'interp'
   'surf', ...
   %'slice', ...
   %'start',[1 1 1],'count',[-1 -1 -1],'stride',[1 1 1], ...
   %'sx',[100]*1e3,'sy',[100]*1e3,'sz',[0]*1e3, ...
   'daspect',[5 5 1], ...
   'caxis',[-0.1 0.1], ...
   'outdir','${pnm_out}, ...
   '${figfmt}', ...
   '${vnm}');
exit
EOF

pnmtmp="fig.tmp"
pnmout="fig.${figfmt}"

mkdir -p $pnmtmp $pnmout

for ((n=${n1};n<=${n2};n=n+${dn})); do
    tfmt=${n}
    while [ ${#tfmt} -lt 5 ]; do tfmt="0${tfmt}"; done
    fnm_skel="${vnm}_ndim${tfmt}_skel.${figfmt}"
    fnm_mask="${vnm}_ndim${tfmt}_mask.${figfmt}"
     fnm_out="${vnm}_ndim${tfmt}.${figfmt}"

if [ -f ${fnm_skel} ]; then
    convert -resize 1024x768 $fnm_skel tmp.${figfmt}
    composite -compose Multiply -gravity center tmp.${figfmt} $fnm_mask $fnm_out
    mv $fnm_skel $fnm_mask $pnmtmp/
    mv $fnm_out $pnmout/
fi

done

rm -rf tmp.${figfmt} $pnmtmp


