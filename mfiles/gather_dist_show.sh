#!/bin/sh

vnm="Vz"
figfmt="png"

matlab -nojvm -nodisplay <<EOF
gather_dist_show('id',1, ...
   'nodisplay', ...
   'light', ...
   'flat', ... %'interp'
   'surf', ...
   %'slice', ...
   %'start',[1 1 1],'count',[-1 -1 -1],'stride',[1 1 1], ...
   %'sx',[100]*1e3,'sy',[100]*1e3,'sz',[0]*1e3, ...
   'daspect',[5 5 1], ...
   'caxis',[-0.1 0.1], ...
   '${figfmt}', ...
   '${vnm}')
exit
EOF

pnmtmp="fig.tmp"
pnmout="fig.${figfmt}"

mkdir -p $pnmtmp $pnmout

fnm_skel="PG${vnm}_skel.${figfmt}"
fnm_mask="PG${vnm}_mask.${figfmt}"
 fnm_out="PG${vnm}.${figfmt}"

if [ -f ${fnm_skel} ]; then
   convert -resize 1024x768 $fnm_skel tmp.${figfmt}
   composite -compose Multiply -gravity center tmp.${figfmt} $fnm_mask $fnm_out
   mv $fnm_skel $fnm_mask $pnmtmp/
   mv $fnm_out $pnmout/
fi

rm -rf tmp.${figfmt} $pnmtmp

