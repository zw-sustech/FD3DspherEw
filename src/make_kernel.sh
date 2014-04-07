#!/bin/bash

fnm_log=make.kernel.log

hr="-"; while [ ${#hr} -lt 70 ]; do hr=${hr}"-"; done

function compile_code {
    make kernel -f Makefile.kernel 2>&1 | tee $fnm_log
    ierr=${PIPESTATUS[0]}
}

time compile_code

echo -e "\n"
echo $hr
if [ $ierr -eq 0 ]; then
   echo "compile kernel succeeded"
    pdir=`cd .. && pwd`
    #echo "move bin/ to $pdir/bin"
    mkdir -p ../bin
    for i in bin/SI_*; do
        echo "cp -f $i $pdir/bin"
        cp -f $i ../bin
    done
else
   echo "compile kernel failed, please check $fnm_log"
fi
echo $hr
