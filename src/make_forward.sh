#!/bin/bash

fnm_log=make.forward.log

hr="-"; while [ ${#hr} -lt 70 ]; do hr=${hr}"-"; done

function compile_code {
    make 2>&1 | tee $fnm_log
    ierr=${PIPESTATUS[0]}
}

time compile_code

echo -e "\n"
echo $hr

if [ $ierr -eq 0 ]; then
    echo "compile forward code succeeded"
    pdir=`cd .. && pwd`
    #echo "move bin/ to $pdir/bin"
    mkdir -p ../bin
    for i in bin/seis*; do
        echo "cp -f $i $pdir/bin"
        cp -f $i ../bin
    done
else
    echo "compile forward code failed, please check $fnm_log"
fi
echo $hr
