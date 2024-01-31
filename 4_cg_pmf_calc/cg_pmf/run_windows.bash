#!/bin/bash
limit=36 
num=0
for w in window_*
do
    echo $w
    cd $w

    ~/Documents/LFC/lmp_serial -in input.pull -var SEED 999 -log log.pull > /dev/null &
    cd ../
    ((num++))
    if [ $num -ge $limit ]
    then
        wait
        num=0
    fi

done

