#!/bin/bash
limit=36
num=0
echo "#path loc_win_min spring" > metadata.dat
i=0
coord=$1
#for i in {0..40}
until [ ! -d window_${i} ]
do
    cd window_${i}

    python ../../process_log_for_wham.py y &
    ((num++))
    if [ $num -ge $limit ]
    then
        wait
        num=0
    fi


    if [ ${coord} == 'x' ]
    then 
    	spring=$( tail -n 1 spring.settings | awk '{print $8, $7}' )
    elif [ ${coord} == 'y' ]
    then 
    	spring=$( tail -n 1 spring.settings | awk '{print $9, $7}' )
    elif [ ${coord} == 'z' ]
    then 
    	spring=$( tail -n 1 spring.settings | awk '{print $10, $7}' )
    fi
    echo "window_${i}/CV.dat $spring" >> ../metadata.dat
    ((i++))
    cd ../
done
wait
