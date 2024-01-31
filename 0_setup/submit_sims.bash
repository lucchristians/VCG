#!/bin/bash
#script=$1
jid=''

for i in {1..7}
do

    if [ $i -eq 1 ]
    then
        script=init_Qwt_v1277K.sbatch
    elif [ $i -eq 2 ]
    then
        script=fstep_Qwt_v1277K.sbatch
    else
        script=run_Qwt_v1277K.sbatch
    fi
    cp ${script} job.md${i}

    # prepare SLURM submission script for each run
    if [ $i -eq 1 ]
    then
	sed -i "s/(INIT)/1/g" job.md${i}
    else
	sed -i "s/(INIT)/0/g" job.md${i}
    fi
    sed -i "s/(SID)/${i}/g" job.md${i}

    # submit each run using SBATCH
    if [ -z $jid ]
    then
	# this is the first run, so there is no job dependency
	jid=$( sbatch --parsable job.md${i} )
	echo $jid
    else
	# these are follow up runs, so there IS job dependency
	jid=$( sbatch --parsable --dependency=afterok:${jid} job.md${i} )
	echo $jid
    fi

done

