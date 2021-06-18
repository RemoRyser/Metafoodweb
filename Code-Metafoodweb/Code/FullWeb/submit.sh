#!/bin/bash
 



#SBATCH -J name
#SBATCH --mail-user=remo.ryser@idiv.de
##SBATCH --mail-type=FAIL     ###don't use until it works
#SBATCH --mem-per-cpu=4G
#SBATCH -t 20-00:00
#SBATCH --cpus-per-task=1
#SBATCH -o /work/%u/%x-%A-%a.out
#SBATCH -e /work/%u/%x-%A-%a.err

SCENARIO=$2
NUT=$SLURM_ARRAY_TASK_ID    ### or do I need %a ?
LANDSCAPESIZE=1
SEED=$4
LANDSCAPE=$4
WEB=$5
NAME='NUT'_${NUT}_'LANDSCAPE'_${LANDSCAPE}_'WEB'_${WEB}
NUTRATE=0.25
F_BASH=0.05
SUPPLY=1

if [[ $SCENARIO -gt 1 ]] ; then
EMIGR=0.05
else
EMIGR=0
fi

if [[ $SCENARIO -gt 2 ]] ; then
KILL=1
else
KILL=0
fi

LOSS=1



module load foss/2019b GSL/1.16


OUTPUTDIR="/work/$USER/$JOB_NAME/$SCENARIO/"
INPUTDIR=$3
mkdir -p "$OUTPUTDIR"

export OUTPUTDIR
export NAME
export LANDSCAPE
export SEED
export WEB
export INPUTDIR
export NUTRATE
export F_BASH
export SUPPLY
export EMIGR
export LANDSCAPESIZE
export NUT
export LOSS
export KILL

./simulation




