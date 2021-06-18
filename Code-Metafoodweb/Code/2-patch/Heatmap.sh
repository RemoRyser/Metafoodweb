
for REP in $(seq 1 1 5)                 ## loop acrooss replicates
do
{

declare -a arr3=($(seq 0 1 20))  ## array for nutrient supply

declare -a arr4=($(seq 0.001 0.05 1.001))  ## array for dispersal loss



for NUT in $(seq 1 ${#arr3[*]})              ## loop acrooss Nutrients
do

for LOSS in $(seq 1 ${#arr4[*]})              ## loop acrooss Dispersal Loss
do
echo 'REP'_${REP}_'NUT'_${arr3[$NUT-1]}_'LOSS'_${arr4[$LOSS-1]}

SEED=$(($NUT + $REP * 500 + $LOSS * 500))
NAME='REP'_${REP}_'NUT'_${arr3[$NUT-1]}_'LOSS'_${arr4[$LOSS-1]}
LANDSCAPE=1
OUTPUTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/Heatmap/"
WEB=1
INPUTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
NUTRATE=0.25
F_BASH=0.05
SUPPLY=${arr3[$NUT-1]}
EMIGR=0.05
LANDSCAPESIZE=${arr4[$LOSS-1]}
LOSS=1
##SUP2=$((30 - $SUPPLY))
SUP2=1
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
export LOSS
export SUP2

./simulation SEED NAME LANDSCAPE OUTPUTDIR WEB INPUTDIR NUTRATE F_BASH SUPPLY EMIGR LANDSCAPESIZE LOSS SUP2
done

done
}&
done
