
for REP in $(seq 1 1 5)
do
{
######## Nutrient Enrichment
declare -a arr1=($(seq 0 0.1 10))  ## array for nutrient supply $(seq 0.1 0.1 10)



for VAR in $(seq 1 ${#arr1[*]})
do
echo ${arr1[$VAR-1]}
SEED=$(($VAR + $REP * 500))
NAME=${VAR}_${REP}
LANDSCAPE=1
OUTPUTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/Nutrients/"
WEB=2
INPUTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
NUTRATE=0.25
F_BASH=0.05
SUPPLY=${arr1[$VAR-1]}
EMIGR=0
LANDSCAPESIZE=1
LOSS=0
SUP2=0

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


########## Emigration Loss


declare -a arr2=($(seq 0 0.001 0.150))  ## array for nutrient supply $(seq 0.1 0.1 10)

## add loop for different scenarios

for VAR in $(seq 1 ${#arr2[*]})
do
echo ${arr2[$VAR-1]}
SEED=$(($VAR + $REP * 500))
NAME=${VAR}_${REP}
LANDSCAPE=1
OUTPUTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/Emigration/"
WEB=2
INPUTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
NUTRATE=0.25
F_BASH=0.05
SUPPLY=10
EMIGR=${arr2[$VAR-1]}
LANDSCAPESIZE=2
LOSS=1
SUP2=0

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
} &
done
