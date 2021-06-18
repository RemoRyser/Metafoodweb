#!/bin/bash

if [[ -z $2 ]] ; then
echo "usage: bash $0 jobname landscape"
exit 1
fi

NAME=$1



INPUTDIR=$(dirname $0)
LIST=""

for SCENARIO in 2
do

if [[ $SCENARIO -gt 1 ]] ; then
LANDSC=$2
else
LANDSC=1
fi


for LANDSCAPE in $(seq 1 $LANDSC)
do

for WEB in 1
do



JOB_ID=$(
sbatch \
--parsable \
-a 1-15,30-45 \
-J $NAME \
$(dirname $0)/submit.sh $NAME $SCENARIO $INPUTDIR $LANDSCAPE $WEB| cut -d. -f1)

LIST+="$JOB_ID,"


done
done

sbatch --dependency=afterany:${LIST%%,} $(dirname $0)/summarize.sh /work/$USER/$NAME/$SCENARIO $SCENARIO $NAME

done

