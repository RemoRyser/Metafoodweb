#!/bin/bash




#SBATCH -J summarize
#SBATCH --mail-user=remo.ryser@idiv.de
##SBATCH --mail-type=FAIL     ###don't use until it works
#SBATCH --mem-per-cpu=1G
#SBATCH -t 0-05:00
#SBATCH --cpus-per-task=1
#SBATCH -o /work/%u/%x-%A-%a.out
#SBATCH -e /work/%u/%x-%A-%a.err




INPUTDIR=$1

NAME=$3

cd "$INPUTDIR" || exit

OUTFILEGLOBAL=SummaryGlobal.out
OUTFILEMASS=SummaryMass.out

mkdir /data/idiv_brose/remo/Patchdistance/Metafoodweb/output/"$NAME""$WEB"

echo 'nut, landscape, web, number.of.spp, number.of.plants, number.of.consumers, number.of.patch, rng.seed, max.emigr.rate.plants, shape.emigr.rate.plants, max.emigr.rate.consumers, shape.emigr.rate.consumers, D_0, theta, eps, mean.patch.dist, sd.patch.dist, mean.nn.dist, sd.nn.dist, mean.con.rgg, sd.con.rgg, mean.con.rgg.plants, sd.con.rgg.plants, mean.con.rgg.consumers, sd.con.rgg.consumers, Ricker' > $OUTFILEGLOBAL
for file in global_*.out ; do


L=$(echo "$file" | grep -Eo 'LANDSCAPE_[0-9]+')
L=$(echo "$L" | grep -Eo '[0-9]+')

W=$(echo "$file" | grep -Eo 'WEB_[0-9]+')
W=$(echo "$W" | grep -Eo '[0-9]+')

I=$(echo "$file" | grep -Eo 'NUT_[0-9]+')
I=$(echo "$I" | grep -Eo '[0-9]+')

ANALYSIS=$I,

LANDSCAPE=$L,
WEB=$W,
dat=$(sed '1d' global_NUT_"${I}"_LANDSCAPE_"${L}"_WEB_"${W}".out)

echo "$ANALYSIS" "$LANDSCAPE" "$WEB" "$dat"  >> $OUTFILEGLOBAL

done





echo 'nut, landscape, web, patch, species, body.mass, int.biomass, mean.biomass,  biomass.variance, tot.mean.biomass, tot.biomass.variance, if.basal.spp, biomass.tend-20k, biomass.tend-10k,biomass.tend, spatial.connectance, Mean.net.growth1, Mean.net.growth2, net.growth.sd'  > $OUTFILEMASS


for file in mass_*.out ; do



L=$(echo "$file" | grep -Eo 'LANDSCAPE_[0-9]+')
L=$(echo "$L" | grep -Eo '[0-9]+')

W=$(echo "$file" | grep -Eo 'WEB_[0-9]+')
W=$(echo "$W" | grep -Eo '[0-9]+')

I=$(echo "$file" | grep -Eo 'NUT_[0-9]+')
I=$(echo "$I" | grep -Eo '[0-9]+')

sed '1d' mass_NUT_"${I}"_LANDSCAPE_"${L}"_WEB_"${W}".out > temp.out

ANALYSIS=$I,

LANDSCAPE=$L,
WEB=$W,
while  read -r p; do


echo "$ANALYSIS" "$LANDSCAPE" "$WEB" "$p" >> $OUTFILEMASS

done < temp.out

rm temp.out




done



