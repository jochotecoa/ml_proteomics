#!/bin/bash


echo $tissue
compound=$1 # "Epirubicin"
compid=$2 # "Epi"  ## The abbreviated form of the compound name to write on all output names
if [ $3 = "CONTROL" ]; then
CONTROL_SAMPLES="true"
else
CONTROL_SAMPLES="false"
fi

echo $CONTROL_SAMPLES
sleep 1

compoundfolder="/ngs-data/data/hecatos/${tissue}/${compound}"
inpfilpre=""  ## Prefix common to all input files (There are some compounds that do not have this prefix)
# startid="1186"  ## The first folder to be opened, which is named after a concrete ID number
# jump="false" ## Set to "true" if NOT all IDsamples are ordered by increasing the ID by 1 (true: (...) 1, 2, 8, 9, (...)) (false: (...) 1, 2, 3, 4, (...))
# jumpstart="918"  ## According to the previous example, this would be "2"
# jumpend="982"  ## According to the previous example, this would be "8"
#inputdir="/ngs-data/data/hecatos/hiseq_data/Total_RNA/${tissue}/${compoundfolder}/raw_files/"
# trimmed="trimmmed"
maindir="/ngs-data/data/Juantxo_tools"
inputdir="${compoundfolder}/TotalRNA/${trimmed_dir}"    ## Directory where all raw data folders are located
if [ "$compid"  ==  "DAU" ]; then
inputdir="${compoundfolder}/TotalRNA/${trimmed_dir}"    ## Directory where all raw data folders are located
fi
outputdir="/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/${tissue}/${compound}"  ## Directory where all output data will be located
index="/share/data/hecatos/juantxo/circRNA/salmon_index"
salmon="/home/jochoteco/miniconda3/bin/salmon"
#####################################################################################################################
# replicates=("1" "2" "3")
# timeTHE=("000" "002" "008" "024" "072" "168" "240" "336") #"000" 
# timeTOX=("002" "008" "024" "072" "168" "240")

#timeTHE=("002" "008" "024" "072" "168" "240" "336")
#timeTOX=("002" "008" "024" "072" "168" "240" "336")
# dose=("The" "Tox")
#####################################################################################################################
echo $compound
#rm -r ${outputdir}/*
iDATE=$(date +%s)

# Generate the sample names

unset outnames

for h in ${dose[@]}; do
if [ "$h"  ==  "Tox" ]; then
time=${timeTOX[@]}
else
time=${timeTHE[@]}
fi
for i in ${time[@]}; do
for j in ${replicates[@]}; do
if ${CONTROL_SAMPLES}; then
# time=${timeTHE[@]}
outnames+=("${inputdir}/${compid}_${i}_${j}");
else
outnames+=("${inputdir}/${compid}_${h}_${i}_${j}");
fi
done;
done;
done
totnumsam=${#outnames[@]}
namepos=0

DATE1=$(date +%s)

# Run Salmon in each paired-end samples

for fn in ${outnames[@]};
do

DATE1=$(date +%s)
samp=`basename ${fn}`
mkdir -p ${outputdir}/${samp}_quant
echo "-----------------------------------------"
echo "Processing sample ${samp}"
echo "-----------------------------------------"
${salmon} quant -i ${index} -l A -1 ${fn}_R1_trimmed_PE.fastq.gz -2 ${fn}_R2_trimmed_PE.fastq.gz -p 24 -o ${outputdir}/${samp}_quant

namepos=$(($namepos+1))
DATEB=$(date +%s)
sec=$(( $DATEB - $DATE1))
h1=$(($sec/3600))
m1=$((($sec-$h1*3600)/60))
s1=$(($sec-$h1*3600-$m1*60))
timeperfolder=$((($DATEB - $iDATE)/$namepos))
remfol=$(($totnumsam - $namepos))
remtime=$(($timeperfolder * $remfol))
h2=$(($remtime/3600))
m2=$((($remtime-$h2*3600)/60))
s2=$(($remtime-$h2*3600-$m2*60))
echo "${thiSample} took ${h1} hours, ${m1} minutes and ${s1} seconds   |  $(($namepos*100/$totnumsam))% completed     |  Time remaining: ${h2} hours, ${m2} minutes and ${s2} seconds";
done

