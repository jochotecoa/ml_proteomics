# complist=("5-FU" "Con_0.1_DMSO" "Con_UNTR" "Doxorubicin" "Epirubicin") ("Celecoxib" "Docetaxel")

##### Insert path where the miRNA files are located

# if [ dir_joa = "JOA" ]; then
# 	ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA_JOA"
# else
# 	ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/"
# fi

ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/fastq/"
ppath="/share/tools/miRge2/bin/miRge2.0"
opath="/ngs-data/analysis/hecatos/juantxo/miRNA/miRge2-2021/${comp}"
newfiles=true
##### Save all files in a variable, get a string where all file paths are connected by commas
DATE1=$(date +%s)

cd ${ipath}
gunzip *.fastq.gz

##### Run the script in the output folder
mkdir -p ${opath}
cd ${opath}
echo "${ppath} annotate -s ${ipath}/*.fastq -d miRBase -pb /share/tools/bowtie-1.1.1/ -lib /share/tools/miRge2/sp -sp human -ad illumina -gff -cpu 3"
# ${ppath} annotate --adapter illumina --species human --diff-isomirs --SampleFiles ${cpath}
${ppath} annotate -s ${ipath}/*.fastq -d miRBase -pb /share/tools/bowtie-1.1.1/ -lib /share/tools/miRge2/sp -sp human -ad illumina -gff -cpu 3
echo "miRge done, lasted ${dur2}s, time to rezip!"
if ${newfiles}; then
  cd ${ipath}
  gzip *.fastq
fi
DATE4=$(date +%s)
dur3=$(($DATE4-$DATE1))
echo "miRge done, lasted ${dur3}s"

cd /share/script/hecatos/juantxo/ml_proteomics
