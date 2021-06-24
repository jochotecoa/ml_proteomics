# bash script/miRge2/miRge2.sh 5-FU JOA
# bash script/miRge2/miRge2.sh Amiodarone JOA
# bash script/miRge2/miRge2.sh Celecoxib JOA
# bash script/miRge2/miRge2.sh Con_0.1_DMSO 
# bash script/miRge2/miRge2.sh Con_Flu_DMSO_Gal
# bash script/miRge2/miRge2.sh Con_UNTR
# bash script/miRge2/miRge2.sh Docetaxel JOA
# bash script/miRge2/miRge2.sh Doxorubicin 
# bash script/miRge2/miRge2.sh Epirubicin  
# bash script/miRge2/miRge2.sh Mitoxantrone JOA
# bash script/miRge2/miRge2.sh Paclitaxel JOA

tissue="Cardiac"
comp="t0_controls_ML"

ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/5FU"
. ./script/miRge2/miRge2.sh 
ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/AMI"
. ./script/miRge2/miRge2.sh 
ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/CEL"
. ./script/miRge2/miRge2.sh 
ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/DAU"
. ./script/miRge2/miRge2.sh 
ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/DF2"
. ./script/miRge2/miRge2.sh 
ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/DMSO_0.1"
. ./script/miRge2/miRge2.sh 
ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/DOC"
. ./script/miRge2/miRge2.sh 
ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/MXT"
. ./script/miRge2/miRge2.sh 
ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/PTX"
. ./script/miRge2/miRge2.sh 
ipath="/ngs-data/data/hecatos/${tissue}/${comp}/miRNA/UNTR"
. ./script/miRge2/miRge2.sh 



# tissue="Hepatic"
# 
# comp="t0_controls"
# . ./script/miRge2/miRge2.sh
# 
# comp="Azathioprine"
# dir_joa = "JOA"
# . ./script/miRge2/miRge2.sh
# 
# comp="Diclofenac"
# . ./script/miRge2/miRge2.sh
# 
# 
# comp="Isoniazid"
# . ./script/miRge2/miRge2.sh
# 
# 
# comp="Methotrexate"
# . ./script/miRge2/miRge2.sh
# 
# 
# 
# comp="Cyclosporin"
# . ./script/miRge2/miRge2.sh
# 
# 
# comp="Isoniazid"
# . ./script/miRge2/miRge2.sh
# 
# 
# comp="Phenytoin"
# . ./script/miRge2/miRge2.sh
# 
# 
# comp="Rifampicin"
# . ./script/miRge2/miRge2.sh
# 
# 
# comp="Valproic_Acid"
# . ./script/miRge2/miRge2.sh
# 
# comp="Con_0.1_DMSO"
# . ./script/miRge2/miRge2.sh
# 
# 
# comp="Acetaminophen"
# . ./script/ciri2/ciri2v6.1_APA.sh
# 
