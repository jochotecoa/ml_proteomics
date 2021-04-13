# bash scripts/miRge2/miRge2.sh 5-FU JOA
# bash scripts/miRge2/miRge2.sh Amiodarone JOA
# bash scripts/miRge2/miRge2.sh Celecoxib JOA
# bash scripts/miRge2/miRge2.sh Con_0.1_DMSO 
# bash scripts/miRge2/miRge2.sh Con_Flu_DMSO_Gal
# bash scripts/miRge2/miRge2.sh Con_UNTR
# bash scripts/miRge2/miRge2.sh Docetaxel JOA
# bash scripts/miRge2/miRge2.sh Doxorubicin 
# bash scripts/miRge2/miRge2.sh Epirubicin  
# bash scripts/miRge2/miRge2.sh Mitoxantrone JOA
# bash scripts/miRge2/miRge2.sh Paclitaxel JOA


tissue="Hepatic"

comp="t0_controls"
. ./script/miRge2/miRge2.sh

comp="Azathioprine"
dir_joa = "JOA"
. ./script/miRge2/miRge2.sh

comp="Diclofenac"
. ./script/miRge2/miRge2.sh


comp="Isoniazid"
. ./script/miRge2/miRge2.sh


comp="Methotrexate"
. ./script/miRge2/miRge2.sh



comp="Cyclosporin"
. ./script/miRge2/miRge2.sh


comp="Isoniazid"
. ./script/miRge2/miRge2.sh


comp="Phenytoin"
. ./script/miRge2/miRge2.sh


comp="Rifampicin"
. ./script/miRge2/miRge2.sh


comp="Valproic_Acid"
. ./script/miRge2/miRge2.sh

comp="Con_0.1_DMSO"
. ./script/miRge2/miRge2.sh


comp="Acetaminophen"
. ./scripts/ciri2/ciri2v6.1_APA.sh

