tissue="Cardiac"       ## set to "Cardiac" or "Hepatic"
trimmed_dir="trimmed_reads"
# dose=("The" "Tox")
dose=("The")

# timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
timeTHE=("002" "008" "024" "072" "168" "240" "336")
# timeTOX=("002" "008" "024" "072" "168" "240")

replicates=("1" "2" "3")

# . ./salmon_APA.sh Acetaminophen APA
# . ./salmon.sh Azathioprine AZA
# . ./salmon.sh Con_0.1_DMSO ConDMSO CONTROL
. ./script/salmon/salmon.sh Con_UNTR ConUNTR CONTROL

timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")

. ./script/salmon/salmon.sh Con_0.1_DMSO ConDMSO CONTROL

# . ./salmon.sh Cyclosporin CYC
# . ./salmon.sh Diclofenac DIC
# . ./salmon.sh Isoniazid ISO
# . ./salmon.sh Methotrexate MTX
# . ./salmon.sh Phenytoin PHE
# . ./salmon.sh Rifampicin RIF
# 
# timeTHE=("000" "002" "008" "024" "072" "168" "240" "336")
# timeTOX=("000" "002" "008" "024" "072" "168" "240")
# 
# . ./salmon.sh Valproic_Acid VPA
