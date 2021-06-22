cd /share/data/hecatos/marcha_juan/transcriptome_salmon/
# Download fasta circrna from circBase3: human_hg19_circRNAs_putative_spliced_sequence.fa
wget http://www.circbase.org/download/human_hg19_circRNAs_putative_spliced_sequence.fa.gz
gunzip human_hg19_circRNAs_putative_spliced_sequence.fa.gz
cat Homo_sapiens.GRCh38.cdna.ncrna.fa human_hg19_circRNAs_putative_spliced_sequence.fa > Homo_sapiens.GRCh38.cdna.ncrna.circbase.fa # /share/data/hecatos/marcha_juan/transcriptome_salmon/human_hg19_circRNAs_putative_spliced_sequence.fa
cp Homo_sapiens.GRCh38.cdna.ncrna.circbase.fa /share/data/hecatos/juantxo/circRNA/Homo_sapiens.GRCh38.cdna.ncrna.circbase.fa
cd /share/data/hecatos/juantxo/circRNA
~/miniconda3/bin/salmon index -t Homo_sapiens.GRCh38.cdna.ncrna.circbase.fa -p 12 -i salmon_index -k 31 # /share/data/hecatos/juantxo/circRNA
