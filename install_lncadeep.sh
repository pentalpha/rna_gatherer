#conda env create -f lncadeep.yml
git clone https://github.com/cyang235/LncADeep.git
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam29.0/Pfam-A.hmm.gz
mv Pfam-A.hmm LncADeep/LncADeep_lncRNA/src
cd LncADeep
chrmod +x configure
./configure