# 1. Porphyromonas gingivalis tree
```sh
cd ~/bin
#setup iqtree3
wget https://github.com/iqtree/iqtree3/releases/download/v3.0.1/iqtree-3.0.1-Linux-intel.tar.gz 
tar -xzvf iqtree-3.0.1-Linux-intel.tar.gz
mv ~/bin/iqtree-3.0.1-Linux-intel/bin/iqtree ./
#setup 3seq
wget https://mol.ax/rxiv/3seq_build_170612.zip
unzip 3seq_build_170612.zip
mv 3seq\ build\ 170612 3seq_build_170612
cd 3seq_build_170612
make
mv 3seq ../
./3seq -g myPvalueTable500 500

# strart making phylogenies
cd ~/rna_dohmain/11-perio/06-phylogenies/prokka
mkdir p_gingivalis && cd p_gingivalis
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
grep Porphyromonas ./ALL_genomes.ffn -A 1 | grep gingivalis -A 1 | grep -vA 1 cangingivalis > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ../../single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
# ls c*align.fa | sed 's/.align.fa//' | while read line; do Phi -f $line.align.fa > $line.rec; done
Rscript ../../../pairwise.R
# 3seq -f c2502.align.fa -ptable ~/bin/myPvalueTable500 -id c2502.align.rec
# 3seq -f c1153.align.fa -id c1153.align.rec
ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
cat recomb | while read line; do mv $line $line.rd; done
# grep Bonferroni *3s.log | sed 's/.3s.*= /.fa\t/' | sort -k2 -n | awk '(NR>1) && ($2 > 0.05 )' > no_rec
# grep Bonferroni *3s.log | sed 's/.3s.*= /.fa\t/' | sort -k2 -n | awk '(NR>1) && ($2 < 0.05 )' > recom

# ls c*align.fa | sed 's/.align.fa//' | while read line; do run_gubbins.py $line.align.fa --threads 130 --model GTRCAT --prefix $line.gubbins > $line.rec; done
# run_gubbins.py c3118.align.fa --threads 60 --model GTRCAT --prefix c3118.gubbins

# find . -type f -name "*.embl" -empty

# combine single copy core
# combine alignments
python3 ../../../combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core -safe
```
# 2. Trepoenma denticola tree
```sh
# strart making phylogenies
cd ~/rna_dohmain/11-perio/06-phylogenies/prokka
mkdir t_dent && cd t_dent
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
grep Treponema ./ALL_genomes.ffn -A 1 | grep denticola -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ../../single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
# 3seq -f c2502.align.fa -ptable ~/bin/myPvalueTable500 -id c2502.align.rec
# 3seq -f c1153.align.fa -id c1153.align.rec
ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
cat recomb | while read line; do mv $line $line.rd; done
# grep Bonferroni *3s.log | sed 's/.3s.*= /.fa\t/' | sort -k2 -n | awk '(NR>1) && ($2 > 0.05 )' > no_rec
# grep Bonferroni *3s.log | sed 's/.3s.*= /.fa\t/' | sort -k2 -n | awk '(NR>1) && ($2 < 0.05 )' > recom

# ls c*align.fa | sed 's/.align.fa//' | while read line; do run_gubbins.py $line.align.fa --threads 130 --model GTRCAT --prefix $line.gubbins > $line.rec; done
# run_gubbins.py c3118.align.fa --threads 60 --model GTRCAT --prefix c3118.gubbins

# find . -type f -name "*.embl" -empty

# combine single copy core
# combine alignments
python3 ../../../combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix tdent.core -safe
```
# 3. Tannerella forsythia tree
```sh
# strart making phylogenies
cd ~/rna_dohmain/11-perio/06-phylogenies/prokka
mkdir t_for && cd t_for
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ALL_genomes.ffn
grep Tannerella ./ALL_genomes.ffn -A 1 | grep forsythia -A 1 > all.fnn
sed -i 's/ .*//' all.fnn
sed -i '/^--$/d' all.fnn
# cluster genes
vsearch --cluster_fast all.fnn --otutabout gene_cluster.tab --uc uc --id 0.5 --threads 190 --clusters c --log vsearchlog #cluster
sed -i 's/#OTU ID/OTU_ID/g' gene_cluster.tab
# find single copy core genes
python3 ../../single_copy.py
grep -w -f single-copy-core-tags c* > single-copy-core-cluster-ids #get the ids
sed -i 's/:>/\t/g' single-copy-core-cluster-ids
cut -f 1 single-copy-core-cluster-ids > single-copy-core-cluster-ids2
#align
mkdir ./align
cp $(cat single-copy-core-cluster-ids2) ./align && cd align
ls c* | grep -v "align\|rec" | while read line; do mafft --thread -1 $line > $line.align.fa; done 
sed -i 's/_.*//g' *align.fa

# test for recombination
# 3seq -f c2502.align.fa -ptable ~/bin/myPvalueTable500 -id c2502.align.rec
# 3seq -f c1153.align.fa -id c1153.align.rec
ls c*align.fa | sed 's/.align.fa//' | while read line; do yes Y | 3seq -f $line.align.fa -id $line.align; done
grep Q_ACCNUM *3s.rec -A 1 | grep SEQ | sed 's/.3s.rec.*/.fa/' | sort > recomb # get genes that are recombinant
cat recomb | while read line; do mv $line $line.rd; done
# grep Bonferroni *3s.log | sed 's/.3s.*= /.fa\t/' | sort -k2 -n | awk '(NR>1) && ($2 > 0.05 )' > no_rec
# grep Bonferroni *3s.log | sed 's/.3s.*= /.fa\t/' | sort -k2 -n | awk '(NR>1) && ($2 < 0.05 )' > recom

# ls c*align.fa | sed 's/.align.fa//' | while read line; do run_gubbins.py $line.align.fa --threads 130 --model GTRCAT --prefix $line.gubbins > $line.rec; done
# run_gubbins.py c3118.align.fa --threads 60 --model GTRCAT --prefix c3118.gubbins

# find . -type f -name "*.embl" -empty

# combine single copy core
# combine alignments
python3 ../../../combine_core.py
grep ">" core_genome.align.fa -c

# make trees
fasttree -nt core_genome.align.fa > core_genome.fast.tre
# iqtree
# mkdir alignments
# cp $(ls c*.align.fa | grep -v core) ./alignments
# iqtree3 -p ./alignments --out-aln core_genome.align.phylip --out-format Raxml -redo # making partition file
# iqtree2 -s core_genome.align.phylip -p core_genome.align.phylip.partitions -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix pging.core_genome -safe
iqtree -s core_genome.align.fa -m MFP+MERGE --seed 332244 -B 1000 -T AUTO --prefix tfor.core -safe
```