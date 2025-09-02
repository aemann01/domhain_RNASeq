# 1. Build bowtie2 db from rpoC from HOMD genome
```sh
cd ~/rna_dohmain/rpoc/database
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ~/rna_dohmain/homd_map/ALL_genomes.ffn > ./ALL_genomes.fnn
grep "DNA-directed RNA polymerase subunit beta'" -A 1 ./ALL_genomes.fnn > rpoC_homd.fa
bowtie2-build rpoC_homd.fa rpoC_homd.db
```
# 2. Map the samples to rpoC homd db
```sh
cd ../
mkdir bowtie2
ls ~/rna_dohmain/rpoc/data/filtN/*R1*.fastq.gz | sed 's/_S.*1.fastq.gz//' | sed 's/.*DM/DM/' | while read line; do
    fq1=$(ls ~/rna_dohmain/rpoc/data/filtN/${line}*R1_001.fastq.gz)
    fq2=$(ls ~/rna_dohmain/rpoc/data/filtN/${line}*R2_001.fastq.gz)
    bowtie2 -x ./database/rpoC_homd.db \
        -1 "$fq1" \
        -2 "$fq2" \
        --end-to-end \
        --qc-filter \
        --no-unal \
        -t \
        --threads 60 \
        --seed 465565 \
        2> ./bowtie2/$line.out \
        1> ./bowtie2/$line.err | samtools view -bS - > ./bowtie2/$line.bam
done
```
# 3. Run feature counts
```sh
cd ~/rna_dohmain/rpoc/bowtie2
grep rpoC ~/rna_dohmain/homd_map/ALL_genomes.gtf > rpoC.gtf
mkdir featurecounts
ls *bam | sed 's/.bam//' | while read line; do featureCounts -f -p -C -B -M -a rpoC.gtf -o featurecounts/$line.out -T 60 $line\.bam -t CDS -g transcript_id; done
cd featurecounts
paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=7;i<=NF;i+=7) printf "%s\t", $i; printf "\n"}' > read_counts.txt
# clean up sample names
sed 's/.bam//g' read_counts.txt | sed 's/.homd//'g > temp
mv temp read_counts.txt