# 1. Download databases
```sh
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.dat.gz
zcat uniprot_trembl_bacteria.dat.gz > uniprot_trembl_bacteria.dat
wget https://raw.githubusercontent.com/tseemann/prokka/master/bin/prokka-uniprot_to_fasta_db
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz
gunzip uniprot_sprot_bacteria.dat.gz
chmod 755 prokka-uniprot_to_fasta_db
```
# 2. Make database
```sh
sed -i '1 i\//' uniprot_sprot_bacteria.dat
cat uniprot_trembl_bacteria.dat uniprot_sprot_bacteria.dat > uniprot_bacteria.dat
#make prokka db
nohup ./prokka-uniprot_to_fasta_db --verbose --evidence=1 uniprot_bacteria.dat > uniprot
mv uniprot ~/bin/miniconda3/envs/2024-HIV_RNASeq/db/kingdom/Bacteria/sprot

cp ~/bin/miniconda3/envs/2024-HIV_RNASeq/db/kingdom/Bacteria/sprot ./sprot_og

./prokka-uniprot_to_fasta_db --verbose --evidence=1 uniprot_bacteria.dat > uniprot_db 2>uniport_log
cp uniprot_db ~/bin/miniconda3/envs/2024-HIV_RNASeq/db/kingdom/Bacteria/sprot

```
# 3. Run prokka
```sh
sudo docker pull staphb/prokka:latest
awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' ~/rna_dohmain/homd_map/ALL_genomes.fna  > ./ALL_genomes.fna
sed -i 's/ .*//' ./ALL_genomes.fna
sed -i 's/\|/_/' ./ALL_genomes.fna

sudo docker run --rm \
  -v $(pwd):$(pwd) -w $(pwd) \
  staphb/prokka \
  prokka --cpus 60 \
	   --force \
	   --prefix prokka_trembl  \
	   --norrna \
	   --notrna \
	   --proteins ~/rna_dohmain/11-perio/05-TrEMBL/uniprot_db \
	   ~/rna_dohmain/11-perio/05-TrEMBL/ALL_genomes.fna 
```
# 4. Star mappign to refseq genomes for red complex
```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

paste -d "/" <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > trep_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > porph_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > tann_query

rm query
cat *query > query
wget -i query &
cat <(cat query | sed 's/.*\///') <(ls GCF*) | sort | uniq -u > missing

paste -d "/" <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > trep_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > porph_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > tann_query

rm query
cat *query > query
wget -i query &
cat <(cat query | sed 's/.*\///') <(ls GCF*.fna.gz) | sort | uniq -u > missing
gzip -d GCF*.gz
ls GCF* | sed 's/\..*//' | sort | uniq -u | while read line; do rm $line*; done

paste -d "/" <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gbff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > trep_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gbff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > porph_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gbff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > tann_query

rm query
cat *query > query
wget -i query &
cat <(cat query | sed 's/.*\///') <(ls GCF*) | sort | uniq -u > missing
```
Combine gff file
```py 
from Bio import SeqIO

def merge_gff_files(file_list, output_file):
    with open(output_file, 'w') as outfile:
        for file in file_list:
            with open(file, 'r') as infile:
                for line in infile:
                    if not line.startswith('#'):  # Skip comment lines
                        outfile.write(line)

# List of GFF files
gff_files = [
    "GCF_000007585.1_ASM758v1_genomic.gff",
    "GCF_002754075.1_ASM275407v1_genomic.gff",
    "GCF_024181485.1_ASM2418148v1_genomic.gff",
    "GCF_000008185.1_ASM818v1_genomic.gff",
    "GCF_002754095.1_ASM275409v1_genomic.gff",
    "GCF_024181505.1_ASM2418150v1_genomic.gff",
    "GCF_000010505.1_ASM1050v1_genomic.gff",
    "GCF_002754115.1_ASM275411v1_genomic.gff",
    "GCF_024181545.1_ASM2418154v1_genomic.gff",
    "GCF_000238215.1_ASM23821v1_genomic.gff",
    "GCF_002754135.1_ASM275413v1_genomic.gff",
    "GCF_024181565.1_ASM2418156v1_genomic.gff",
    "GCF_000270225.1_ASM27022v1_genomic.gff",
    "GCF_002754155.1_ASM275415v1_genomic.gff",
    "GCF_024181605.1_ASM2418160v1_genomic.gff",
    "GCF_001263815.1_ASM126381v1_genomic.gff",
    "GCF_002892555.1_ASM289255v1_genomic.gff",
    "GCF_024181625.1_ASM2418162v1_genomic.gff",
    "GCF_001274615.1_ASM127461v1_genomic.gff",
    "GCF_002892575.1_ASM289257v1_genomic.gff",
    "GCF_024181645.1_ASM2418164v1_genomic.gff",
    "GCF_001314265.1_ASM131426v1_genomic.gff",
    "GCF_002892595.1_ASM289259v1_genomic.gff",
    "GCF_024400535.1_ASM2440053v1_genomic.gff",
    "GCF_001444325.1_ASM144432v1_genomic.gff",
    "GCF_018141745.1_ASM1814174v1_genomic.gff",
    "GCF_024400725.1_ASM2440072v1_genomic.gff",
    "GCF_001547855.1_ASM154785v1_genomic.gff",
    "GCF_018141765.1_ASM1814176v1_genomic.gff",
    "GCF_025905525.1_ASM2590552v1_genomic.gff",
    "GCF_001547875.1_ASM154787v1_genomic.gff",
    "GCF_023822425.1_ASM2382242v1_genomic.gff",
    "GCF_028335085.1_ASM2833508v1_genomic.gff",
    "GCF_002753935.1_ASM275393v1_genomic.gff",
    "GCF_023822445.1_ASM2382244v1_genomic.gff",
    "GCF_028335105.1_ASM2833510v1_genomic.gff",
    "GCF_002753955.1_ASM275395v1_genomic.gff",
    "GCF_023822465.1_ASM2382246v1_genomic.gff",
    "GCF_028335125.1_ASM2833512v1_genomic.gff",
    "GCF_002753975.1_ASM275397v1_genomic.gff",
    "GCF_024181405.1_ASM2418140v1_genomic.gff",
    "GCF_030252365.1_ASM3025236v1_genomic.gff",
    "GCF_002754015.1_ASM275401v1_genomic.gff",
    "GCF_024181425.1_ASM2418142v1_genomic.gff",
    "GCF_030440475.1_ASM3044047v1_genomic.gff",
    "GCF_002754035.1_ASM275403v1_genomic.gff",
    "GCF_024181445.1_ASM2418144v1_genomic.gff",
    "GCF_030440495.1_ASM3044049v1_genomic.gff",
    "GCF_002754055.1_ASM275405v1_genomic.gff",
    "GCF_024181465.1_ASM2418146v1_genomic.gff",
    "GCF_900638305.1_57043_C01_genomic.gff"
]

output_file = "combined.gff"
merge_gff_files(gff_files, output_file)
```
Make db
```sh
cat GCF*fna | sed 's/ .*//' > combined.fna

gffread combined.gff -Tv -o combined.gtf
STAR --runMode genomeGenerate \
	--genomeFastaFiles combined.fna  \
	--runThreadN 7 \
	--sjdbGTFfile combined.gtf \
	--genomeChrBinNbits 10 \
	--limitGenomeGenerateRAM 1493556046773 \
	--sjdbGTFfeatureExon CDS 2> db.log
```
Move db onto palmetto
```sh
#GenomeDir
rsync -a ~/rna_dohmain/filtered scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/

rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/combined.fna scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/combined.gtf scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
#prokka
rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/GenomeDir scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
```
Run star mapping
```sh
#!/bin/bash

#SBATCH --job-name redmap-SAMPLE
#SBATCH --nodes 1
#SBATCH --tasks-per-node 40
#SBATCH --cpus-per-task 1
#SBATCH --mem 750gb
#SBATCH --time 72:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# load module
module add biocontainers star/2.7.10b

# move into scratch
cd /scratch/scrull/hiv_rnaseq/filtered

# run command
STAR --runThreadN 40 \
   --genomeDir /scratch/scrull/hiv_rnaseq/red-mapping/GenomeDir \
   --readFilesIn <(gunzip -c SAMPLE.1.fastq.gz) <(gunzip -c SAMPLE.2.fastq.gz) \
   --outFileNamePrefix /scratch/scrull/hiv_rnaseq/red-mapping/SAMPLE.red \
   --outSAMtype BAM Unsorted \
   --outReadsUnmapped /scratch/scrull/hiv_rnaseq/red-mapping/SAMPLE.red.unmapped.bam \
   --quantMode TranscriptomeSAM GeneCounts \
   --alignIntronMax 1 \
   --chimOutType SeparateSAMold
```
Creapt scripts and submit them
```sh
ls /scratch/scrull/hiv_rnaseq/filtered | sed 's/\..*//' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done

ls /home/scrull/rna_scripts/02-operon-mapping/DM*.sh | sed 's/\..*//' | sort | uniq | sed 's/.*\///' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done
# submit
ls D*sh | while read line; do sbatch $line; done
```
Make feauture counts
```sh
#!/bin/bash

#SBATCH --job-name redFeature
#SBATCH --nodes 1
#SBATCH --tasks-per-node 60
#SBATCH --cpus-per-task 1
#SBATCH --mem 500gb
#SBATCH --time 24:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module add subread/2.0.3
cd /scratch/scrull/hiv_rnaseq/red-mapping
mkdir featurecounts
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -C -B -M -a ./combined.gtf -o featurecounts/$line.out -T 60 $line\Aligned.out.bam -t CDS -g transcript_id; done

cd featurecounts
paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=7;i<=NF;i+=7) printf "%s\t", $i; printf "\n"}' > read_counts.txt
# clean up sample names
sed 's/Aligned.out.bam//g' read_counts.txt | sed 's/.homd//'g > temp
mv temp read_counts.txt
```
Move feature counts to stella
```sh
scp  scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/featurecounts/read_counts.txt ~/rna_dohmain/11-perio/05-TrEMBL
```
Make gene map key
```sh
parallel -a <(awk '{print $1}' read_counts.txt | sed '1d') -j 7 -k " grep CDS ./combined.gff | grep -wm 1 '{}'" > gene.txt
#get the contig to taxa
grep gene= gene.txt > gene2.txt
paste -d "\t" <(awk '{print $1}' gene2.txt) <(sed 's/.*Parent=//' gene2.txt | sed 's/;.*gene=/\t/' | sed 's/;.*product=/\t/' | sed 's/;.*//' | sed 's/ /_/g') > annots.txt
paste -d "/" <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > trep_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > porph_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > tann_query

rm query
cat *query > query
wget -i query &
cat <(cat query | sed 's/.*\///') <(ls GCF*.fna.gz) | sort | uniq -u > missing
gzip -d GCF*
#get taxonomy
parallel -a <(awk '{print $1}' annots.txt) -j 7 -k '
    grep -wm 1 "{}" ./*genomic.fna > /dev/null && \
    grep -wm 1 "{}" ./*genomic.fna || \
    echo "{} not found"
' > taxa.txt
sed -i 's/not found/Treponema denticola ATCC 35405, complete genome/' taxa.txt

paste -d "\t" annots.txt <(awk -F " " '{print $2}' taxa.txt) <(awk -F " " '{print $2, $3}' taxa.txt | sed 's/ /_/') <(awk -F " " '{print $2, $3, $4, $5}' taxa.txt | sed 's/ /_/g' | sed 's/,_.*//') > red_annots.txt
grep -v cangingivalis red_annots.txt > temp
mv temp red_annots.txt
sed -i '1i\contig\ttag\tgene\tproduct\tgenus\tspecies\tstrain' red_annots.txt
cd ../06-red-complex
#give uniq identifier
awk '
{
    # Initialize or increment the counter for each unique gene ID
    if (!seen[$1]) {
        seen[$1] = 1
    } else {
        seen[$1]++
    }
    
    # Print the updated line with the new suffix
    printf("%s_%d", $1, seen[$1])
    for (i = 2; i <= NF; i++) {
        printf("\t%s", $i)
    }
    print ""
}
' ../05-TrEMBL/read_counts.txt | sed 's/Geneid_1/Geneid/' > read_counts.txt

awk '
BEGIN { FS = OFS = "\t" }  # Set input and output field separators to tab
NR == 1 { print; next }    # Print the header line unchanged
{
    if (!seen[$2]) {
        seen[$2] = 1
    } else {
        seen[$2]++
    }
    
    # Print the modified line with the incremented suffix in the tag column
    $2 = $2 "_" seen[$2]
    print
}
' ../05-TrEMBL/red_annots.txt > red_annots.txt

parallel -a <(awk '{print $2}' ./red_annots.txt | sed 's/tag/Geneid/') -j 7 -k " grep -w '{}' ./read_counts.txt" > red_counts.txt
```
Move onto 06-red-complex
