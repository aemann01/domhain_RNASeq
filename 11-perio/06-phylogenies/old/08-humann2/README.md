# 1. Download databases
```sh
sudo docker pull biobakery/humann

sudo docker run -it --rm \
    -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    humann_databases --download chocophlan full /data/humann_dbs --update-config yes

sudo docker run -it --rm \
    -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    humann_databases --download uniref uniref90_diamond /data/humann_dbs --update-config yes

sudo docker run -it --rm \
    -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    humann_databases --download utility_mapping full /data/humann_dbs --update-config yes
# Make files to run docker
ls ~/rna_dohmain/filtered//*fastq.gz | sed 's/\..*//' | sed 's/.*\///' | sort | uniq | while read line; do zcat ~/rna_dohmain/filtered/$line.1.fastq.gz ~/rna_dohmain/filtered/$line.2.fastq.gz | gzip > $line.fastq.gz; done 
```
Running docker
```sh
sudo docker run -it --rm \
	-v ~/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    metaphlan --install --bowtie2db /data/metaphlan_db 
cd /home/suzanne/rna_dohmain/11-perio/08-humann2
```
Running metaphlan
```sh
ls DM*.gz | sed 's/\..*//' | sort | uniq -u| while read -r line; do
    docker run \
        -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
        biobakery/humann:latest \
        metaphlan "/data/${line}.fastq.gz" --input_type fastq \
        --nproc 60 --bowtie2db /data/metaphlan_db \
        -o "/data/${line}.metaphlan.txt" || {
            echo "Error processing ${line}.fastq.gz"
        }
done
```
Script for running hummann after metaphlan
```sh
#!/bin/bash

# List and process files
ls DM*.gz | sed 's/\..*//' | sort | uniq -u | while IFS= read -r line; do
    echo "Processing ${line}..."
    if ! docker run --rm \
    -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    humann --input /data/${line}.fastq.gz \
    --output /data/humann \
    --threads 60 \
    --nucleotide-database /data/humann_dbs/chocophlan \
    --protein-database /data/humann_dbs/uniref \
    --taxonomic-profile /data/${line}.metaphlan.txt; then
        echo "Error processing ${line}"
    fi
done

sudo docker run --rm \
    -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    humann --input /data/DM00103V1PQ16.fastq.gz \
    --output /data/humann \
    --threads 60 \
    --nucleotide-database /data/humann_dbs/chocophlan \
    --protein-database /data/humann_dbs/uniref \
    --taxonomic-profile /data/DM00103V1PQ16.metaphlan.txt
```
```sh
chmod 755 humann3.sh
sudo nohup ./humann3.sh
singularity build humann3.sif docker://biobakery/humann:latest

```
Running on palmetto
```sh
#!/bin/bash

#SBATCH --job-name humann-SAMPLE
#SBATCH --nodes 1
#SBATCH --tasks-per-node 40
#SBATCH --cpus-per-task 1
#SBATCH --mem 650gb
#SBATCH --time 72:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

cd /scratch/scrull/hiv_rnaseq/08-humann2
apptainer exec humann3.sif humann \
    --input /scratch/scrull/hiv_rnaseq/08-humann2/SAMPLE.fastq.gz \
    --output /scratch/scrull/hiv_rnaseq/08-humann2/humann \
    --threads 40 \
    --nucleotide-database /scratch/scrull/hiv_rnaseq/08-humann2/humann_dbs/chocophlan \
    --protein-database /scratch/scrull/hiv_rnaseq/08-humann2/humann_dbs/uniref \
    --taxonomic-profile /scratch/scrull/hiv_rnaseq/08-humann2/SAMPLE.metaphlan.txt

rm -r ./humann/SAMPLE_humann_temp
rm ./SAMPLE.fastq.gz
rm ./SAMPLE.fastq.gz.bowtie2out.txt
rm ./SAMPLE.metaphlan.txt
```

Creapt scripts and submit them
```sh
ls /home/scrull/rna_scripts/02-operon-mapping/DM*.sh | sed 's/\..*//' | sed 's/.*\///' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done

ls /home/scrull/rna_scripts/02-operon-mapping/DM*.sh | sed 's/\..*//' | sed 's/.*\///' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done
# submit
ls D*sh | while read line; do sbatch $line; done
ls /scratch/scrull/hiv_rnaseq/08-humann2/DM*fastq.gz | sed 's/.*DM/DM/' | sed 's/.fastq.gz/.sh/' | while read line; do sbatch $line; done
```
# Combine
```sh
sudo docker run --rm \
    -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    humann_join_tables -i /data/humann \
    -o /data/hmp_subset_genefamilies.tsv \
    --file_name genefamilies

#normalize
sudo docker run --rm \
    -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    humann_renorm_table -i /data/hmp_subset_genefamilies.tsv \
    -o /data/hmp_subset_genefamilies-cpm.tsv \
    --units cpm


sudo docker run --rm \
    -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    humann_join_tables -i /data/humann \
    -o /data/hmp_subset_pathabundance.tsv \
    --file_name pathabundance

sudo docker run --rm \
    -v /home/suzanne/rna_dohmain/11-perio/08-humann2:/data \
    biobakery/humann:latest \
    humann_renorm_table -i /data/hmp_subset_pathabundance.tsv \
    -o /data/hmp_subset_pathabundance-cpm.tsv \
    --units cpm
python3 add_hiv.py
sed -i 's/# Pathway/FEATURE \\ SAMPLE/' updated_abundance.tsv

humann_barplot --input updated_abundance.tsv --focal-metadata hiv_status --last-metadata hiv_status \
    --output plot1.png --focal-feature METSYN-PWY 
humann_barplot --input updated_abundance.tsv --focal-metadata hiv_status --last-metadata hiv_status \
    --output plot2.png --focal-feature METSYN-PWY --sort sum
humann_barplot --input updated_abundance.tsv --focal-metadata hiv_status --last-metadata hiv_status \
    --output plot3.png --focal-feature METSYN-PWY --sort sum metadata --scaling logstack


merge_metaphlan_tables.py *.metaphlan.txt > ./merged_abundance_table.txt
sed -i '1d' merged_abundance_table.txt
python3 transpose.py merged_abundance_table.txt merged_abundance_trans.txt
sed -i 's/.metaphlan//'  merged_abundance_trans.txt

sed 's/_Abundance//g' hmp_subset_pathabundance-cpm.tsv > pathabundance_rel.tsv
```
# Maaslyn2
```R
library(Maaslin2)

#data
df_input_data = read.table(file = "merged_abundance_trans.txt", header = TRUE, sep = "\t", row.names = 1,stringsAsFactors = FALSE)
colnames(df_input_data) <- df_input_data[1, ]
df_input_data <- df_input_data[-1, ]
df_input_data[1:5, 1:5]
#metadata
df_input_metadata = read.table(file = "~/rna_dohmain/metadata.txt", header = TRUE, sep = "\t", row.names = 1,stringsAsFactors = FALSE)
df_input_metadata$sample <- rownames(df_input_metadata)
df_input_metadata <- subset(df_input_metadata, !sample %in% c("DM00103V1PQ16", "DM00035V2PQ16"))
df_input_metadata[1:5, ]
#pathways
df_input_path = read.csv("./pathabundance_rel.tsv", sep = "\t", stringsAsFactors = FALSE, row.names = 1)
df_input_path[1:5, 1:5]

#running Maaslin2
fit_data2 = Maaslin2(input_data     = df_input_data, 
                     input_metadata = df_input_metadata, 
                     min_prevalence = 0,
                     normalization  = "NONE",
                     output         = "demo_output2", 
                     fixed_effects  = c("hiv_status"),
                     reference      = c("hiv_status,HUU"),
                     transform      = "NONE",
                     correction = "BH")


df_input_metadata$hiv_status <- factor(df_input_metadata$hiv_status)
levels(df_input_metadata$hiv_status)

fit_func = Maaslin2(input_data     = df_input_path, 
                    input_metadata = df_input_metadata, 
                    output         = "hiv_functional", 
                    fixed_effects  = c("hiv_status"),
                    reference      = c("hiv_status,HUU"),
                    transform      = "NONE",
                    correction = "BH",)

# Filter the results, keeping rows where pval <= 0.05
filtered_results <- fit_func$results[fit_func$results$pval <= 0.05, , drop = FALSE]
filtered_results_prorymonas <- filtered_results[grep("Porphyromonas_gingivalis", filtered_results$feature), , drop = FALSE]
filtered_results_tannerella <- filtered_results[grep("Tannerella_forsythia", filtered_results$feature), , drop = FALSE]
filtered_results_treponema <- filtered_results[grep("Treponema_denticola", filtered_results$feature), , drop = FALSE]
filtered_results_prorymonas
filtered_results_tannerella
filtered_results_treponema
```

# DRAM
```sh
cd ~/rna_dohmain/11-perio/05-TrEMBL
rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/GCA*genomic.fna scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/05-TrEMBL/d

sudo docker pull jinlongru/dram
#setup databases
mkdir dram
sudo docker run -d \
  -v /home/suzanne/rna_dohmain/11-perio/05-TrEMBL:/mnt/data \
  -w /mnt/data/dram \
  --name dram-container \
  jinlongru/dram \
  DRAM-setup.py prepare_databases \
  --output_dir /mnt/data/dram \
  --threads 40 &

sudo docker logs -f dram-container #for monitoring
sudo docker start dram-container

sudo docker exec -w /mnt/data/dram -it dram-container \
  DRAM-setup.py print_config


sudo docker exec -w /mnt/data/ -it dram-container \
  DRAM.py annotate -i '/mnt/data/*_genomic.fna' -o '/mnt/data/dramannots' --threads 50 1> /mnt/data/dram.log 2> /mnt/data/dram.err


DRAM-setup.py print_config
#run
DRAM.py annotate -i '*_genomic.fna' -o dramannots --threads 50 1> dram.log 2> dram.err 

```
```sh
#!/bin/bash

#SBATCH --job-name dram-install
#SBATCH --nodes 1
#SBATCH --tasks-per-node 20
#SBATCH --cpus-per-task 1
#SBATCH --mem 20gb
#SBATCH --time 48:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module add dram/1.4.6
source /project/cugbf/software/gbf/DRAM/1.4.6.sh
cd /scratch/scrull/hiv_rnaseq/05-TrEMBL
cp /project/cugbf/software/gbf/DRAM/DRAM_CUGBF_config.txt ./
DRAM-setup.py import_config --config_loc DRAM_CUGBF_config.txt
DRAM-setup.py prepare_databases --output_dir /scratch/scrull/hiv_rnaseq/05-TrEMBL
```
```sh
sbatch dram_setup.sh 
```









Palmetto
Make metaphlan db

```sh
#!/bin/bash

#SBATCH --job-name metaphlan-db
#SBATCH --nodes 1
#SBATCH --tasks-per-node 5
#SBATCH --cpus-per-task 1
#SBATCH --mem 150gb
#SBATCH --time 10:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# load module
module add humann/3.9

# move into scratch
cd /scratch/scrull/hiv_rnaseq/filtered

metaphlan --install --bowtie2db /scratch/scrull/hiv_rnaseq/humann/metaphlan_db

humann_databases --download chocophlan full /scratch/scrull/hiv_rnaseq/humann/humann_dbs --update-config yes

humann_databases --download uniref uniref90_diamond /scratch/scrull/hiv_rnaseq/humann/humann_dbs --update-config yes

humann_databases --download utility_mapping full /scratch/scrull/hiv_rnaseq/humann/humann_dbs --update-config yes
```
Run metaphlan
```sh
#!/bin/bash

#SBATCH --job-name metaphlan-SAMPLE
#SBATCH --nodes 1
#SBATCH --tasks-per-node 40
#SBATCH --cpus-per-task 1
#SBATCH --mem 150gb
#SBATCH --time 72:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# load module
module add humann/3.9

# move into scratch
cd /scratch/scrull/hiv_rnaseq/filtered

# run command
metaphlan <(cat SAMPLE.1.fastq.gz SAMPLE.fastq.gz) --input_type fastq \
	--nproc 40 --bowtie2db /scratch/scrull/hiv_rnaseq/humann/metaphlan_db \
	-o ./SAMPLE.metaphlan.txt
```
Submit jobs
```sh
ls /scratch/scrull/hiv_rnaseq/filtered | grep DM | sed 's/\..*//' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done
ls D*sh | while read line; do sbatch $line; done
```






# 2. Cat together the paired in
```sh
#!/bin/bash

#SBATCH --job-name cat_fastq
#SBATCH --nodes 1
#SBATCH --tasks-per-node 10
#SBATCH --cpus-per-task 1
#SBATCH --mem 50gb
#SBATCH --time 8:00:00

#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

cd /scratch/scrull/hiv_rnaseq/humann

ls /scratch/scrull/hiv_rnaseq/filtered/*fastq.gz | sed 's/\..*//' | sed 's/.*\///' | sort | uniq | while read line; do zcat /scratch/scrull/hiv_rnaseq/filtered/$line/ scratch/scrull/hiv_rnaseq/filtered/$line.2.fastq.gz | gzip > $line.fastq.gz; done
echo "Done"
```
# 3. Run HUMANN3
```sh
humann --input demo.fastq.gz --output demo_fastq -v --diamond /home/suzanne/rna_dohmain/11-perio/08-humann2 --threads 30


humann --input DM00008V1PQ16-2.fastq.gz --output DM00008V1PQ16-2.humann2.fastq.gz -v --threads 50 

ls *fastq* | sed 's/\..*//' | grep DM | while read line; do humann -i $line.fastq.gz -o humann3 -v --diamond /home/suzanne/rna_dohmain/11-perio/08-humann2 --threads 30; done 1>> humann3.log 2>> humann3.err
```