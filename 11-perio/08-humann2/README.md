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
#SBATCH --mem 750gb
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

rm SAMPLE_humann_temp
```

Creapt scripts and submit them
```sh
ls /home/scrull/rna_scripts/02-operon-mapping/DM*.sh | sed 's/\..*//' | sed 's/.*\///' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done

ls /home/scrull/rna_scripts/02-operon-mapping/DM*.sh | sed 's/\..*//' | sed 's/.*\///' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done
# submit
ls D*sh | while read line; do sbatch $line; done
```














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