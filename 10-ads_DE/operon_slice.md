# Operon identification script

The first thing we need to do is limit our analysis to genes that are present in a full operon in our genomes -- to be considered "full" they need to include arcA, arcB, and arcC in syntenty with the caveat that some genomes may have intermediate genes that are dispersed between the genes. We can get the information that we need from the gff file. 

Pull gene information by EC number:
arcA -- Arginine deiminase: 3.5.3.6
arcB -- Ornithine carbamoyltransferase: 2.1.3.3
arcC -- Carbamate kinase: 2.7.2.2

```bash
conda activate 2024-HIV_RNASeq
mamba install matplotlib
mamba install seaborn
```

```bash
cd ~/domhain_RNAseq/07-ads_expression/
# pull records for all arc genes
grep -E "eC_number=3.5.3.6|eC_number=2.1.3.3|eC_number=2.7.2.2" /home/allie/domhain_RNAseq/03-star_map/homd_db/ALL_genomes.gff > arc_genes.gff
# clean up to only include the data that we need and format for parsing
awk -F"\t" '{print $1 "\t" $4 "\t" $5 "\t" $9}' arc_genes.gff | sed 's/;/\t/g' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > arc_genes.clean.txt
# add in header line
sed -i '1i homdID\tstart\tend\tlocustag\tEC' arc_genes.clean.txt
```

Next I want to get a file that has each genome followed by the coordinates for each gene in the operon to filter out truncated or otherwise problematic genomes (i.e., incomplete, poor assembly, etc.)

```python
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

arcgenes = pd.read_csv("arc_genes.clean.txt", header=0, sep="\t")
# get length of each concatenated gene
arcgenes["length"] = abs(arcgenes["end"] - arcgenes["start"])
# add gene names
arcgenes.loc[arcgenes['EC'] == "eC_number=3.5.3.6", 'EC'] = 'arcA'
arcgenes.loc[arcgenes['EC'] == "eC_number=2.1.3.3", 'EC'] = 'arcB'
arcgenes.loc[arcgenes['EC'] == "eC_number=2.7.2.2", 'EC'] = 'arcC'

# first need to check that each homdID has at lest the three genes required for the ADS pathway
ads_genes = {'arcA', 'arcB', 'arcC'}

def check_genes(homdID):
    return ads_genes.issubset(set(homdID['EC']))

genespres = arcgenes.groupby('homdID').filter(check_genes)

# check to see if the leptos are still present
genespres["homdID"].str.contains('SEQF2444').any() # should return true

# slice out locus tag number and add to dataframe
genespres["locusnum"] = genespres.locustag.str.split('_', expand=True)[1].astype("int")

# define the allowable deviation from synteny
allowdev = 2

def find_consecutive_groups(sub_df):
    sub_df = sub_df.sort_values('locusnum').reset_index(drop=True)
    sub_df['diff'] = sub_df['locusnum'].diff().fillna(0)
    sub_df['group'] = (sub_df['diff'] > allowdev).cumsum()
    consecutive_groups = sub_df.groupby('group')['locustag'].apply(list).tolist()
    return [group for group in consecutive_groups if len(group) >= 3]

# group by genome ID
result = genespres.groupby('homdID').apply(find_consecutive_groups).tolist()

# flatten result
flatres = [group for sublist in result for group in sublist]

# make a new dataframe and merge with previous data
result_df = pd.DataFrame({'syntenous_loci': flatres})

# explode lists into indivudal rows
exploded_df = result_df.explode('syntenous_loci').reset_index(drop=True)

# check to make sure leptos still present
# exploded_df["syntenous_loci"].str.contains('SEQF2444').any() 

# merge with previous data
filtarcgenes = pd.merge(arcgenes, exploded_df, left_on="locustag", right_on="syntenous_loci")

# save this to file for troubleshooting purposes
filtarcgenes.to_csv('filtered_arc_genes.txt', index=False, sep="\t")

# now group by homdID and get min/max coordinates for the operon
filt_df = filtarcgenes.groupby("homdID")["length"].sum().reset_index().merge(filtarcgenes.groupby("homdID")["start"].min().reset_index()).merge(filtarcgenes.groupby("homdID")["end"].max().reset_index())
# get actual operon length from start and end coordinates
filt_df["operon_len"] = filt_df["end"] - filt_df["start"]

# first remove extreme outliers (longer than 1.5x the median)
filt_df = filt_df[(filt_df['operon_len'] <= (filt_df['operon_len'].median() * 1.5))]

# calculate summary statistics across operon length column after removing outliers
summary_stats = {
    'mean': filt_df['operon_len'].mean(),
    'median': filt_df['operon_len'].median(),
    'std': filt_df['operon_len'].std(),
    'min': filt_df['operon_len'].min(),
    'max': filt_df['operon_len'].max(),
    '25%': filt_df['operon_len'].quantile(0.25),
    '50%': filt_df['operon_len'].quantile(0.50),  # same as median
    '75%': filt_df['operon_len'].quantile(0.75),
    'count': filt_df['operon_len'].count(),
}
print(summary_stats)
# {'mean': 3569.771777003484, 'median': 3406.0, 'std': 535.5225522927637, 'min': 1954, 'max': 5114, '25%': 3280.0, '50%': 3406.0, '75%': 3460.75, 'count': 1722}

# since data seems to be generally normally distributed (mean ~= median), will use IQR to define length distribution threshold 
iqr = (filt_df["operon_len"].quantile(0.75) - filt_df["operon_len"].quantile(0.25)) * 2

# finally, filter the dataframe based on the 25th and 75th percentiles to remove outliers (plus minus iqr)
q25 = filt_df["operon_len"].quantile(0.25) - iqr
q75 = filt_df["operon_len"].quantile(0.75) + iqr

lenfilt_df = filt_df[(filt_df['operon_len'] >= q25) & (filt_df['operon_len'] <= q75)]

# check that you didn't remove leptos
# lenfilt_df['homdID'].str.contains('SEQF2444').any()

# finally, save to file
lenfilt_df.to_csv('filtered_ads_operon.txt', index=False, sep="\t")
```

Now that we have a filtered operon table, can slice these out of the original fasta file using seqtk

```bash
conda activte 2024-HIV_RNASeq
# mamba install seqtk
# get coordinates in bed like format
awk -F"\t" '{print $1 "\t" $3 "\t" $4}' filtered_ads_operon.txt > filtcoord
# filter operons from homd database
seqtk subseq ~/domhain_RNAseq/03-star_map/homd_db/ALL_genomes.fna filtcoord > ads_operons.fna
```

Generate alignment and clean up to find poorly resolved operons

```bash
# mamba install trimal
# mamba install mafft
# mamba install raxml
# dereplicate
vsearch --derep_fulllength ads_operons.fna \
   --output ads_operons.uniq.fna
# align
mafft --thread 40 \
   --adjustdirectionaccurately \
   ads_operons.uniq.fna > ads_operons.align.fna
# clean up alignment 
trimal -in ads_operons.align.fna \
	-out ads_operons.trim.fna \
	-htmlout ads_operons.trim.html \
	-gt 0.5 \
	-resoverlap 0.5 \
	-seqoverlap 50
# build a tree with our trimmed alignment to make sure it looks ok
rm *tre
raxmlHPC-PTHREADS-SSE3 -T 40 \
	-m GTRCAT \
	-c 25 \
	-e 0.001 \
	-p 31514 \
	-f a \
	-N 100 \
	-x 02938 \
	-n tre -s ads_operons.trim.fna
# get annotation file from sequence identifiers
grep ">" ads_operons.trim.fna | sed 's/>//' > seqids
awk -F"|" '{print $1}' seqids | sed 's/_R_//' | while read line; do grep -m1 -w $line ../03-star_map/homd_db/SEQID_info.txt | awk -F"\t" '{print $3"_"$4}'; done > species
paste seqids species > annotations
# add in header
sed -i "1i seqid\ttaxonomy" annotations
```

Tree looks really good, go ahead and pull the genes from this analysis to run differential expression 

```bash
# get arc genes of interest
awk -F"\t" '{print $4}' filtered_arc_genes.txt | sed 's/ID=//' > arcGene.ids
# pull these from the read counts table from homd
cat arcGene.ids | while read line; do grep $line -w -m1 ~/domhain_RNAseq/03-star_map/02-HOMD_map/featurecounts/read_counts.txt; done > arcGene_read_counts.txt
# add back in header
sed "1i $(head -n1 ~/domhain_RNAseq/03-star_map/02-HOMD_map/featurecounts/read_counts.txt)" arcGene_read_counts.txt > temp
mv temp arcGene_read_counts.txt
```