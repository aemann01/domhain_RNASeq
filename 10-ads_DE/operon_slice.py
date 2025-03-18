#!/usr/bin/env python3

'''TODO: DESCRIPTION
'''

# import libraries
import pandas as pd


# arguments
import argparse
parser = argparse.ArgumentParser()
requireparser = parser.add_argument_group('required arguments')
requireparser.add_argument('-i', '--input', help='Input file. See example for formatting requirements.')
requireparser.add_argument('-g', '--genes', nargs='+', type=str)

args = parser.parse_args()

# input variables
inputf = pd.read_csv(args.input, header=0, sep="\t")
req_genes = {args.accumulate(args.genes)}


# check that all required genes are present for contig to pass to next filtering step




# get length of each concatenated gene
arcgenes["length"] = abs(arcgenes["end"] - arcgenes["start"])
# add gene names
arcgenes.loc[arcgenes['EC'] == "eC_number=3.5.3.6", 'EC'] = 'arcA'
arcgenes.loc[arcgenes['EC'] == "eC_number=2.1.3.3", 'EC'] = 'arcB'
arcgenes.loc[arcgenes['EC'] == "eC_number=2.7.2.2", 'EC'] = 'arcC'


def check_genes(homdID):
    return ads_genes.issubset(set(homdID['EC']))

genespres = arcgenes.groupby('homdID').filter(check_genes)

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

# merge with previous data
filtarcgenes = pd.merge(arcgenes, exploded_df, left_on="locustag", right_on="syntenous_loci")

# save this to file for troubleshooting purposes
filtarcgenes.to_csv('filtered_arc_genes.txt', index=False, sep="\t")

# now group by homdID and get min/max coordinates for the operon
filt_df = filtarcgenes.groupby("homdID")["length"].sum().reset_index().merge(filtarcgenes.groupby("homdID")["start"].min().reset_index()).merge(filtarcgenes.groupby("homdID")["end"].max().reset_index())
# get actual operon length from start and end coordinates
filt_df["operon_len"] = filt_df["end"] - filt_df["start"]

# calculate summary statistics across operon length column
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
# {'mean': 10640.963326446281, 'median': 3410.0, 'std': 146077.62077091116, 'min': 1954, 'max': 4071749, '25%': 3298.0, '50%': 3410.0, '75%': 3728.0, 'count': 1936}

# since the data is not normally distributed (mean is much higher than median) -- will use interquartile range to define our length variation

# finally, filter the dataframe based on the 25th and 75th percentiles to remove outliers
q25 = filt_df["operon_len"].quantile(0.25)
q75 = filt_df["operon_len"].quantile(0.75)

lenfilt_df = filt_df[(filt_df['operon_len'] >= q25) & (filt_df['operon_len'] <= q75)]

# finally, save to file
lenfilt_df.to_csv('filtered_ads_operon.txt', index=False, sep="\t")