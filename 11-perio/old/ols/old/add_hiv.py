import pandas as pd

metadata_file = "../../metadata.txt"  
metadata_df = pd.read_csv(metadata_file, sep='\t')

hiv_status_dict = dict(zip(metadata_df['sample_id'], metadata_df['hiv_status']))

abundance_file = "hmp_subset_pathabundance-cpm.tsv"  
abundance_df = pd.read_csv(abundance_file, sep='\t')

sample_ids = [col.split('_')[0] for col in abundance_df.columns[1:]]  

hiv_status_row = ['hiv_status3  '] 
for sample_id in sample_ids:
    hiv_status_row.append(hiv_status_dict.get(sample_id, 'Unknown')) 

# Add hiv
abundance_df.loc[-1] = hiv_status_row  
abundance_df.index = abundance_df.index + 1  
abundance_df.sort_index(inplace=True)  

abundance_df.to_csv("updated_abundance.tsv", index=False, sep='\t')

print("HIV status has been added successfully!")
