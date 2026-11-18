import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap

# List of PCOS-related genes
pcos_genes = [
    'ALOX15B', 'LHCGR', 'AMH', 'INSR', 'FSHR', 'ESR2', 'AR', 'CYP11A1', 
    'CYP17A1', 'CYP19A1', 'SHBG', 'THADA', 'DENND1A', 'YAP1', 'HMGA2', 
    'RAB5B', 'SUOX', 'KHDRBS3', 'GATA4', 'NEIL2', 'IRS1', 'IRS2', 'PIK3CA', 
    'FSHB', 'LHB', 'AMHR2', 'INHA', 'INHB', 'IGF1', 'IGF2', 'IGFBP1', 
    'IGFBP2', 'IGFBP3', 'LEP', 'LEPR', 'POMC', 'MC4R', 'FTO', 'PPARG', 'ADIPOQ',
    'FST', 'STAR', 'HSD3B2', 'SRD5A1', 'SRD5A2', 'INPPL1'
]

# Find all TSV files
files = glob.glob('*_annotated.vep_annotated_missense.tsv')

data_list = []
for file in files:
    sample_name = file.split('_')[0]
    df = pd.read_csv(file, sep='\t', usecols=['CHROM', 'POS', 'REF', 'ALT', 'Gene', 'Protein_Position', 'AA_Change'])
    df_pcos = df[df['Gene'].isin(pcos_genes)]
    
    # We want to count unique variants. VEP annotations have multiple transcripts per variant.
    # Group by variant to get a single representative AA change per variant
    def get_aa_str(row):
        try:
            # handle cases where AA_Change might be missing or not in X/Y format
            pos = str(row['Protein_Position']).split('-')[0] # in case of ranges
            aa = str(row['AA_Change'])
            if '/' in aa:
                ref_aa, alt_aa = aa.split('/')
                return f"{ref_aa}{pos}{alt_aa}"
            return f"p.{pos}{aa}"
        except:
            return ""
    
    # Drop duplicates to keep one annotation per variant per gene
    df_unique = df_pcos.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT', 'Gene'])
    
    df_unique = df_unique.copy()
    df_unique['AA_Str'] = df_unique.apply(get_aa_str, axis=1)
    
    df_unique['Sample'] = sample_name
    data_list.append(df_unique)

if not data_list:
    print("No data found or no PCOS genes found in samples.")
    exit()

all_data = pd.concat(data_list, ignore_index=True)

# Count total unique variants per gene to identify top genes (same logic as before)
gene_counts = all_data.groupby('Gene').size().sort_values(ascending=False)
top_genes = gene_counts.head(20).index

plot_data = all_data[all_data['Gene'].isin(top_genes)]

# Now we need to create two pivot tables:
# 1. counts (for coloring the heatmap)
# 2. annotations (the strings of AA changes for the cells)

# Group by Sample and Gene
grouped = plot_data.groupby(['Sample', 'Gene'])

# 1. Counts
counts_df = grouped.size().reset_index(name='Count')
pivot_counts = counts_df.pivot(index='Gene', columns='Sample', values='Count').fillna(0)

# 2. Annotations (Combined Count + AAs)
def create_combined_annot(group):
    count = len(group)
    aas = [x for x in group['AA_Str'] if x]
    # Join AA changes
    joined_aas = ", ".join(aas)
    # Wrap AA changes text
    wrapped_aas = "\n".join(textwrap.wrap(joined_aas, width=35))
    # Combine Count and AAs (Count on top)
    if count > 0:
        return f"{count}\n({wrapped_aas})" if wrapped_aas else f"{count}"
    return ""

annots_df = grouped.apply(create_combined_annot).reset_index(name='Combined_Annot')
pivot_annots = annots_df.pivot(index='Gene', columns='Sample', values='Combined_Annot').fillna("")

# Ensure alignment of rows and columns
pivot_counts = pivot_counts.reindex(index=top_genes).fillna(0)
pivot_annots = pivot_annots.reindex(index=top_genes, columns=pivot_counts.columns).fillna("")

sns.set_theme(style="whitegrid")

plt.figure(figsize=(38, 34)) # Adjusted to requested size
# Heatmap uses pivot_counts for color, pivot_annots for text
ax = sns.heatmap(pivot_counts, cmap='YlOrRd', annot=pivot_annots, fmt='s', linewidths=.5, 
                 cbar_kws={'label': 'Mutation Count'},
                 annot_kws={"size": 24, "weight": "bold", "va": "center"}) 

# Set colorbar label font size and weight
cbar = ax.collections[0].colorbar
cbar.set_label('Mutation Count', fontsize=38, weight='bold')
cbar.ax.tick_params(labelsize=30)

plt.title('Heatmap of Mutated PCOS Genes: Counts and AA Changes', fontsize=52, weight='bold', pad=35)
plt.xlabel('Sample', fontsize=44, weight='bold', labelpad=30)
plt.ylabel('Gene', fontsize=44, weight='bold', labelpad=30)

# Customize tick label sizes and weight
plt.xticks(rotation=0, ha='center', fontsize=34, weight='bold')
plt.yticks(rotation=0, fontsize=34, weight='bold')

plt.tight_layout()
plt.savefig('pcos_gene_aa_mutation_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved pcos_gene_aa_mutation_heatmap.png")
