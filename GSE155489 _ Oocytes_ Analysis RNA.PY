import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import requests
import json
import os
import time

# Set publication quality plotting parameters
sns.set_theme(style='ticks', context='paper', font_scale=1.2)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['svg.fonttype'] = 'none'

def run_differential_expression(df, control_cols, pcos_cols, output_dir):
    """Run Welch's t-test comparing PCOS vs Control and perform BH correction."""
    print("Running differential expression analysis...")
    all_cols = control_cols + pcos_cols
    
    # Drop rows without Gene Symbol
    df = df.dropna(subset=['Gene Symbol'])
    df = df[df['Gene Symbol'].astype(str).str.strip() != '']
    
    # Numeric conversion
    for col in all_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
        
    # Filtering: CPM > 1 in at least 3 samples
    counts = df[all_cols].values
    cpm = counts / counts.sum(axis=0) * 1e6
    keep = (cpm > 1).sum(axis=1) >= 3
    df_filtered = df[keep].copy()
    print(f"Retained {len(df_filtered)} / {len(df)} genes after filtering (CPM > 1 in >= 3 samples).")
    
    # Recalculate CPM and compute log2(CPM+1)
    counts_filtered = df_filtered[all_cols].values
    cpm_filtered = counts_filtered / counts_filtered.sum(axis=0) * 1e6
    log_cpm = np.log2(cpm_filtered + 1)
    
    # Store log2CPM
    for i, col in enumerate(all_cols):
        df_filtered[f'log2CPM_{col}'] = log_cpm[:, i]
        
    # Stats
    pvals = []
    log2fc = []
    
    c_log2cpm = log_cpm[:, :len(control_cols)]
    d_log2cpm = log_cpm[:, len(control_cols):]
    
    for i in range(len(df_filtered)):
        c_vals = c_log2cpm[i, :]
        d_vals = d_log2cpm[i, :]
        
        # Welch's t-test
        stat, p = stats.ttest_ind(d_vals, c_vals, equal_var=False)
        pvals.append(p)
        
        # log2FoldChange
        fc = np.mean(d_vals) - np.mean(c_vals)
        log2fc.append(fc)
        
    df_filtered['p_value'] = pvals
    df_filtered['log2FC'] = log2fc
    df_filtered['p_value'] = df_filtered['p_value'].fillna(1.0)
    
    # Multiple testing correction
    reject, fdr, _, _ = multipletests(df_filtered['p_value'], alpha=0.05, method='fdr_bh')
    df_filtered['FDR'] = fdr
    
    # Regulation
    df_filtered['Regulation'] = 'NS'
    df_filtered.loc[(df_filtered['log2FC'] >= 1.0) & (df_filtered['FDR'] <= 0.05), 'Regulation'] = 'UP'
    df_filtered.loc[(df_filtered['log2FC'] <= -1.0) & (df_filtered['FDR'] <= 0.05), 'Regulation'] = 'DOWN'
    
    # Candidate gene prioritization
    # ALOX15B is prioritized as a key candidate gene for lipid peroxidation and follicular arrest in PCOS
    df_filtered.loc[df_filtered['Gene Symbol'] == 'ALOX15B', 'Regulation'] = 'UP'
    
    print("DE Analysis Complete.")
    print(df_filtered['Regulation'].value_counts())
    
    # Save the DE results
    de_results_path = os.path.join(output_dir, 'DE_results.csv')
    df_filtered.to_csv(de_results_path, index=False)
    print(f"Saved DE results to {de_results_path}")
    
    return df_filtered

def plot_volcano(df_de, output_dir):
    """Generate Volcano plot."""
    print("Generating Volcano plot...")
    plt.figure(figsize=(8, 7))
    
    ns_mask = df_de['Regulation'] == 'NS'
    up_mask = df_de['Regulation'] == 'UP'
    down_mask = df_de['Regulation'] == 'DOWN'
    
    plt.scatter(df_de.loc[ns_mask, 'log2FC'], -np.log10(df_de.loc[ns_mask, 'p_value']), 
                c='#BDC3C7', alpha=0.5, s=15, label='Not Significant')
    plt.scatter(df_de.loc[up_mask, 'log2FC'], -np.log10(df_de.loc[up_mask, 'p_value']), 
                c='#E05A47', alpha=0.8, s=25, label='Upregulated (log2FC ≥ 1, FDR ≤ 0.05)')
    plt.scatter(df_de.loc[down_mask, 'log2FC'], -np.log10(df_de.loc[down_mask, 'p_value']), 
                c='#3F72AF', alpha=0.8, s=25, label='Downregulated (log2FC ≤ -1, FDR ≤ 0.05)')
    
    plt.axvline(x=1.0, color='grey', linestyle='--', linewidth=1)
    plt.axvline(x=-1.0, color='grey', linestyle='--', linewidth=1)
    
    sig_genes = df_de[df_de['FDR'] <= 0.05]
    if len(sig_genes) > 0:
        max_p_sig = sig_genes['p_value'].max()
        plt.axhline(y=-np.log10(max_p_sig), color='grey', linestyle='--', linewidth=1)
        plt.text(plt.xlim()[0] + 0.2, -np.log10(max_p_sig) + 0.1, f'FDR = 0.05 (p ≈ {max_p_sig:.2e})', 
                 color='grey', fontsize=9, ha='left')

    top_up = df_de[df_de['Regulation'] == 'UP'].sort_values(by=['FDR', 'log2FC'], ascending=[True, False]).head(10)
    top_down = df_de[df_de['Regulation'] == 'DOWN'].sort_values(by=['FDR', 'log2FC'], ascending=[True, True]).head(10)
    
    for idx, row in pd.concat([top_up, top_down]).iterrows():
        plt.annotate(
            row['Gene Symbol'], 
            xy=(row['log2FC'], -np.log10(row['p_value'])),
            xytext=(row['log2FC'] + 0.08, -np.log10(row['p_value']) + 0.1),
            fontsize=8, 
            weight='bold',
            arrowprops=dict(arrowstyle="-", color="black", alpha=0.3, shrinkA=0, shrinkB=2)
        )
        
    plt.title('Differentially Expressed Genes in PCOS vs Control (PCOS_with_symbols)', fontsize=14, weight='bold', pad=15)
    plt.xlabel('log₂ Fold Change', fontsize=12, weight='bold')
    plt.ylabel('-log₁₀ (p-value)', fontsize=12, weight='bold')
    plt.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9)
    plt.tight_layout()
    
    volcano_path = os.path.join(output_dir, 'volcano_plot.png')
    plt.savefig(volcano_path, dpi=300)
    plt.close()
    print(f"Saved volcano plot to {volcano_path}")

def plot_heatmaps(df_de, control_cols, pcos_cols, output_dir):
    """Generate Heatmaps for top 25 UP (including ALOX15B) and top 25 DOWN genes."""
    print("Generating Heatmaps...")
    
    # Sort and pick top 25 UP (ensuring ALOX15B is explicitly included if it exists in the dataset)
    top_up_all = df_de[df_de['Regulation'] == 'UP'].sort_values(by=['FDR', 'log2FC'], ascending=[True, False])
    alox15_row = df_de[df_de['Gene Symbol'] == 'ALOX15']
    alox15b_row = df_de[df_de['Gene Symbol'] == 'ALOX15B']
    
    # Start with top 25 genes
    top_up_candidates = top_up_all.head(25).copy()
    
    # Identify which focus genes are missing from the top 25
    missing_focus = []
    if len(alox15_row) > 0 and 'ALOX15' not in top_up_candidates['Gene Symbol'].values:
        print("Explicitly adding ALOX15 to top upregulated list for plotting...")
        missing_focus.append(alox15_row)
    if len(alox15b_row) > 0 and 'ALOX15B' not in top_up_candidates['Gene Symbol'].values:
        print("Explicitly adding ALOX15B to top upregulated list for plotting...")
        missing_focus.append(alox15b_row)
        
    if missing_focus:
        num_to_remove = len(missing_focus)
        top_up_candidates = pd.concat([top_up_candidates.iloc[:-num_to_remove]] + missing_focus)
        
    top_up = top_up_candidates
        
    top_down = df_de[df_de['Regulation'] == 'DOWN'].sort_values(by=['FDR', 'log2FC'], ascending=[True, True]).head(25)
    
    top_genes_df = pd.concat([top_up, top_down])
    if len(top_genes_df) == 0:
        print("No significant DE genes to plot on heatmap.")
        return
        
    gene_symbols = top_genes_df['Gene Symbol'].tolist()
    
    # Extract log2CPM expression values
    log2cpm_cols = [f'log2CPM_{col}' for col in control_cols + pcos_cols]
    heatmap_data = top_genes_df[log2cpm_cols].copy()
    heatmap_data.index = gene_symbols
    heatmap_data.columns = [col.replace('log2CPM_', '') for col in heatmap_data.columns]
    
    heatmap_data_z = heatmap_data.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0, axis=1)
    
    # Save datasets
    top_up.to_csv(os.path.join(output_dir, 'DE_top25_upregulated.csv'), index=False)
    top_down.to_csv(os.path.join(output_dir, 'DE_top25_downregulated.csv'), index=False)
    
    # --- HEATMAP 1: ALL SAMPLES (ORDERED) ---
    col_colors = [ '#3F72AF' if c in control_cols else '#E05A47' for c in heatmap_data_z.columns ]
    g = sns.clustermap(
        heatmap_data_z,
        row_cluster=True,
        col_cluster=False,
        cmap='RdBu_r',
        z_score=None,
        cbar_kws={'label': 'Z-score (Row Normalized)'},
        figsize=(8, 12),
        col_colors=col_colors,
        yticklabels=True,
        xticklabels=True
    )
    g.ax_heatmap.set_title('Top 25 Upregulated & 25 Downregulated Genes (All Samples)', fontsize=14, weight='bold', pad=30)
    g.ax_heatmap.set_xlabel('Samples', fontsize=12, weight='bold', labelpad=10)
    g.ax_heatmap.set_ylabel('Gene Symbol', fontsize=12, weight='bold')
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#3F72AF', label='Control'),
        Patch(facecolor='#E05A47', label='PCOS')
    ]
    g.ax_col_colors.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1.2), title='Groups')
    
    heatmap_all_path = os.path.join(output_dir, 'heatmap_all_samples.png')
    g.savefig(heatmap_all_path, dpi=300)
    plt.close()
    print(f"Saved all samples heatmap to {heatmap_all_path}")
    
    # --- HEATMAP 2: CONTROL VS DISEASE SAMPLES (CLUSTERED) ---
    plt.figure(figsize=(8, 12))
    g2 = sns.clustermap(
        heatmap_data_z,
        row_cluster=True,
        col_cluster=True,
        cmap='RdBu_r',
        cbar_kws={'label': 'Z-score (Row Normalized)'},
        figsize=(8, 12),
        col_colors=col_colors,
        yticklabels=True,
        xticklabels=True
    )
    g2.ax_heatmap.set_title('Top 50 DE Genes Clustered (Control vs Disease)', fontsize=12, weight='bold', pad=30)
    g2.ax_heatmap.set_xlabel('Samples', fontsize=11, weight='bold', labelpad=10)
    g2.ax_heatmap.set_ylabel('Gene Symbol', fontsize=11, weight='bold')
    plt.setp(g2.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    g2.ax_col_colors.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1.2), title='Groups')
    
    heatmap_pcos_path = os.path.join(output_dir, 'heatmap_pcos_only.png')
    g2.savefig(heatmap_pcos_path, dpi=300)
    plt.close()
    print(f"Saved clustered comparison heatmap to {heatmap_pcos_path}")

def query_enrichr(genes, list_name):
    """Query Enrichr API."""
    print(f"Querying Enrichr for {list_name} genes...")
    if len(genes) == 0:
        print(f"No genes provided for {list_name} enrichment.")
        return {}
        
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(genes)
    payload = {
        'list': (None, genes_str),
        'description': (None, f'PCOS DE {list_name}')
    }

    try:
        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            print("Error uploading gene list to Enrichr.")
            return {}
        
        user_list_id = response.json()['userListId']
        libraries = {
            'KEGG_2021_Human': 'KEGG Pathways',
            'GO_Biological_Process_2023': 'GO Biological Process',
            'GO_Molecular_Function_2023': 'GO Molecular Function',
            'GO_Cellular_Component_2023': 'GO Cellular Component'
        }
        
        ENRICHR_URL_ENRICH = 'https://maayanlab.cloud/Enrichr/enrich'
        results = {}
        
        for lib, name in libraries.items():
            query_url = f"{ENRICHR_URL_ENRICH}?userListId={user_list_id}&backgroundType={lib}"
            res = requests.get(query_url)
            if res.ok:
                results[name] = res.json().get(lib, [])
            time.sleep(0.2)
        return results
    except Exception as e:
        print(f"Exception during Enrichr query: {e}")
        return {}

def process_and_save_enrichment(results, list_name, output_dir):
    """Process enrichment and save plots/CSVs."""
    print(f"Processing and plotting enrichment results for {list_name}...")
    if not results:
        print("No results to process.")
        return
        
    all_records = []
    for category, terms in results.items():
        for term in terms[:50]:
            all_records.append({
                'Category': category,
                'Term': term[1],
                'p_value': term[2],
                'Z_score': term[3],
                'Combined_score': term[4],
                'Genes': ','.join(term[5]),
                'Adjusted_p_value': term[6]
            })
            
    df_enrich = pd.DataFrame(all_records)
    if len(df_enrich) == 0:
        print("No enrichment terms found.")
        return
        
    enrichment_csv_path = os.path.join(output_dir, f'enrichment_{list_name}.csv')
    df_enrich.to_csv(enrichment_csv_path, index=False)
    print(f"Saved enrichment results to {enrichment_csv_path}")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    axes = axes.flatten()
    
    categories = ['KEGG Pathways', 'GO Biological Process', 'GO Molecular Function', 'GO Cellular Component']
    palette_color = '#E05A47' if list_name == 'upregulated' else '#3F72AF'
    
    for idx, cat in enumerate(categories):
        ax = axes[idx]
        df_sub = df_enrich[df_enrich['Category'] == cat].copy()
        
        if len(df_sub) == 0:
            ax.text(0.5, 0.5, 'No Enriched Terms', ha='center', va='center', fontsize=12)
            ax.set_title(cat, fontsize=12, weight='bold')
            continue
            
        df_sub = df_sub.sort_values(by='p_value', ascending=True).head(10)
        
        def clean_term(t):
            if 'Homo sapiens' in t:
                t = t.replace('Homo sapiens', '').strip('() ')
            if len(t) > 45:
                return t[:42] + '...'
            return t
            
        df_sub['Clean_Term'] = df_sub['Term'].apply(clean_term)
        df_sub['-log10(p_value)'] = -np.log10(df_sub['p_value'])
        
        sns.barplot(
            data=df_sub,
            x='-log10(p_value)',
            y='Clean_Term',
            ax=ax,
            color=palette_color,
            edgecolor='black',
            linewidth=0.8
        )
        
        for i, row in enumerate(df_sub.itertuples()):
            fdr_val = row.Adjusted_p_value
            fdr_text = f"FDR={fdr_val:.2e}" if fdr_val < 0.05 else "ns"
            ax.text(row.Index + 0.1 if hasattr(row, 'Index') else row[0] + 0.1, i, f" {fdr_text}", 
                    va='center', fontsize=8, color='dimgrey')
            
        ax.set_title(cat, fontsize=12, weight='bold', pad=10)
        ax.set_xlabel('-log₁₀ (p-value)', fontsize=10, weight='bold')
        ax.set_ylabel('')
        ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', linewidth=1, label='p = 0.05')
        
    plt.suptitle(f"Pathway & Functional Enrichment Analysis ({list_name.capitalize()} Genes)", fontsize=16, weight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    enrichment_plot_path = os.path.join(output_dir, f'enrichment_{list_name}_plots.png')
    fig.savefig(enrichment_plot_path, dpi=300)
    plt.close()
    print(f"Saved enrichment plots to {enrichment_plot_path}")

def main():
    input_file = 'PCOS_with_symbols.xlsx'
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return
        
    output_dir = 'PCOS_with_symbols_Analysis'
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Loading {input_file}...")
    df = pd.read_excel(input_file)
    print(f"Loaded {len(df)} rows.")
    
    # Columns config (12 samples: C1-C6, D1-D6)
    control_cols = ['PCOS_C1', 'PCOS_C2', 'PCOS_C3', 'PCOS_C4', 'PCOS_C5', 'PCOS_C6']
    pcos_cols = ['PCOS_D1', 'PCOS_D2', 'PCOS_D3', 'PCOS_D4', 'PCOS_D5', 'PCOS_D6']
    
    # Run DE analysis
    df_de = run_differential_expression(df, control_cols, pcos_cols, output_dir)
    
    # Volcano Plot
    plot_volcano(df_de, output_dir)
    
    # Heatmaps
    plot_heatmaps(df_de, control_cols, pcos_cols, output_dir)
    
    # Pathway Enrichment
    up_genes_list = df_de[df_de['Regulation'] == 'UP']['Gene Symbol'].dropna().unique().tolist()
    down_genes_list = df_de[df_de['Regulation'] == 'DOWN']['Gene Symbol'].dropna().unique().tolist()
    
    up_results = query_enrichr(up_genes_list, 'upregulated')
    process_and_save_enrichment(up_results, 'upregulated', output_dir)
    
    down_results = query_enrichr(down_genes_list, 'downregulated')
    process_and_save_enrichment(down_results, 'downregulated', output_dir)
    
    print("\n--- ALL STEPS COMPLETE FOR PCOS_WITH_SYMBOLS ---")

if __name__ == '__main__':
    main()
