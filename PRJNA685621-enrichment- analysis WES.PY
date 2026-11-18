import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

# Exclude terms containing any of these keywords (infections, cancers, unrelated autoimmune)
exclude_keywords = [
    'infection', 'influenza', 'leishmaniasis', 'tuberculosis', 'amoebiasis', 
    'malaria', 'myocarditis', 'pertussis', 'graft-versus-host', 'lupus', 
    'arthritis', 'bowel disease', 'hepatitis', 'measles', 'hiv', 'covid', 
    'herpes', 'cancer', 'carcinoma', 'melanoma', 'glioma', 'leukemia', 
    'lymphoma', 'toxoplasmosis', 'salmonella', 'shigellosis', 'chagas', 
    'trypanosomiasis', 'vibrio', 'pneumonia', 'colitis', 'asthma', 
    'alzheimer', 'huntington', 'parkinson', 'prion', 'coccal'
]

# Set publication plotting parameters
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

def clean_term_name(t):
    # Shorten long names and remove Homo sapiens suffix
    if 'Homo sapiens' in t:
        t = t.replace('Homo sapiens', '').strip('() ')
    # Remove GO ID if present
    if ' (GO:' in t:
        t = t.split(' (GO:')[0]
    if len(t) > 55:
        return t[:52] + '...'
    return t

def generate_plots_for_dataset(csv_path, output_dir):
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found.")
        return
        
    df = pd.read_csv(csv_path)
    # Calculate Gene Count
    df['Count'] = df['Genes'].apply(lambda x: len(str(x).split(',')) if pd.notnull(x) else 0)
    
    categories = {
        'KEGG Pathways': ('kegg', 'PCOS-Related KEGG Pathway Enrichment'),
        'GO Biological Process': ('go_bp', 'PCOS-Related GO Biological Process (BP) Enrichment'),
        'GO Molecular Function': ('go_mf', 'PCOS-Related GO Molecular Function (MF) Enrichment'),
        'GO Cellular Component': ('go_cc', 'PCOS-Related GO Cellular Component (CC) Enrichment')
    }
    
    os.makedirs(output_dir, exist_ok=True)
    
    for cat_name, (file_prefix, title_prefix) in categories.items():
        df_cat = df[df['Category'] == cat_name].copy()
        
        # Filter terms
        mask = df_cat['Term'].apply(lambda x: not any(kw in x.lower() for kw in exclude_keywords))
        df_filtered = df_cat[mask].copy()
        
        if len(df_filtered) == 0:
            print(f"Warning: No terms remain for category '{cat_name}' in {csv_path} after filtering.")
            continue
            
        # Get top 15 terms (or all if fewer)
        df_sub = df_filtered.sort_values(by='p_value', ascending=True).head(15).copy()
        df_sub['Clean_Term'] = df_sub['Term'].apply(clean_term_name)
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 6.5))
        
        # Color gradient: Red (most significant, low Adjusted_p_value) to Blue (less significant, high Adjusted_p_value)
        fdr_min = df_sub['Adjusted_p_value'].min()
        fdr_max = df_sub['Adjusted_p_value'].max()
        
        # Prevent division by zero if all values are identical
        if fdr_min == fdr_max:
            norm = mcolors.Normalize(vmin=fdr_min - 0.001, vmax=fdr_min + 0.001)
        else:
            norm = mcolors.Normalize(vmin=fdr_min, vmax=fdr_max)
            
        cmap = mcolors.LinearSegmentedColormap.from_list("red_to_blue", ["#E05A47", "#3F72AF"])
        colors = [cmap(norm(val)) for val in df_sub['Adjusted_p_value']]
        
        bars = ax.barh(df_sub['Clean_Term'], df_sub['Count'], color=colors, edgecolor='black', linewidth=0.6)
        ax.invert_yaxis()  # Most significant term at the top
        
        ax.set_xlabel('Gene Count', fontsize=11, weight='bold')
        ax.set_title(f"{title_prefix}\n(Upregulated Genes)", fontsize=12, weight='bold', pad=15)
        ax.grid(axis='x', linestyle='--', alpha=0.5)
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02, shrink=0.75)
        cbar.set_label('p.adjust (FDR)', fontsize=10, weight='bold')
        
        # Fine-tune tick parameters
        ax.tick_params(axis='y', labelsize=9.5)
        
        plt.tight_layout()
        
        # Save locally in dataset analysis directory
        local_filename = f"{file_prefix}_upregulated.png"
        local_path = os.path.join(output_dir, local_filename)
        plt.savefig(local_path, dpi=300)
        plt.close()
        
        print(f"Saved: {local_path}")

def resolve_path(folder_name):
    """Find analysis directory relative to script or working directory."""
    candidates = [
        folder_name,
        os.path.join("..", folder_name),
        os.path.join(os.path.dirname(__file__), folder_name),
        os.path.join(os.path.dirname(__file__), "..", folder_name),
    ]
    for c in candidates:
        if os.path.exists(c):
            return os.path.abspath(c)
    return folder_name

def main():
    geneid_dir = resolve_path("GeneID_Analysis")
    symbols_dir = resolve_path("PCOS_with_symbols_Analysis")
    
    # 1. Dataset 1: Gene-ID (8-sample dataset)
    print("\n--- Plotting for Gene-ID (8-sample dataset) ---")
    generate_plots_for_dataset(
        csv_path=os.path.join(geneid_dir, "enrichment_upregulated.csv"),
        output_dir=geneid_dir
    )
    
    # 2. Dataset 2: PCOS with symbols (12-sample dataset)
    print("\n--- Plotting for PCOS with symbols (12-sample dataset) ---")
    generate_plots_for_dataset(
        csv_path=os.path.join(symbols_dir, "enrichment_upregulated.csv"),
        output_dir=symbols_dir
    )
    print("\nAll enrichment plots successfully generated.")

if __name__ == '__main__':
    main()
