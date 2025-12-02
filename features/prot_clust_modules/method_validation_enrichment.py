import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")

"""
Module: validation_analysis.py
Description: 
    Iterates over clusters, performs ORA using curated Pathway Databases only 
    (Reactome, WikiPathways, KEGG) to ensure specific, non-redundant results.
    
    Adjustments for strictness:
    1. Sources: REAC, WP, KEGG (No GO:BP).
    2. Correction: Bonferroni (Stricter than FDR).
    3. Filters: Term size 5-300 + Deduplication.
"""

import pandas as pd
import pickle
import os
import matplotlib.pyplot as plt

# Try importing gprofiler-official
try:
    from gprofiler import GProfiler
    GPROFILER_AVAILABLE = True
except ImportError:
    GPROFILER_AVAILABLE = False

def run_validation_analysis(pickle_path, output_dir, clust_algo_used="unknown"):
    """
    Validates clustering quality using curated pathway databases.
    """

    if not GPROFILER_AVAILABLE:
        print("❌ Error: 'gprofiler-official' library not installed.")
        return None

    print(f"\n--- 🧪 Validation Analysis (Curated Pathways): {clust_algo_used} ---")
    print("--- ⚙️  Sources: Reactome (REAC), WikiPathways (WP), KEGG ---")
    print("--- ⚙️  Correction: Bonferroni (Strict) | P < 0.05 ---")

    # --- 1. Load Data ---
    if not os.path.exists(pickle_path):
        print(f"❌ Error: Pickle file not found at {pickle_path}")
        return None

    try:
        with open(pickle_path, "rb") as f:
            data = pickle.load(f)
        
        communities = data['communities']
        G = data['graph']
        background_gene_list = list(G.nodes())
        print(f"[INFO] Background universe size: {len(background_gene_list)} genes.")

    except Exception as e:
        print(f"❌ Error loading pickle: {e}")
        return None

    # --- 2. Initialize g:Profiler ---
    gp = GProfiler(return_dataframe=True)
    validation_results = []

    print("[INFO] Starting iteration over clusters...")

    # --- 3. Iteration Loop ---
    for i, cluster_set in enumerate(communities):
        cluster_gene_list = list(cluster_set)
        n_total_proteins = len(cluster_gene_list)

        # Skip tiny clusters
        if n_total_proteins < 3:
            continue

        try:
            # ---------------------------------------------------------
            # CHANGE 1: Sources & Correction Method
            # ---------------------------------------------------------
            results = gp.profile(
                organism='hsapiens',
                query=cluster_gene_list,
                background=background_gene_list,
                # ONLY use curated pathway DBs. Removing GO:BP reduces noise significantly.
                sources=['REAC'], 
                user_threshold=0.05, 
                # Bonferroni is stricter than FDR. 
                # It divides alpha by the number of tests. Good for reducing count.
                significance_threshold_method='bonferroni', 
                no_evidences=False 
            )

            n_sig_pathways = 0
            unique_sig_proteins = set()

            if not results.empty:
                # ---------------------------------------------------------
                # CHANGE 2: Filter by Term Size (Size 5 to 300)
                # ---------------------------------------------------------
                # Pathways > 300 genes are usually too broad (e.g. "Metabolism")
                filtered_results = results[
                    (results['term_size'] >= 5) & 
                    (results['term_size'] <= 100)
                ].copy()

                col_name = 'intersections' if 'intersections' in filtered_results.columns else 'intersection'
                
                # ---------------------------------------------------------
                # CHANGE 3: Deduplication
                # ---------------------------------------------------------
                if not filtered_results.empty and col_name in filtered_results.columns:
                    def get_gene_signature(val):
                        if isinstance(val, list): return tuple(sorted(val))
                        elif isinstance(val, dict): return tuple(sorted(val.keys()))
                        return str(val)

                    filtered_results['gene_signature'] = filtered_results[col_name].apply(get_gene_signature)
                    # Keep only the single most significant pathway for any specific set of genes
                    filtered_results = filtered_results.sort_values('p_value').drop_duplicates(subset='gene_signature')

                n_sig_pathways = len(filtered_results)

                # Extract genes from filtered results
                if n_sig_pathways > 0 and col_name in filtered_results.columns:
                    for item in filtered_results[col_name]:
                        if isinstance(item, list): unique_sig_proteins.update(item)
                        elif isinstance(item, dict): unique_sig_proteins.update(item.keys())
                        elif isinstance(item, str): unique_sig_proteins.add(item)

            # Calculate Stats
            percentage = (len(unique_sig_proteins) / n_total_proteins) * 100

            validation_results.append({
                'Cluster_ID': i,
                'Total_Proteins': n_total_proteins,
                'Significant_Pathways': n_sig_pathways,
                'Unique_Significant_Proteins': len(unique_sig_proteins),
                'Percentage_Involved': percentage
            })
            
            if (i + 1) % 10 == 0:
                print(f"   ... Processed {i + 1}/{len(communities)} clusters")

        except Exception as e:
            print(f"⚠️ Cluster {i} Error: {e}")

    # --- 4. Save & Plot ---
    if not validation_results:
        print("❌ No significant results found with these strict settings.")
        return None

    df_validation = pd.DataFrame(validation_results)
    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, f"validation_curated_{clust_algo_used}.csv")
    df_validation.to_csv(csv_path, index=False)
    print(f"\n✅ Results saved: {csv_path}")

    # Stats
    avg_pct = df_validation['Percentage_Involved'].mean()
    avg_path = df_validation['Significant_Pathways'].mean()
    print(f"\n--- 📊 Global Stats ({clust_algo_used}) ---")
    print(f"Avg % Proteins in Pathways: {avg_pct:.2f}%")
    print(f"Avg # Pathways per Cluster: {avg_path:.2f}")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.hist(df_validation['Percentage_Involved'], bins=20, color='#69b3a2', edgecolor='black', alpha=0.8)
    plt.axvline(avg_pct, color='red', linestyle='dashed', label=f'Mean: {avg_pct:.1f}%')
    plt.title(f'Cluster Validation (Curated Pathways Only)\n{clust_algo_used} | Bonferroni | REAC+WP+KEGG')
    plt.xlabel('% Proteins in Significant Pathways')
    plt.ylabel('Number of Clusters')
    plt.legend()
    
    plot_path = os.path.join(output_dir, f"validation_plot_curated_{clust_algo_used}.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    return df_validation