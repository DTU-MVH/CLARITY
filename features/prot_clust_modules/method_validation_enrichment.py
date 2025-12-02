import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")

"""
Module: validation_analysis.py
Description: 
    Iterates over ALL clusters in a clustering result, performs ORA using 
    gprofiler-official, and calculates quality metrics (pathway coverage).
    Generates a summary CSV and a histogram of protein coverage.
"""
import pandas as pd
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np

# Try importing gprofiler-official
try:
    from gprofiler import GProfiler
    GPROFILER_AVAILABLE = True
except ImportError:
    GPROFILER_AVAILABLE = False

def run_validation_analysis(pickle_path, output_dir, clust_algo_used="unknown"):
    """
    Validates clustering quality by running ORA on all clusters.

    Parameters:
        pickle_path (str): Path to the .pkl file containing 'communities' and 'graph'.
        output_dir (str): Directory where results and plots will be saved.
        clust_algo_used (str): Name of the algorithm (for file naming).

    Returns:
        pd.DataFrame: Summary dataframe of validation results.
    """

    if not GPROFILER_AVAILABLE:
        print("❌ Error: 'gprofiler-official' library not installed.")
        print("   Run: pip install gprofiler-official")
        return None

    print(f"\n--- 🧪 Validation Analysis (g:Profiler): {clust_algo_used} ---")

    # --- 1. Load Data ---
    if not os.path.exists(pickle_path):
        print(f"❌ Error: Pickle file not found at {pickle_path}")
        return None

    try:
        with open(pickle_path, "rb") as f:
            data = pickle.load(f)
        
        communities = data['communities']
        G = data['graph']
        
        # Define Background (Universe)
        background_gene_list = list(G.nodes())
        print(f"[INFO] Loaded {len(communities)} clusters.")
        print(f"[INFO] Background universe size: {len(background_gene_list)} genes.")

    except Exception as e:
        print(f"❌ Error loading pickle: {e}")
        return None

    # --- 2. Initialize g:Profiler ---
    gp = GProfiler(return_dataframe=True)
    validation_results = []

    print("[INFO] Starting iteration over clusters (this may take time)...")

    # --- 3. Iteration Loop ---
    for i, cluster_set in enumerate(communities):
        # Convert to list
        cluster_gene_list = list(cluster_set)
        n_total_proteins = len(cluster_gene_list)

        # Skip tiny clusters (optional, but good for noise reduction)
        if n_total_proteins < 3:
            continue

        try:
            # Run ORA
            # We use no_evidences=False to ensure we get the 'intersections' column
            results = gp.profile(
                organism='hsapiens',
                query=cluster_gene_list,
                background=background_gene_list,
                sources=['GO:BP', 'REAC'], 
                user_threshold=0.05,
                significance_threshold_method='fdr',
                no_evidences=False 
            )

            n_sig_pathways = 0
            unique_sig_proteins = set()

            if not results.empty:
                # 'p_value' is already the adjusted p-value in g:Profiler outputs
                # Just to be safe, filter again
                sig_results = results[results['p_value'] < 0.05]
                n_sig_pathways = len(sig_results)

                if n_sig_pathways > 0:
                    # Extract genes involved in significant pathways.
                    # g:Profiler returns 'intersections' as a list of gene IDs (strings) 
                    # OR sometimes a dict if evidence codes are mixed.
                    
                    col_name = 'intersections' if 'intersections' in sig_results.columns else 'intersection'
                    
                    if col_name in sig_results.columns:
                        for item in sig_results[col_name]:
                            if isinstance(item, list):
                                unique_sig_proteins.update(item)
                            elif isinstance(item, dict):
                                unique_sig_proteins.update(item.keys())
                            elif isinstance(item, str):
                                unique_sig_proteins.add(item)

            # Calculate Stats
            n_unique_sig_proteins = len(unique_sig_proteins)
            percentage = (n_unique_sig_proteins / n_total_proteins) * 100

            # Store results
            validation_results.append({
                'Cluster_ID': i,
                'Total_Proteins': n_total_proteins,
                'Significant_Pathways': n_sig_pathways,
                'Unique_Significant_Proteins': n_unique_sig_proteins,
                'Percentage_Involved': percentage
            })
            
            # Simple progress indicator
            if (i + 1) % 10 == 0:
                print(f"   ... Processed {i + 1}/{len(communities)} clusters")

        except Exception as e:
            print(f"⚠️ Error processing Cluster {i}: {e}")

    # --- 4. Analysis & Visualization ---
    if not validation_results:
        print("❌ No validation results generated.")
        return None

    df_validation = pd.DataFrame(validation_results)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Save to CSV
    csv_filename = f"validation_summary_{clust_algo_used}.csv"
    csv_path = os.path.join(output_dir, csv_filename)
    df_validation.to_csv(csv_path, index=False)
    print(f"\n✅ Validation results saved to: {csv_path}")

    # Calculate Global Averages
    avg_percentage = df_validation['Percentage_Involved'].mean()
    avg_pathways = df_validation['Significant_Pathways'].mean()

    print(f"\n--- 📊 Global Stats ({clust_algo_used}) ---")
    print(f"Average % Proteins in Significant Pathways: {avg_percentage:.2f}%")
    print(f"Average # Significant Pathways per Cluster: {avg_pathways:.2f}")

    # --- 5. Plotting ---
    plt.figure(figsize=(10, 6))
    
    # Histogram
    plt.hist(df_validation['Percentage_Involved'], bins=20, color='#69b3a2', edgecolor='black', alpha=0.8)
    
    # Add vertical line for mean
    plt.axvline(avg_percentage, color='red', linestyle='dashed', linewidth=1.5, label=f'Mean: {avg_percentage:.1f}%')
    
    plt.title(f'Cluster Validation: Protein Coverage by Pathways\n({clust_algo_used}, Sources: GO:BP, REAC)', fontsize=14)
    plt.xlabel('Percentage of Proteins in Cluster Involved in Sig. Pathways (%)', fontsize=12)
    plt.ylabel('Number of Clusters', fontsize=12)
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    
    # Save Plot
    plot_filename = f"validation_coverage_plot_{clust_algo_used}.png"
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path, dpi=300)
    print(f"✅ Plot saved to: {plot_path}")
    plt.close() # Close to prevent display issues in non-notebook environments

    return df_validation

# --- Example Usage (if run as script) ---
if __name__ == "__main__":
    # Define paths (Edit these as needed)
    PICKLE_FILE = "results/louvain_results.pkl" # Example path
    OUT_DIR = "results/validation"
    ALGO = "louvain"
    
    if os.path.exists(PICKLE_FILE):
        run_validation_analysis(PICKLE_FILE, OUT_DIR, ALGO)
    else:
        print(f"File {PICKLE_FILE} not found. Please import this module and run 'run_validation_analysis'.")