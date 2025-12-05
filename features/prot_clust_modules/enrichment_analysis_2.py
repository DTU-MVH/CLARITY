import warnings
# Filter specific warning if it persists, though less likely with gprofiler
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")

"""
Module: enrichment_analysis.py
Description: 
    Performs Over-Representation Analysis (ORA) on a specific cluster 
    using the gprofiler-official library (replacing gseapy).
    Sources used: GO:BP (Biological Process) and REAC (Reactome).
"""

import pandas as pd
import pickle
import os
import sys
from IPython.display import display

# Try importing gprofiler-official
try:
    from gprofiler import GProfiler
    GPROFILER_AVAILABLE = True
except ImportError:
    GPROFILER_AVAILABLE = False

def run_enrichment_analysis(target_protein, pickle_path, output_dir, clust_algo_used=None):
    """
    Runs ORA on the cluster containing the target protein.

    Parameters:
        target_protein (str): The name of the protein to analyze (e.g., "SNCAIP").
        pickle_path (str): Path to the .pkl file (Louvain or MCL data).
        output_dir (str): Directory where results will be saved.
        clust_algo_used (str, optional): Name of the clustering algorithm used.

    Returns:
        pd.DataFrame: The enrichment results dataframe (or None if failed).
    """

    if not GPROFILER_AVAILABLE:
        print(" Error: 'gprofiler-official' library not installed.")
        print("   Run: pip install gprofiler-official")
        return None

    print(f"\n---  Enrichment Analysis (g:Profiler): {target_protein} ---")

    # --- 1. Load Data ---
    if not os.path.exists(pickle_path):
        print(f" Error: Pickle file not found at {pickle_path}")
        return None

    try:
        with open(pickle_path, "rb") as f:
            data = pickle.load(f)
        
        # We need communities list and graph nodes for background
        communities = data['communities']
        G = data['graph']
        
        # Rebuild node map if needed
        if 'node2cluster' in data:
            node2cluster = data['node2cluster']
        else:
            node2cluster = {}
            for cid, comm in enumerate(communities):
                for node in comm:
                    node2cluster[node] = cid

    except Exception as e:
        print(f" Error loading pickle: {e}")
        return None

    # --- 2. Identify Target Cluster ---
    if target_protein not in node2cluster:
        print(f" Error: Protein '{target_protein}' not found in the graph.")
        return None

    cluster_id = node2cluster[target_protein]
    target_gene_list = list(communities[cluster_id])
    
    print(f"[INFO] Target found in Cluster #{cluster_id}")
    print(f"[INFO] Analyzing {len(target_gene_list)} genes against GO:BP and Reactome...")

    # --- 3. Define Background ---
    background_gene_list = list(G.nodes())
    print(f"[INFO] Background size: {len(background_gene_list)} genes")

    # --- 4. Run ORA with g:Profiler ---
    try:
        print("[INFO] Initializing g:Profiler...")
        gp = GProfiler(return_dataframe=True)
        
        print("[INFO] Querying g:Profiler API...")
        
        try:
            # Run the profile query
            # Sources: GO:BP (Biological Process), REAC (Reactome)
            results = gp.profile(
                organism='hsapiens',
                query=target_gene_list,
                background=background_gene_list,
                sources=['REAC', "kegg"], 
                user_threshold=0.05,
                significance_threshold_method='fdr', # False Discovery Rate (Benjamini-Hochberg)
                no_evidences=False  # Keep True if you don't need evidence codes, False speeds it up slightly
            )
        except KeyError as ke:
            print(f" KeyError in g:Profiler response: {ke}")
            print("[INFO] This often happens due to missing 'intersection' column in the API response.")
            print("[INFO] Retrying with 'no_evidences=True' parameter...")
            
            # Retry without evidences which might simplify response
            results = gp.profile(
                organism='hsapiens',
                query=target_gene_list,
                background=background_gene_list,
                sources=['GO:BP', 'REAC'], 
                user_threshold=0.05,
                significance_threshold_method='fdr',
                no_evidences=True  # Simplify response
            )

        if results.empty:
            print("[RESULT] No statistically significant pathways found.")
            return None
        
        # --- 5. Format Results to match previous structure ---
        # g:Profiler output cols: source, native, name, p_value, significant, description, term_size, query_size, intersection_size, effective_domain_size, intersections, evidences, etc.
        
        results['Overlap'] = results['intersection_size'].astype(str) + "/" + results['term_size'].astype(str)
        
        # g:Profiler 'p_value' is already corrected/adjusted based on the threshold method chosen above
        results.rename(columns={
            'name': 'Term',
            'p_value': 'Adjusted P-value',
            'source': 'Source'
        }, inplace=True)

        # Convert list of intersection genes to string (e.g., "GENE1;GENE2") if it's a list
        # g:Profiler returns 'intersections' (plural) as a dict with gene IDs as keys
        # Handle case where 'intersections' or 'intersection' column might not exist
        if 'intersections' in results.columns:
            results['Genes'] = results['intersections'].apply(lambda x: ';'.join(x.keys()) if isinstance(x, dict) else str(x))
        elif 'intersection' in results.columns:
            results['Genes'] = results['intersection'].apply(lambda x: ';'.join(x) if isinstance(x, list) else str(x))
        else:
            # Fallback: use empty string if intersection column doesn't exist
            results['Genes'] = ""

        # We don't have a raw P-value from g:Profiler output easily, so we duplicate Adj P-value or leave blank
        results['P-value'] = results['Adjusted P-value'] 

        # Select and reorder columns
        cols_to_keep = ['Source', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Genes']
        results_df = results[cols_to_keep].sort_values('Adjusted P-value')
        
        # --- 6. Save Results ---
        algo_tag = clust_algo_used if clust_algo_used else "unspecified"
        output_file_name = f"enrichment_cluster_{cluster_id}_{algo_tag}.csv"

        os.makedirs(output_dir, exist_ok=True)
        csv_path = os.path.join(output_dir, output_file_name)
        results_df.to_csv(csv_path, index=False)
        print(f"   Full results saved to: {csv_path}")

        # --- 7. Display ---
        print("\n--- Top Enriched Pathways (FDR < 0.05) ---")
        
        if not results_df.empty:
            # Format for display
            display_df = results_df[['Source', 'Term', 'Overlap', 'Adjusted P-value', 'Genes']].head(10).style.format({
                'Adjusted P-value': '{:.2e}'
            })
            display(display_df)
        
        return results_df

    except Exception as e:
        print(f"Analysis failed: {e}")
        # Print full traceback for debugging
        import traceback
        print("\n[DEBUG] Full traceback:")
        traceback.print_exc()
        # Helpful debugging for g:profiler specific errors
        if "404" in str(e) or "Connection" in str(e):
            print("   (Check your internet connection)")
        return None