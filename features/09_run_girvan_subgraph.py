import pandas as pd
import networkx as nx
import time
import os
import pickle
import sys
import random
from networkx.algorithms.community import girvan_newman as nx_girvan_newman
from networkx.algorithms.community.quality import modularity as calculate_modularity

# --- CRITICAL DEPENDENCY CHECK ---
# Check for required modules before starting
try:
    import pandas as pd
    import networkx as nx
    from networkx.algorithms.community import girvan_newman as nx_girvan_newman
except ImportError as e:
    print(f"\nDEPENDENCY ERROR: Critical module is missing. Please ensure your environment has all packages installed.")
    print(f"Missing module: {e}")
    sys.exit(1)
# -----------------------------------


# ==========================================
# CONFIGURATION
# ==========================================
INPUT_CSV_PATH = "data/cleaned_data.csv"
OUTPUT_DIR = "data/output"
# MODIFIED: Output file now indicates SUBGRAPH and FIRST PARTITION run
OUTPUT_PICKLE_NAME = "girvan_newman_clust_FIRST_PARTITION.pkl" 
# This variable is no longer used in the execution logic but kept for consistency
MAX_PARTITIONS_TO_CHECK = 50 
RANDOM_SEED = 42
SUBGRAPH_RADIUS = 1 # Defines the size of the network chunk

# --- PATH OVERRIDE ---
FULL_PATH_OVERRIDE = ""
# ---------------------

# Ensure random operations are reproducible
random.seed(RANDOM_SEED)

def make_graph(df):
    """Builds the weighted NetworkX graph from a DataFrame."""
    G = nx.Graph()
    for p1, p2, score in zip(df['protein1'], 
                             df['protein2'], 
                             df['combined_score']):
        G.add_edge(p1, p2, combined_score=score)
    return G

def get_first_girvan_newman_partition(G, weight_attr=None):
    """
    MODIFIED: Executes GN and returns only the FIRST partition generated, 
    ignoring modularity maximization (for documentation purposes).
    """
    print("[INFO] Starting Girvan-Newman (First Partition Only)...")
    
    t0 = time.time()
    
    # GN generator initialized
    communities_generator = nx_girvan_newman(G)

    try:
        # Get the partition after the first edge removal (first item from iterator)
        first_community = next(communities_generator)
    except StopIteration:
        # Handle case where graph is too small/empty
        first_community = [list(G.nodes())]
        
    # Calculate modularity on the extracted partition
    modularity_score = calculate_modularity(G, first_community, weight=weight_attr)

    runtime = time.time() - t0

    return first_community, modularity_score, runtime

def main():
    # Set up output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    G_full = None
    
    # --- Data Path Resolution ---
    if FULL_PATH_OVERRIDE:
        resolved_input_path = FULL_PATH_OVERRIDE
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        resolved_input_path = os.path.join(script_dir, '..', INPUT_CSV_PATH)
    # ---------------------------------------

    # --- Attempt to Load Data ---
    try:
        print(f"[INFO] Loading data from {resolved_input_path}...")
        df_interactions = pd.read_csv(resolved_input_path)
        G_full = make_graph(df_interactions)
        print(f"[STATS] Full Graph Nodes: {G_full.number_of_nodes()}, Edges: {G_full.number_of_edges()}")

    except FileNotFoundError:
        print(f"[FATAL ERROR] Input file not found at {resolved_input_path}. Cannot proceed.")
        sys.stdout.flush()
        sys.exit(1)

    # --- 2. Define Run Graph (The Chunk) ---
    all_nodes = list(G_full.nodes())
    random_start_node = random.choice(all_nodes) 
    G_run = nx.ego_graph(G_full, random_start_node, radius=SUBGRAPH_RADIUS)
    
    print(f"[STATS] Running GN on Subgraph (Ego of {random_start_node})")
    print(f"[STATS] Graph Nodes: {G_run.number_of_nodes()}, Edges: {G_run.number_of_edges()}")
        
    # --- 3. Run Girvan-Newman Clustering (FIRST PARTITION) ---
    print(f"[WARNING] Results represent the partition after the *first* edge removal (non-optimal).")
    
    communities, mod_score, runtime = get_first_girvan_newman_partition(
        G_run, 
        weight_attr=None 
    )
        
    # --- 4. Calculate and Print Final Outputs (Easily Accessible) ---
    num_communities = len(communities)
    cluster_sizes = sorted([len(c) for c in communities], reverse=True)
    
    print("\n--- Girvan-Newman Results Summary (FIRST PARTITION) ---")
    print(f"[RESULT] Run time on chunk: {runtime:.2f}s")
    print(f"[RESULT] Detected {num_communities} communities")
    print(f"[RESULT] Modularity Score (Unweighted): {mod_score:.4f}")
    print(f"[RESULT] Cluster Sizes (Top 5): {cluster_sizes[:5]}...")

    # --- 5. Save to Pickle ---
    data_to_save = {
        "graph": G_run, # Saving the SUBGRAPH that was actually clustered
        "communities": communities,
        "modularity": mod_score,
        "parameters": {
            "algorithm": "Girvan-Newman (First Partition Only)",
            "max_partitions_checked": 1,
            "seed_used": RANDOM_SEED,
            "runtime_s": runtime,
            "subgraph_root": random_start_node,
            "subgraph_radius": SUBGRAPH_RADIUS
        }
    }

    output_path = os.path.join(OUTPUT_DIR, OUTPUT_PICKLE_NAME)
    print(f"[INFO] Saving results to {output_path}...")
    with open(output_path, "wb") as f:
        pickle.dump(data_to_save, f)
    
    print(f"[DONE] Script completed successfully.")
    sys.stdout.flush() 
    
if __name__ == "__main__":
    main()