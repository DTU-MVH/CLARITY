"""
Method validation using the project's enrichment (g:Profiler) pipeline.

Loads a clustering pickle (saved by `01_run_louvain.py`) and runs ORA for
each cluster using `gprofiler-official` (same approach as
`enrichment_analysis_2.py`). Produces a CSV with per-cluster metrics and a
histogram showing the percentage of proteins in each cluster involved in
significant pathways.

Usage:
    python features/prot_clust_modules/method_validation_enrichment.py

Adjust `PICKLE_PATH` and `OUTPUT_DIR` variables below as needed.
"""

import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from gprofiler import GProfiler


def run_validation(pickle_path="data/output/louvain_clust_julle.pkl", output_dir="data/output",
                   pval_thresh=0.05, sources=["GO:BP", "REAC"]):
    if not os.path.exists(pickle_path):
        raise FileNotFoundError(f"Pickle file not found: {pickle_path}")

    with open(pickle_path, "rb") as f:
        data = pickle.load(f)

    communities = data.get("communities")
    G = data.get("graph")

    if communities is None or G is None:
        raise ValueError("Pickle does not contain 'communities' and 'graph'.")

    background_gene_list = list(G.nodes())

    print(f"Background universe size: {len(background_gene_list)} genes.")

    gp = GProfiler()

    validation_results = []

    for cid, comm in enumerate(communities):
        cluster_gene_list = list(comm)
        n_total_proteins = len(cluster_gene_list)

        if n_total_proteins == 0:
            continue

        try:
            results = gp.profile(
                organism="hsapiens",
                query=cluster_gene_list,
                background=background_gene_list,
                sources=sources,
                user_threshold=pval_thresh,
                significance_threshold_method='fdr',
                no_evidences=True,
                as_dataframe=True
            )

            # If results empty or no significant rows, handle gracefully
            if results is None or (hasattr(results, 'empty') and results.empty):
                n_sig_pathways = 0
                unique_sig_proteins = set()
            else:
                # Debug: show what columns we received
                try:
                    print(f"[DEBUG] Cluster {cid} -- result columns: {list(results.columns)}")
                    # print a tiny sample
                    try:
                        print(results.head(3))
                    except Exception:
                        pass
                except Exception:
                    pass

                # Determine p-value column (robust to different wrapper versions)
                pval_col = None
                for candidate in ['p_value', 'Adjusted P-value', 'pvalue', 'p_val', 'pval']:
                    if candidate in results.columns:
                        pval_col = candidate
                        break

                if pval_col is None:
                    # No p-value column found; treat as no significant results
                    n_sig_pathways = 0
                    unique_sig_proteins = set()
                else:
                    # Extract genes participating in each enriched term
                    if 'intersections' in results.columns:
                        results['Genes'] = results['intersections'].apply(
                            lambda x: ';'.join(x.keys()) if isinstance(x, dict) else str(x)
                        )
                    elif 'intersection' in results.columns:
                        results['Genes'] = results['intersection'].apply(
                            lambda x: ';'.join(x) if isinstance(x, list) else str(x)
                        )
                    else:
                        # fallback: try known gene columns
                        results['Genes'] = results.get('intersections', results.get('intersection', ''))

                    # Use the determined p-value column for filtering
                    sig_results = results[results[pval_col] < pval_thresh]
                    n_sig_pathways = len(sig_results)

                    unique_sig_proteins = set()
                    if n_sig_pathways > 0:
                        all_genes_string = ";".join(sig_results['Genes'].astype(str).tolist())
                        split_genes = [g for g in all_genes_string.split(';') if g]
                        unique_sig_proteins = set(split_genes)

        except Exception as e:
            print(f"Error processing cluster {cid}: {e}")
            n_sig_pathways = 0
            unique_sig_proteins = set()

        n_unique_sig_proteins = len(unique_sig_proteins)
        percentage = (n_unique_sig_proteins / n_total_proteins) * 100

        validation_results.append({
            'Cluster_ID': cid,
            'Total_Proteins': n_total_proteins,
            'Significant_Pathways': n_sig_pathways,
            'Unique_Significant_Proteins': n_unique_sig_proteins,
            'Percentage_Involved': percentage
        })

    # Save and report
    df_validation = pd.DataFrame(validation_results)
    os.makedirs(output_dir, exist_ok=True)
    output_filename = os.path.join(output_dir, "method_validation_results.csv")
    df_validation.to_csv(output_filename, index=False)

    print("\n--- Validation Results Summary (first 10 rows) ---")
    print(df_validation.head(10))
    print(f"\nResults saved to {output_filename}")

    avg_percentage = df_validation['Percentage_Involved'].mean()
    avg_pathways = df_validation['Significant_Pathways'].mean()

    print(f"\nAverage Percentage of Proteins involved in Significant Pathways: {avg_percentage:.2f}%")
    print(f"Average number of Significant Pathways per Cluster: {avg_pathways:.2f}")

    # Plot histogram
    plt.figure(figsize=(10, 6))
    plt.hist(df_validation['Percentage_Involved'].dropna(), bins=20, color='skyblue', edgecolor='black')
    plt.title('Method Validation: Protein Coverage by Significant Pathways')
    plt.xlabel('Percentage of Proteins in Cluster Involved in Sig. Pathways (%)')
    plt.ylabel('Number of Clusters')
    plt.axvline(avg_percentage, color='red', linestyle='dashed', linewidth=1, label=f'Mean: {avg_percentage:.1f}%')
    plt.legend()
    plt.grid(axis='y', alpha=0.5)
    plt.show()

    return df_validation


__all__ = ["run_validation"]

# Usage from other modules:
# from features.prot_clust_modules.method_validation_enrichment import run_validation
# df = run_validation(pickle_path="data/output/louvain_clust_julle.pkl", output_dir="data/output", pval_thresh=0.05)
