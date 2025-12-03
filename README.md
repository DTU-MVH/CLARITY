# CLARITY: Protein-Protein Interaction Clustering Analysis Pipeline

## Overview

**CLARITY** is a comprehensive bioinformatics tool designed to analyze protein-protein interaction (PPI) networks using multiple graph clustering algorithms. The pipeline implements two complementary clustering approaches—**Louvain** and **Markov Clustering (MCL)**—to identify functional modules and protein complexes from STRING database interaction data. This dual-method approach enables researchers to validate findings across different algorithmic paradigms and gain deeper biological insights into protein network organization.

### Key Features

- **Data Import & Preprocessing**: Automatic download and cleaning of protein interaction data from STRING database
- **Multiple Clustering Algorithms**: Louvain (hierarchical community detection) and MCL (flow-based clustering)
- **Network Analysis**: Comprehensive cluster statistics, density measurements, and topological analysis
- **Bottleneck Detection**: Identifies critical hub proteins using betweenness centrality, degree centrality, and clustering coefficients
- **Enrichment Analysis**: Gene Ontology (GO) and Reactome pathway enrichment using g:Profiler
- **Reproducibility Assessment**: Validates clustering consistency across different random seeds
- **Interactive Visualization**: Interactive HTML visualizations using PyVis for network exploration
- **Method Validation**: Cross-validation between clustering algorithms to assess robustness

---

## 📊 Clustering Approaches Explained

### A. Louvain Algorithm
**Type**: Hierarchical community detection based on modularity optimization

**How it works**: The Louvain algorithm is a greedy heuristic that iteratively groups nodes to maximize modularity—the difference between observed edge density within communities and expected density in a random null model.

**Biological Interpretation**: View Louvain clusters as broad "functional neighborhoods" in the PPI network. For example, if the network were a city, these would represent distinct districts (e.g., the "Metabolism District", "Signal Transduction District"). This approach excels at identifying:
- Large-scale biological processes
- Global network topology
- Hierarchical functional organization

**Strengths**:
- Fast and scalable
- Excellent for large networks
- Resolves broad, functionally coherent groups

**Limitations**:
- Suffers from **resolution limit**: small but distinct protein complexes may be merged into larger super-modules if they share even a few connections
- May miss tightly-knit small complexes

---

### B. Markov Clustering (MCL)
**Type**: Flow-based algorithm simulating random walks

**How it works**: MCL simulates random walks through the network, with the intuition that a random walker is more likely to remain inside a dense cluster than to exit it. The algorithm alternates between:
- **Expansion**: Spreading flow through the network
- **Inflation**: Pruning weak connections to reinforce strong clusters

**Biological Interpretation**: View MCL clusters as tight-knit "families" or "molecular machines"—identifying specific protein complexes rather than broad functional categories. This approach is excellent for:
- Detecting specific molecular complexes (e.g., the Ribosome)
- Finding tightly coordinated protein interactions
- Handling noise in interaction data

**Strengths**:
- Often considered the gold standard for defining protein complexes
- Handles noise well
- Identifies tight, functionally coherent modules

**Limitations**:
- Sensitive to the **inflation parameter**:
  - Too high: fragments valid pathways into small pieces
  - Too low: lumps distinct complexes together
- Slower on very large networks

---

### Comparison: When to Use Each Method

| Aspect | Louvain | MCL |
|--------|---------|-----|
| **Best for** | Broad functional pathways | Specific protein complexes |
| **Scale** | Large networks (10k+ nodes) | Medium networks (1k-15k nodes) |
| **Resolution** | Low (broad communities) | High (tight complexes) |
| **Speed** | Fast | Slower |
| **Interpretation** | Biological processes | Molecular machines |

---

## 📁 Project File Structure

```
CLARITY/
├── main.ipynb                           # Main Jupyter notebook orchestrating the pipeline
├── requirements.txt                     # Python package dependencies
├── 06_data_import.py                   # Script to download and organize data from Google Drive
├── .gitignore                           # Git ignore file
│
├── _raw/                               # Raw data directory (downloaded from STRING)
│                                        # [Contents generated at runtime - not tracked in git]
│
├── data/                               # Processed data directory
│                                        # [Contents generated at runtime - not tracked in git]
│
└── features/                           # Core analysis modules
    ├── preprocessing.py               # Data cleaning pipeline
    ├── 01_run_louvain.py             # Louvain clustering execution
    ├── 02_analyze_clusters.py        # Cluster statistics & characterization
    ├── 03_reproducibility.py         # Seed-based reproducibility validation
    ├── 04_target_bottleneck.py       # Hub protein identification
    ├── 05_markov_clust_highperformance_object_out.py  # MCL clustering
    ├── analyze_results_louvain.py    # Post-clustering analysis
    ├── plotting/                      # Visualization modules
    └── prot_clust_modules/            # Reusable analysis modules
```

---

## 🚀 Getting Started: Installation & Setup

### Prerequisites
- **Python**: 3.8 or higher
- **Operating System**: macOS, Linux, or Windows
- **Internet Connection**: Required for initial data download

### Step 1: Clone or Download the Repository

```bash
cd /path/to/CLARITY
```

### Step 2: Create a Virtual Environment

We recommend using a virtual environment to avoid package conflicts:

**Option A: Using VS Code (Recommended)**
1. Press `Cmd + Shift + P`
2. Type `Python: Create Environment`
3. Select `Venv` and choose `Quick Create venv`
4. VS Code automatically installs all dependencies from `requirements.txt`

**Option B: Using Terminal**
```bash
# Create virtual environment
python3 -m venv venv

# Activate it (macOS/Linux)
source venv/bin/activate

# Activate it (Windows)
venv\Scripts\activate
```

### Step 3: Install Dependencies

```bash
pip install -r requirements.txt
```

**Required packages**:
- `pandas` – Data manipulation
- `networkx` – Graph algorithms
- `python-louvain` – Louvain clustering
- `markov-clustering` – MCL clustering
- `scikit-learn` – Machine learning utilities
- `scipy` – Scientific computing
- `matplotlib`, `seaborn` – Visualization
- `pyvis` – Interactive network visualization
- `gprofiler-official` – Functional enrichment
- `gdown` – Google Drive downloads
- `ipykernel` – Jupyter kernel support

---

## 📖 How to Use: Complete Guide

### Overview: Main Workflow

The pipeline runs sequentially through these stages:

```
Data Import → Preprocessing → Clustering → Analysis → Enrichment → Validation
```

All workflows are orchestrated through `main.ipynb`.

---

### Running the Complete Pipeline

#### 1. **Start Jupyter Notebook**

```bash
jupyter notebook main.ipynb
```

The notebook opens in your browser and guides you through each step.

---

#### 2. **Environment Setup** (First-time only)

The first cell sets up a helper function to execute standalone Python scripts:

```python
def run_script(script_path):
    """Executes a Python script in a clean subprocess."""
    # Implementation provided in notebook
```

No action needed—this is automatic.

---

#### 3. **Data Import**

**Run Cell**: "3. Importing data"

```python
run_script("06_data_import.py")
```

This downloads and organizes data:
- Downloads protein metadata and interaction files from Google Drive
- Extracts raw data to `_raw/`
- Pre-cleaned data is saved to `data/cleaned_data.csv`
- Pre-computed clustering results saved to `data/output/`

**Output**:
- `_raw/9606.protein.info.v12.0.txt` – Protein metadata
- `_raw/9606.protein.links.v12.0.min400.onlyAB.tsv` – Interaction edges
- `data/cleaned_data.csv` – Ready-to-use interaction data

**Time**: ~2-5 minutes depending on connection speed

---

#### 4. **Data Preprocessing** (Optional)

**Run Cell**: "4. Data Preprocessing" (if you modify parameters)

```python
run_script("features/preprocessing.py")
```

This script:
- Loads raw STRING data files
- Converts protein IDs to human-readable names (gene symbols)
- Removes nodes with only 1 edge (singletons)
- Constructs weighted graph where weights = interaction confidence scores
- Saves cleaned data to `data/protein_link_cleaned.tsv`

**Skip this if** `data/cleaned_data.csv` already exists—preprocessing is already done.

**Configuration** (in `preprocessing.py`):
```python
info_file = '_raw/9606.protein.info.v12.0.txt'
interaction_file = '_raw/9606.protein.links.v12.0.min400.onlyAB.tsv'
output_file = 'data/protein_link_cleaned.tsv'
```

---

#### 5. **Clustering: Louvain Algorithm**

**Run Cell**: "5.1 Louvain clustering analysis"

```python
run_script("features/01_run_louvain.py")
```

This executes Louvain community detection:

**What it does**:
1. Loads interaction data from `data/cleaned_data.csv`
2. Constructs weighted NetworkX graph
3. Runs Louvain algorithm with specified resolution
4. Calculates modularity score (quality metric)
5. Saves results to pickle file

**Configuration** (in `features/01_run_louvain.py`):
```python
INPUT_CSV_PATH = "data/cleaned_data.csv"
OUTPUT_DIR = "data/output"
OUTPUT_PICKLE_NAME = "louvain_clust_julle.pkl"
LOUVAIN_RESOLUTION = 1.0        # ← Adjust for finer/coarser clustering
RANDOM_SEED = 42                 # ← Set for reproducibility
```

**Parameters explained**:
- `LOUVAIN_RESOLUTION`: Controls cluster granularity
  - Lower values (0.5) → Fewer, larger clusters
  - Higher values (1.5) → More, smaller clusters
  - Default (1.0) → Balanced

**Output**:
- `data/output/louvain_clust_julle.pkl` – Graph + communities + modularity

**Time**: 30 seconds to 5 minutes (depending on network size)

---

#### 6. **Clustering: Markov Clustering (MCL)**

**Run Cell**: "5.2 MCL clustering analysis" (if available in your version)

```python
run_script("features/05_markov_clust_highperformance_object_out.py")
```

This executes MCL clustering:

**What it does**:
1. Loads same interaction data
2. Constructs graph with MCL-optimized parameters
3. Runs MCL algorithm with inflation parameter
4. Saves results to pickle file

**Configuration** (typical settings):
```python
inflation_parameter = 1.5  # ← Adjust cluster tightness
# Lower (1.2) → looser clusters
# Higher (2.0) → tighter clusters
```

**Output**:
- `data/output/mcl_data.pkl` – Similar structure to Louvain results

---

#### 7. **Analyze Clusters**

**Run Cell**: "6. Cluster Analysis"

```python
run_script("features/02_analyze_clusters.py")
```

Generates detailed cluster statistics:

**What it does**:
1. Loads clustering results from pickle
2. Calculates per-cluster metrics:
   - **Size**: Number of nodes in cluster
   - **Density**: Ratio of actual edges to possible edges
   - **Clustering Coefficient**: Local interconnectedness
   - **Modularity**: Quality of partition

3. Saves comprehensive report to CSV

**Output**:
- `data/output/cluster_statistics.csv`

**Example output table**:
```
Cluster ID | Size | Edges | Density | Avg Clustering
0          | 1205 | 8903  | 0.0122  | 0.342
1          | 432  | 2145  | 0.0231  | 0.456
2          | 156  | 1032  | 0.0848  | 0.623
...
```

---

#### 8. **Set Target Protein**

**Run Cell**: "5. Clustering Algorithms" section

```python
TARGET_PROTEIN = "AGTRAP"  # ← Change to your protein of interest
```

All downstream analyses (bottleneck, enrichment) use this target.

Common proteins to analyze:
- `SNCAIP` – α-Synuclein-associated protein (neuroscience)
- `CHMP2B` – Charged multivesicular body protein (membrane biology)
- `TP53` – Tumor suppressor (cancer biology)
- `MYC` – Proto-oncogene (cell proliferation)

---

#### 9. **Bottleneck Analysis**

**Run Cell**: "7. Target Bottleneck Analysis"

```python
run_script("features/04_target_bottleneck.py")
```

Identifies critical hub proteins in the target cluster:

**What it does**:
1. Locates the cluster containing `TARGET_PROTEIN`
2. Calculates three centrality metrics:
   - **Betweenness Centrality (BC)**: How often a protein bridges other pairs
   - **Degree Centrality**: Number of direct interactions
   - **Clustering Coefficient**: Local network density (want LOW values)

3. Combines rankings to identify "bottleneck" proteins (high BC + high degree + low clustering)
4. Generates visualizations

**How bottlenecks work**:
- High BC + High Degree → Network hub
- Low Clustering Coeff → Not in tight local cluster
- **Result**: Protein that bridges different functional groups

**Output files**:
- `data/output/bottlenecks_cluster_*.csv` – Ranked proteins
- `data/output/images/cluster_*_interactive.html` – Interactive visualization
- `data/output/target_cluster_nodes.txt` – Cluster membership

**Configuration** (in `features/04_target_bottleneck.py`):
```python
TARGET_PROTEIN = "SNCAIP"    # ← Set in main notebook
TOP_N_RANKING = 40            # ← How many top candidates to report
```

---

#### 10. **Enrichment Analysis**

**Run Cell**: "8. Functional Enrichment"

```python
run_script("features/prot_clust_modules/enrichment_analysis_2.py")
```

Performs Gene Ontology and pathway enrichment:

**What it does**:
1. Takes proteins in target cluster as "query" genes
2. Uses g:Profiler to find significantly enriched:
   - GO Biological Processes (GO:BP)
   - GO Molecular Functions (GO:MF)
   - Reactome Pathways
   
3. Calculates p-values and adjusts for multiple testing
4. Saves enrichment report

**Interpretation**:
- High-scoring terms represent biological functions overrepresented in the cluster
- Example: If analyzing AGTRAP cluster, you might find enrichment in "G protein-coupled receptor signaling"

**Output**:
- `data/output/enrichment_cluster_*.csv` – GO terms + p-values + gene counts

**Prerequisites**:
- Requires internet connection (queries live g:Profiler API)
- Requires `gprofiler-official` package

---

#### 11. **Reproducibility Check**

**Run Cell**: "9. Method Reproducibility"

```python
run_script("features/03_reproducibility.py")
```

Tests how consistent the clustering is across random seeds:

**What it does**:
1. Runs Louvain clustering 10 times with different random seeds
2. Compares results using:
   - **Adjusted Rand Index (ARI)**: Measures partition similarity (-1 to 1, higher = more similar)
   - **Normalized Mutual Information (NMI)**: Information-theoretic similarity (0 to 1)

3. Plots cluster size distributions to visualize stability

**Example output**:
```
Reproducibility Results (10 runs):
Mean ARI: 0.9523 ± 0.0145
Mean NMI: 0.9834 ± 0.0089
```

High scores (>0.95) indicate robust, reproducible results.

**Output**:
- Console output with ARI/NMI statistics
- `data/output/plots/reproducibility_*.png` – Distribution comparison

---

#### 12. **Method Validation**

**Run Cell**: "10. Cross-Algorithm Validation"

```python
run_script("features/prot_clust_modules/method_validation_enrichment.py")
```

Compares Louvain and MCL clustering results:

**What it does**:
1. Loads both clustering results
2. Compares community assignments using:
   - ARI (Adjusted Rand Index)
   - NMI (Normalized Mutual Information)
   - Overlap statistics

3. Tests enrichment results consistency
4. Generates validation report

**Interpretation**:
- High agreement (ARI > 0.7) → Methods find similar structure
- Low agreement → Methods detect different aspects of network
- Both results are valid—just answering different biological questions

**Output**:
- `data/output/method_validation_results.csv` – Comparison metrics

---

## ⚙️ Configuration Reference

### Key Configuration Parameters

**File**: `features/01_run_louvain.py`
```python
LOUVAIN_RESOLUTION = 1.0        # [0.5-2.0] Granularity of clustering
RANDOM_SEED = 42                # [int] For reproducibility
```

**File**: `features/04_target_bottleneck.py`
```python
TARGET_PROTEIN = "SNCAIP"       # [str] Protein to analyze
TOP_N_RANKING = 40              # [int] Top bottleneck candidates
```

**File**: `features/06_data_import.py`
```python
DRIVE_FOLDER_URL = "..."        # Google Drive folder containing data
files_to_raw = [...]            # Files going to _raw/
files_to_data = [...]           # Files going to data/
files_to_output = [...]         # Files going to data/output/
```

**File**: `requirements.txt`
All Python dependencies—modify to add new packages or change versions.

---

## 📊 Expected Output Examples

### 1. Cluster Statistics
```
Cluster ID | Size | Edges | Density | Avg Clustering | Description
0          | 1205 | 8903  | 0.0122  | 0.342          | Metabolism-related
1          | 432  | 2145  | 0.0231  | 0.456          | Signal transduction
23         | 156  | 1032  | 0.0848  | 0.623          | Target cluster
```

### 2. Bottleneck Analysis
```
Protein  | BC_Rank | Deg_Rank | Clust_Rank | Total_Score | Interpretation
AGTRAP   | 2       | 1        | 5          | 8           | Strong hub
GPCR3    | 1       | 3        | 2          | 6           | Excellent bottleneck
PKA      | 4       | 2        | 1          | 7           | Good hub
```

### 3. Enrichment Results
```
Term                              | p-value | Gene_Count | Database
G protein-coupled receptor signal | 1.2e-15 | 24         | GO:BP
MAPK cascade                      | 3.4e-12 | 18         | GO:BP
Phosphorylation                   | 2.1e-08 | 42         | Reactome
```

---

## 🔧 Troubleshooting

### Issue: "Module not found" errors
**Solution**: 
```bash
pip install -r requirements.txt
```

### Issue: Pickle file not found
**Solution**: 
- Run data import script first: `run_script("06_data_import.py")`
- Or run clustering: `run_script("features/01_run_louvain.py")`

### Issue: Out of memory on large networks
**Solution**:
- Use MCL with approximate centrality: Set `use_approximation=True` in bottleneck analysis
- Reduce network size by filtering low-confidence interactions

### Issue: g:Profiler enrichment fails
**Solution**:
```bash
pip install --upgrade gprofiler-official
```

---

## 📚 Understanding the Output Files

| File | Purpose | Format |
|------|---------|--------|
| `louvain_clust_julle.pkl` | Louvain clustering results | Python pickle (graph + communities) |
| `cluster_statistics.csv` | Per-cluster metrics | CSV table |
| `bottlenecks_cluster_*.csv` | Hub protein rankings | CSV table |
| `enrichment_cluster_*.csv` | GO/pathway enrichment | CSV table |
| `method_validation_results.csv` | Algorithm comparison | CSV table |
| `*.html` | Interactive visualizations | HTML (open in browser) |

---

## 🎯 Example Analysis Workflow

**Scenario**: Analyze a novel Alzheimer's-related protein

```python
# Step 1: Set target
TARGET_PROTEIN = "PSEN1"  # Presenilin-1 (Alzheimer's-associated)

# Step 2: Run through pipeline
run_script("features/01_run_louvain.py")       # Cluster
run_script("features/02_analyze_clusters.py")  # Get stats
run_script("features/04_target_bottleneck.py") # Find hubs
run_script("features/prot_clust_modules/enrichment_analysis_2.py")  # Get function
run_script("features/03_reproducibility.py")   # Validate

# Step 3: Examine outputs
# - bottlenecks_cluster_*.csv → Which proteins are crucial?
# - enrichment_cluster_*.csv → What functions are enriched?
# - cluster_23_interactive.html → Visualize the network
```

---

## 📝 Citation & References

**Data Source**: STRING Database v12.0
- Human proteins (Homo sapiens, NCBI taxonomy ID: 9606)
- Minimum combined score: 400

**Algorithms**:
- **Louvain**: Blondel et al., *J. Stat. Mech.* (2008)
- **MCL**: Enright et al., *Nucl. Acids Res.* (2002)

**Enrichment**: g:Profiler (Raudvere et al., *Nucl. Acids Res.* 2019)

**Validation Metrics**: 
- ARI: Hubert & Arabie (1985)
- NMI: Danon et al., *J. Stat. Mech.* (2005)

---

## 📄 License & Contact

This project is part of DTU's Computational Tools for Data Science course.

**Repository**: github.com/DTU-MVH/CLARITY

**Questions?** See the inline comments in `main.ipynb` or individual feature scripts.

---

## 🚀 Advanced: Customization Guide

### Modify Clustering Resolution

Edit `features/01_run_louvain.py`:
```python
LOUVAIN_RESOLUTION = 0.5  # Coarser clustering
# or
LOUVAIN_RESOLUTION = 2.0  # Finer clustering
```

Lower resolution merges communities, higher resolution splits them.

### Analyze Multiple Target Proteins

Create a loop in the notebook:
```python
targets = ["AGTRAP", "CHMP2B", "TP53"]
for target in targets:
    TARGET_PROTEIN = target
    run_script("features/04_target_bottleneck.py")
```

### Filter Interactions by Confidence

In `preprocessing.py`, add threshold:
```python
df_interactions = df_interactions[df_interactions['combined_score'] >= 700]
```

---

## 📞 Support

For issues or questions:
1. Check the troubleshooting section above
2. Review inline comments in the relevant Python script
3. Run with `print()` debugging enabled
4. Contact project maintainer

---

**Last Updated**: December 2025
**Version**: 2.0
**Status**: Production Ready