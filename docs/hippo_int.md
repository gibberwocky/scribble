# Hippocampal foetal organoid snRNA-seq

Here we take 10X snRNA-seq data of developing hippocampal tissue and of an organoid, and use `scribble` to pre-process both datasets, integrate and cluster for downstream annotation.


## Generate velocyto loom files

Velocyto requires an older version of Python (3.8) to install cleanly. So we create an env to support this.

```bash
mamba create --name velocyto python=3.8
mamba activate velocyto
mamba install numpy scipy cython numba matplotlib scikit-learn h5py click pysam scvelo
pip install velocyto
```

Next, we run velocyto to generate spliced and unspliced count matrices from the 10X BAM output. The `velocyto.sh` script reference below is provided in the `scribble/scripts` floder. In the below example, the HPC431 sample is processed, the `BAM` and `BARCODES` file paths should be updated to process the K2HO120 sample.

```bash
mamba activate velocyto
PROJECT=/uoa/home/s14dw4/sharedscratch/KangLab/hippocampus
BAM=${PROJECT}/cellranger/HPC431/outs/possorted_genome_bam.bam
BARCODES=${PROJECT}/cellranger/HPC431/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
GENES=/uoa/home/s14dw4/sharedscratch/software/cellranger-10.0.0/refdata-gex-GRCh38-2024-A/genes/genes.gtf
RPT=/uoa/home/s14dw4/sharedscratch/software/cellranger-10.0.0/refdata-gex-GRCh38-2024-A/GRCh38.rpt.gtf
mkdir -p ${PROJECT}/velocyto

sbatch --partition uoa-compute \
    -o ${PROJECT}/logs/velo.%j.out \
    -e ${PROJECT}/logs/velo.%j.err \
    /uoa/home/s14dw4/sharedscratch/scripts/velocyto.sh \
        --barcodes ${BARCODES} \
        --out ${PROJECT}/velocyto \
        --repeats ${RPT} \
        --bam ${BAM} \
        --genes ${GENES}
```

## Pre-process data with scribble

<br>

### Import data

I have `scribble` installed in an env which uses `python=3.12`, and in which `torch` was installed whilst connected to a GPU node. The import tool expects an `xlsx` file containing the worksheet `meta`, which includes the column `sample` whose values match the sample names provided to `--samples`, which should also be directories in the cellranger and velocyto directories. After importing the data and appending the metadata, a `combined.h5ad` file is written to the `scribble/adata` directory in `--project_dir`.

```bash
mamba activate scribble

SCRATCH=/uoa/home/s14dw4/sharedscratch/KangLab/hippocampus

# Scribble: import data
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 1 --mem 24G --time=4:00:00 \
    -o ${SCRATCH}/logs/sc_import.%j.out -e ${SCRATCH}/logs/sc_import.%j.err \
    scribble import \
    --project_dir ${SCRATCH} \
    --cellranger_dir ${SCRATCH}/cellranger \
    --velocyto_dir ${SCRATCH}/velocyto \
    --samples HPC431 K2HO120 \
    --metadata_file ${SCRATCH}/samples.xlsx
```

<br>
<table>
  <tr>
    <td><img src="../img/hippo_int/HPC431_qc_panel.png" alt="HPC431 QC panel"></td>
  </tr>
  <tr>
    <td><img src="../img/hippo_int/K2HO120_qc_panel.png" alt="K2HO120 QC panel"></td>
  </tr>
</table>

<br>

### Identify MT outliers

After generating the `combined.h5ad` file, the next step is to annotate it with mitchondrial (MT) metrics. The `--nmads` parameters sets the number of median absolute deviations as a threshold for which to label cells as MT outliers. This outputs a new `h5ad` file in `scribble/adata`.

```bash
# Scribble: MT QC
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 1 --mem 4G --time=2:00:00 \
    -o ${SCRATCH}/logs/sc_mt.%j.out -e ${SCRATCH}/logs/sc_mt.%j.err \
    scribble mt \
    --project_dir ${SCRATCH} \
    --input ${SCRATCH}/scribble/adata/combined.h5ad \
    --nmads 8
```

<br>
<table>
  <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8.png" alt="MT outliers"></td>
  </tr>
</table>

<br>

### Identify doublets

Net we need to label doublets. Here we run `dbl` in `hybrid` mode to apply both quantile and `scrublet` methods, with an expected doublet fraction of `0.07` and minimum cell count of 200 for a sample to be processed with scrublet. This outputs a new `h5ad` file in `scribble/adata`.

```bash
# Scribble: Doublet QC
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 1 --mem 32G --time=2:00:00 \
    -o ${SCRATCH}/logs/sc_dbl.%j.out -e ${SCRATCH}/logs/sc_dbl.%j.err \
    scribble dbl \
    --project_dir ${SCRATCH} \
    --input ${SCRATCH}/scribble/adata/combined_mtqc_nMADs-8.h5ad \
    --expected 0.07 \
    --mode hybrid \
    --min_cells 200
```

<br>
<table>
  <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07.png" alt="Doublets summary"></td>
  </tr>
</table>
<table>
  <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_HPC431_doublet_hist.png" alt="HPC431 doublets"></td>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_K2HO120_doublet_hist.png" alt="K2HO120 doublets"></td>
  </tr>
</table>

<br>

### Visually evaluate QC effects

After annotating MT and doublets, we generate a PCA based on a subset of highly variable genes `hvgs` to visually evaluate the QC effects in PCA space. This applies filtering to remove cells labeled as MT outliers or doublets, and applies min `mingens` and max `maxgenes` thresholds to n_genes_by_counts - the number of genes where count > 0 in a cell. A low n_genes_by_counts value indicates a low quality-cell or empty droplet, whilst a very high n_genes_by_counts value can be indicative of a doublet. It returns before and after PCA plots showing log10 counts, doubelt score and %MT.

``` bash
# Scribble: PCA before/after
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 1 --mem 16G --time=2:00:00 \
    -o ${SCRATCH}/logs/sc_pca.%j.out -e ${SCRATCH}/logs/sc_pca.%j.err \
    scribble pca \
    --project_dir ${SCRATCH} \
    --input ${SCRATCH}/scribble/adata/combined_mtqc_nMADs-8_dblqc_exp-0.07.h5ad \
    --mingenes 100 \
    --maxgenes 9000 \
    --hvgs 3000 \
    --vmax 0.99
```

<br>
<table>
  <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_pca.png" alt="PCA"></td>
  </tr>
</table>

<br>

### Apply filtering

Having reviewed the results, we next filter the data. Here we pass an `xlsx` file which has a `filters` workseet containing fields for `sample`, `min_genes` and `max_genes`. This enables sample-specific filtering thresholds. Alternatively, fixed thresholds could be applied across all samples by setting `--mingenes` and `--maxgenes`. This outputs a new `h5ad` file in `scribble/adata`.

```bash
# Filtering
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 1 --mem 4G --time=2:00:00 \
    -o ${SCRATCH}/logs/sc_filter.%j.out -e ${SCRATCH}/logs/sc_filter.%j.err \
    scribble filter \
    --project_dir ${SCRATCH} \
    --input ${SCRATCH}/scribble/adata/combined_mtqc_nMADs-8_dblqc_exp-0.07.h5ad \
    --filter_xlsx ${SCRATCH}/samples.xlsx
```

<bt>

### Pre-integration processing

Prior to batch interration, we run the `preintegration` tool to pre-process the data. This step preserves raw counts and metadata to avoid later loss during any transformations. It emoves genes expressed in too few cells (n = 3) to reduce noise and sparsity. It calculates HVGs, performs normalisation and log transformation and then subsets the data to the HVGs. It can optionally perform regression to remove effets of covariates, e.g. depth and %MT, and scales data to standardise gene expression (use `--no-scale` to disable). It then performs PCA, generates a KNN graph, and plots UMAP(s) cololured by the specified variables `vars`. This outputs a new `h5ad` file in `scribble/adata`.

```bash
# Pre-integration checks
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 1 --mem 4G --time=2:00:00 \
    -o ${SCRATCH}/logs/sc_preintegration.%j.out -e ${SCRATCH}/logs/sc_preintegration.%j.err \
    scribble preintegration \
    --project_dir ${SCRATCH} \
    --input ${SCRATCH}/scribble/adata/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered.h5ad \
    --min_cells_per_gene 3 \
    --hvgs 3000 \
    --npcs 50 \
    --neighbors 15 \
    --batch sample \
    --vars sample \
    --regress total_counts pct_counts_mt
```

<br>
<table>
  <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_pca_counts.png" alt="PCA counts"></td>
  </tr>
</table>
<table>
  <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_pca_vars.png" alt="PCA sample"></td>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_umap.png" alt="Pre-integration UMAP"></td>
  </tr>
</table>

<br>

### Batch integration

The data is then ready for batch-interation. Currently this is achieved using `Harmony`. To determine the optimal theta for a given dataset, it is necessary to repeat this process for a range of `theta` values. This outputs a new `h5ad` file in `scribble/adata`.

```bash
thetas=(3 6 9 12)
for theta in ${thetas[@]}
do
    # Integration
    sbatch -p uoa-compute --ntasks 1 --cpus-per-task 1 --mem 4G --time=2:00:00 \
        -o ${SCRATCH}/logs/sc_harmony.%j.out -e ${SCRATCH}/logs/sc_harmony_%j.err \
        scribble harmony \
        --project_dir ${SCRATCH} \
        --input ${SCRATCH}/scribble/adata/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration.h5ad \
        --npcs 50 \
        --neighbors 15 \
        --theta ${theta} \
        --batch sample \
        --vars sample
done
```

<br>
<table>
  <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-3_umap.png" alt="theta = 3"></td>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-6_umap.png" alt="theta = 6"></td>
  </tr>
  <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-9_umap.png" alt="theta = 9"></td>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-12_umap.png" alt="theta = 12"></td>
  </tr>
</table>

<br>

### Clustering

Once the integration runs have completed, they are processed to perform Leiden clustering. The optimal resolution can be determined using the `--auto_resolution` method. This requires specifying the lower `--res_min` and upper `--res_max` bounds for the resolution, and the number of resolutions `--res_steps` to test within that range. A silhouette score is calculated from each run and the optimal coarse resolution determined by comparing these scores. That resolution is then used as an anchor to refine resolution based on a `--fine_width`, e.g. if the optimal coarse resolution is 1.0 and `--fine_width 0.2` with `--res_steps 10` then the fine resolution search will have a lower bound of `1.0-0.2`, an upper bound of `1.0+0.2`, and test `10` resolutions within that range. As with the coarse resolution run, silhouette score are calculated for each resolution and the optimal identified and used for clustering. Clustering is performed for `n_repeats`, and the cell to cluster stabiility recorded along with cluster entropy (sample mixing). After clustering, cluster makers are identified for the top `--nmarkers` based on Wilcoxon P value, % expression difference, and log fold change, thereby prioritising significance and specificity of the markers.

```bash
# Once integration compete, cluster
for theta in ${thetas[@]}
do
    # Clustering on Harmony-integrated data
    sbatch -p uoa-compute --ntasks 1 --cpus-per-task 1 --mem 16G --time=2:00:00 \
        -o ${SCRATCH}/logs/sc_cluster.%j.out -e ${SCRATCH}/logs/sc_cluster_%j.err \
        scribble cluster \
        --project_dir ${SCRATCH} \
        --input ${SCRATCH}/scribble/adata/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-${theta}.h5ad \
        --embedding X_pca_harmony \
        --neighbors 15 \
        --auto_resolution \
        --res_min 0.2 \
        --res_max 2.0 \
        --res_steps 10 \
        --fine_width 0.3 \
        --vars sample cluster_stability \
        --n_repeats 10 \
        --nmarkers 100
done
```

<br>

### Evaluate clustering

The resulting `cluster_summary.tsv` outputs are then imported to the `evaluate` tool. This tool evaluates clustering quality, suggests clusters to merge or subset, and outputs a decision table `cluster_summary_decisions.tsv`. Results are scored based on mean stability, mean entropy, low stabiility fraction, and a cluster penalty (if too few or too many clusters). The best score favours high stability, good mixing, few unstable clusters, and a reasonable cluster number.

```bash
# Evaluate clustering
INPUTS=(${SCRATCH}/scribble/tables/*harmony*cluster_summary.tsv)
scribble evaluate \
    --project_dir ${SCRATCH} \
    --input ${INPUTS[@]} \
    --min_cells 200 \
    --large_cells 800 \
    --low_stability 0.75 \
    --high_stability 0.95 \
    --low_entropy 0.5 \
    --merge_size_ratio 2.5 \
    --merge_stability_tol 0.1 \
    --merge_entropy_tol 0.2
```

These results indicate theta = 9 to be optimal.

```tsv
file	score	mean_stability	mean_entropy	low_stability_fraction	n_clusters
theta-9_cluster_summary.tsv	1.2143174350522337	0.9510374553290778	0.526559959446312	0.0	15
theta-12_cluster_summary.tsv	1.1313952975256716	0.9591771675775617	0.5549625756856938	0.05263157894736842	19
theta-6_cluster_summary.tsv	1.000175978412794	0.9482019704201696	0.314474331774722	0.05263157894736842	19
theta-3_cluster_summary.tsv	0.20253667864620173	0.8503816779565261	0.30431000137935127	0.0	2
```

The below plots from `scribble cluster` are from the theta = 9 Harmony run. Note, the cluster labels start at 0.

<br>
<table>
  <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-9_resolution_optimisation.png" alt="Resolution optimisation"></td>
  </tr>
</table>
<table>
    <tr>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-9_clusters.png" alt="Leiden clusters"></td>
    <td><img src="../img/hippo_int/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-9_stability.png" alt="UMAP and stability"></td>
  </tr>
</table>
<br>

The `cluster_summary.tsv` from `scribble cluster` at theta = 9, provided below, indicates for each cluster a range of metrics. These include the number of cells, mean stability, emedian stability, fraction of total cells, sample cell counts, entropy, whether the cluster has low stability or low sample mixing, the clustering resolution, embedding used, and number of clustering iterations undertaken.

```tsv
cluster	n_cells	mean_stability	median_stability	fraction	HPC431	K2HO120	sample_entropy	low_stability	low_mixing	resolution	embedding	n_repeats
1	3989	0.9740035096515418	1.0	0.25993744298188454	3504	485	0.37006987143163084	False	True	0.22777777777777777	X_pca_harmony	10
2	2258	1.0	1.0	0.147139319692428	2037	221	0.3203869033059639	False	True	0.22777777777777777	X_pca_harmony	10
4	2116	0.8223062381852552	0.7	0.1378860940961814	1170	946	0.6875334863432727	False	False	0.22777777777777777	X_pca_harmony	10
6	1830	0.9865027322404372	1.0	0.11924931578261436	1105	725	0.6714301921816759	False	False	0.22777777777777777	X_pca_harmony	10
10	1458	0.9742112482853225	1.0	0.0950084712628698	489	969	0.63793088082967	False	False	0.22777777777777777	X_pca_harmony	10
3	1228	0.9178338762214983	1.0	0.08002085233937183	758	470	0.6653877334401116	False	False	0.22777777777777777	X_pca_harmony	10
0	771	1.0	1.0	0.05024110517398671	635	136	0.46588170850890376	False	True	0.22777777777777777	X_pca_harmony	10
11	421	1.0	1.0	0.027433858986054997	252	169	0.673585297162481	False	False	0.22777777777777777	X_pca_harmony	10
8	394	1.0	1.0	0.02567444285155741	326	68	0.4599674821976243	False	True	0.22777777777777777	X_pca_harmony	10
5	355	0.7907042253521126	0.8	0.02313306399061645	269	86	0.5536636244312005	False	False	0.22777777777777777	X_pca_harmony	10
12	168	1.0	1.0	0.01094747817020722	75	93	0.687396352150425	False	False	0.22777777777777777	X_pca_harmony	10
13	163	0.9000000000000001	0.9	0.010621660367522481	100	63	0.6671581372394269	False	False	0.22777777777777777	X_pca_harmony	10
7	103	0.9	0.9	0.006711846735305617	90	13	0.37912498845561804	False	True	0.22777777777777777	X_pca_harmony	10
9	73	1.0	1.0	0.004756939919197185	27	46	0.6588827340166754	False	False	0.22777777777777777	X_pca_harmony	10
14	19	1.0	1.0	0.001238107650202007	19	0	0.0	False	True	0.22777777777777777	X_pca_harmony	10
```

The `summary_decisions.tsv` output from `scribble evaluate` applied to the theta = 9 results is provided below. This step classifies clusters based on the number of cells, cluster stability, and entropy. In some datasets a recommendation will be made to subset or merge clusetrs.

```tsv
cluster	action	reason	detail	priority	merge_group
1	flag_bias	sample_specific_cluster	n=3989; stability=0.97; entropy=0.37	medium
2	flag_bias	sample_specific_cluster	n=2258; stability=1.00; entropy=0.32	medium
4	subset	heterogeneous_large_cluster	n=2116; stability=0.82; entropy=0.69	high
6	keep	well_defined_cluster	n=1830; stability=0.99; entropy=0.67	low
10	keep	well_defined_cluster	n=1458; stability=0.97; entropy=0.64	low
3	subset	heterogeneous_large_cluster	n=1228; stability=0.92; entropy=0.67	high
0	flag_bias	sample_specific_cluster	n=771; stability=1.00; entropy=0.47	medium
11	keep	well_defined_cluster	n=421; stability=1.00; entropy=0.67	low
8	flag_bias	sample_specific_cluster	n=394; stability=1.00; entropy=0.46	medium
5	keep	well_defined_cluster	n=355; stability=0.79; entropy=0.55	low
12	keep	small_cluster	n=168; stability=1.00; entropy=0.69	low
13	keep	small_cluster	n=163; stability=0.90; entropy=0.67	low
7	keep	small_cluster	n=103; stability=0.90; entropy=0.38	low
9	keep	small_cluster	n=73; stability=1.00; entropy=0.66	low
14	keep	small_cluster	n=19; stability=1.00; entropy=0.00	low
```

<br>

### Refine clustering

The `cluster_summary_decisions.tsv` file can then be parsed to refine the clustering. This fully re-processes the cells from clusters that are tagged for subset or merging. Briefly, it subsets the cells, identifies HVGs, runs normalisation, etc.., re-integration with Harmony, clustering, builds hierarchical cluster labels and extracts refined markers (within lineage and global).

```bash
# Refine clusters
theta=9 # optimal
sbatch -p uoa-compute --ntasks 1 --cpus-per-task 1 --mem 16G --time=2:00:00 \
    -o ${SCRATCH}/logs/sc_refine.%j.out -e ${SCRATCH}/logs/sc_refine_%j.err \
    scribble refine \
        --project_dir ${SCRATCH} \
        --decisions ${SCRATCH}/scribble/tables/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-${theta}_cluster_summary_decisions.tsv \
        --input ${SCRATCH}/scribble/adata/combined_mtqc_nMADs-8_dblqc_exp-0.07_filtered_preintegration_harmony_theta-${theta}_clustered.h5ad \
        --min_cells_per_gene 3 \
        --batch sample \
        --hvgs 3000 \
        --npcs 50 \
        --neighbors 15 \
        --theta 7 \
        --auto_resolution \
        --res_min 0.1 \
        --res_max 1.5 \
        --res_steps 10 \
        --fine_width 0.2 \
        --min_cells_per_group 500 \
        --n_repeats 10 \
        --nmarkers 100 \
        --max_refine_depth 5 \
        --stability_threshold 0.9 \
        --min_cells_per_cluster 50 \
        --marker_strength_threshold 1.0
```
