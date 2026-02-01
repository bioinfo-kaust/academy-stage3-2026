#!/usr/bin/env python3
# =============================================================================
# 10_visualization.py - Create Plots
# =============================================================================
# Run with: python3 10_visualization.py
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

print("==============================================")
print("Step 10: Visualization")
print("==============================================")

# -----------------------------------------------------------------------------
# Create output directory
# -----------------------------------------------------------------------------
os.makedirs(os.path.expanduser("~/genomics/results/figures"), exist_ok=True)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
print("Loading data...")

de_results = pd.read_csv(os.path.expanduser("~/genomics/results/tables/deseq2_all_results.csv"))
norm_counts = pd.read_csv(os.path.expanduser("~/genomics/results/tables/normalized_counts.csv"), index_col="gene_id")

samples = ["SRR975551", "SRR975552", "SRR975553", "SRR975554", "SRR975555", "SRR975556"]
conditions = ["Normal", "Normal", "Normal", "Tumor", "Tumor", "Tumor"]

print(f"  Loaded {len(de_results)} genes")

# -----------------------------------------------------------------------------
# 1. Volcano Plot
# -----------------------------------------------------------------------------
print("Creating volcano plot...")

fig, ax = plt.subplots(figsize=(10, 8))

de_results["-log10_padj"] = -np.log10(de_results["padj"].fillna(1))
de_results["significant"] = (de_results["padj"] < 0.05) & (abs(de_results["log2FoldChange"]) > 1)

colors = ["gray"] * len(de_results)
for i in range(len(de_results)):
    if de_results["significant"].iloc[i]:
        if de_results["log2FoldChange"].iloc[i] > 0:
            colors[i] = "red"
        else:
            colors[i] = "blue"

ax.scatter(de_results["log2FoldChange"], de_results["-log10_padj"], c=colors, alpha=0.5, s=10)
ax.axhline(y=-np.log10(0.05), color="black", linestyle="--", linewidth=0.5)
ax.axvline(x=1, color="black", linestyle="--", linewidth=0.5)
ax.axvline(x=-1, color="black", linestyle="--", linewidth=0.5)
ax.set_xlabel("log2 Fold Change (Tumor vs Normal)")
ax.set_ylabel("-log10(adjusted p-value)")
ax.set_title("Volcano Plot")

plt.tight_layout()
plt.savefig(os.path.expanduser("~/genomics/results/figures/volcano_plot.png"), dpi=150)
plt.close()
print("  Saved volcano_plot.png")

# -----------------------------------------------------------------------------
# 2. MA Plot
# -----------------------------------------------------------------------------
print("Creating MA plot...")

fig, ax = plt.subplots(figsize=(10, 8))

ax.scatter(np.log10(de_results["baseMean"] + 1), de_results["log2FoldChange"], c=colors, alpha=0.5, s=10)
ax.axhline(y=0, color="black", linestyle="-", linewidth=0.5)
ax.axhline(y=1, color="black", linestyle="--", linewidth=0.5)
ax.axhline(y=-1, color="black", linestyle="--", linewidth=0.5)
ax.set_xlabel("log10(Mean Expression)")
ax.set_ylabel("log2 Fold Change")
ax.set_title("MA Plot")

plt.tight_layout()
plt.savefig(os.path.expanduser("~/genomics/results/figures/ma_plot.png"), dpi=150)
plt.close()
print("  Saved ma_plot.png")

# -----------------------------------------------------------------------------
# 3. Heatmap of Top Genes
# -----------------------------------------------------------------------------
print("Creating heatmap...")

sig_genes = de_results[de_results["significant"]].nsmallest(50, "padj")

if len(sig_genes) > 0:
    top_genes = sig_genes["gene_id"].tolist()
    heatmap_data = norm_counts.loc[norm_counts.index.isin(top_genes), samples]

    # Log transform and z-score
    heatmap_data = np.log2(heatmap_data + 1)
    heatmap_data = (heatmap_data.T - heatmap_data.mean(axis=1)).T
    heatmap_data = (heatmap_data.T / heatmap_data.std(axis=1)).T

    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(heatmap_data, cmap="RdBu_r", center=0, ax=ax)
    ax.set_title(f"Top {len(heatmap_data)} Differentially Expressed Genes")

    plt.tight_layout()
    plt.savefig(os.path.expanduser("~/genomics/results/figures/heatmap.png"), dpi=150)
    plt.close()
    print("  Saved heatmap.png")

# -----------------------------------------------------------------------------
# 4. PCA Plot
# -----------------------------------------------------------------------------
print("Creating PCA plot...")

from sklearn.decomposition import PCA

pca_data = np.log2(norm_counts[samples] + 1).T
pca = PCA(n_components=2)
pca_result = pca.fit_transform(pca_data)

fig, ax = plt.subplots(figsize=(10, 8))

color_map = {"Normal": "#3498db", "Tumor": "#e74c3c"}
for i, (sample, condition) in enumerate(zip(samples, conditions)):
    ax.scatter(pca_result[i, 0], pca_result[i, 1], c=color_map[condition], s=100)
    ax.annotate(sample, (pca_result[i, 0], pca_result[i, 1]), xytext=(5, 5), textcoords="offset points")

ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("PCA of Gene Expression")

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor="#3498db", label="Normal"),
                   Patch(facecolor="#e74c3c", label="Tumor")]
ax.legend(handles=legend_elements)

plt.tight_layout()
plt.savefig(os.path.expanduser("~/genomics/results/figures/pca_plot.png"), dpi=150)
plt.close()
print("  Saved pca_plot.png")

# -----------------------------------------------------------------------------
# 5. Bar Charts for Top 10 DE Genes
# -----------------------------------------------------------------------------
print("Creating bar charts for top 10 DE genes...")

# Ensure 'padj' column is correctly handled for nsmallest
# It's already done for sig_genes, so we can reuse that logic or re-filter
# We need to make sure we're using the 'de_results' DataFrame after filling NaNs in 'padj'
# and calculating 'significant' status.
top10_de_genes = de_results[de_results["significant"]].nsmallest(10, "padj")

if not top10_de_genes.empty:
    for _, row in top10_de_genes.iterrows():
        gene_id = row["gene_id"]
        log2fc = row["log2FoldChange"]
        padj = row["padj"]

        # Get normalized counts for this gene
        gene_counts = norm_counts.loc[gene_id, samples]

        # Create a DataFrame for seaborn
        plot_data = pd.DataFrame({
            "Sample": samples,
            "Condition": conditions,
            "Normalized Count": gene_counts.values
        })

        fig, ax = plt.subplots(figsize=(8, 6))
        sns.barplot(x="Sample", y="Normalized Count", hue="Condition", data=plot_data, ax=ax, palette={"Normal": "#3498db", "Tumor": "#e74c3c"})
        ax.set_title(f"Gene: {gene_id}\nlog2FC: {log2fc:.2f}, adj. p-value: {padj:.2e}")
        ax.set_ylabel("Normalized Count")
        ax.set_xlabel("Sample")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(os.path.expanduser(f"~/genomics/results/figures/barchart_{gene_id}.png"), dpi=150)
        plt.close()
        print(f"  Saved barchart_{gene_id}.png")
else:
    print("  No significant genes found to create bar charts for.")

print("==============================================")
print("Visualization complete!")
print("Figures in: ~/genomics/results/figures/")
print("==============================================")
