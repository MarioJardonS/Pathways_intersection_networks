# Scripts/03_compute_unifrac_for_network_pairs.py

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skbio import TreeNode
from skbio.diversity import beta_diversity


TABLE_FILE = Path("Data/pathway_sample_unifrac_abundance_table.tsv")
TREE_FILE = Path("Data/pathway_sample_tree_pruned.nwk")
METADATA_FILE = Path("Data/metadata.csv")
NETWORKS_DIR = Path("Networks")

OUT_DIR = Path("Data/data_samples/network_unifrac")
PLOT_DIR = Path("Plots/network_unifrac")

OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

N_PER_GROUP = 'none' #5
GROUPS_TO_USE = ["Obese", "Healthy"]
TOP_N = 20


def pathway_id(x):
    return str(x).split(":", 1)[0].strip()


def compute_pairwise_unifrac(array_x, array_y, taxa, tree, metric):
    counts = np.vstack([
        np.asarray(array_x, dtype=float),
        np.asarray(array_y, dtype=float)
    ])

    dm = beta_diversity(
        metric=metric,
        counts=counts,
        ids=["x", "y"],
        taxa=taxa,
        tree=tree,
    )

    return dm["x", "y"]

def plot_top_nonzero_unweighted(results, sample_id, diagnosis, plot_dir, top_n=20):
    df = results.dropna(subset=["Unweighted_UniFrac"]).copy()
    df = df[df["Unweighted_UniFrac"] > 0].copy()

    if df.empty:
        print(f"No nonzero Unweighted_UniFrac values to plot for {sample_id}.")
        return

    df["PairID"] = df["PathwayA_ID"] + " vs " + df["PathwayB_ID"]

    closest = df.sort_values("Unweighted_UniFrac", ascending=True).head(top_n).copy()
    distant = df.sort_values("Unweighted_UniFrac", ascending=False).head(top_n).copy()

    prefix = f"{sample_id}_{diagnosis}"

    closest.to_csv(plot_dir / f"{prefix}_unweighted_closest_pairs_top{top_n}.tsv", sep="\t", index=False)
    distant.to_csv(plot_dir / f"{prefix}_unweighted_most_distant_pairs_top{top_n}.tsv", sep="\t", index=False)

    if not closest.empty:
        plot_df = closest.sort_values("Unweighted_UniFrac", ascending=False)

        plt.figure(figsize=(10, 8))
        plt.barh(plot_df["PairID"], plot_df["Unweighted_UniFrac"])
        plt.xlabel("Unweighted UniFrac distance")
        plt.ylabel("Pathway pair")
        plt.title(f"{sample_id} ({diagnosis}): closest nonzero pathway pairs")
        plt.tight_layout()
        plt.savefig(plot_dir / f"{prefix}_unweighted_closest_pairs_top{top_n}.png", dpi=300)
        plt.close()

    if not distant.empty:
        plot_df = distant.sort_values("Unweighted_UniFrac", ascending=True)

        plt.figure(figsize=(10, 8))
        plt.barh(plot_df["PairID"], plot_df["Unweighted_UniFrac"])
        plt.xlabel("Unweighted UniFrac distance")
        plt.ylabel("Pathway pair")
        plt.title(f"{sample_id} ({diagnosis}): most distant pathway pairs")
        plt.tight_layout()
        plt.savefig(plot_dir / f"{prefix}_unweighted_most_distant_pairs_top{top_n}.png", dpi=300)
        plt.close()


def compute_unifrac_for_network_file(network_file, diagnosis, table, tree, otu_cols):
    sample_id = Path(network_file).stem

    sample_table = table[table["SampleID"] == sample_id].copy()

    if sample_table.empty:
        print(f"Skipping {sample_id}: sample not found in abundance table.")
        return None

    sample_table = sample_table.set_index("PathwayID", drop=False)
    available_pathways = set(sample_table["PathwayID"])

    net = pd.read_csv(network_file, sep="\t")

    pair_cols = net.columns[:2]
    pairs = net[list(pair_cols)].copy()
    pairs.columns = ["PathwayA", "PathwayB"]

    pairs = pairs[
        ~pairs["PathwayA"].astype(str).isin(["UNINTEGRATED", "UNMAPPED"]) &
        ~pairs["PathwayB"].astype(str).isin(["UNINTEGRATED", "UNMAPPED"])
    ].copy()

    pairs["PathwayA_ID"] = pairs["PathwayA"].map(pathway_id)
    pairs["PathwayB_ID"] = pairs["PathwayB"].map(pathway_id)

    pairs = pairs[pairs["PathwayA_ID"] != pairs["PathwayB_ID"]].copy()

    pairs["pair_key"] = pairs.apply(
        lambda r: "__".join(sorted([r["PathwayA_ID"], r["PathwayB_ID"]])),
        axis=1
    )

    pairs = pairs.drop_duplicates("pair_key").copy()

    results = []

    for _, row in pairs.iterrows():
        pathway_a_id = row["PathwayA_ID"]
        pathway_b_id = row["PathwayB_ID"]

        pathway_a_available = pathway_a_id in available_pathways
        pathway_b_available = pathway_b_id in available_pathways

        weighted = np.nan
        unweighted = np.nan

        pathway_a_total = np.nan
        pathway_b_total = np.nan
        pathway_a_nonzero = np.nan
        pathway_b_nonzero = np.nan
        shared_otus = np.nan
        union_otus = np.nan

        pathway_a_name = row["PathwayA"]
        pathway_b_name = row["PathwayB"]

        if pathway_a_available and pathway_b_available:
            row_a = sample_table.loc[pathway_a_id]
            row_b = sample_table.loc[pathway_b_id]

            array_a = row_a[otu_cols].values
            array_b = row_b[otu_cols].values

            pathway_a_name = row_a["Pathway"]
            pathway_b_name = row_b["Pathway"]

            weighted = compute_pairwise_unifrac(
                array_a, array_b,
                taxa=otu_cols,
                tree=tree,
                metric="weighted_unifrac"
            )

            unweighted = compute_pairwise_unifrac(
                array_a, array_b,
                taxa=otu_cols,
                tree=tree,
                metric="unweighted_unifrac"
            )

            pathway_a_total = np.sum(array_a)
            pathway_b_total = np.sum(array_b)
            pathway_a_nonzero = np.sum(array_a > 0)
            pathway_b_nonzero = np.sum(array_b > 0)
            shared_otus = np.sum((array_a > 0) & (array_b > 0))
            union_otus = np.sum((array_a > 0) | (array_b > 0))

        results.append({
            "SampleID": sample_id,
            "Diagnosis": diagnosis,
            "PathwayA_ID": pathway_a_id,
            "PathwayB_ID": pathway_b_id,
            "PairKey": row["pair_key"],
            "PathwayA": pathway_a_name,
            "PathwayB": pathway_b_name,
            "Weighted_UniFrac": weighted,
            "Unweighted_UniFrac": unweighted,
            "PathwayA_available": pathway_a_available,
            "PathwayB_available": pathway_b_available,
            "PathwayA_total_abundance": pathway_a_total,
            "PathwayB_total_abundance": pathway_b_total,
            "PathwayA_nonzero_otus": pathway_a_nonzero,
            "PathwayB_nonzero_otus": pathway_b_nonzero,
            "Shared_nonzero_otus": shared_otus,
            "Union_nonzero_otus": union_otus,
        })

    results = pd.DataFrame(results)

    out_file = OUT_DIR / f"{sample_id}_{diagnosis}_network_unifrac.tsv"
    results.to_csv(out_file, sep="\t", index=False)

    print("\nSample:", sample_id, diagnosis)
    print("Network pairs:", len(pairs))
    print("Computable pairs:", results["Unweighted_UniFrac"].notna().sum())
    print("Pairs with NaN:", results["Unweighted_UniFrac"].isna().sum())
    print("Output:", out_file)

    plot_top_nonzero_unweighted(
        results=results,
        sample_id=sample_id,
        diagnosis=diagnosis,
        plot_dir=PLOT_DIR,
        top_n=TOP_N
    )

    return results


def plot_common_pairs_by_condition(common_df):
    df = common_df.dropna(subset=["Unweighted_UniFrac"]).copy()
    df = df[df["Unweighted_UniFrac"] > 0].copy()

    if df.empty:
        print("No nonzero common-pair values to plot.")
        return

    summary = (
        df.groupby(["Diagnosis", "PairKey"])["Unweighted_UniFrac"]
        .mean()
        .reset_index()
    )

    pivot = summary.pivot(
        index="PairKey",
        columns="Diagnosis",
        values="Unweighted_UniFrac"
    )

    if not all(g in pivot.columns for g in GROUPS_TO_USE):
        print("Not all groups are present in common-pair pivot.")
        return

    pivot["Delta_Obese_minus_Healthy"] = pivot["Obese"] - pivot["Healthy"]
    pivot = pivot.dropna()

    out_summary = OUT_DIR / "common_pairs_mean_unweighted_by_condition.tsv"
    pivot.to_csv(out_summary, sep="\t")

    top_delta = pivot.reindex(
        pivot["Delta_Obese_minus_Healthy"].abs()
        .sort_values(ascending=False)
        .head(TOP_N)
        .index
    ).copy()

    plt.figure(figsize=(10, 8))
    plt.barh(
        top_delta.index,
        top_delta["Delta_Obese_minus_Healthy"]
    )
    plt.xlabel("Mean Unweighted UniFrac difference: Obese - Healthy")
    plt.ylabel("Pathway pair")
    plt.title("Common pathway pairs with largest condition difference")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "common_pairs_obese_vs_healthy_delta_top20.png", dpi=300)
    plt.close()

    print("\nCommon-pair condition summary saved:")
    print(out_summary)
    print(PLOT_DIR / "common_pairs_obese_vs_healthy_delta_top20.png")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

table = pd.read_csv(TABLE_FILE, sep="\t")
tree = TreeNode.read(str(TREE_FILE))
meta = pd.read_csv(METADATA_FILE, sep=",")

otu_cols = table.columns[3:].tolist()

table_samples = set(table["SampleID"].astype(str))

network_files = sorted(NETWORKS_DIR.glob("*.tsv"))
network_files = [
    p for p in network_files
    if p.name not in ["edges.tsv", "network.tsv"]
]

network_sample_ids = set(p.stem for p in network_files)

meta["SampleID"] = meta["SampleID"].astype(str)

valid_meta = meta[
    meta["SampleID"].isin(table_samples) &
    meta["SampleID"].isin(network_sample_ids) &
    meta["Diagnosis"].isin(GROUPS_TO_USE)
].copy()

# # Selects N files per group
# selected_meta = (
#     valid_meta
#     .groupby("Diagnosis", group_keys=False)
#     .head(N_PER_GROUP)
#     .copy()
# )

selected_meta = valid_meta[
    valid_meta["Diagnosis"].isin(GROUPS_TO_USE)
].copy()

selected_meta = valid_meta.copy()

if N_PER_GROUP is not None:
    selected_meta = (
        selected_meta
        .groupby("Diagnosis", group_keys=False)
        .head(N_PER_GROUP)
        .copy()
    )

print("\nSelected samples:")
print(selected_meta[["SampleID", "Diagnosis"]].to_string(index=False))

sample_to_diagnosis = dict(zip(selected_meta["SampleID"], selected_meta["Diagnosis"]))

selected_network_files = [
    p for p in network_files
    if p.stem in sample_to_diagnosis
]

all_results = []

for network_file in selected_network_files:
    sample_id = network_file.stem
    diagnosis = sample_to_diagnosis[sample_id]

    res = compute_unifrac_for_network_file(
        network_file=network_file,
        diagnosis=diagnosis,
        table=table,
        tree=tree,
        otu_cols=otu_cols
    )

    if res is not None and not res.empty:
        all_results.append(res)

if len(all_results) == 0:
    raise SystemExit("No results generated.")

all_results = pd.concat(all_results, ignore_index=True)

all_results = all_results.merge(
    selected_meta[
        ["SampleID", "Diagnosis", "Project", "GMHI", "hiPCA", "Shannon_entropy"]
    ],
    on=["SampleID", "Diagnosis"],
    how="left"
)

# out_all = OUT_DIR / "pilot_5_obese_5_healthy_network_unifrac_all_pairs.tsv"
out_all = OUT_DIR / "obese_healthy_network_unifrac_all_pairs.tsv"
all_results.to_csv(out_all, sep="\t", index=False)

print("\nCombined output saved:")
print(out_all)

print("\nSamples per diagnosis:")
print(
    all_results[["SampleID", "Diagnosis"]]
    .drop_duplicates()["Diagnosis"]
    .value_counts()
)

print("\nUnweighted UniFrac by diagnosis:")
print(
    all_results
    .dropna(subset=["Unweighted_UniFrac"])
    .groupby("Diagnosis")["Unweighted_UniFrac"]
    .describe()
)

# ---------------------------------------------------------------------
# Per-sample summary
# ---------------------------------------------------------------------

sample_summary = (
    all_results
    .dropna(subset=["Unweighted_UniFrac"])
    .groupby(["SampleID", "Diagnosis"])["Unweighted_UniFrac"]
    .agg(["mean", "median", "std", "count"])
    .reset_index()
)

out_sample_summary = OUT_DIR / "pilot_5_obese_5_healthy_sample_summary.tsv"
sample_summary.to_csv(out_sample_summary, sep="\t", index=False)

print("\nPer-sample Unweighted UniFrac summary:")
print(sample_summary.to_string(index=False))

print("\nPer-sample summary by diagnosis:")
print(
    sample_summary
    .groupby("Diagnosis")[["mean", "median"]]
    .describe()
)

print("\nSample summary saved:")
print(out_sample_summary)



# Keep only pathway pairs represented in both conditions
pair_condition_counts = (
    all_results
    .dropna(subset=["Unweighted_UniFrac"])
    .groupby("PairKey")["Diagnosis"]
    .nunique()
)

common_pair_keys = pair_condition_counts[pair_condition_counts == len(GROUPS_TO_USE)].index

common_results = all_results[
    all_results["PairKey"].isin(common_pair_keys)
].copy()

out_common = OUT_DIR / "pilot_5_obese_5_healthy_common_pathway_pairs.tsv"
common_results.to_csv(out_common, sep="\t", index=False)

print("\nCommon pairs output saved:")
print(out_common)
print("Common pathway pairs:", common_results["PairKey"].nunique())

plot_common_pairs_by_condition(common_results)