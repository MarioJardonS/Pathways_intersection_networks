# 02_prepare_unifrac_inputs.py

from pathlib import Path
import pandas as pd
from skbio import TreeNode

print("RUNNING 02_prepare_unifrac_inputs.py")

RAW_FILE    = Path("Data/relab_table_joined.tsv")
MAPPED_FILE = Path("Data/unifrac/mapped_exact.tsv")
GTDB_TREE   = Path("Data/unifrac/bac120.tree.gz")

PATHWAY_TYPES_FILE = Path("Data/Pathways_Types_complete.txt")

OUT_TABLE = Path("Data/pathway_sample_unifrac_abundance_table.tsv")
OUT_TREE  = Path("Data/pathway_sample_tree_pruned.nwk")

FEATURE_COL = "# Pathway"
EXAMPLE_SAMPLE = "ERR209069"


def normalize_accession(x):
    return str(x).replace("_", " ").strip()


# 1. Load relative abundance table
df = pd.read_csv(RAW_FILE, sep="\t")

print("Input table:", df.shape)



# ---------------------------------------------------------------------
# CHECK 0: abundance scale in raw relab_table_joined.tsv before filtering
# ---------------------------------------------------------------------

raw_original_sample_cols = df.columns[1:]

raw_original_sums = df[raw_original_sample_cols].sum(axis=0)

print("\nRaw abundance scale check before filtering")
print("Column sums across all rows:")
print(raw_original_sums.describe())

print("\nTop 10 raw sample sums:")
print(raw_original_sums.sort_values(ascending=False).head(10).to_string())

print("\nBottom 10 raw sample sums:")
print(raw_original_sums.sort_values(ascending=True).head(10).to_string())






# 2. Split first column: Pathway|taxon
split = df[FEATURE_COL].astype(str).str.split("|", n=1, expand=True)

df["Pathway"] = split[0]
df["original_taxon"] = split[1]

# PathwayID = part before ":"
df["PathwayID"] = (
    df["Pathway"]
    .astype(str)
    .str.split(":", n=1)
    .str[0]
    .str.strip()
)

# 3. Remove rows without taxon and non-biological aggregate categories
df = df[df["original_taxon"].notna()].copy()
df = df[~df["PathwayID"].isin(["UNINTEGRATED", "UNMAPPED"])].copy()

print("Rows with Pathway|taxon after removing UNINTEGRATED/UNMAPPED:", len(df))

# 4. Load GTDB mapping
mapped = pd.read_csv(MAPPED_FILE, sep="\t")

otu_to_acc = dict(
    zip(
        mapped["original_taxon"],
        mapped["gtdb_accession"].map(normalize_accession)
    )
)

# 5. Keep only taxa mapped exactly to GTDB
df = df[df["original_taxon"].isin(otu_to_acc)].copy()

# 6. Replace original taxon names by GTDB accessions
df["accession"] = df["original_taxon"].map(otu_to_acc)

print("Rows after GTDB exact mapping:", len(df))
print("Mapped taxa/accessions:", df["accession"].nunique())
print("Mapped original taxa:", df["original_taxon"].nunique())
print("Pathway IDs:", df["PathwayID"].nunique())
print("Pathway names:", df["Pathway"].nunique())

# 7. Sample columns are all original abundance columns
sample_cols = [
    c for c in df.columns
    if c not in [FEATURE_COL, "Pathway", "PathwayID", "original_taxon", "accession"]
]

print("Samples before cleaning names:", len(sample_cols))


# ---------------------------------------------------------------------
# CHECK 2: abundance scale in relab_table_joined.tsv
# ---------------------------------------------------------------------

raw_sample_sums = df[sample_cols].sum(axis=0)
raw_sample_nonzero = (df[sample_cols] > 0).sum(axis=0)

print("\nAbundance scale check after filtering/mapping")
print("Column sums across remaining Pathway|taxon rows:")
print(raw_sample_sums.describe())

print("\nTop 10 samples by total abundance after filtering/mapping:")
print(raw_sample_sums.sort_values(ascending=False).head(10).to_string())

print("\nBottom 10 samples by total abundance after filtering/mapping:")
print(raw_sample_sums.sort_values(ascending=True).head(10).to_string())

print("\nNonzero Pathway|taxon entries per sample:")
print(raw_sample_nonzero.describe())
# ---------------------------------------------------------------------





# 8. Wide to long format
long_df = df.melt(
    id_vars=["PathwayID", "Pathway", "accession"],
    value_vars=sample_cols,
    var_name="SampleID",
    value_name="abundance"
)

# Clean sample IDs: ERR209069_concat_Abundance -> ERR209069
long_df["SampleID"] = (
    long_df["SampleID"]
    .astype(str)
    .str.replace("_concat_Abundance", "", regex=False)
)

print("Samples after cleaning names:", long_df["SampleID"].nunique())



# ---------------------------------------------------------------------
# CHECK 1: sample names before vs after cleaning
# ---------------------------------------------------------------------

sample_name_check = pd.DataFrame({
    "original_sample_col": sample_cols,
    "clean_sample_id": (
        pd.Series(sample_cols)
        .astype(str)
        .str.replace("_concat_Abundance", "", regex=False)
    )
})

duplicated_clean_ids = sample_name_check[
    sample_name_check["clean_sample_id"].duplicated(keep=False)
].sort_values("clean_sample_id")

print("\nSample name cleaning check")
print("Original sample columns:", sample_name_check["original_sample_col"].nunique())
print("Clean sample IDs:", sample_name_check["clean_sample_id"].nunique())
print("Duplicated clean sample IDs:", duplicated_clean_ids["clean_sample_id"].nunique())

if len(duplicated_clean_ids) > 0:
    print("\nDuplicated sample IDs after cleaning:")
    print(duplicated_clean_ids.to_string(index=False))

# ---------------------------------------------------------------------



# 9. Build final table:
# rows    = PathwayID x Pathway x SampleID
# columns = GTDB accessions
# values  = relative abundance
# IMPORTANT: no normalization is applied.
pathway_sample_table = long_df.pivot_table(
    index=["PathwayID", "Pathway", "SampleID"],
    columns="accession",
    values="abundance",
    aggfunc="sum",
    fill_value=0
)

# 10. Remove rows where all taxa are zero
row_sums = pathway_sample_table.sum(axis=1)
pathway_sample_table = pathway_sample_table.loc[row_sums > 0].copy()

# 11. Convert MultiIndex back to columns
pathway_sample_table = pathway_sample_table.reset_index()

# 12. Save abundance table
pathway_sample_table.to_csv(OUT_TABLE, sep="\t", index=False)

print("\nOutput table:", OUT_TABLE)
print("Shape:", pathway_sample_table.shape)
print("Pathway IDs:", pathway_sample_table["PathwayID"].nunique())
print("Pathway names:", pathway_sample_table["Pathway"].nunique())
print("Samples:", pathway_sample_table["SampleID"].nunique())
print("Pathway-Sample rows:", len(pathway_sample_table))
print("Taxa columns:", pathway_sample_table.shape[1] - 3)

# 13. OTU columns
otu_cols = pathway_sample_table.columns[3:]

row_sums_check = pathway_sample_table[otu_cols].sum(axis=1)

print("\nRow sums across OTUs:")
print(row_sums_check.describe())

print("\nFirst rows:")
print(pathway_sample_table.iloc[:5, :9].to_string(index=False))

# 14. Load and prune tree
tree = TreeNode.read(str(GTDB_TREE))

tips_needed = set(otu_cols)

pruned_tree = tree.shear(tips_needed)

with open(OUT_TREE, "w") as f:
    pruned_tree.write(f)

print("\nPruned tree saved:", OUT_TREE)
print("Tips needed:", len(tips_needed))
print("Tips in pruned tree:", len([tip.name for tip in pruned_tree.tips()]))

# ---------------------------------------------------------------------
# SANITY CHECKS BEFORE UNIFRAC
# ---------------------------------------------------------------------

print("\n" + "=" * 70)
print("SANITY CHECKS BEFORE UNIFRAC")
print("=" * 70)

# Basic numbers
n_pathways = pathway_sample_table["PathwayID"].nunique()
n_samples = pathway_sample_table["SampleID"].nunique()
n_rows = len(pathway_sample_table)
n_otus = pathway_sample_table.shape[1] - 3

print("\nBasic counts")
print("Number of pathway IDs:", n_pathways)
print("Number of samples:", n_samples)
print("Number of Pathway-Sample rows:", n_rows)
print("Number of OTUs:", n_otus)

# Mapped taxa
print("\nMapped taxa")
print("Number of mapped taxa/accessions:", df["accession"].nunique())
print("Number of mapped original taxa:", df["original_taxon"].nunique())

# Load official pathway list
pathway_types = pd.read_csv(
    PATHWAY_TYPES_FILE,
    sep=r"\s+|\t",
    engine="python",
    header=None,
    dtype=str
)

official_pathways = (
    pathway_types.iloc[:, 0]
    .dropna()
    .astype(str)
    .str.strip()
)

official_set = set(official_pathways)

represented_set = set(
    pathway_sample_table["PathwayID"]
    .dropna()
    .astype(str)
    .str.strip()
)

represented_official = represented_set.intersection(official_set)
excluded_official = official_set - represented_set
extra_represented = represented_set - official_set

print("\nPathways_Types_complete.txt check")
print("Official pathways:", len(official_set))
print("Official pathways represented:", len(represented_official))
print("Official pathways excluded:", len(excluded_official))
print("Represented pathway IDs not in official list:", len(extra_represented))

print("\nFirst 20 represented official pathways:")
print(pd.Series(sorted(represented_official)).head(20).to_string(index=False))

print("\nFirst 20 excluded official pathways:")
print(pd.Series(sorted(excluded_official)).head(20).to_string(index=False))

print("\nFirst 20 represented pathway IDs not in official list:")
print(pd.Series(sorted(extra_represented)).head(20).to_string(index=False))

# Row sums
row_sums = pathway_sample_table[otu_cols].sum(axis=1)

print("\nRow sums across OTUs")
print(row_sums.describe())

print("\nTop 10 Pathway-Sample rows by total abundance:")
top_rows = pathway_sample_table[["PathwayID", "Pathway", "SampleID"]].copy()
top_rows["row_sum"] = row_sums

print(
    top_rows
    .sort_values("row_sum", ascending=False)
    .head(10)
    .to_string(index=False)
)

print("\nBottom 10 nonzero Pathway-Sample rows by total abundance:")
print(
    top_rows
    .sort_values("row_sum", ascending=True)
    .head(10)
    .to_string(index=False)
)









# ---------------------------------------------------------------------
# CHECK 3: find valid sample files in Networks/
# ---------------------------------------------------------------------

NETWORKS_DIR = Path("Networks")

network_files = sorted(NETWORKS_DIR.glob("*.tsv"))
network_sample_ids = [p.stem for p in network_files]

table_sample_ids = set(pathway_sample_table["SampleID"].astype(str))

valid_network_samples = sorted(
    set(network_sample_ids).intersection(table_sample_ids)
)

print("\nNetwork sample file check")
print("Network .tsv files:", len(network_sample_ids))
print("Samples in pathway_sample_table:", len(table_sample_ids))
print("Valid samples present in both:", len(valid_network_samples))

print("\nFirst 20 valid samples:")
print(valid_network_samples[:20])

print("\nFirst 20 network samples not in table:")
print(sorted(set(network_sample_ids) - table_sample_ids)[:20])

print("\nFirst 20 table samples not in Networks:")
print(sorted(table_sample_ids - set(network_sample_ids))[:20])

# ELEGIR SAMPLE VALIDO
if len(valid_network_samples) > 0:
    EXAMPLE_SAMPLE = valid_network_samples[0]
    print("\nUsing valid example sample:", EXAMPLE_SAMPLE)
else:
    raise SystemExit("No valid sample found between Networks/ and pathway_sample_table.")


# Example sample 
available_samples = pathway_sample_table["SampleID"].unique().tolist()

if EXAMPLE_SAMPLE not in available_samples:
    matches = [
        s for s in available_samples
        if EXAMPLE_SAMPLE in str(s)
    ]

    print(f"\nSample {EXAMPLE_SAMPLE} was not found.")
    print("Possible matches:")
    print(matches[:20])

    if len(matches) > 0:
        EXAMPLE_SAMPLE = matches[0]
        print("Using first match:", EXAMPLE_SAMPLE)
    else:
        EXAMPLE_SAMPLE = None

if EXAMPLE_SAMPLE is not None:
    example_df = pathway_sample_table[
        pathway_sample_table["SampleID"] == EXAMPLE_SAMPLE
    ].copy()

    example_sums = example_df[otu_cols].sum(axis=1)

    example_summary = example_df[["PathwayID", "Pathway", "SampleID"]].copy()
    example_summary["row_sum"] = example_sums
    example_summary["n_nonzero_otus"] = (example_df[otu_cols] > 0).sum(axis=1)

    print(f"\nExample sample: {EXAMPLE_SAMPLE}")
    print("Rows/pathways in this sample:", len(example_summary))
    print("Pathways with row_sum > 0:", (example_summary["row_sum"] > 0).sum())

    print("\nTop 20 pathways in example sample by total abundance:")
    print(
        example_summary
        .sort_values("row_sum", ascending=False)
        .head(20)
        .to_string(index=False)
    )

    print("\nTop OTUs for the most abundant pathway in example sample:")
    top_pathway_id = (
        example_summary
        .sort_values("row_sum", ascending=False)
        .iloc[0]["PathwayID"]
    )

    top_pathway_row = example_df[
        example_df["PathwayID"] == top_pathway_id
    ].iloc[0]

    top_otus = (
        top_pathway_row[otu_cols]
        .sort_values(ascending=False)
        .head(20)
    )

    print("PathwayID:", top_pathway_id)
    print("Pathway:", top_pathway_row["Pathway"])
    print(top_otus.to_string())


