from pathlib import Path 
import re
import pandas as pd
from skbio import TreeNode

RAW_FILE = Path("Data/pathways_species2.tsv")
GTDB_TAXONOMY = Path("Data/gtdb/bac120_taxonomy.tsv.gz")
GTDB_TREE = Path("Data/gtdb/bac120.tree.gz")
GTDB_METADATA = Path("Data/gtdb/bac120_metadata.tsv.gz")

OUT_MAPPING_ALL = Path("Data/gtdb_mapping_all.tsv")
OUT_MAPPED_EXACT = Path("Data/mapped_exact.tsv")
OUT_MAPPED_GENUS_ONLY = Path("Data/mapped_genus_only.tsv")
OUT_UNMAPPED = Path("Data/unmapped.tsv")
OUT_REVIEW = Path("Data/gtdb_mapping_review.tsv")


BACTERIA_COL = 1

def read_raw_no_header(path):
    return pd.read_csv(path, sep=r"\s+", header=None, dtype=str, engine="python")

def parse_original_taxon(taxon):
    original_taxon = taxon

    issues = []

    if taxon is None:
        return None, {
            "original_taxon": original_taxon,
            "issue": "missing_value",
            "details": "Taxon is None"
        }

    taxon = str(taxon)

    if taxon != taxon.strip():
        issues.append("leading_or_trailing_spaces")

    taxon_clean = taxon.strip()

    if re.search(r"\s", taxon_clean):
        issues.append("internal_whitespace")

    genus_matches = re.findall(r"g__([^.;\s]+)", taxon_clean)
    species_matches = re.findall(r"s__([^.;\s]+)", taxon_clean)

    if len(genus_matches) == 0:
        return None, {
            "original_taxon": original_taxon,
            "issue": "missing_genus",
            "details": "No g__ pattern found"
        }

    if len(genus_matches) > 1:
        issues.append("multiple_genus_patterns")

    if len(species_matches) > 1:
        issues.append("multiple_species_patterns")

    genus = genus_matches[0]

    if len(species_matches) == 0:
        parsed = {
            "original_taxon": taxon_clean,
            "genus": genus,
            "species_name": None,
            "gtdb_genus": f"g__{genus}",
            "gtdb_species": None,
        }

        return parsed, {
            "original_taxon": original_taxon,
            "issue": "missing_species",
            "details": "No s__ pattern found"
        }

    species_raw = species_matches[0]

    if species_raw == "":
        return None, {
            "original_taxon": original_taxon,
            "issue": "empty_species",
            "details": "s__ was found but species name is empty"
        }

    species_name = species_raw.replace("_", " ")

    if not species_name.startswith(genus):
        issues.append("species_does_not_start_with_genus")
        species_name = f"{genus} {species_name}"

    parsed = {
        "original_taxon": taxon_clean,
        "genus": genus,
        "species_name": species_name,
        "gtdb_genus": f"g__{genus}",
        "gtdb_species": f"s__{species_name}",
    }

    if len(issues) > 0:
        return parsed, {
            "original_taxon": original_taxon,
            "issue": ";".join(issues),
            "details": f"Parsed but suspicious: {taxon_clean}"
        }

    return parsed, None

# #

def get_rank(taxonomy, prefix):
    for part in str(taxonomy).split(";"):
        part = part.strip()
        if part.startswith(prefix):
            return part
    return None


def main():
    raw = read_raw_no_header(RAW_FILE)

    taxa = (
        raw.iloc[:, BACTERIA_COL]
        .dropna()
        .astype(str)
        .str.strip()
    )

    taxa = sorted(taxa[taxa.str.contains("g__", regex=False)].unique())

    wanted = [] #taxa parseables
    weird_taxa = [] #taxa con problemas o advertencias

    for taxon in taxa:
        parsed, issue = parse_original_taxon(taxon)

        if parsed is not None:
            wanted.append(parsed)

        if issue is not None:
            weird_taxa.append(issue)

    wanted     = pd.DataFrame(wanted)
    weird_taxa = pd.DataFrame(weird_taxa)

    wanted.to_csv("Data/parsed_taxa.tsv", sep="\t", index=False)
    weird_taxa.to_csv("Data/weird_taxa.tsv", sep="\t", index=False)

    # print(f"Taxa parseados: {len(wanted)}")
    # print(f"Taxa raros o problemáticos: {len(weird_taxa)}")

    tree      = TreeNode.read(str(GTDB_TREE))
    tree_tips = {tip.name for tip in tree.tips()}

    def normalize_accession(x):
        return str(x).replace("_", " ").strip()

    tree_tips = {
        normalize_accession(x)
        for x in tree_tips
    }

    # print("\nÁrbol")
    # print("tree.count():", tree.count())
    # print("len(tree_tips):", len(tree_tips))

    gtdb = pd.read_csv(
        GTDB_TAXONOMY,
        sep="\t",
        header=None,
        names=["gtdb_accession", "taxonomy"],
        dtype=str,
    )


# META ----------------------------------------------------------------
    meta = pd.read_csv(
        GTDB_METADATA,
        sep="\t",
        dtype=str,
    )
    meta["ncbi_species"] = meta["ncbi_taxonomy"].apply(
        lambda x: get_rank(x, "s__")
    )    
    meta.loc[meta["ncbi_species"] == "s__", "ncbi_species"] = None

    meta["accession_norm"] = (
        meta["accession"]
        .astype(str)
        .str.replace("_", " ", regex=False)
        .str.strip()
    )

    meta = meta[meta["accession_norm"].isin(tree_tips)].copy()

    meta["metadata_gtdb_genus"] = meta["gtdb_taxonomy"].apply(lambda x: get_rank(x, "g__"))
    meta["metadata_gtdb_species"] = meta["gtdb_taxonomy"].apply(lambda x: get_rank(x, "s__"))

    meta["ncbi_organism_species"] = (
        "s__" + meta["ncbi_organism_name"].fillna("").str.strip()
    )

    # print("\nMetadata columns:")
    # print(meta.columns.tolist())
    # print("\nPrimeras especies NCBI")
    # print(meta["ncbi_species"].dropna().head(20).tolist())
    # print("\nPrimeras taxonomías NCBI")
    # print(meta["ncbi_taxonomy"].dropna().head(10).tolist())


# META ----------------------------------------------------------------

    gtdb["accession_norm"] = (
        gtdb["gtdb_accession"]
        .astype(str)
        .str.replace("_", " ", regex=False)
        .str.strip()
    )

    # print("\nPrueba de intersección")

    taxonomy_acc = set(gtdb["accession_norm"])

    # print("taxonomy_acc:", len(taxonomy_acc))
    # print("tree_tips:", len(tree_tips))

    inter = taxonomy_acc.intersection(tree_tips)

    # print("intersección:", len(inter))

    # if len(inter) > 0:
    #     print("ejemplos:", list(inter)[:10])

    # print("\nPrimer accession normalizado")
    # print(repr(gtdb["accession_norm"].iloc[0]))

    # print("\nPrimer tip")
    first_tip = next(iter(tree_tips))
    # print(repr(first_tip))
    # print(type(first_tip))

    # El árbol de referencia sólo contiene representantes.
    # Por eso filtramos la taxonomía a accessions presentes en el árbol.
    #gtdb = gtdb[gtdb["gtdb_accession"].isin(tree_tips)].copy()
    gtdb = gtdb[gtdb["accession_norm"].isin(tree_tips)].copy()

    gtdb["gtdb_domain"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "d__"))
    gtdb["gtdb_phylum"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "p__"))
    gtdb["gtdb_class"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "c__"))
    gtdb["gtdb_order"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "o__"))
    gtdb["gtdb_family"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "f__"))
    gtdb["gtdb_genus"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "g__"))
    gtdb["gtdb_species"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "s__"))

    target = "g__Bacteroides"

    # print(
    #     gtdb[
    #         gtdb["gtdb_genus"] == target
    #     ][["gtdb_species"]]
    #     .head(50)
    # )

    # print("\nGTDB shape")
    # print("GTDB después de filtrar:", gtdb.shape)
    # print("\nPrimeras taxonomías")
    # print(gtdb["taxonomy"].head())
    # print("\nPrimeros géneros")
    # print(gtdb["gtdb_genus"].dropna().head(20).tolist())
    # print("\nPrimeras especies")
    # print(gtdb["gtdb_species"].dropna().head(20).tolist())
    # print("\nNúmero de géneros únicos")
    # print(gtdb["gtdb_genus"].nunique())
    # print("\nNúmero de especies únicas")
    # print(gtdb["gtdb_species"].nunique())

    accepted = []
    review = []

    mapping = []
    review = []

    for _, row in wanted.iterrows():                      #---------- LOOP
        original = row["original_taxon"]
        target_species = row["gtdb_species"]
        target_genus = row["gtdb_genus"]

        base = {
            "original_taxon": original,
            "input_genus": row["gtdb_genus"],
            "input_species": row["gtdb_species"],
            "gtdb_accession": None,
            "match_status": None,
            "gtdb_species": None,
            "gtdb_genus": None,
            "gtdb_family": None,
            "gtdb_order": None,
            "gtdb_class": None,
            "gtdb_phylum": None,
            "gtdb_taxonomy": None,
        }

        if pd.isna(target_species) or target_species is None:
            base["match_status"] = "no_species_in_original_name"
            mapping.append(base)

            review.append({
                "original_taxon": original,
                "reason": "no_species_in_original_name",
                "candidate_accession": None,
                "candidate_species": None,
                "candidate_taxonomy": None,
            })
            continue

        # 1. GTDB species exact
        exact = gtdb[gtdb["gtdb_species"] == target_species]

        if len(exact) == 1:
            hit = exact.iloc[0]

            base.update({
                "gtdb_accession": hit["gtdb_accession"],
                "match_status": "exact_species_match",
                "gtdb_species": hit["gtdb_species"],
                "gtdb_genus": hit["gtdb_genus"],
                "gtdb_family": hit["gtdb_family"],
                "gtdb_order": hit["gtdb_order"],
                "gtdb_class": hit["gtdb_class"],
                "gtdb_phylum": hit["gtdb_phylum"],
                "gtdb_taxonomy": hit["taxonomy"],
            })

            mapping.append(base)
            continue


        # 2. GTDB species prefix
        prefix = gtdb[
            gtdb["gtdb_species"].fillna("").str.startswith(target_species + "_")
        ]

        if len(prefix) == 1:
            hit = prefix.iloc[0]

            base.update({
                "gtdb_accession": hit["gtdb_accession"],
                "match_status": "unique_species_prefix",
                "gtdb_species": hit["gtdb_species"],
                "gtdb_genus": hit["gtdb_genus"],
                "gtdb_family": hit["gtdb_family"],
                "gtdb_order": hit["gtdb_order"],
                "gtdb_class": hit["gtdb_class"],
                "gtdb_phylum": hit["gtdb_phylum"],
                "gtdb_taxonomy": hit["taxonomy"],
            })

            mapping.append(base)
            continue


        # 3. NCBI organism name exact
        ncbi_org = meta[meta["ncbi_organism_species"] == target_species]

        if len(ncbi_org) == 1:
            hit = ncbi_org.iloc[0]

            base.update({
                "gtdb_accession": hit["accession"],
                "match_status": "ncbi_organism_name_match",
                "gtdb_species": hit["metadata_gtdb_species"],
                "gtdb_genus": hit["metadata_gtdb_genus"],
                "gtdb_family": get_rank(hit["gtdb_taxonomy"], "f__"),
                "gtdb_order": get_rank(hit["gtdb_taxonomy"], "o__"),
                "gtdb_class": get_rank(hit["gtdb_taxonomy"], "c__"),
                "gtdb_phylum": get_rank(hit["gtdb_taxonomy"], "p__"),
                "gtdb_taxonomy": hit["gtdb_taxonomy"],
            })

            mapping.append(base)
            continue

        elif len(ncbi_org) > 1:
            review.append({
                "original_taxon": original,
                "reason": "ambiguous_ncbi_organism_name_match",
                "candidate_accession": None,
                "candidate_species": target_species,
                "candidate_taxonomy": None,
            })


        # 4. NCBI taxonomy species exact
        ncbi_tax = meta[meta["ncbi_species"] == target_species]

        if len(ncbi_tax) == 1:
            hit = ncbi_tax.iloc[0]

            base.update({
                "gtdb_accession": hit["accession"],
                "match_status": "ncbi_taxonomy_species_match",
                "gtdb_species": hit["metadata_gtdb_species"],
                "gtdb_genus": hit["metadata_gtdb_genus"],
                "gtdb_family": get_rank(hit["gtdb_taxonomy"], "f__"),
                "gtdb_order": get_rank(hit["gtdb_taxonomy"], "o__"),
                "gtdb_class": get_rank(hit["gtdb_taxonomy"], "c__"),
                "gtdb_phylum": get_rank(hit["gtdb_taxonomy"], "p__"),
                "gtdb_taxonomy": hit["gtdb_taxonomy"],
            })

            mapping.append(base)
            continue

        elif len(ncbi_tax) > 1:
            review.append({
                "original_taxon": original,
                "reason": "ambiguous_ncbi_taxonomy_species_match",
                "candidate_accession": None,
                "candidate_species": target_species,
                "candidate_taxonomy": None,
            })


        # 5. Genus only
        genus_hits = gtdb[gtdb["gtdb_genus"] == target_genus].copy()

        if len(genus_hits) > 0:
            base["match_status"] = "genus_match_only"
            mapping.append(base)

            for _, cand in genus_hits.head(20).iterrows():
                review.append({
                    "original_taxon": original,
                    "reason": "genus_match_only",
                    "candidate_accession": cand["gtdb_accession"],
                    "candidate_species": cand["gtdb_species"],
                    "candidate_taxonomy": cand["taxonomy"],
                })

            continue

        # 6. Unmapped
        base["match_status"] = "no_species_or_genus_match"
        mapping.append(base)

        review.append({
            "original_taxon": original,
            "reason": "no_species_or_genus_match",
            "candidate_accession": None,
            "candidate_species": None,
            "candidate_taxonomy": None,
        })



# ---------------------------------------------------------------------------


    mapping = pd.DataFrame(mapping)
    review = pd.DataFrame(review)

    mapped_exact = mapping[
        mapping["match_status"].isin([
            "exact_species_match",
            "unique_species_prefix",
            "ncbi_organism_name_match",
            "ncbi_taxonomy_species_match",
        ])
    ].copy()

    mapped_genus_only = mapping[
        mapping["match_status"] == "genus_match_only"
    ].copy()

    unmapped = mapping[
        mapping["match_status"].isin([
            "no_species_or_genus_match",
            "no_species_in_original_name",
        ])
    ].copy()

    mapping.to_csv(OUT_MAPPING_ALL, sep="\t", index=False)
    mapped_exact.to_csv(OUT_MAPPED_EXACT, sep="\t", index=False)
    mapped_genus_only.to_csv(OUT_MAPPED_GENUS_ONLY, sep="\t", index=False)
    unmapped.to_csv(OUT_UNMAPPED, sep="\t", index=False)
    review.to_csv(OUT_REVIEW, sep="\t", index=False)

    print(f"Mapeos exactos/especie: {len(mapped_exact)}")
    print(f"Mapeos solo por género: {len(mapped_genus_only)}")
    print(f"No mapeados: {len(unmapped)}")
    print(f"Para revisión manual: {review['original_taxon'].nunique() if len(review) else 0}")

    print(f"Archivo completo: {OUT_MAPPING_ALL}")
    print(f"Archivo exactos: {OUT_MAPPED_EXACT}")
    print(f"Archivo solo género: {OUT_MAPPED_GENUS_ONLY}")
    print(f"Archivo no mapeados: {OUT_UNMAPPED}")
    print(f"Archivo revisión: {OUT_REVIEW}")


# SANITY -----------------------------------------------------------------------

    mapped = pd.read_csv("Data/mapped_exact.tsv", sep="\t")
    print("Filas:", len(mapped))
    print("Taxa originales:", mapped["original_taxon"].nunique())
    print("Accessions únicos:", mapped["gtdb_accession"].nunique())
    print("Especies GTDB únicas:", mapped["gtdb_species"].nunique())
    # deberia obtener:
    # Filas: 275 Taxa originales: 275 Accessions únicos: 275 Especies GTDB únicas: 275

    tips = set(tree_tips)
    mapped["accession_norm"] = (
        mapped["gtdb_accession"]
        .str.replace("_", " ")
    )
    present = mapped["accession_norm"].isin(tips)
    print(present.sum(), "/", len(mapped))
    # debo obtener 275/275
# SANITY -----------------------------------------------------------------------



# PRUNE TREE --------------------------------------------------------------------
    # pruned_tree = tree.shear(accessions_275)

    # save
    # Data/tree_275.nwk

    

if __name__ == "__main__":
    main()