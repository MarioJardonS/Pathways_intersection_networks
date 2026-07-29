# Scripts/03_compute_unifrac_for_network_pairs.py

# ---------------------------------------------------------------------
# Objetivo:
# Calcular las distancias UniFrac (weighted y unweighted) entre pares de
# pathways conectados en la red de cada muestra.
#
# Entradas:
# TABLE_FILE: Tabla de abundancias donde cada fila corresponde a una
# combinación PathwayID x SampleID. Las tres primeras columnas contienen
# PathwayID, Pathway y SampleID; las columnas restantes corresponden a
# los taxones de GTDB y sus abundancias.
#
# TREE_FILE: Árbol filogenético de GTDB podado a los taxones presentes
# en TABLE_FILE.
#
# METADATA_FILE: Tabla de metadatos con la información de cada muestra
# (SampleID, Diagnosis, Project, GMHI, hiPCA y Shannon_entropy).
#
# NETWORKS_DIR: Directorio con una red por muestra (<SampleID>.tsv).
# Las dos primeras columnas contienen los pares de pathways conectados.
#
# Salida:
# MASTER_OUT_FILE: Tabla maestra con una fila por combinación
# SampleID x par de pathways, incluyendo las distancias UniFrac,
# métricas de abundancia y los metadatos de la muestra.

# NOTA: Cada pathway se representa como un vector de abundancias sobre los mismos taxones del árbol
#                Taxón1  Taxón2  Taxón3  Taxón4 ...
# PWY1 Sample1     0.10    0.30    0.00    0.05
# PWY2 Sample1     0.00    0.15    0.20    0.01

# y ambos vectores utilizan exactamente:
# el mismo conjunto de taxones (otu_cols);
# el mismo árbol (TREE_FILE).


# Mildred SM Julio 2026
# ---------------------------------------------------------------------

from pathlib import Path
import numpy as np
import pandas as pd
from skbio import TreeNode
from skbio.diversity import beta_diversity

# tabla con abundancias de pathways por muestra y taxón
TABLE_FILE    = Path("Data/unifrac/pathway_sample_unifrac_abundance_table.tsv") 
TREE_FILE     = Path("Data/unifrac/pathway_sample_tree_pruned.nwk") # árbol filogenético de GTDB podado
METADATA_FILE = Path("Data/metadata.csv") # tabla de metadatos
NETWORKS_DIR  = Path("Networks") # directorio con archivos de red por muestra

OUT_DIR = Path("Data/data_samples/network_unifrac") # directorio de salida para los resultados
MASTER_OUT_FILE = OUT_DIR / "network_pathway_pair_unifrac_master.tsv" # archivo maestro de salida con todas las combinaciones SampleID x par de pathways
OUT_DIR.mkdir(parents=True, exist_ok=True) # crear directorio de salida si no existe

N_PER_GROUP = None #5 # Número de muestras por grupo a seleccionar (None para usar todas las muestras disponibles)
GROUPS_TO_USE = ["Obese", "Healthy"] # Grupos de diagnóstico a incluir en el análisis

# ---------------------------------------------------------------------
# funcion para extraer el ID del pathway de la columna PathwayID (formato: "PathwayID:taxon")
def pathway_id(x):
    return str(x).split(":", 1)[0].strip()

# funcion que calcula la distancia UniFrac entre dos arrays de abundancias (array_x y array_y) para un par de pathways
def compute_pairwise_unifrac(array_x, array_y, taxa, tree, metric):
    counts = np.vstack([
        np.asarray(array_x, dtype=float),
        np.asarray(array_y, dtype=float)
    ])

    # esta función permite puede calcular Bray-Curtis, 
    # Jaccard, UniFrac, etc. El tipo de distancia depende del argumento metric
    # counts es la matriz de abundancia, cada fila representa un pathway y cada columna un taxón,
    # tree es el árbol filogenético de los taxones
    # taxa es la lista de taxones correspondientes a las columnas, tiene el mismo orden que las columnas en counts
    # devuelve un objeto DistanceMatrix con las distancias entre los pathways (filas) y sus IDs (x, y)
    # Si el input tiene más filas, calcula todas las distancias posibles entre todos los pares 

    dm = beta_diversity(
        metric=metric,
        counts=counts,
        ids=["x", "y"],
        taxa=taxa,
        tree=tree,
    )
    return dm["x", "y"]

# funcion que calcula las distancias UniFrac para todos los pares de pathways conectados en la red de una muestra
def compute_unifrac_for_network_file(network_file, diagnosis, table, tree, otu_cols):
    sample_id    = Path(network_file).stem
    sample_table = table[table["SampleID"] == sample_id].copy()
    if sample_table.empty:
        print(f"Skipping {sample_id}: sample not found in abundance table.")
        return None
    
    # sample_table es una tabla con las abundancias de pathways para la muestra actual, filtrada por SampleID
    sample_table       = sample_table.set_index("PathwayID", drop=False)

    # available_pathways es un conjunto de PathwayID presentes en la tabla de abundancias para la muestra actual
    available_pathways = set(sample_table["PathwayID"])

    # Lee la network de la muestra actual y filtra los pares de pathways conectados
    net = pd.read_csv(network_file, sep="\t")

    # Filtra los pares de pathways que no sean "UNINTEGRATED" o "UNMAPPED" y que tengan PathwayID diferentes
    pair_cols = net.columns[:2]
    pairs = net[list(pair_cols)].copy()
    pairs.columns = ["PathwayA", "PathwayB"]
    pairs = pairs[
        ~pairs["PathwayA"].astype(str).isin(["UNINTEGRATED", "UNMAPPED"]) &
        ~pairs["PathwayB"].astype(str).isin(["UNINTEGRATED", "UNMAPPED"])
    ].copy()

    # Agrega columnas con los PathwayID extraídos de las columnas PathwayA y PathwayB
    pairs["PathwayA_ID"] = pairs["PathwayA"].map(pathway_id)
    pairs["PathwayB_ID"] = pairs["PathwayB"].map(pathway_id)

    # Filtra los pares de pathways que tengan PathwayID diferentes
    pairs = pairs[pairs["PathwayA_ID"] != pairs["PathwayB_ID"]].copy()

    # Crea una columna "pair_key" que combina los PathwayID de cada par 
    # de pathways en orden alfabético, para evitar duplicados
    pairs["pair_key"] = pairs.apply(
        lambda r: "__".join(sorted([r["PathwayA_ID"], r["PathwayB_ID"]])),
        axis=1
    )
    # Elimina los pares duplicados basados en la columna "pair_key"
    pairs = pairs.drop_duplicates("pair_key").copy()

    results = []

    # Itera sobre cada par de pathways conectados en la network de la muestra actual
    for _, row in pairs.iterrows():
        pathway_a_id = row["PathwayA_ID"]
        pathway_b_id = row["PathwayB_ID"]

        pathway_a_available = pathway_a_id in available_pathways
        pathway_b_available = pathway_b_id in available_pathways

        weighted   = np.nan
        unweighted = np.nan

        pathway_a_total   = np.nan
        pathway_b_total   = np.nan
        pathway_a_nonzero = np.nan
        pathway_b_nonzero = np.nan
        shared_otus       = np.nan
        union_otus        = np.nan

        pathway_a_name = row["PathwayA"]
        pathway_b_name = row["PathwayB"]

        # Si ambos pathways están disponibles en la tabla de abundancias, 
        # calcula las distancias UniFrac y las métricas de abundancia
        if pathway_a_available and pathway_b_available:
            row_a = sample_table.loc[pathway_a_id]
            row_b = sample_table.loc[pathway_b_id]

            array_a = row_a[otu_cols].values
            array_b = row_b[otu_cols].values

            pathway_a_name = row_a["Pathway"]
            pathway_b_name = row_b["Pathway"]

            # Las ramas del árbol reciben un peso proporcional a la abundancia.
            # responde qué tan diferente es la composición filogenética 
            # de los microorganismos que participan en dos pathways distintos
            weighted = compute_pairwise_unifrac(
                array_a, array_b,
                taxa=otu_cols,
                tree=tree,
                metric="weighted_unifrac"
            )

            # calcula qué fracción de la longitud del árbol pertenece a linajes
            # presentes en un pathway pero ausentes en el otro
            unweighted = compute_pairwise_unifrac(
                array_a, array_b,
                taxa=otu_cols,
                tree=tree,
                metric="unweighted_unifrac"
            )

            pathway_a_total   = np.sum(array_a)
            pathway_b_total   = np.sum(array_b)
            pathway_a_nonzero = np.sum(array_a > 0)
            pathway_b_nonzero = np.sum(array_b > 0)
            shared_otus       = np.sum((array_a > 0) & (array_b > 0))
            union_otus        = np.sum((array_a > 0) | (array_b > 0))

        results.append({
            "ResultKey":                f"{sample_id}__{row['pair_key']}",
            "SampleID":                 sample_id,
            "Diagnosis":                diagnosis,
            "PathwayA_ID":              pathway_a_id,
            "PathwayB_ID":              pathway_b_id,
            "PairKey":                  row["pair_key"],
            "PathwayA":                 pathway_a_name,
            "PathwayB":                 pathway_b_name,
            "Weighted_UniFrac":         weighted,
            "Unweighted_UniFrac":       unweighted,
            "PathwayA_available":       pathway_a_available,
            "PathwayB_available":       pathway_b_available,
            "PathwayA_total_abundance": pathway_a_total,
            "PathwayB_total_abundance": pathway_b_total,
            "PathwayA_nonzero_otus":    pathway_a_nonzero,
            "PathwayB_nonzero_otus":    pathway_b_nonzero,
            "Shared_nonzero_otus":      shared_otus,
            "Union_nonzero_otus": union_otus,
        })

    results = pd.DataFrame(results)

    # print("\nSample:", sample_id, diagnosis)
    # print("Network pairs:", len(pairs))
    # print("Computable pairs:", results["Unweighted_UniFrac"].notna().sum())
    # print("Pairs with NaN:", results["Unweighted_UniFrac"].isna().sum())

    return results




# ---------------------------------------------------------------------
# Main
table = pd.read_csv(TABLE_FILE, sep="\t")
tree  = TreeNode.read(str(TREE_FILE))
meta  = pd.read_csv(METADATA_FILE, sep=",")

# lista de los nombres de las columnas que contienen 
# las abundancias de los taxones en la tabla de abundancias
otu_cols = table.columns[3:].tolist()

# conjunto de los SampleID presentes en la tabla de abundancias
table_samples = set(table["SampleID"].astype(str))

# crea lista de los archivos de las networks en el directorio 
# NETWORKS_DIR, excluyendo "edges.tsv" y "network.tsv"
network_files = sorted(NETWORKS_DIR.glob("*.tsv"))
network_files = [
    p for p in network_files
    if p.name not in ["edges.tsv", "network.tsv"]
]
# conjunto de los SampleID presentes en los archivos de network
network_sample_ids = set(p.stem for p in network_files)

# Filtra la tabla de metadatos para incluir solo las muestras 
# que están presentes en la tabla de abundancias, 
# en los archivos de network y que pertenecen a los grupos de diagnóstico especificados
meta["SampleID"] = meta["SampleID"].astype(str)
valid_meta = meta[
    meta["SampleID"].isin(table_samples) &
    meta["SampleID"].isin(network_sample_ids) &
    meta["Diagnosis"].isin(GROUPS_TO_USE)
].copy()

selected_meta = valid_meta.copy()

if N_PER_GROUP is not None:
    if not isinstance(N_PER_GROUP, int) or N_PER_GROUP <= 0:
        raise ValueError("N_PER_GROUP must be None or a positive integer.")

    selected_meta = (
        selected_meta
        .groupby("Diagnosis", group_keys=False)
        .head(N_PER_GROUP)
        .copy()
    )

# print("\nSelected samples:")
# print(selected_meta[["SampleID", "Diagnosis"]].to_string(index=False))

# Crea un diccionario que mapea cada SampleID a su Diagnosis correspondiente
sample_to_diagnosis = dict(zip(selected_meta["SampleID"], selected_meta["Diagnosis"]))

# Filtra la lista de archivos de network para incluir solo aquellos que 
# corresponden a los SampleID seleccionados
selected_network_files = [
    p for p in network_files
    if p.stem in sample_to_diagnosis
]

all_results = []

# Itera sobre cada archivo de network seleccionado 
# y calcula las distancias UniFrac para los pares de pathways conectados
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



# ---------------------------------------------------------------------
# crea tabla de resultados maestra con todas las combinaciones SampleID x par de pathways
if MASTER_OUT_FILE.exists():
    old_results = pd.read_csv(MASTER_OUT_FILE, sep="\t")

    combined = pd.concat(
        [old_results, all_results],
        ignore_index=True
    )

    combined = (
        combined
        .drop_duplicates(subset=["ResultKey"], keep="last")
        .sort_values(["Diagnosis", "SampleID", "PairKey"])
        .reset_index(drop=True)
    )
else:
    combined = (
        all_results
        .sort_values(["Diagnosis", "SampleID", "PairKey"])
        .reset_index(drop=True)
    )

combined.to_csv(MASTER_OUT_FILE, sep="\t", index=False)

# print("\nMaster output saved/updated:")
# print(MASTER_OUT_FILE)

# print("\nTotal rows in master table:")
# print(len(combined))

# print("\nSamples per diagnosis in master table:")
# print(
#     combined[["SampleID", "Diagnosis"]]
#     .drop_duplicates()["Diagnosis"]
#     .value_counts()
# )