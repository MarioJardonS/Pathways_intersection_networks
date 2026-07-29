# El script toma nombres bacterianos de pathways_species2.tsv 
# e intenta asociar cada taxón con un genoma representante 
# del árbol filogenético de GTDB.

# El procedimiento sigue una estrategia progresiva:
# 1. Buscar una coincidencia exacta con la especie GTDB.
# 2. Coincidencia con una variante nomenclatural de la especie GTDB
#    (sufijos como _A, _B, etc.).
# 3. Coincidencia usando el nombre del organismo en NCBI.
# 4. Coincidencia usando la especie de la taxonomía NCBI.
# 5. Si no existe una coincidencia de especie, recuperar candidatos del
#    mismo género para revisión manual.
# 6. Si no existe ninguna coincidencia, registrar el taxón como no mapeado.

# Después genera archivos separados para:
# mapeos a nivel de especie;
# coincidencias solo por género;
# taxones no mapeados;
# casos que requieren revisión manual

# Mildred S-M. Version Julio 2026
# ---------------------------------------------------------------------------

from curses import raw
from pathlib import Path 
import re
import pandas as pd
from skbio import TreeNode

RAW_FILE      = Path("Data/pathways_species2.tsv")

# estos archivos contienen la taxonomía y el árbol de 
# referencia de GTDB, así como metadatos adicionales
GTDB_TAXONOMY = Path("Data/gtdb/bac120_taxonomy.tsv.gz")
GTDB_TREE     = Path("Data/gtdb/bac120.tree.gz")
GTDB_METADATA = Path("Data/gtdb/bac120_metadata.tsv.gz")

# 
OUT_MAPPING_ALL       = Path("Data/gtdb_mapping_all.tsv") # mapeos completos, incluyendo coincidencias exactas, solo género y no mapeados
OUT_MAPPED_EXACT      = Path("Data/mapped_exact.tsv") # coincidencias exactas
OUT_MAPPED_GENUS_ONLY = Path("Data/mapped_genus_only.tsv") # coincidencias solo por género
OUT_UNMAPPED          = Path("Data/unmapped.tsv") # taxones no mapeados
OUT_REVIEW            = Path("Data/gtdb_mapping_review.tsv") # casos que requieren revisión manual

BACTERIA_COL = 1 # identificador de la columna que contiene los nombres bacterianos en el archivo pathways_species2.tsv

# leemos el archivo de entrada sin encabezado, usando espacios como separadores
def read_raw_no_header(path):
    return pd.read_csv(path, sep=r"\s+", header=0, dtype=str, engine="python")

# esta función intenta extraer el género y la especie de un nombre con formato: 
# g__Bacteroides.s__Bacteroides_fragilis
# devuelve: 
# parsed: información taxonómica extraída.
# issue: advertencia o problema detectado.
def parse_original_taxon(taxon):
    original_taxon = taxon

    issues = []

    # lee el taxón y verifica si es nulo
    if taxon is None:
        return None, {
            "original_taxon": original_taxon,
            "issue": "missing_value",
            "details": "Taxon is None"
        }

    taxon = str(taxon)

    # verifica si hay espacios al inicio o al final del taxón
    if taxon != taxon.strip():
        issues.append("leading_or_trailing_spaces")

    # verifica si hay espacios internos en el taxón
    taxon_clean = taxon.strip()

    if re.search(r"\s", taxon_clean):
        issues.append("internal_whitespace")

    # si tiene espacios al inicio o al final, se registra:
    # leading_or_trailing_spaces
    # despues se eliminan estos espacios

    # verifica si hay caracteres no alfanuméricos en el taxón
    genus_matches = re.findall(r"g__([^.;\s]+)", taxon_clean)
    species_matches = re.findall(r"s__([^.;\s]+)", taxon_clean)

    # taxon sin genero : sin un patrón g__, el taxón no se incorpora a la lista de taxones que serán mapeados
    if len(genus_matches) == 0:
        return None, {
            "original_taxon": original_taxon,
            "issue": "missing_genus",
            "details": "No g__ pattern found"
        }

    # si aparecen varios géneros o especies, se registra una advertencia
    # más adelante solo se usa el primer resultado los demás se ignoran
    if len(genus_matches) > 1:
        issues.append("multiple_genus_patterns")
    if len(species_matches) > 1:
        issues.append("multiple_species_patterns")

    genus = genus_matches[0]

    # axon con género pero sin especie
    # el taxon se conserva para mapeo pero no tiene especie
    # se reporta como "missing_species"
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

    # lo adaptamos para que el nombre de la especie tenga espacios en lugar de guiones bajos
    # y verificamos si el nombre de la especie comienza con el género
    species_raw = species_matches[0]

    if species_raw == "":
        return None, {
            "original_taxon": original_taxon,
            "issue": "empty_species",
            "details": "s__ was found but species name is empty"
        }

    species_name = species_raw.replace("_", " ")

    # verificacion del genero dentro del nombre de la especie 
    if not species_name.startswith(genus):
        issues.append("species_does_not_start_with_genus")
        species_name = f"{genus} {species_name}"

    # construye el resultado
    # ejemplo:
    # g__Bacteroides.s__Bacteroides_fragilis
    # se convierte en 
    # {
    # "original_taxon": "g__Bacteroides.s__Bacteroides_fragilis",
    # "genus": "Bacteroides",
    # "species_name": "Bacteroides fragilis",
    # "gtdb_genus": "g__Bacteroides",
    # "gtdb_species": "s__Bacteroides fragilis"
    # }
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








# Esta función recibe una taxonomía completa como cadena 
# de texto y un prefijo (por ejemplo, "g__" para género o "s__" para especie) 
# y devuelve la parte de la taxonomía que coincide con ese prefijo. 
# Si no se encuentra ninguna coincidencia, devuelve None.
def get_rank(taxonomy, prefix):
    for part in str(taxonomy).split(";"):
        part = part.strip()
        if part.startswith(prefix):
            return part
    return None


# funcion que contiene el flujo princial del script, 
# que realiza la lectura de datos, el mapeo de taxones y la generación de archivos de salida.
def main():
    # carga archivo
    raw = read_raw_no_header(RAW_FILE) 

    # selecciona columna de taxones bacterianos, elimina valores nulos, convierte a cadena 
    # y elimina espacios al inicio o al final
    taxa = (
    raw["OTU"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    # filtra y ordena los taxones que contienen el prefijo "g__" (género)
    taxa = sorted(taxa[taxa.str.contains("g__", regex=False)].unique())

    wanted     = [] # almacenara los taxones parseados correctamente
    weird_taxa = [] # almacenara los taxones que presentan problemas o advertencias durante el parseo

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


    # CARGA DE DATOS ----------------------------------------------------------------
    # Carga del árbol filogenético de GTDB y extracción de los nombres de las hojas 
    # (tips) del árbol. Se normalizan los nombres de los tips para eliminar guiones bajos
    tree      = TreeNode.read(str(GTDB_TREE))
    tree_tips = {tip.name for tip in tree.tips()} # extrae los nombres de las hojas del árbol en un set

    # convierte accesion a formato normalizado (sin guiones bajos y sin espacios al inicio o al final)
    # los hace compatibles con los de las tablas (GB_GCA_000123 -> GB GCA 000123)
    def normalize_accession(x):
        return str(x).replace("_", " ").strip()

    tree_tips = {
        normalize_accession(x)
        for x in tree_tips
    }

    # lectura de la taxonomía de GTDB, que contiene los accesiones y sus taxonomías correspondientes
    gtdb = pd.read_csv(
        GTDB_TAXONOMY,
        sep="\t",
        header=None,
        names=["gtdb_accession", "taxonomy"],
        dtype=str,
    )


# META ----------------------------------------------------------------
# lectura de los metadatos de GTDB, que incluyen información adicional 
# como el nombre del organismo en NCBI y la taxonomía NCBI.
    meta = pd.read_csv(
        GTDB_METADATA,
        sep="\t",
        dtype=str,
    )
    meta["ncbi_species"] = meta["ncbi_taxonomy"].apply(
        lambda x: get_rank(x, "s__")
    )    
    meta.loc[meta["ncbi_species"] == "s__", "ncbi_species"] = None

    # se aplica la misma transformacion usada con los tips del arbol
    meta["accession_norm"] = (
        meta["accession"]
        .astype(str)
        .str.replace("_", " ", regex=False)
        .str.strip()
    )
    # de la taxonomia NCBI completa se extrae el rango de especie
    # si la taxonomia solo contiene el prefijo vacio s__, se considera que la especie esta ausente

    meta = meta[meta["accession_norm"].isin(tree_tips)].copy()
    # se eliminan los genomas que no estan presentes como hojas del arbola

    # estas columnas indican la clasificacion GTDB del genome accession
    meta["metadata_gtdb_genus"] = meta["gtdb_taxonomy"].apply(lambda x: get_rank(x, "g__"))
    meta["metadata_gtdb_species"] = meta["gtdb_taxonomy"].apply(lambda x: get_rank(x, "s__"))

    # construccion del  nombre de especie NCBI a partir del nombre del organismo NCBI, 
    # agregando el prefijo s__ y eliminando espacios al inicio o al final
    # ojo: algunos nombres de organismos NCBI no contienen el nombre de la especie, 
    # por lo que esta columna puede contener valores nulos
    meta["ncbi_organism_species"] = (
        "s__" + meta["ncbi_organism_name"].fillna("").str.strip()
    )

    # -----------------------------------------------------------------------------------
    # preparacion de la tabla de taxones GTDB para el mapeo, 
    # normalizando los accesiones y filtrando solo aquellos presentes en el arbol filogenetico
    gtdb["accession_norm"] = (
        gtdb["gtdb_accession"]
        .astype(str)
        .str.replace("_", " ", regex=False)
        .str.strip()
    )

    # se calcula la intersección entre los accesiones de la taxonomía GTDB y los tips del árbol filogenético
    # esto permite identificar qué accesiones de GTDB están presentes en el árbol y cuáles no
    # la intersección se almacena en la variable 'inter'
    taxonomy_acc = set(gtdb["accession_norm"])

    inter = taxonomy_acc.intersection(tree_tips)

    # first_tip = next(iter(tree_tips))

    # La tabla taxonómica queda restringida a genomas que aparecen como tips en el árbol
    gtdb = gtdb[gtdb["accession_norm"].isin(tree_tips)].copy()

    # a partir de aqui cualquier mapeo exacto apunta a un accession que debería poder usarse para podar el árbol


    # separacion de la taxonomía GTDB en rangos taxonómicos individuales 
    # (dominio, filo, clase, orden, familia, género y especie)
    gtdb["gtdb_domain"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "d__"))
    gtdb["gtdb_phylum"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "p__"))
    gtdb["gtdb_class"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "c__"))
    gtdb["gtdb_order"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "o__"))
    gtdb["gtdb_family"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "f__"))
    gtdb["gtdb_genus"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "g__"))
    gtdb["gtdb_species"] = gtdb["taxonomy"].apply(lambda x: get_rank(x, "s__"))

    #Esta variable se usaba para inspeccionar especies del género Bacteroides, 
    # pero no influye en el resultado
    target = "g__Bacteroides"

    accepted = []
    review   = []
    mapping  = []

    # El script procesa cada taxón interpretado
    for _, row in wanted.iterrows():                     
        original = row["original_taxon"]
        target_species = row["gtdb_species"]
        target_genus = row["gtdb_genus"]

        # crea un diccionario base con la información del taxón original y los campos de mapeo inicializados a None
        # esta estructura almaecena la información de mapeo para cada taxón y se actualizará a medida que se realicen 
        # las búsquedas en GTDB y NCBI
        # conserva el género y especie de entrada
        # deja vacíos el accession, la taxonomía GTDB y el estado del mapeo.
        # sespués, dependiendo del tipo de coincidencia, se completan los campos.
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

        # si no se encuentra especie en el nombre original, 
        # se marca como "no_species_in_original_name" y se agrega a la lista de revisión para inspección manual
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
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # OJO: aqui quiero corregir algo porque aunque el género exista en GTDB, 
        # el código no busca sus especies candidatas. Podriamos guardar especies candidatas para revisión manual, 
        # pero no se hace actualmente.
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # 1. GTDB species exact
        # Busca genomas representantes cuya especie GTDB sea exactamente igual a la especie de entrada.
        exact = gtdb[gtdb["gtdb_species"] == target_species]

        # Solo se acepta automáticamente cuando existe exactamente una fila
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
        #Busca nombres que empiecen con la especie original seguida de _.
        # GTDB utiliza sufijos como _A, _B, etc., para distinguir grupos 
        # que históricamente compartían un mismo nombre
        prefix = gtdb[
            gtdb["gtdb_species"].fillna("").str.startswith(target_species + "_")
        ]

        # Si solo existe un candidato, se acepta con: "unique_species_prefix"
        # Si hay varios candidatos con diferentes sufijos, no se elige ninguno y se continúa
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
        # Busca el nombre original dentro de ncbi_organism_name
        # útil cuando el nombre NCBI corresponde al nombre de entrada, 
        # pero GTDB reclasificó el organismo con otro género o especie
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
            # Si hay exactamente una coincidencia, se acepta y se agrega a la lista de mapeos.

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
        # busca la especie dentro de la taxonomía NCBI formal almacenada en los metadatos
        # nombre incluye: identificadores de cepa; texto adicional; nombres informales; subespecies
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
        # Si no hubo coincidencia a nivel de especie, se buscan todos los representantes del mismo género GTDB
        genus_hits = gtdb[gtdb["gtdb_genus"] == target_genus].copy()

        if len(genus_hits) > 0:
            base["match_status"] = "genus_match_only"
            mapping.append(base)

            # se guardan candidatos para revisión manual, limitando a los primeros 20 resultados
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
        # Si tampoco existe el género se agrega una fila a revision
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
# convertimos resultados a tablas y generamos archivos de salida

    mapping = pd.DataFrame(mapping) # contiene una fila por cada taxón procesado
    review  = pd.DataFrame(review)
    # review no tiene necesariamente una fila por cada taxon, 
    # ya que algunos taxones pueden generar múltiples candidatos para revisión manual.

    # separacion de mapeos por categoria:

    # mapeo a nivel de especie. Incluye: coincidencias exactas, 
    # coincidencias con sufijo único, 
    # coincidencias con nombre de organismo NCBI y 
    # coincidencias con especie de taxonomía NCBI
    mapped_exact = mapping[
        mapping["match_status"].isin([
            "exact_species_match",
            "unique_species_prefix",
            "ncbi_organism_name_match",
            "ncbi_taxonomy_species_match",
        ])
    ].copy()

    # solo genero
    mapped_genus_only = mapping[
        mapping["match_status"] == "genus_match_only"
    ].copy()

    # no mapeados (Los taxones sin especie se incluyen como no mapeados, aunque su género podría existir en GTDB)
    unmapped = mapping[
        mapping["match_status"].isin([
            "no_species_or_genus_match",
            "no_species_in_original_name",
        ])
    ].copy()

    mapping.to_csv(OUT_MAPPING_ALL, sep="\t", index=False) # Contiene todos los resultados
    mapped_exact.to_csv(OUT_MAPPED_EXACT, sep="\t", index=False) # contiene solo los mapeos exactos
    mapped_genus_only.to_csv(OUT_MAPPED_GENUS_ONLY, sep="\t", index=False) # contiene solo los mapeos por género
    unmapped.to_csv(OUT_UNMAPPED, sep="\t", index=False) # contiene solo los taxones no mapeados
    review.to_csv(OUT_REVIEW, sep="\t", index=False)# contiene los casos que requieren revisión manual

    # Número de taxones aceptados a nivel de especie
    print(f"Mapeos exactos/especie: {len(mapped_exact)}") 

    # Número de taxones que solo pudieron relacionarse con un género
    print(f"Mapeos solo por género: {len(mapped_genus_only)}")

    # Número de taxones que no pudieron mapearse a GTDB
    print(f"No mapeados: {len(unmapped)}")

    # Número de taxones que requieren revisión manual
    print(f"Para revisión manual: {review['original_taxon'].nunique() if len(review) else 0}")

    print(f"Archivo completo: {OUT_MAPPING_ALL}")
    print(f"Archivo exactos: {OUT_MAPPED_EXACT}")
    print(f"Archivo solo género: {OUT_MAPPED_GENUS_ONLY}")
    print(f"Archivo no mapeados: {OUT_UNMAPPED}")
    print(f"Archivo revisión: {OUT_REVIEW}")

# SANITY -----------------------------------------------------------------------
    # vuelve a leer el archivo que acaba de guardar
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
    # Esto comprueba que todos los accessions seleccionados existen en el árbol
    #  y, por tanto, podrían usarse para podarlo
# SANITY -----------------------------------------------------------------------


if __name__ == "__main__":
    main()







# ---------------------------------------------------------------------------
# TODO 1:
# Decidir cómo tratar los taxones que contienen género pero no especie.
#
# Líneas pertinentes:
# - Editar líneas 519-532:
#   bloque que clasifica inmediatamente el taxón como
#   "no_species_in_original_name" y ejecuta continue.
#
# - Reutilizar o convertir en función las líneas 660-677:
#   bloque que busca candidatos mediante target_genus y guarda hasta
#   20 candidatos en review.
#
# - Revisar líneas 722-728:
#   clasificación final de "no_species_in_original_name" dentro de unmapped.
#
# Pendiente:
# - En las líneas 519-532, buscar target_genus en gtdb antes de ejecutar continue.
# - Si el género existe, usar "genus_match_only" o un estado más específico
#   como "genus_match_without_input_species".
# - Guardar las especies candidatas en review.
# - Si el género no existe, conservar "no_species_in_original_name"
#   o usar "genus_not_found".
# - Ajustar las líneas 722-728 de acuerdo con los nuevos estados.

# TODO 2:
# Detectar y registrar coincidencias ambiguas en las búsquedas GTDB.
#
# Líneas pertinentes:
# - Líneas 541-560:
#   búsqueda y aceptación de coincidencia exacta GTDB.
#   Falta manejar len(exact) > 1 después de la línea 560.
#
# - Líneas 567-589:
#   búsqueda y aceptación por prefijo.
#   Falta manejar len(prefix) > 1 después de la línea 589.
#
# - Líneas 617-624:
#   ejemplo actual de cómo se registra una coincidencia NCBI ambigua.
#
# - Líneas 650-657:
#   segundo ejemplo actual de ambigüedad NCBI.
#
# - Líneas 504-517:
#   diccionario base donde podría guardarse el estado o conteo de ambigüedades.
#
# Pendiente:
# - Agregar un bloque elif len(exact) > 1 después de la línea 560.
# - Agregar un bloque elif len(prefix) > 1 después de la línea 589.
# - Guardar en review una fila por candidato, no solamente una advertencia
#   general sin accession.
# - Registrar en base que se encontraron coincidencias ambiguas.
# - Decidir si el proceso debe continuar hacia NCBI o detenerse.


# TODO 3:
# Cambiar el nombre mapped_exact por uno que describa correctamente su contenido.
#
# Líneas pertinentes:
# - Línea 199:
#   renombrar OUT_MAPPED_EXACT y cambiar "mapped_exact.tsv".
#
# - Líneas 704-715:
#   renombrar mapped_exact a mapped_species_level.
#
# - Línea 731:
#   actualizar la escritura del archivo.
#
# - Línea 737:
#   actualizar el nombre de la variable en el resumen impreso.
#
# - Línea 749:
#   actualizar el texto y la constante mostrada.
#
# - Línea 756:
#   cambiar la ruta leída por la sección SANITY.
#
# Pendiente:
# - Cambiar OUT_MAPPED_EXACT por OUT_MAPPED_SPECIES.
# - Cambiar "Data/mapped_exact.tsv" por "Data/mapped_species.tsv".
# - Cambiar mapped_exact por mapped_species_level.
# - Cambiar el comentario "coincidencias exactas".
# - Actualizar todas las referencias posteriores.

# TODO 4:
# Evitar que review mezcle advertencias provisionales con casos que finalmente
# sí fueron resueltos.
#
# Líneas pertinentes:
# - Líneas 617-624:
#   agrega inmediatamente la ambigüedad de ncbi_organism_name a review.
#
# - Líneas 630-648:
#   el mismo taxón todavía puede resolverse posteriormente mediante
#   ncbi_taxonomy_species_match.
#
# - Líneas 650-657:
#   agrega inmediatamente otra ambigüedad a review.
#
# - Líneas 660-690:
#   determinan el resultado final cuando no se resuelve a nivel de especie.
#
# - Líneas 697-700:
#   conversión final de review a DataFrame.
#
# - Línea 734:
#   escritura del archivo de revisión.
#
# - Línea 746:
#   conteo de taxones supuestamente pendientes de revisión.
#
# Pendiente:
# - Decidir si review guarda todos los eventos de diagnóstico o únicamente
#   casos cuyo resultado final requiere revisión.
# - No agregar directamente a review las advertencias provisionales de las
#   líneas 617-624 y 650-657.
# - Guardarlas temporalmente dentro de la iteración.
# - Añadirlas a review solo si el taxón no se resuelve después.
# - Como alternativa, crear diagnostics y review como listas separadas.
# - Ajustar la escritura y los conteos de las líneas 697-746.


# TODO 6:
# Eliminar o reactivar variables que actualmente no participan en el resultado.
#
# Líneas pertinentes:
# - Líneas 459-464:
#   taxonomy_acc e inter se calculan, pero no se utilizan.
#
# - Línea 466:
#   first_tip ya está comentada y puede eliminarse definitivamente.
#
# - Líneas 484-486:
#   target se define, pero no se utiliza.
#
# - Línea 488:
#   accepted se inicializa, pero no se utiliza.
#
# - Líneas 489-490:
#   review y mapping sí se utilizan y deben conservarse.
#
# Pendiente:
# - Eliminar las líneas 462-464 si ya se validó la intersección.
# - Eliminar la línea 466.
# - Eliminar las líneas 484-486.
# - Eliminar la línea 488.
# - Conservar únicamente review y mapping en las líneas 489-490.

# TODO 7:
# Revisar las pruebas SANITY que esperan una correspondencia uno a uno.
#
# Líneas pertinentes:
# - Línea 756:
#   lectura del archivo de mapeos a nivel de especie.
#
# - Líneas 757-760:
#   conteo de filas, taxones, accessions y especies.
#
# - Líneas 761-762:
#   expectativa fija de 275 para todas las categorías.
#
# - Líneas 764-770:
#   comprobación de que los accessions están presentes en el árbol.
#
# - Línea 771:
#   expectativa fija de 275/275.
#
# Pendiente:
# - Mantener la comprobación de una sola fila por original_taxon.
# - No exigir que todos los accessions sean únicos.
# - No exigir que todas las especies GTDB sean únicas.
# - Identificar y mostrar accessions compartidos.
# - Identificar y mostrar especies GTDB compartidas.
# - Sustituir las expectativas fijas de las líneas 761-762 y 771 por
#   comprobaciones basadas en el tamaño real del DataFrame.
# - Usar normalize_accession() en las líneas 765-768 para mantener una sola
#   regla de normalización.

# TODO 8:
# Guardar en mapping un resumen de las coincidencias ambiguas encontradas
# durante el proceso.
#
# Líneas pertinentes:
# - Líneas 504-517:
#   definición del diccionario base. Aquí deben agregarse las columnas
#   diagnósticas.
#
# - Línea 541:
#   después de crear exact, registrar len(exact).
#
# - Línea 567:
#   después de crear prefix, registrar len(prefix).
#
# - Línea 596:
#   después de crear ncbi_org, registrar len(ncbi_org).
#
# - Línea 630:
#   después de crear ncbi_tax, registrar len(ncbi_tax).
#
# - Líneas 617-624 y 650-657:
#   lugares donde ya se detectan ambigüedades NCBI.
#
# - Líneas 697 y 730:
#   conversión y escritura de mapping, que conservarán automáticamente
#   las nuevas columnas.
#
# Pendiente:
# - Agregar a base columnas para los conteos de candidatos.
# - Agregar una columna ambiguous_methods.
# - Actualizar los conteos inmediatamente después de cada búsqueda.
# - Acumular los nombres de los métodos ambiguos.
# - Guardar el resumen antes de ejecutar mapping.append(base).


# ---------------------------------------------------------------------------
