# Pipeline de Análisis UniFrac para Pathways Metabólicos

Este repositorio contiene los scripts necesarios para mapear taxas bacterianas contra la base de datos GTDB, preparar matrices de abundancia, calcular distancias UniFrac (ponderadas y no ponderadas) entre pathways metabólicos y visualizar los resultados mediante un análisis de coordenadas principales (PCoA).

---

## Estructura de los Scripts

### 📁 01_map_bacteria_to_gtdb.py
[cite_start]Este script toma los nombres bacterianos originales y los intenta mapear contra los existentes en el árbol filogenético descargado de GTDB[cite: 1].

* **INPUT:**
  * [cite_start]`Data/pathways_species2.tsv` 
  * [cite_start]**Bases GTDB:** * `Data/gtdb/bac120_taxonomy.tsv.gz` [cite: 2]
    * [cite_start]`Data/gtdb/bac120.tree.gz` [cite: 2]
    * [cite_start]`Data/gtdb/bac120_metadata.tsv.gz` [cite: 2]
      Bajar archivos de aqui: https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/
      

* **Metodología:**
  1. [cite_start]Lee la tabla original sin encabezado en `pathways_species2.tsv` y usa la columna 1 como la columna de bacterias[cite: 2]. Extrae los taxa que contienen `g__` y limpia/parsea buscando patrones tipo `g__Generos__Genero_especie`[cite: 3].
  2. [cite_start]Carga el árbol GTDB `bac120`, toma sus tips y normaliza los accessions reemplazando `_` por espacio (crucial debido a diferencias de formato entre el árbol y las tablas)[cite: 4, 5].
  3. [cite_start]Carga la taxonomía y metadata de GTDB, filtrando para conservar solo los accessions presentes en el árbol[cite: 6]. [cite_start]Esto evita mapear bacterias a accessions inexistentes en la filogenia[cite: 7].
  4. Realiza el mapeo en el siguiente orden de prioridad:
     * [cite_start]Match exacto contra especie GTDB[cite: 8].
     * Match por prefijo de especie GTDB[cite: 8].
     * [cite_start]Match exacto contra `ncbi_organism_name` en la metadata[cite: 9].
     * [cite_start]Match exacto contra `ncbi_taxonomy`[cite: 9].
     * [cite_start]Si no encuentra especie, intenta match por género[cite: 9]. [cite_start]Si tampoco encuentra, queda como **no mapeada**[cite: 10].

* **OUTPUT:**
  * [cite_start]`Data/gtdb_mapping_all.tsv`: Todos los intentos de mapeo[cite: 10].
  * `Data/mapped_exact.tsv`: Taxa mapeados a nivel de especie (utilizados para el análisis UniFrac)[cite: 10, 11].
  * [cite_start]`Data/mapped_genus_only.tsv`: Taxa mapeados solo por género (no se usan automáticamente para UniFrac de manera estricta)[cite: 11, 12].
  * [cite_start]`Data/unmapped.tsv`: Taxa no mapeados[cite: 12].
  * `Data/gtdb_mapping_review.tsv`: Candidatos ambiguos que requieren revisión manual[cite: 13].

> 🔍 **Nota de control:** Al final realiza *sanity checks* esperando 275 filas/taxa/accessions únicos presentes en el árbol[cite: 13, 14]. Estos números dependen del refinamiento; se espera obtener un número mayor en etapas posteriores (hasta 615)[cite: 15].

---

### 📁 02_prepare_unifrac_inputs.py
Prepara los dos objetos necesarios para comparar pathways mediante UniFrac: una tabla de abundancia pathway × taxón y un árbol filogenético GTDB podado[cite: 16].

* **INPUT:**
  * `Data/relab_table_joined.tsv` [cite: 17]
  * [cite_start]`Data/gtdb/mapped_exact.tsv` [cite: 17]
  * `Data/gtdb/bac120.tree.gz` [cite: 17]

* **Metodología:**
  1. [cite_start]Carga la tabla de abundancia relativa por muestra (`Pathway | taxón`)[cite: 18].
  2. [cite_start]Separa los nombres de los pathways de los nombres originales de los taxa[cite: 19].
  3. [cite_start]Reemplaza los nombres originales por sus respectivos accessions GTDB usando el archivo `mapped_exact.tsv`[cite: 20, 22].
  4. [cite_start]Convierte los valores de las muestras a formato numérico (valores inválidos/ausentes se transforman en cero)[cite: 23].
  5. [cite_start]Agrupa los datos y construye una nueva matriz con la estructura: **Filas = pathways**, **Columnas = taxa (accessions GTDB)**, **Valores = abundancia**[cite: 25].
  6. [cite_start]Normaliza cada pathway para que la suma de sus abundancias sea igual a 1, permitiendo su interpretación como una distribución taxonómica comparable[cite: 26].
  7. [cite_start]Poda el árbol filogenético completo de GTDB conservando únicamente los taxa presentes en esta nueva tabla[cite: 27].

* **OUTPUT:**
  * [cite_start]**Tabla de abundancia:** `Data/pathway_unifrac_abundance_table.tsv` [cite: 27]
  * **Árbol podado:** `Data/pathway_tree_pruned.nwk` [cite: 27]

> 🧬 **Objetivo biológico:** A diferencia del uso clásico (muestras biológicas), este análisis utiliza a cada **pathway como una "unidad ecológica"**[cite: 28]. Busca responder qué pathways están asociados a conjuntos filogenéticamente distintos de bacterias[cite: 30].

---

### 📁 03_compute_unifrac.py
Calcula las matrices de distancia UniFrac entre los pathways[cite: 31].

* **INPUT:**
  * `Data/pathway_unifrac_abundance_table.tsv` [cite: 32]
  * [cite_start]`Data/pathway_tree_pruned.nwk` [cite: 32]

* **Metodología:**
  * [cite_start]Calcula la matriz usando **Weighted UniFrac** (incorpora abundancia relativa; distancias pequeñas implican linajes similares con distribuciones parecidas)[cite: 34].
  * Calcula la matriz usando **Unweighted UniFrac** (usa presencia/ausencia; distancias pequeñas implican que comparten los mismos linajes independientemente de su abundancia)[cite: 36].

* **OUTPUT:**
  * `Data/pathway_weighted_unifrac_distance_matrix.tsv` [cite: 38]
  * [cite_start]`Data/pathway_unweighted_unifrac_distance_matrix.tsv` [cite: 38]

---

### 📁 identify_common_pathway_pairs.py
[cite_start]Identifica pathways comunes entre la tabla original de asociaciones y la matriz UniFrac, generando todos los pares posibles ordenados por su distancia no ponderada[cite: 44, 45, 52].

* **INPUT:**
  * [cite_start]`Data/pathways_species2.tsv` [cite: 46]
  * [cite_start]`Data/pathway_unweighted_unifrac_distance_matrix.tsv` [cite: 46]

* **Metodología:**
  1. Normaliza y corrige diferencias de formato en los nombres de los pathways (ej. convierte `:_` a `: ` y guiones bajos a espacios)[cite: 48].
  2. Identifica intersecciones, genera todos los pares posibles y recupera sus distancias[cite: 49, 50, 51].
  3. Ordena los pares de menor a mayor distancia e imprime ejemplos ilustrativos[cite: 52, 53].

* **OUTPUT:**
  * `Data/common_pathways_normalized.tsv` [cite: 54]
  * [cite_start]`Data/example_pathway_pairs_from_unifrac.tsv` [cite: 54]

---

### 📁 04_plot_unweighted_unifrac_pcoa.py
[cite_start]Visualiza las relaciones filogenéticas entre pathways proyectando las distancias en un espacio bidimensional mediante **PCoA (Principal Coordinates Analysis)**[cite: 58, 59].

* **INPUT:**
  * [cite_start]`Data/pathway_unweighted_unifrac_distance_matrix.tsv` [cite: 60]

* **Metodología:**
  1. [cite_start]Destaca un par cercano seleccionado manualmente (ej. `METHANOGENESIS-PWY` y `PWY-8112`, asociados a arqueas metanógenas con distancia cercana a cero)[cite: 61].
  2. [cite_start]Encuentra y destaca automáticamente el par más distante de toda la matriz[cite: 62].
  3. [cite_start]Realiza el PCoA (vía `scikit-bio`) y grafica la estructura global usando baja opacidad para los puntos de fondo[cite: 66, 67, 69].
  4. [cite_start]Resalta los pares clave con marcadores diferenciados (círculos para cercanos, cuadrados para distantes) y añade etiquetas ajustadas para evitar sobrelapamientos[cite: 65, 70, 71].

* **OUTPUT:**
  * [cite_start]**Figura PCoA:** `Plots/unweighted_unifrac_pcoa_pathways_selected_pairs.png` [cite: 74]

---

## 📌 Guía de Interpretación de Resultados

* [cite_start]**Distancia cercana a 0:** Los dos pathways están asociados a comunidades microbianas filogenéticamente muy parecidas (funciones ejecutadas por organismos evolutivamente similares)[cite: 43, 78].
* **Distancia cercana a 1:** Los pathways están asociados a conjuntos de linajes muy diferentes e independientes[cite: 43, 79].

⚠️ **Punto Crítico del Pipeline:** Aunque el pipeline calcula tanto la métrica *Weighted* como la *Unweighted*, actualmente **solo se guarda y grafica la versión Unweighted (presencia/ausencia)**[cite: 96]. Si el objetivo biológico requiere evaluar los cambios basados en la abundancia bacteriana, se deberá modificar el script de graficado para explotar la matriz ponderada[cite: 96].

---

## 📚 Marco Teórico y Referencias

### Lozupone & Knight (2007)
UniFrac cuantifica cuánta historia evolutiva es exclusiva de cada comunidad medido como la fracción de longitud de ramas que lleva a descendientes presentes en una comunidad o en la otra, pero no en ambas[cite: 81, 82, 83]. En su estudio masivo determinaron que la salinidad es el factor que más fuertemente separa a las comunidades bacterianas[cite: 84, 86].

### Métrica UniFrac (Definición General)
A diferencia de índices clásicos como Bray-Curtis, UniFrac incorpora las relaciones filogenéticas entre organismos[cite: 91].
$$\text{Unweighted UniFrac} = \frac{\text{Longitud de ramas no compartidas}}{\text{Longitud total del árbol}}$$
