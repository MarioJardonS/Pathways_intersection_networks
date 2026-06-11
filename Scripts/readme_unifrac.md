# Pipeline de Análisis UniFrac para Pathways Metabólicos

Este repositorio contiene los scripts necesarios para mapear taxas bacterianas contra la base de datos GTDB, preparar matrices de abundancia, calcular distancias UniFrac (ponderadas y no ponderadas) entre pathways metabólicos y visualizar los resultados mediante un análisis de coordenadas principales (PCoA).

---

## Estructura de los Scripts

### 📁 01_map_bacteria_to_gtdb.py
Este script toma los nombres bacterianos originales y los intenta mapear contra los existentes en el árbol filogenético descargado de GTDB.

* **INPUT:**
  * `Data/pathways_species2.tsv`
  * **Bases GTDB:** * `Data/gtdb/bac120_taxonomy.tsv.gz`
    * `Data/gtdb/bac120.tree.gz`
    * `Data/gtdb/bac120_metadata.tsv.gz`
    * *Nota: Puedes bajar los archivos de GTDB desde aquí: https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/*

* **Metodología:**
  1. Lee la tabla original sin encabezado en `pathways_species2.tsv` y usa la columna 1 como la columna de bacterias. Extrae los taxa que contienen `g__` y limpia/parsea buscando patrones tipo `g__Genero__especie`.
  2. Carga el árbol GTDB `bac120`, toma sus tips y normaliza los accessions reemplazando `_` por espacio (crucial debido a diferencias de formato entre el árbol y las tablas).
  3. Carga la taxonomía y metadata de GTDB, filtrando para conservar solo los accessions presentes en el árbol. Esto evita mapear bacterias a accessions inexistentes en la filogenia.
  4. Realiza el mapeo en el siguiente orden de prioridad:
     * Match exacto contra especie GTDB.
     * Match por prefijo de especie GTDB.
     * Match exacto contra `ncbi_organism_name` en la metadata.
     * Match exacto contra `ncbi_taxonomy`.
     * Si no encuentra especie, intenta match por género. Si tampoco encuentra, queda como **no mapeada**.

* **OUTPUT:**
  * `Data/gtdb_mapping_all.tsv`: Todos los intentos de mapeo.
  * `Data/mapped_exact.tsv`: Taxa mapeados a nivel de especie (utilizados para el análisis UniFrac).
  * `Data/mapped_genus_only.tsv`: Taxa mapeados solo por género (no se usan automáticamente para UniFrac de manera estricta).
  * `Data/unmapped.tsv`: Taxa no mapeados.
  * `Data/gtdb_mapping_review.tsv`: Candidatos ambiguos que requieren revisión manual.

> 🔍 **Nota de control:** Al final realiza *sanity checks* esperando 275 filas/taxa/accessions únicos presentes en el árbol. Estos números dependen del refinamiento; se espera obtener un número mayor en etapas posteriores (hasta 615).

---

### 📁 02_prepare_unifrac_inputs.py
Prepara los dos objetos necesarios para comparar pathways mediante UniFrac: una tabla de abundancia pathway × taxón y un árbol filogenético GTDB podado.

* **INPUT:**
  * `Data/relab_table_joined.tsv`
  * `Data/gtdb/mapped_exact.tsv`
  * `Data/gtdb/bac120.tree.gz`

* **Metodología:**
  1. Carga la tabla de abundancia relativa por muestra (`Pathway | taxón`).
  2. Separa los nombres de los pathways de los nombres originales de los taxa.
  3. Reemplaza los nombres originales por sus respectivos accessions GTDB usando el archivo `mapped_exact.tsv`.
  4. Convierte los valores de las muestras a formato numérico (valores inválidos/ausentes se transforman en cero).
  5. Agrupa los datos y construye una nueva matriz con la estructura: **Filas = pathways**, **Columnas = taxa (accessions GTDB)**, **Valores = abundancia**.
  6. Normaliza cada pathway para que la suma de sus abundancias sea igual a 1, permitiendo su interpretación como una distribución taxonómica comparable.
  7. Poda el árbol filogenético completo de GTDB conservando únicamente los taxa presentes en esta nueva tabla.

* **OUTPUT:**
  * **Tabla de abundancia:** `Data/pathway_unifrac_abundance_table.tsv`
  * **Árbol podado:** `Data/pathway_tree_pruned.nwk`

> 🧬 **Objetivo biológico:** A diferencia del uso clásico (muestras biológicas), este análisis utiliza a cada **pathway como una "unidad ecológica"**. Busca responder qué pathways están asociados a conjuntos filogenéticamente distintos de bacterias.

---

### 📁 03_compute_unifrac.py
Calcula las matrices de distancia UniFrac entre los pathways.

* **INPUT:**
  * `Data/pathway_unifrac_abundance_table.tsv`
  * `Data/pathway_tree_pruned.nwk`

* **Metodología:**
  * Calcula la matriz usando **Weighted UniFrac** (incorpora abundancia relativa; distancias pequeñas implican linajes similares con distribuciones parecidas).
  * Calcula la matriz usando **Unweighted UniFrac** (usa presencia/ausencia; distancias pequeñas implican que comparten los mismos linajes independientemente de su abundancia).

* **OUTPUT:**
  * `Data/pathway_weighted_unifrac_distance_matrix.tsv`
  * `Data/pathway_unweighted_unifrac_distance_matrix.tsv`

---

### 📁 identify_common_pathway_pairs.py
Identifica pathways comunes entre la tabla original de asociaciones y la matriz UniFrac, generando todos los pares posibles ordenados por su distancia no ponderada.

* **INPUT:**
  * `Data/pathways_species2.tsv`
  * `Data/pathway_unweighted_unifrac_distance_matrix.tsv`

* **Metodología:**
  1. Normaliza y corrige diferencias de formato en los nombres de los pathways (ej. convierte `:_` a `: ` y guiones bajos a espacios).
  2. Identifica intersecciones, genera todos los pares posibles y recupera sus distancias.
  3. Ordena los pares de menor a mayor distancia e imprime ejemplos ilustrativos.

* **OUTPUT:**
  * `Data/common_pathways_normalized.tsv`
  * `Data/example_pathway_pairs_from_unifrac.tsv`

---

### 📁 04_plot_unweighted_unifrac_pcoa.py
Visualiza las relaciones filogenéticas entre pathways proyectando las distancias en un espacio bidimensional mediante **PCoA (Principal Coordinates Analysis)**.

* **INPUT:**
  * `Data/pathway_unweighted_unifrac_distance_matrix.tsv`

* **Metodología:**
  1. Destaca un par cercano seleccionado manualmente (ej. `METHANOGENESIS-PWY` y `PWY-8112`, asociados a arqueas metanógenas con distancia cercana a cero).
  2. Encuentra y destaca automáticamente el par más distante de toda la matriz.
  3. Realiza el PCoA (vía `scikit-bio`) y grafica la estructura global usando baja opacidad para los puntos de fondo.
  4. Resalta los pares clave con marcadores diferenciados (círculos para cercanos, cuadrados para distantes) y añade etiquetas ajustadas para evitar sobrelapamientos.

* **OUTPUT:**
  * **Figura PCoA:** `Plots/unweighted_unifrac_pcoa_pathways_selected_pairs.png`

---

## 📌 Guía de Interpretación de Resultados

* **Distancia cercana a 0:** Los dos pathways están asociados a comunidades microbianas filogenéticamente muy parecidas (funciones ejecutadas por organismos evolutivamente similares).
* **Distancia cercana a 1:** Los pathways están asociados a conjuntos de linajes muy diferentes e independientes.

⚠️ **Punto Crítico del Pipeline:** Aunque el pipeline calcula tanto la métrica *Weighted* como la *Unweighted*, actualmente **solo se guarda y grafica la versión Unweighted (presencia/ausencia)**. Si el objetivo biológico requiere evaluar los cambios basados en la abundancia bacteriana, se deberá modificar el script de graficado para explotar la matriz ponderada.

---

## 📚 Marco Teórico y Referencias

### Lozupone & Knight (2007)
UniFrac cuantifica cuánta historia evolutiva es exclusiva de cada comunidad medido como la fracción de longitud de ramas que lleva a descendientes presentes en una comunidad o en la otra, pero no en ambas. En su estudio masivo determinaron que la salinidad es el factor que más fuertemente separa a las comunidades bacterianas.

### Métrica UniFrac (Definición General)
A diferencia de índices clásicos como Bray-Curtis, UniFrac incorpora las relaciones filogenéticas entre organismos.

$$\text{Unweighted UniFrac} = \frac{\text{Longitud de ramas no compartidas}}{\text{Longitud total del árbol}}$$
