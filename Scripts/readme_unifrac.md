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
      OJO: estos archivos .gz no estan en el repositorio debido a que son muy pesados, se pueden descargar de aqui https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/
      se usaron los datos del 2026-03-24

* **Metodología:**
  1. Lee la tabla original sin encabezado en `pathways_species2.tsv` y usa la columna 1 como la columna de bacterias. Extrae los taxa que contienen `g__` y limpia/parsea buscando patrones tipo `g__Generos__Genero_especie`.
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
  4. Convierte los valores de las muestras a formato numérico (valores invál
