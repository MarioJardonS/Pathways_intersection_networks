#!/bin/bash

INPUT_DIR="Networks"
OUTPUT_DIR="Networks_Separadas"

mkdir -p "$OUTPUT_DIR"

echo "Procesando archivos .tsv (Solo Object ID)..."

for file in "$INPUT_DIR"/*.tsv; do
    filename=$(basename "$file")
    output_file="$OUTPUT_DIR/$filename"
    
    awk 'BEGIN {FS="\t"; OFS="\t"} 
    {
        # --- Para la Columna 1 ---
        id1 = $1; 
        if ($1 ~ /:_/) { sub(/:_.*/, "", id1); }
        
        # --- Para la Columna 2 ---
        id2 = $2; 
        if ($2 ~ /:_/) { sub(/:_.*/, "", id2); }
        
        # --- Imprimir la nueva fila ---
        # Imprimimos únicamente el ID de la col 1 y el ID de la col 2
        printf "%s\t%s", id1, id2;
        
        # Imprimir el resto de las columnas numéricas (OTUs, pesos, etc.)
        for(i=3; i<=NF; i++) {
            printf "\t%s", $i;
        }
        
        printf "\n";
        
    }' "$file" > "$output_file"

done

echo "¡Completado exitosamente! Archivos listos en $OUTPUT_DIR/"
