# 🧬 Sec61-Target-Scout

Pipeline automático para la extracción, curación y deduplicación de secuencias de **Péptidos Señal (Signal Peptides)** desde UniProt, orientado a la síntesis química y ensayos de inhibición del translocón Sec61.

---

## 📋 Tabla de Contenidos

- [Objetivo](#objetivo)
- [Contexto Biológico](#contexto-biológico)
- [Instalación](#instalación)
- [Uso Rápido](#uso-rápido)
- [Presets Disponibles](#presets-disponibles)
- [Argumentos CLI](#argumentos-cli)
- [Columnas del Output](#columnas-del-output)
- [Lógica de Parsing por Grupo](#lógica-de-parsing-por-grupo)
- [Evidencia y Calidad](#evidencia-y-calidad)
- [Deduplicación](#deduplicación)
- [Features Computados](#features-computados)
- [Estructura del Proyecto](#estructura-del-proyecto)
- [Extensibilidad](#extensibilidad)

---

## Objetivo

Generar un listado curado de péptidos señal para:

1. **Síntesis química** (wet lab)
2. **Ensayos de sensibilidad a inhibición de Sec61**

La precisión en la secuencia de aminoácidos y la distinción entre tipos de señales es crítica.

---

## Contexto Biológico

### Signal Peptides Type 1 (Clivables)

Péptidos señal clásicos N-terminales que dirigen la proteína al retículo endoplásmico (RE) vía Sec61 y son clivados por la señal peptidasa. Estructura tripartita:

- **n-region**: N-terminal, carga positiva
- **h-region**: Núcleo hidrofóbico (7-15 aa)
- **c-region**: C-terminal con motivo de corte (residuos pequeños en -1 y -3)

### Signal Anchors Type 2 (No clivables)

Dominios transmembrana N-terminales que funcionan como señal de translocación pero NO son clivados. Funcionan como ancla permanente en la membrana del RE. Ejemplo: Neuraminidasa de Influenza.

### Señales Virales

**Flavivirus** (Dengue, Zika, WNV, etc.):
La poliproteína viral contiene múltiples señales internas:

- **anchC** (C-terminal de Capsid) → señal para translocación de prM
- **prM C-term TM** → señal para translocación de E
- **E C-term TM2** → señal para translocación de NS1
- Todas son **Sec61-dependientes**

**Alphavirus** (Chikungunya, Sindbis, etc.):
Polipéptido estructural: C → E3 → E2 → 6K → E1

- **Capsid (C)**: Se auto-cliva en el citoplasma, **no pasa por Sec61**
- **E3**: Señal para translocación de E2 (precursor p62). **Sec61-dependiente** ✅
- **6K**: Señal para translocación de E1. **Sec61-dependiente** ✅

---

## Instalación

```bash
# Clonar o navegar al directorio del proyecto
cd sps_pipeline

# Crear entorno virtual (recomendado)
python -m venv venv
venv\Scripts\activate       # Windows
# source venv/bin/activate  # Linux/Mac

# Instalar dependencias
pip install -r requirements.txt
```

**Dependencias**: `requests`, `pandas`, `tqdm`, `pyyaml`

---

## Uso Rápido

```bash
# Un solo preset (SPs humanos Type 1, solo Swiss-Prot)
python -m src.main --preset human_sp1 --output csv

# Múltiples presets en una sola corrida
python -m src.main --preset human_sp1 human_sa2 influenza_ha --output csv

# Todos los presets
python -m src.main --preset all --output csv fasta

# Query personalizada
python -m src.main --query "(gene:INS) AND (organism_id:9606)" --output csv

# Solo evidencia experimental
python -m src.main --preset human_sp1 --evidence experimental

# Incluir TrEMBL (entradas no revisadas)
python -m src.main --preset human_sp1 --include-unreviewed

# Testing rápido (limitar a 50 resultados)
python -m src.main --preset human_sp1 --max-results 50

# Sin deduplicación
python -m src.main --preset human_sp1 --no-dedup

# Múltiples formatos
python -m src.main --preset human_sp1 --output csv fasta tsv

# Nombrar el archivo de salida con un alias
python -m src.main --preset human_sp1 --name mi_experimento --output csv
# → genera mi_experimento.csv (sin alias: human_sp1_20260217_115522.csv)

# Listar presets disponibles
python -m src.main --list-presets
```

Los archivos se generan en la carpeta `output/` del proyecto.

---

## Presets Disponibles

| Preset                    | Descripción                                  | Tipo SP | Base de datos       |
| ------------------------- | -------------------------------------------- | ------- | ------------------- |
| **Humanos**               |                                              |         |                     |
| `human_sp1`               | SPs clivables humanos                        | Type 1  | Swiss-Prot          |
| `human_sp1_all`           | SPs clivables humanos                        | Type 1  | Swiss-Prot + TrEMBL |
| `human_sa2`               | Signal Anchors humanos                       | Type 2  | Swiss-Prot          |
| `human_sa2_all`           | Signal Anchors humanos                       | Type 2  | Swiss-Prot + TrEMBL |
| **Influenza**             |                                              |         |                     |
| `influenza_ha`            | HA de Influenza A (todos los subtipos)       | Type 1  | Todos               |
| `influenza_na`            | NA de Influenza A (todos los subtipos)       | Type 2  | Todos               |
| **Flavivirus/Alphavirus** |                                              |         |                     |
| `flavivirus`              | Señales de poliproteína (Dengue, Zika, etc.) | Mixto   | Swiss-Prot          |
| `alphavirus`              | E3 + 6K — Género completo                    | Mixto   | Todos               |
| `veev`                    | VEEV (Venezuelan Equine Encephalitis)        | Mixto   | Todos               |
| `sinv`                    | SINV (Sindbis Virus)                         | Mixto   | Todos               |
| `eilv`                    | EILV (Eilat Virus, insect-specific)          | Mixto   | Todos (TrEMBL)      |
| **Orthobunyavirus**       |                                              |         |                     |
| `orthobunyavirus`         | GPC signal — Género completo                 | Type 1  | Todos               |
| `orov`                    | OROV (Oropouche Virus)                       | Type 1  | Todos (TrEMBL)      |
| `lacv`                    | LACV (La Crosse Virus)                       | Type 1  | Todos               |
| `sbv`                     | SBV (Schmallenberg Virus)                    | Type 1  | Todos               |
| **Orthohantavirus**       |                                              |         |                     |
| `orthohantavirus`         | GPC signal — Género completo                 | Type 1  | Todos               |
| `andv`                    | ANDV (Andes Virus)                           | Type 1  | Todos               |
| `snv`                     | SNV (Sin Nombre Virus)                       | Type 1  | Todos               |
| `bccv`                    | BCCV (Black Creek Canal Virus)               | Type 1  | Todos               |

Cuando se usan **múltiples presets**, los resultados se concatenan en un único archivo con una columna `query_group` que identifica el origen de cada registro, y luego se deduplican globalmente.

---

## Argumentos CLI

| Argumento                                     | Descripción                                 | Default         |
| --------------------------------------------- | ------------------------------------------- | --------------- |
| `--preset NAME [NAME ...]`                    | Preset(s) a ejecutar. Usar `all` para todos | Requerido\*     |
| `--query STRING`                              | Query personalizada para UniProt            | Requerido\*     |
| `--evidence {all,experimental,by_similarity}` | Filtrar por nivel de evidencia              | `all`           |
| `--include-unreviewed`                        | Incluir TrEMBL (no revisadas)               | Solo Swiss-Prot |
| `--max-results N`                             | Límite de resultados (solo testing)         | Sin límite      |
| `--output {csv,tsv,fasta} [...]`              | Formato(s) de salida                        | `csv`           |
| `--output-dir PATH`                           | Directorio de salida                        | `output/`       |
| `--name STRING`                               | Alias para el archivo de salida             | Auto-generado   |
| `--no-dedup`                                  | Desactivar deduplicación                    | Dedup activa    |
| `--verbose`                                   | Logging detallado                           | Off             |
| `--list-presets`                              | Listar presets y salir                      | —               |

\* `--preset` y `--query` son mutuamente excluyentes.

---

## Columnas del Output

| Columna                              | Descripción                                                                  |
| ------------------------------------ | ---------------------------------------------------------------------------- |
| `accession`                          | Accession de UniProt (ej: P01308)                                            |
| `gene`                               | Nombre del gen                                                               |
| `organism`                           | Especie (nombre científico)                                                  |
| `taxonomy_id`                        | ID de taxonomía NCBI                                                         |
| `query_group`                        | Preset de origen (ej: `human_sp1`, `influenza_ha`)                           |
| `sp_type`                            | Tipo: `SIGNAL_PEPTIDE_TYPE1`, `SIGNAL_ANCHOR_TYPE2`, `VIRAL_INTERNAL_SIGNAL` |
| `sp_subtype`                         | Para virales: `E3→p62/E2 signal`, `6K→E1 signal`, etc.                       |
| `evidence_category`                  | `EXPERIMENTAL`, `BY_SIMILARITY`, `PREDICTED`                                 |
| `evidence_codes`                     | Códigos ECO crudos (pipe-separated)                                          |
| `sp_sequence_aa`                     | Secuencia de aminoácidos del péptido señal                                   |
| `sp_length`                          | Largo en residuos                                                            |
| `start_pos` / `end_pos`              | Posición en la proteína precursora (1-based)                                 |
| `cleavage_site_motif`                | Motivo: `...AFA\|DPVV...` o `N/A` para anchors                               |
| `hydrophobicity_mean`                | Media Kyte-Doolittle                                                         |
| `net_charge_ph7`                     | Carga neta estimada a pH 7                                                   |
| `n_region` / `h_region` / `c_region` | Secuencias de las 3 regiones del SP                                          |
| `full_sequence`                      | Secuencia completa de la proteína                                            |
| `duplicate_count`                    | Nro de entradas con SP idéntico                                              |
| `all_accessions`                     | Todos los accessions con mismo SP                                            |
| `source`                             | `Swiss-Prot` o `TrEMBL`                                                      |
| `reviewed`                           | `true` / `false`                                                             |

---

## Lógica de Parsing por Grupo

### Humanos Type 1 (`sp1`)

- Busca feature `type: "Signal"` en el JSON de UniProt
- Extrae secuencia del SP (1-based a Python slice)
- Computa motivo de sitio de corte con ±5 residuos de contexto

### Humanos Type 2 / Signal Anchors (`sa2`)

- **Estrategia 1**: Feature `type: "Signal"` en entradas con keyword KW-0735
- **Estrategia 2**: Feature `type: "Transmembrane"` que empiece antes de la posición 60 (candidato a signal anchor)
- Solo toma el primer anchor encontrado por proteína

### Influenza HA (`sp1`)

- Idéntico al parser Type 1. HA tiene un SP clivable clásico.

### Influenza NA (`sa2`)

- Idéntico al parser Signal Anchor. NA tiene un dominio TM N-terminal no clivable.

### Flavivirus (`viral_flavi`)

- Parsea features `Signal` de la poliproteína
- Usa features `Chain` para mapear cada señal a su proteína target
- Clasifica como `SIGNAL_PEPTIDE_TYPE1` (si está en pos ≤5) o `VIRAL_INTERNAL_SIGNAL`
- Identifica roles: `anchC→prM`, `E→NS1`, etc.

### Alphavirus (`viral_alpha`)

- Busca features `Signal` y los asigna a E3 o 6K según posición y Chain features
- Fallback: si no hay `Signal`, busca `Transmembrane` cerca del N-terminal
- Todos se clasifican como `VIRAL_INTERNAL_SIGNAL`

---

## Evidencia y Calidad

La clasificación sigue los códigos ECO de UniProt:

| Categoría         | Códigos ECO                                        | Significado                         |
| ----------------- | -------------------------------------------------- | ----------------------------------- |
| **EXPERIMENTAL**  | ECO:0000269, ECO:0000303, ECO:0000305              | Evidencia directa de publicación    |
| **BY_SIMILARITY** | ECO:0000250, ECO:0000255                           | Inferido por similitud de secuencia |
| **PREDICTED**     | ECO:0000256, ECO:0007829, ECO:0000259, ECO:0000312 | Anotación automática                |

El filtro `--evidence experimental` conserva solo registros con evidencia directa. Para virus, la mayoría son predicciones — usar `--evidence all` (default).

---

## Deduplicación

- Agrupa por **secuencia exacta** (`sp_sequence_aa`) + **tipo** (`sp_type`)
- Conserva el registro con **mejor evidencia** (EXPERIMENTAL > BY_SIMILARITY > PREDICTED)
- Fusiona metadata: todos los accessions compartiendo ese SP se listan en `all_accessions`
- El campo `duplicate_count` indica cuántas entradas comparten el mismo SP
- Desactivar con `--no-dedup`

---

## Features Computados

Calculados para cada SP sin herramientas externas:

- **Hidrofobicidad media** (escala Kyte-Doolittle): índice de la capacidad del SP para insertarse en la membrana
- **Carga neta a pH 7**: basada en residuos cargados (K, R, H+, D-, E-)
- **Regiones n/h/c**: estimación heurística basada en perfil de hidrofobicidad con ventana deslizante

---

## Estructura del Proyecto

```
sps_pipeline/
├── README.md                          ← Este archivo
├── requirements.txt                   ← Dependencias Python
├── .gitignore
├── config/
│   └── presets.yaml                   ← Presets de queries
├── src/
│   ├── __init__.py
│   ├── __main__.py                    ← python -m src
│   ├── main.py                        ← CLI y pipeline principal
│   ├── api/
│   │   └── uniprot_client.py          ← Cliente REST UniProt con paginación
│   ├── parsers/
│   │   ├── signal_peptide.py          ← Parser SP Type 1
│   │   ├── signal_anchor.py           ← Parser Signal Anchor Type 2
│   │   └── viral_signals.py           ← Parser señales virales
│   ├── processing/
│   │   ├── deduplication.py           ← Dedup por secuencia exacta
│   │   ├── evidence.py                ← Clasificación ECO
│   │   └── features.py                ← Features computados
│   └── output/
│       ├── writers.py                 ← CSV, TSV, FASTA
│       └── run_report.py             ← Reporte JSON de reproducibilidad
├── output/                            ← Archivos generados (gitignored)
└── logs/                              ← Logs de ejecución (gitignored)
```

---

## Extensibilidad

### Agregar un nuevo preset

Editar `config/presets.yaml` y agregar una entrada con:

```yaml
mi_preset:
  description: "Descripción"
  query: "(taxonomy_id:XXXX) AND (keyword:KW-YYYY)"
  parse_mode: "sp1" # sp1, sa2, viral_flavi, viral_alpha, auto
  biology_notes: "Notas biológicas"
```

### Agregar un nuevo parser

1. Crear archivo en `src/parsers/`
2. Implementar función que reciba `(entry: dict, query_group: str) -> list[dict]`
3. Registrar el parse_mode en `src/main.py` → `PARSE_FUNCTIONS`

### Tail-Anchored Proteins (futuro)

Las TA proteins tienen un dominio TM en el C-terminal (últimos ~50 aa) y usan la vía GET/TRC40. Para soportarlas:

1. Nuevo parser `tail_anchor.py` que busque `Transmembrane` en los últimos 50 residuos
2. Nuevo preset con keyword KW-0812
3. Nuevo `sp_type`: `TAIL_ANCHOR`
