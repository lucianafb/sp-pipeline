"""
Sec61-Target-Scout: CLI Entry Point.

Pipeline para extracción y curación de péptidos señal desde UniProt.

Uso:
    python -m src.main --preset human_sp1 --output csv
    python -m src.main --preset human_sp1 human_sa2 influenza_ha --output csv
    python -m src.main --preset all --output csv fasta
    python -m src.main --query "(gene:INS) AND (organism_id:9606)" --output csv
"""

import os
import sys
import argparse
import logging
from datetime import datetime
from pathlib import Path

import yaml
import pandas as pd

from .api.uniprot_client import fetch_entries
from .parsers.signal_peptide import parse_signal_peptides
from .parsers.signal_anchor import parse_signal_anchors
from .parsers.viral_signals import parse_flavivirus_signals, parse_alphavirus_signals
from .parsers.bunyavirus_signals import parse_bunyavirus_signals
from .processing.deduplication import deduplicate_by_sequence
from .processing.evidence import get_evidence_priority
from .output.writers import write_output
from .output.run_report import RunReport

# ---------------------------------------------------------------------------
# Configuración de logging
# ---------------------------------------------------------------------------

def setup_logging(log_dir: str, verbose: bool = False):
    """Configura logging a archivo y consola."""
    os.makedirs(log_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"run_{timestamp}.log")

    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.FileHandler(log_file, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )

    return log_file


# ---------------------------------------------------------------------------
# Carga de presets
# ---------------------------------------------------------------------------

def load_presets(config_path: str) -> dict:
    """Carga presets de queries desde YAML."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_project_root() -> Path:
    """Retorna la raíz del proyecto (donde está src/)."""
    return Path(__file__).parent.parent


# ---------------------------------------------------------------------------
# Mapeo de parse_mode a funciones
# ---------------------------------------------------------------------------

PARSE_FUNCTIONS = {
    "sp1": parse_signal_peptides,
    "sa2": parse_signal_anchors,
    "viral_flavi": parse_flavivirus_signals,
    "viral_alpha": parse_alphavirus_signals,
    "viral_bunya": parse_bunyavirus_signals,
}


def parse_entry(entry: dict, parse_mode: str, query_group: str) -> list[dict]:
    """
    Parsea una entrada según el modo indicado.

    Para mode 'auto', prueba todos los parsers.
    """
    if parse_mode == "auto":
        # Probar SP1, SA2, y ambos virales
        results = []
        results.extend(parse_signal_peptides(entry, query_group))
        if not results:
            results.extend(parse_signal_anchors(entry, query_group))
        return results

    parser_fn = PARSE_FUNCTIONS.get(parse_mode)
    if not parser_fn:
        logging.getLogger(__name__).warning(f"parse_mode desconocido: {parse_mode}")
        return []

    return parser_fn(entry, query_group)


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def run_pipeline(
    presets_config: dict,
    preset_names: list[str],
    custom_query: str = None,
    evidence_filter: str = "all",
    include_unreviewed: bool = False,
    max_results: int = None,
    no_dedup: bool = False,
    output_formats: list[str] = None,
    output_dir: str = "output",
    output_name: str = None,
):
    """
    Ejecuta el pipeline completo.

    1. Para cada preset (o custom query), consulta UniProt
    2. Parsea los resultados según el parse_mode
    3. Concatena todos los resultados
    4. Filtra por evidencia si corresponde
    5. Deduplica
    6. Genera outputs
    """
    if output_formats is None:
        output_formats = ["csv"]

    report = RunReport()
    report.presets_used = preset_names if not custom_query else ["custom"]
    report.evidence_filter = evidence_filter
    report.include_unreviewed = include_unreviewed
    report.deduplication_enabled = not no_dedup

    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    all_records = []

    # -------------------------------------------------------------------------
    # Fase 1: Queries y parsing
    # -------------------------------------------------------------------------

    queries_to_run = []

    if custom_query:
        queries_to_run.append({
            "name": "custom",
            "query": custom_query,
            "parse_mode": "auto",
        })
    else:
        for preset_name in preset_names:
            if preset_name not in presets_config:
                print(f"⚠️  Preset desconocido: '{preset_name}'. Disponibles: {list(presets_config.keys())}")
                report.add_error(f"Preset desconocido: {preset_name}")
                continue

            preset = presets_config[preset_name]
            query = preset["query"]

            # Si no incluye unreviewed y la query no tiene reviewed:
            # Presets con include_unreviewed: true no se filtran por reviewed
            preset_includes_unreviewed = preset.get("include_unreviewed", False)
            if not include_unreviewed and not preset_includes_unreviewed and "reviewed:" not in query:
                query += " AND (reviewed:true)"

            queries_to_run.append({
                "name": preset_name,
                "query": query,
                "parse_mode": preset.get("parse_mode", "auto"),
            })

    for query_info in queries_to_run:
        preset_name = query_info["name"]
        query = query_info["query"]
        parse_mode = query_info["parse_mode"]

        print(f"\n🔬 Ejecutando preset: {preset_name}")
        print(f"   Query: {query}")
        print(f"   Parse mode: {parse_mode}")

        api_count = 0
        parsed_count = 0

        try:
            for entry in fetch_entries(query, max_results=max_results):
                api_count += 1
                records = parse_entry(entry, parse_mode, preset_name)
                parsed_count += len(records)
                all_records.extend(records)

        except Exception as e:
            error_msg = f"Error en preset '{preset_name}': {e}"
            print(f"❌ {error_msg}")
            report.add_error(error_msg)
            logging.getLogger(__name__).exception(error_msg)

        report.add_query(preset_name, query, api_count, parsed_count)
        print(f"   ✅ API: {api_count} entradas | Parseados: {parsed_count} SPs")

    # -------------------------------------------------------------------------
    # Fase 2: Crear DataFrame
    # -------------------------------------------------------------------------

    if not all_records:
        print("\n⚠️  No se encontraron resultados.")
        report.print_summary()
        return

    df = pd.DataFrame(all_records)
    print(f"\n📊 Total registros pre-procesamiento: {len(df)}")

    # -------------------------------------------------------------------------
    # Fase 3: Filtro de evidencia
    # -------------------------------------------------------------------------

    if evidence_filter != "all":
        before = len(df)
        if evidence_filter == "experimental":
            df = df[df["evidence_category"] == "EXPERIMENTAL"]
        elif evidence_filter == "by_similarity":
            df = df[df["evidence_category"].isin(["EXPERIMENTAL", "BY_SIMILARITY"])]

        print(f"🔍 Filtro de evidencia '{evidence_filter}': {before} → {len(df)}")

    # -------------------------------------------------------------------------
    # Fase 4: Deduplicación
    # -------------------------------------------------------------------------

    if not no_dedup and not df.empty:
        df = deduplicate_by_sequence(df)
        report.set_dedup_results(len(df))
    else:
        report.set_dedup_results(len(df))

    # -------------------------------------------------------------------------
    # Fase 5: Generar outputs
    # -------------------------------------------------------------------------

    if df.empty:
        print("\n⚠️  No hay resultados después del procesamiento.")
        report.print_summary()
        return

    # Nombre base del archivo
    if output_name:
        base_name = output_name
    elif custom_query:
        base_name = f"custom_{timestamp}"
    elif len(preset_names) == 1:
        base_name = f"{preset_names[0]}_{timestamp}"
    else:
        base_name = f"multi_preset_{timestamp}"

    base_path = os.path.join(output_dir, base_name)

    generated_files = write_output(df, base_path, output_formats)

    for f in generated_files:
        report.add_output_file(f)

    # Guardar reporte
    report_path = os.path.join(output_dir, f"{base_name}_report.json")
    report.save(report_path)
    report.add_output_file(report_path)

    report.print_summary()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Construye el parser de argumentos CLI."""
    parser = argparse.ArgumentParser(
        prog="sec61-target-scout",
        description="🧬 Sec61-Target-Scout: Pipeline de extracción de péptidos señal desde UniProt",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  # Un preset
  python -m src.main --preset human_sp1 --output csv

  # Múltiples presets
  python -m src.main --preset human_sp1 human_sa2 influenza_ha --output csv

  # Todos los presets
  python -m src.main --preset all --output csv fasta

  # Query personalizada
  python -m src.main --query "(gene:INS) AND (organism_id:9606)" --output csv

  # Solo evidencia experimental
  python -m src.main --preset human_sp1 --evidence experimental

  # Limitar resultados (para testing)
  python -m src.main --preset human_sp1 --max-results 100

Presets disponibles:
  human_sp1       - SPs clivables humanos (Swiss-Prot)
  human_sp1_all   - SPs clivables humanos (Swiss-Prot + TrEMBL)
  human_sa2       - Signal Anchors humanos (Swiss-Prot)
  human_sa2_all   - Signal Anchors humanos (Swiss-Prot + TrEMBL)
  influenza_ha    - SP de Hemaglutinina (Influenza A)
  influenza_na    - Signal Anchor de Neuraminidasa (Influenza A)
  flavivirus      - Señales de poliproteína (Dengue, Zika, etc.)
  alphavirus      - Señales E3/6K (Chikungunya, Sindbis, etc.)
""",
    )

    # Fuente de datos (not required so --list-presets works standalone)
    source = parser.add_mutually_exclusive_group(required=False)
    source.add_argument(
        "--preset",
        nargs="+",
        metavar="NAME",
        help="Preset(s) a ejecutar. Usar 'all' para todos. Se pueden combinar múltiples.",
    )
    source.add_argument(
        "--query",
        type=str,
        help="Query personalizada para UniProt.",
    )

    # Filtros
    parser.add_argument(
        "--evidence",
        type=str,
        choices=["all", "experimental", "by_similarity"],
        default="all",
        help="Filtrar por categoría de evidencia (default: all).",
    )
    parser.add_argument(
        "--include-unreviewed",
        action="store_true",
        help="Incluir entradas TrEMBL (no revisadas). Por defecto solo Swiss-Prot.",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=None,
        help="Limitar resultados por preset (para testing). Sin límite por defecto.",
    )

    # Output
    parser.add_argument(
        "--output",
        nargs="+",
        choices=["csv", "tsv", "fasta"],
        default=["csv"],
        help="Formato(s) de salida (default: csv).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Directorio de salida (default: output/).",
    )
    parser.add_argument(
        "--name",
        type=str,
        default=None,
        help="Alias para el archivo de salida (ej: --name mi_experimento → mi_experimento.csv). Si no se indica, se genera automáticamente.",
    )

    # Procesamiento
    parser.add_argument(
        "--no-dedup",
        action="store_true",
        help="Desactivar deduplicación por secuencia.",
    )

    # Misc
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Activar logging detallado (DEBUG).",
    )
    parser.add_argument(
        "--list-presets",
        action="store_true",
        help="Listar presets disponibles y salir.",
    )

    return parser


def main():
    """Punto de entrada principal."""
    parser = build_parser()
    args = parser.parse_args()

    # Resolver paths
    project_root = get_project_root()
    config_path = project_root / "config" / "presets.yaml"

    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = str(project_root / "output")

    log_dir = str(project_root / "logs")

    # Setup logging
    log_file = setup_logging(log_dir, verbose=args.verbose)

    # Cargar presets
    if not config_path.exists():
        print(f"❌ No se encontró config de presets: {config_path}")
        sys.exit(1)

    presets_config = load_presets(str(config_path))

    # Listar presets
    if args.list_presets:
        print("\n📋 Presets disponibles:\n")
        for name, preset in presets_config.items():
            desc = preset.get("description", "Sin descripción")
            print(f"  {name:20s} → {desc}")
        print()
        sys.exit(0)

    # Validate that --preset or --query was provided
    if not args.preset and not args.query:
        parser.error("one of the arguments --preset --query is required")

    # Resolver presets
    preset_names = []
    if args.preset:
        if "all" in args.preset:
            preset_names = list(presets_config.keys())
            print(f"🔬 Ejecutando TODOS los presets: {preset_names}")
        else:
            preset_names = args.preset

    # Banner
    print("\n" + "=" * 60)
    print("🧬 Sec61-Target-Scout v1.0.0")
    print("   Pipeline de Extracción de Péptidos Señal")
    print("=" * 60)

    if args.preset:
        print(f"  Presets:           {', '.join(preset_names)}")
    else:
        print(f"  Custom query:      {args.query}")
    print(f"  Evidencia:         {args.evidence}")
    print(f"  Include TrEMBL:    {args.include_unreviewed}")
    print(f"  Deduplicación:     {'No' if args.no_dedup else 'Sí'}")
    print(f"  Formatos output:   {', '.join(args.output)}")
    print(f"  Output dir:        {output_dir}")
    print(f"  Log:               {log_file}")

    if args.max_results:
        print(f"  Max results:       {args.max_results} (modo testing)")

    print("=" * 60)

    # Ejecutar pipeline
    run_pipeline(
        presets_config=presets_config,
        preset_names=preset_names,
        custom_query=args.query,
        evidence_filter=args.evidence,
        include_unreviewed=args.include_unreviewed,
        max_results=args.max_results,
        no_dedup=args.no_dedup,
        output_formats=args.output,
        output_dir=output_dir,
        output_name=args.name,
    )


if __name__ == "__main__":
    main()
