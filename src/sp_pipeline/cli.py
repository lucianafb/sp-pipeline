"""Command-line interface for SP-Pipeline.

Usage:
    sp-pipeline query --preset human_type1 --output results.csv
    sp-pipeline query --preset human_type1 human_type2 --evidence experimental
    sp-pipeline query --preset influenza_HA --mode exhaustive
    sp-pipeline custom --organism "Homo sapiens" --sp-type type1_cleaved
    sp-pipeline presets
    sp-pipeline check
"""

from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from sp_pipeline.pipeline import SPPipeline

app = typer.Typer(
    name="sp-pipeline",
    help="Automated pipeline for querying signal peptide data from curated databases.",
    add_completion=False,
)
console = Console()


@app.command()
def query(
    preset: list[str] = typer.Option(
        ...,
        "--preset",
        "-p",
        help="Preset query name(s). Use 'sp-pipeline presets' to see available options.",
    ),
    evidence: Optional[str] = typer.Option(
        None,
        "--evidence",
        "-e",
        help="Evidence filter: experimental, predicted, curated_annotation, all",
    ),
    mode: Optional[str] = typer.Option(
        None,
        "--mode",
        "-m",
        help="Query mode: exhaustive or representative (for viral queries)",
    ),
    include_predictions: bool = typer.Option(
        False,
        "--include-predictions",
        help="Run SignalP/Phobius de novo predictions",
    ),
    output: str = typer.Option(
        "output/results.csv",
        "--output",
        "-o",
        help="Output CSV file path",
    ),
    append: bool = typer.Option(
        False,
        "--append",
        "-a",
        help="Append to existing CSV without duplicating",
    ),
    fasta: bool = typer.Option(
        False,
        "--fasta",
        help="Also export signal peptide sequences as FASTA",
    ),
    config: Optional[str] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to custom config file",
    ),
):
    """Run a query using one or more presets."""
    console.print(f"[bold green]SP-Pipeline[/bold green] - Running presets: {', '.join(preset)}")

    pipeline = SPPipeline(config_path=config)
    output_path = pipeline.run(
        presets=preset,
        evidence=evidence,
        mode=mode,
        include_predictions=include_predictions,
        output=output,
        append=append,
        include_fasta=fasta,
    )

    if output_path:
        console.print(f"\n[bold green]✓[/bold green] Results saved to: {output_path}")
    else:
        console.print("\n[bold red]✗[/bold red] No results found.")


@app.command()
def custom(
    organism: Optional[str] = typer.Option(
        None,
        "--organism",
        help="Organism name (e.g., 'Homo sapiens')",
    ),
    taxonomy_id: Optional[int] = typer.Option(
        None,
        "--taxonomy-id",
        help="NCBI Taxonomy ID",
    ),
    sp_type: str = typer.Option(
        "all",
        "--sp-type",
        help="SP type: type1_cleaved, type2_signal_anchor, internal_signal, all",
    ),
    evidence: str = typer.Option(
        "all",
        "--evidence",
        "-e",
        help="Evidence filter: experimental, predicted, curated_annotation, all",
    ),
    sources: Optional[str] = typer.Option(
        None,
        "--sources",
        "-s",
        help="Comma-separated sources: uniprot,ncbi",
    ),
    output: str = typer.Option(
        "output/custom_results.csv",
        "--output",
        "-o",
        help="Output CSV file path",
    ),
    config: Optional[str] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to custom config file",
    ),
):
    """Run a custom query with specific parameters."""
    console.print("[bold green]SP-Pipeline[/bold green] - Custom query")

    source_list = sources.split(",") if sources else None
    pipeline = SPPipeline(config_path=config)
    output_path = pipeline.run_custom(
        organism=organism,
        taxonomy_id=taxonomy_id,
        sp_type=sp_type,
        evidence=evidence,
        sources=source_list,
        output=output,
    )

    if output_path:
        console.print(f"\n[bold green]✓[/bold green] Results saved to: {output_path}")
    else:
        console.print("\n[bold red]✗[/bold red] No results found.")


@app.command()
def presets(
    config: Optional[str] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to custom config file",
    ),
):
    """List all available query presets."""
    pipeline = SPPipeline(config_path=config)
    available = pipeline.list_available_presets()

    table = Table(title="Available Presets")
    table.add_column("Name", style="cyan", no_wrap=True)
    table.add_column("Description", style="white")

    for p in available:
        table.add_row(p["name"], p["description"])

    console.print(table)


@app.command()
def check(
    config: Optional[str] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to custom config file",
    ),
):
    """Check availability of data sources."""
    console.print("[bold]Checking data sources...[/bold]\n")

    pipeline = SPPipeline(config_path=config)
    status = pipeline.check_sources()

    for source, available in status.items():
        icon = "[green]✓[/green]" if available else "[red]✗[/red]"
        console.print(f"  {icon} {source}")

    console.print()


@app.command()
def clear_cache(
    config: Optional[str] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to custom config file",
    ),
):
    """Clear the query cache."""
    pipeline = SPPipeline(config_path=config)
    pipeline.cache.clear()
    console.print("[bold green]✓[/bold green] Cache cleared.")


def main():
    app()


if __name__ == "__main__":
    main()
