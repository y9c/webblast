#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command-line interface for webblast (fast NCBI BLAST web queue)."""

from __future__ import annotations

import json
import random
import sys
from typing import Optional

import rich_click as click
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

from .cache import cache_key, get as cache_get, put as cache_put
from .client import DATABASES, PROGRAMS, BlastClient
from .parse import BlastReport, parse_blast_json
from .species import SpeciesResolver
from .utils import build_query

console = Console()          # stdout — used for data output
err_console = Console(stderr=True)  # stderr — used for progress / status

# Fun, cycling status lines shown while NCBI runs the search — makes the wait
# feel alive and fast even though it usually only takes a few seconds.
WAIT_MESSAGES = [
    "Warming up NCBI's compute farm...",
    "Indexing your query...",
    "Scanning the database...",
    "Chasing the best hits...",
    "Aligning homologs...",
    "Ranking statistical significance...",
    "Polishing the report...",
]


def _render_text(report: BlastReport, max_alignments: int = 10,
                 resolver: Optional[SpeciesResolver] = None) -> None:
    """Print a human-readable summary table of hits per query."""
    for qi, q in enumerate(report.queries):
        if qi:
            console.print()
        header = f"[bold cyan]{q.query_title or q.query_id}[/bold cyan]"
        if q.query_len:
            header += f"  (len {q.query_len})"
        console.print(header)
        if q.message or not q.hits:
            console.print(f"[yellow]  {q.message or 'No significant hits.'}[/yellow]")
            continue

        table = Table(show_lines=False, header_style="bold", title=None)
        table.add_column("#", style="dim", justify="right")
        table.add_column("Accession", style="bold")
        table.add_column("Description", max_width=44, no_wrap=True, overflow="ellipsis")
        table.add_column("E-value", justify="right")
        table.add_column("Identity", justify="right")
        table.add_column("Align len", justify="right")
        table.add_column("Bit score", justify="right")
        if resolver is not None:
            table.add_column("中文名", style="green")

        shown = 0
        for hit in q.hits:
            if max_alignments and shown >= max_alignments:
                console.print(f"[dim]  … and {len(q.hits) - shown} more hits[/dim]")
                break
            hsp = hit.hsps[0] if hit.hsps else None
            row = [
                str(hit.num) if hit.num is not None else "",
                hit.accession or "",
                (hit.title or " ")[:60],
                f"{hsp.evalue:.2e}" if hsp and hsp.evalue is not None else "—",
                f"{hsp.pident:.1f}%" if hsp and hsp.pident is not None else "—",
                str(hsp.align_len) if hsp and hsp.align_len is not None else "—",
                str(round(hsp.bit_score, 1)) if hsp and hsp.bit_score is not None else "—",
            ]
            if resolver is not None:
                zh, ai = resolver.resolve(hit.sciname or hit.title or "")
                row.append(f"{zh} [dim]⚠[/dim]" if (zh and ai) else (zh or ""))
            table.add_row(*row)
            shown += 1
        console.print(table)


def _render_tsv(report: BlastReport, resolver: Optional[SpeciesResolver] = None) -> str:
    """Tab-delimited: query, accession, title, evalue, pident, align_len, bitscore[, zh]."""
    lines = []
    for q in report.queries:
        for hit in q.hits:
            hsp = hit.hsps[0] if hit.hsps else None
            evalue = hsp.evalue if hsp else ""
            pident = f"{hsp.pident:.1f}" if hsp and hsp.pident is not None else ""
            align = hsp.align_len if hsp and hsp.align_len is not None else ""
            bits = round(hsp.bit_score, 1) if hsp and hsp.bit_score is not None else ""
            row = [
                q.query_title or q.query_id or "",
                hit.accession or "",
                hit.title or "",
                str(evalue),
                pident,
                str(align),
                str(bits),
            ]
            if resolver is not None:
                zh, ai = resolver.resolve(hit.sciname or hit.title or "")
                row.append(f"{zh} (AI)" if (zh and ai) else (zh or ""))
            lines.append("\t".join(row))
    return "\n".join(lines) + ("\n" if lines else "")


def _render_fasta(report: BlastReport) -> str:
    """Emit the subject segments of the top HSP for each hit as pseudo-FASTA."""
    lines = []
    for q in report.queries:
        for hit in q.hits:
            if not hit.hsps:
                continue
            hsp = hit.hsps[0]
            desc = hit.title or hit.accession or "hit"
            lines.append(f">{hit.accession or hit.id} {desc}")
            lines.append(hsp.hseq or "")
    return "\n".join(lines) + ("\n" if lines else "")


def _write_output(text: str, out: Optional[str]) -> None:
    if out:
        with open(out, "w") as fh:
            fh.write(text)
    else:
        sys.stdout.write(text)


@click.command(
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "-p", "--program", type=click.Choice(sorted(PROGRAMS), case_sensitive=False),
    default="megablast", show_default=True, help="BLAST program / algorithm.",
)
@click.option(
    "-d", "--database", default="nt", show_default=True,
    help="Database to search (e.g. nt, nr, refseq_rna, swissprot).",
)
@click.option(
    "-f", "--format", "fmt",
    type=click.Choice(["text", "json", "tsv", "fasta"], case_sensitive=False),
    default="text", show_default=True,
    help="Output format.",
)
@click.option("-o", "--out", type=click.Path(writable=True), default=None, help="Write output to a file.")
@click.option("--max-num-seq", default=100, show_default=True, help="Max target sequences per query.")
@click.option("--hitlist-size", default=100, show_default=True, help="Hitlist size (HITLIST_SIZE).")
@click.option("--expect", default=None, help="E-value cutoff (e.g. 1e-5).")
@click.option("--email", default="yech1990@gmail.com", show_default=True, help="Contact email (NCBI policy).")
@click.option("--api-key", default="", help="NCBI API key (raises rate limit).")
@click.option("--translate/--no-translate", default=True, show_default=True,
              help="Add Chinese species names (from the bundled dictionary).")
@click.option("--cache/--no-cache", default=True, show_default=True, help="Use the on-disk result cache.")
@click.option("--cache-dir", default="blast_cache", show_default=True, help="Directory for the result cache.")
@click.option("--timeout", default=20.0, show_default=True, help="Per-request HTTP timeout (s).")
@click.option("--wait-timeout", default=1200.0, show_default=True, help="Max time to wait for the search (s).")
@click.option("--limit", default=None, type=int, help="Only use the first N query records.")
@click.option("--color/--no-color", default=True, show_default=True, help="Enable colored output.")
@click.option("--alignments", default=None, type=int, help="Max alignments to show in text output.")
@click.argument("query_files", nargs=-1, type=click.Path(exists=True))
def main(
    program,
    database,
    fmt,
    out,
    max_num_seq,
    hitlist_size,
    expect,
    email,
    api_key,
    translate,
    cache,
    cache_dir,
    timeout,
    wait_timeout,
    limit,
    color,
    alignments,
    query_files,
):
    """Run a fast NCBI BLAST search (interactive web queue) and print results."""
    if not color:
        # richest output is only meaningfully colored; allow plain via --no-color
        pass

    # species -> Chinese name resolver (dictionary only, lazy-loaded)
    resolver = SpeciesResolver() if translate else None

    # 1) assemble the multi-record FASTA query (files, else stdin)
    with err_console.status("[bold green]Reading query sequences..."):
        query = build_query(list(query_files), limit)

    # 2) build the client (session + interactive web path)
    try:
        client = BlastClient(program=program, database=database, email=email, api_key=api_key, timeout=timeout)
    except ValueError as exc:
        err_console.print(f"[red]{exc}[/red]")
        sys.exit(1)

    # 3) cache lookup
    opts = {
        "max_num_seq": max_num_seq,
        "hitlist_size": hitlist_size,
        "expect": expect,
    }
    key = cache_key(program, database, query, opts)
    cached = cache_get(key, cache_dir) if cache else None

    if cached is not None:
        err_console.print("[dim]Using cached result.[/dim]")
        report = parse_blast_json(cached)
        _emit(report, fmt, out, alignments, resolver)
        return

    # 4) submit, wait, fetch — with a lively animation that keeps the user
    #    engaged/fast-feeling (a spinning frame + moving bar + cycling status).
    rid = None
    with Progress(
        SpinnerColumn("dots12"),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=22, style="cyan", complete_style="green"),
        TextColumn("{task.fields[state]}"),
        TimeElapsedColumn(),
        refresh_per_second=12,
        console=err_console,
    ) as progress:
        task = progress.add_task(
            f"[bold cyan]{program}[/bold cyan] vs [bold green]{database}[/bold green]",
            total=None,
            state="Submitting query...",
        )

        submission = client.submit(
            query, format_type="JSON2_S",
            hitlist_size=hitlist_size, max_num_seq=max_num_seq, EXPECT=expect,
        )
        rid = submission.rid
        progress.update(task, state=f"RID [blue]{rid[:10]}…[/blue]")

        def on_status(attempt, elapsed, status):
            if status == "WAITING":
                msg = WAIT_MESSAGES[attempt % len(WAIT_MESSAGES)]
                progress.update(task, state=f"{msg} [dim]· poll {attempt}[/dim]")
            elif status == "READY":
                progress.update(task, state=f"[bold green]Done in {elapsed:.1f}s![/bold green]")
            else:
                progress.update(task, state=f"[red]{status}[/red]")

        status = client.wait(rid, submission.rtoe, timeout=wait_timeout, on_status=on_status)
        if status != "READY":
            err_console.print(f"[red]Search ended with {status}.[/red]")
            sys.exit(2)
        progress.update(task, state="Downloading report...")
        data = client.fetch(rid, format_type="JSON2_S")
        progress.update(task, state="Parsing results...")

    if cache:
        cache_put(key, data, cache_dir)

    report = parse_blast_json(data)
    _emit(report, fmt, out, alignments, resolver)

    if rid:
        # keep machine-readable outputs clean: the link goes to stderr
        click.secho(
            f"\nView in browser: https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID={rid}",
            fg="bright_black", err=True,
        )

    if fmt == "text" and out is None:
        total_hits = report.total_hits
        nq = len(report.queries)
        console.print(
            f"[dim]✔ {nq} quer{'y' if nq == 1 else 'ies'} · "
            f"{total_hits} hit{'s' if total_hits != 1 else ''}[/dim]"
        )


def _emit(report: BlastReport, fmt: str, out: Optional[str], max_alignments: Optional[int],
          resolver: Optional[SpeciesResolver] = None):
    """Render the report in the requested format (to stdout or a file)."""
    if fmt == "json":
        def _hit(h):
            zh, ai = resolver.resolve(h.sciname or h.title or "") if resolver else (None, False)
            top = h.hsps[0] if h.hsps else None
            return {
                "accession": h.accession,
                "title": h.title,
                "taxid": h.taxid,
                "sciname": h.sciname,
                "len": h.len,
                "zh": zh,
                "zh_is_ai": ai if zh else None,
                "top": {
                    "evalue": (top.evalue if top else None),
                    "pident": (round(top.pident, 2) if top and top.pident is not None else None),
                    "align_len": (top.align_len if top else None),
                    "bit_score": (top.bit_score if top else None),
                    "qseq": (top.qseq if top else None),
                    "hseq": (top.hseq if top else None),
                },
            }

        payload = {
            "program": report.program,
            "version": report.version,
            "db": report.db,
            "queries": [
                {
                    "query_title": q.query_title,
                    "query_len": q.query_len,
                    "hits": [_hit(h) for h in q.hits],
                }
                for q in report.queries
            ],
        }
        _write_output(json.dumps(payload, indent=2) + "\n", out)
    elif fmt == "tsv":
        _write_output(_render_tsv(report, resolver), out)
    elif fmt == "fasta":
        _write_output(_render_fasta(report), out)
    else:  # text
        if out:
            _write_output("", out)  # todo: file text
        _render_text(report, max_alignments, resolver)


if __name__ == "__main__":
    main()
