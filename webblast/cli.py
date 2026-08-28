#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command-line interface for webblast (fast NCBI BLAST web queue)."""

from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional, Tuple

import rich_click as click
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

from .cache import cache_key, get as cache_get, put as cache_put
from .client import PROGRAMS, BlastClient
from .parse import BlastReport, parse_blast_json
from .species import SpeciesResolver
from .utils import chunk_records, join_records, read_records

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


def _fmt_evalue(e: Optional[float]) -> str:
    if e is None:
        return "—"
    if e == 0:
        return "0"
    if e < 1e-200:
        return "<1e-200"
    if e < 1e-3 or e >= 1e4:
        return f"{e:.1e}"
    return f"{e:.3g}"


def _ident_style(p: Optional[float]) -> str:
    if p is None:
        return "dim"
    if p >= 90:
        return "bold green"
    if p >= 70:
        return "bold yellow"
    return "red"


def _zh(hit, resolver):
    if resolver is None:
        return None, False
    return resolver.resolve(hit.sciname or hit.title or "")


def _badge_zh(hit, resolver) -> str:
    zh, ai = _zh(hit, resolver)
    if not zh:
        return ""
    return f"[green]{zh}[/green]" + (" [dim]⚠[/dim]" if ai else "")


def _best_hit(hit, hsp, resolver) -> None:
    """Highlight the top match: accession, species, 中文名, description, stats."""
    zh, ai = _zh(hit, resolver)
    label = f"[bold green]◆ {hit.accession or '?'}[/bold green]"
    if hit.sciname:
        label += f"  [bold]{hit.sciname}[/bold]"
    if zh:
        label += f"  [green]{zh}[/green]" + (" [dim]⚠[/dim]" if ai else "")
    console.print(label)
    if hit.title:
        console.print(f"   [dim]{hit.title}[/dim]")
    if hsp:
        pident = hsp.pident
        imp = (f"[{_ident_style(pident)}]{pident:.1f}%[/] identity"
               if pident is not None else "identity —")
        console.print(
            f"   [bold]E-value[/bold] {_fmt_evalue(hsp.evalue)}"
            f" · [bold]bits[/bold] {hsp.bit_score:.1f}"
            f" · {imp}"
            f" · [bold]{hsp.align_len} aa[/bold]"
            f" · query {hsp.query_from}–{hsp.query_to}"
            f" · subject {hsp.hit_from}–{hsp.hit_to}"
        )


def _render_alignment(hit, hsp, resolver) -> None:
    """Print the pairwise alignment (query / mid / subject) for an HSP."""
    if not hsp or not hsp.qseq:
        return
    zh, _ = _zh(hit, resolver)
    title = f"   alignment vs [bold]{hit.accession or '?'}[/bold]"
    if zh:
        title += f" [green]({zh})[/green]"
    console.print(f"[dim]{title}[/dim]")
    width = 60
    q, m, s = hsp.qseq, hsp.midline, hsp.hseq
    qstart = hsp.query_from or 1
    sstart = hsp.hit_from or 1
    for i in range(0, len(q), width):
        cq = q[i : i + width]
        cm = m[i : i + width]
        cs = s[i : i + width]
        qpos = qstart + i
        spos = sstart + i
        console.print(f"   [cyan]{cq}[/cyan] [dim]{qpos}[/dim]")
        console.print(f"   {cm}")
        console.print(f"   [green]{cs}[/green] [dim]{spos}[/dim]")


def _render_text(report: BlastReport, top: int = 15, show_align: bool = False,
                 resolver: Optional[SpeciesResolver] = None) -> None:
    """Rich, human-friendly per-query output: best-hit highlight + compact table."""
    for qi, q in enumerate(report.queries):
        if qi:
            console.print()
        header = f"[bold cyan]Query:[/bold cyan] [bold]{q.query_title or q.query_id}[/bold]"
        if q.query_len:
            header += f" [dim]({q.query_len} residues)[/dim]"
        console.print(header)

        if q.message or not q.hits:
            console.print(f"[yellow]  {q.message or 'No significant hits.'}[/yellow]")
            continue

        hits = q.hits
        total = len(hits)
        top = max(top, 1)
        shown = hits[:top]

        # 1) best hit highlight
        _best_hit(shown[0], shown[0].hsps[0] if shown[0].hsps else None, resolver)

        # 2) compact table for the rest of the shown hits
        rest = shown[1:]
        if rest:
            table = Table(show_lines=False, header_style="bold", box=None,
                          pad_edge=False)
            table.add_column("#", style="dim", justify="right")
            table.add_column("Accession", style="bold")
            table.add_column("E-value", justify="right")
            table.add_column("Ident", justify="right")
            table.add_column("Bits", justify="right")
            table.add_column("Len", justify="right")
            if resolver is not None:
                table.add_column("中文名", style="green")
            for hit in rest:
                hsp = hit.hsps[0] if hit.hsps else None
                row = [
                    str(hit.num) if hit.num is not None else "",
                    hit.accession or "",
                    _fmt_evalue(hsp.evalue if hsp else None),
                    f"[{_ident_style(hsp.pident if hsp else None)}]"
                    f"{hsp.pident:.1f}%[/]" if hsp and hsp.pident is not None else "—",
                    f"{hsp.bit_score:.1f}" if hsp and hsp.bit_score is not None else "—",
                    str(hsp.align_len) if hsp and hsp.align_len is not None else "—",
                ]
                if resolver is not None:
                    row.append(_badge_zh(hit, resolver))
                table.add_row(*row)
            console.print(table)

        # 3) footer / more
        if total > top:
            console.print(f"[dim]  … and {total - top} more hits "
                          f"(--top {total} to show all)[/dim]")
        if show_align:
            for hit in shown:
                hsp = hit.hsps[0] if hit.hsps else None
                _render_alignment(hit, hsp, resolver)
        console.print(f"[dim]  {total} hits[/dim]")


def _render_tabular(report: BlastReport) -> str:
    """Classic BLAST 'tabular' output (outfmt 6): qseqid sseqid pident length
    mismatch gapopen qstart qend sstart send evalue bitscore."""
    HEADER = ("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\t"
              "sstart\tsend\tevalue\tbitscore")
    lines = [HEADER]
    for q in report.queries:
        for hit in q.hits:
            hsp = hit.hsps[0] if hit.hsps else None
            if hsp is None:
                continue
            identity = hsp.identity or 0
            gaps = hsp.gaps or 0
            align_len = hsp.align_len or 0
            mismatch = max(0, align_len - identity - gaps)
            lines.append("\t".join([
                q.query_title or q.query_id or "",
                hit.accession or "",
                f"{hsp.pident:.2f}" if hsp.pident is not None else "",
                str(align_len),
                str(round(mismatch)),
                str(round(gaps)),
                str(hsp.query_from or ""),
                str(hsp.query_to or ""),
                str(hsp.hit_from or ""),
                str(hsp.hit_to or ""),
                f"{hsp.evalue:.3g}" if hsp.evalue is not None else "",
                f"{hsp.bit_score:.1f}" if hsp.bit_score is not None else "",
            ]))
    return "\n".join(lines) + "\n"


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


# --------------------------------------------------------------------------- #
# Batch / parallel searches
# --------------------------------------------------------------------------- #
def _search_one(
    program, database, query_text, email, api_key, max_num_seq, hitlist_size,
    expect, cache, cache_dir, timeout, wait_timeout,
) -> Tuple[BlastReport, Optional[str]]:
    """Run a single search: cache lookup -> submit -> wait -> fetch -> parse."""
    client = BlastClient(program=program, database=database, email=email,
                         api_key=api_key, timeout=timeout)
    opts = {"max_num_seq": max_num_seq, "hitlist_size": hitlist_size, "expect": expect}
    key = cache_key(program, database, query_text, opts)
    cached = cache_get(key, cache_dir) if cache else None
    if cached is not None:
        return parse_blast_json(cached), None
    sub = client.submit(query_text, format_type="JSON2_S",
                        hitlist_size=hitlist_size, max_num_seq=max_num_seq, EXPECT=expect)
    status = client.wait(sub.rid, sub.rtoe, timeout=wait_timeout)
    if status != "READY":
        raise RuntimeError(f"BLAST {sub.rid} ended with status {status}")
    data = client.fetch(sub.rid, format_type="JSON2_S")
    if cache:
        cache_put(key, data, cache_dir)
    return parse_blast_json(data), sub.rid


def _merge_reports(reports: List[BlastReport]) -> BlastReport:
    merged = BlastReport()
    if not reports:
        return merged
    merged.program = reports[0].program
    merged.version = reports[0].version
    merged.db = reports[0].db
    merged.params = reports[0].params
    for rep in reports:
        merged.queries.extend(rep.queries)
    return merged


def _run_parallel(program, database, records, jobs, email, api_key, max_num_seq,
                  hitlist_size, expect, cache, cache_dir, timeout, wait_timeout,
                  fmt, out, top, show_align, resolver) -> None:
    """Split records across ``jobs`` threads, run them concurrently, merge."""
    chunks = chunk_records(records, jobs)
    err_console.print(f"[dim]Running {len(chunks)} parallel searches "
                      f"({len(records)} records)...[/dim]")
    rids: List[str] = []
    with ThreadPoolExecutor(max_workers=len(chunks)) as ex:
        futures = [
            ex.submit(_search_one, program, database, join_records(chunk), email,
                      api_key, max_num_seq, hitlist_size, expect, cache, cache_dir,
                      timeout, wait_timeout)
            for chunk in chunks
        ]
        reports = []
        for f in futures:
            try:
                report, rid = f.result()
            except Exception as exc:
                err_console.print(f"[red]{exc}[/red]")
                sys.exit(2)
            reports.append(report)
            if rid:
                rids.append(rid)

    merged = _merge_reports(reports)
    _emit(merged, fmt, out, top, show_align, resolver)
    for rid in rids:
        click.secho(
            f"View in browser: https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID={rid}",
            fg="bright_black", err=True,
        )


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
    type=click.Choice(["text", "json", "tsv", "tabular", "fasta"], case_sensitive=False),
    default="text", show_default=True,
    help="Output format (tabular = classic BLAST outfmt 6).",
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
@click.option("--jobs", default=1, show_default=True, type=int,
              help="Run up to N searches in parallel (splits the input records). "
                   "Big speedup for batch queries.")
@click.option("--color/--no-color", default=True, show_default=True, help="Enable colored output.")
@click.option("--top", default=15, show_default=True, type=int,
              help="Max hits shown in the text view (search still returns --max-num-seq).")
@click.option("--show-align", is_flag=True, default=False,
              help="Also print the pairwise alignment for the shown hits.")
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
    top,
    show_align,
    jobs,
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
        records = read_records(list(query_files), limit)
    query = join_records(records)

    # 1b) parallel batch mode — split records across --jobs concurrent searches
    if jobs > 1 and len(records) > 1:
        _run_parallel(
            program, database, records, jobs, email, api_key, max_num_seq,
            hitlist_size, expect, cache, cache_dir, timeout, wait_timeout,
            fmt, out, top, show_align, resolver,
        )
        return

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
        _emit(report, fmt, out, top, show_align, resolver)
        return

    # 4) submit, wait, fetch — with a lively animation that keeps the user
    #    engaged/fast-feeling (a spinning frame + moving bar + cycling status).
    rid = None
    data = None
    try:
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
    except KeyboardInterrupt:
        err_console.print("\n[bold red]Interrupted.[/bold red] [dim]No result was cached.[/dim]")
        sys.exit(130)

    if cache and data is not None:
        cache_put(key, data, cache_dir)

    report = parse_blast_json(data)
    _emit(report, fmt, out, top, show_align, resolver)

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


def _emit(report: BlastReport, fmt: str, out: Optional[str], top: int = 15,
          show_align: bool = False, resolver: Optional[SpeciesResolver] = None):
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
    elif fmt == "tabular":
        _write_output(_render_tabular(report), out)
    elif fmt == "fasta":
        _write_output(_render_fasta(report), out)
    else:  # text
        if out:
            _write_output("", out)  # todo: file text
        _render_text(report, top, show_align, resolver)


if __name__ == "__main__":
    main()
