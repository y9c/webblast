#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-06-29 01:08


import os
import sys
import time
import urllib.parse
from hashlib import md5

import requests
import rich_click as click
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

# Define available BLAST programs and databases
BLAST_PROGRAMS = {
    "blastn": "BLAST nucleotide",
    "blastp": "BLAST protein",
    "blastx": "Translated nucleotide to protein",
    "tblastn": "Protein query to translated database",
    "tblastx": "Translated nucleotide to translated database",
    "megablast": "Optimized BLAST nucleotide",
    "rpsblast": "Reverse Position Specific BLAST",
}

BLAST_DATABASES = [
    "nt",
    "nr",
    "refseq_rna",
    "refseq_protein",
    "swissprot",
    "pdbaa",
    "pdbnt",
]

CACHE_DIR = "blast_cache"


def get_cache_path(input_hash):
    return os.path.join(CACHE_DIR, f"{input_hash}.cache")


def check_cache(input_hash):
    """Check if the result is already cached and return path if available."""
    cache_path = get_cache_path(input_hash)
    if os.path.exists(cache_path):
        return cache_path
    return None


def write_cache(input_hash, data):
    """Write data to cache."""
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)
    cache_path = get_cache_path(input_hash)
    with open(cache_path, "w") as file:
        file.write(data)


@click.command()
@click.option(
    "--program",
    "-p",
    type=click.Choice(list(BLAST_PROGRAMS.keys()), case_sensitive=False),
    default="blastn",
    show_default=True,
    help="BLAST program to use.",
)
@click.option(
    "--database",
    "-d",
    type=click.Choice(BLAST_DATABASES, case_sensitive=False),
    default="nt",
    show_default=True,
    help="Database to search against.",
)
@click.option(
    "--color/--no-color",
    default=True,
    show_default=True,
    help="Enable or disable colorized output.",
)
@click.option(
    "--cache/--no-cache", default=True, show_default=True, help="Enable caching."
)
@click.argument("query_files", type=click.Path(exists=True), required=False, nargs=-1)
def main(program, database, cache, color, query_files):
    console = Console()

    # Handle query input
    encoded_query = ""
    if not query_files:  # Read from stdin if no files are provided
        console.print("Reading from standard input... (Press Ctrl-D to finish)")
        encoded_query = urllib.parse.quote(sys.stdin.read())
    else:
        for query_file in query_files:
            with open(query_file, "r") as file:
                encoded_query += urllib.parse.quote(file.read())

    if cache:
        console.print("Caching enabled.", style="bold green")
        # Generate hash key for the current query to check cache
        input_hash = md5(
            (program + database + encoded_query).encode("utf-8")
        ).hexdigest()
        cache_file = check_cache(input_hash)
    else:
        console.print("Caching disabled.", style="bold yellow")
        cache_file = None

    if cache_file:
        console.print(f"Using cached results from {cache_file}", style="bold green")
        with open(cache_file, "r") as file:
            if not color:
                for line in file:
                    print(line, end="")
                return
            with console.pager(styles=True):
                console.print(file.read())
                return

    # If no cache, proceed with BLAST request

    if program == "megablast":
        args = f"CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE={database}&QUERY={encoded_query}"
    elif program == "rpsblast":
        args = f"CMD=Put&PROGRAM=blastp&SERVICE=rpsblast&DATABASE={database}&QUERY={encoded_query}"
    else:
        args = f"CMD=Put&PROGRAM={program}&DATABASE={database}&QUERY={encoded_query}"
    url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    response = requests.post(url, headers=headers, data=args)

    if response.status_code != 200:
        click.secho("Failed to submit BLAST request.", fg="red")
        sys.exit(1)

    # Parsing RID and RTOE from response
    rid = None
    rtoe = None
    for line in response.text.splitlines():
        if line.startswith("    RID = "):
            rid = line.split("=")[1].strip()
        if line.startswith("    RTOE = "):
            rtoe = int(line.split("=")[1].strip())

    if not rid or rtoe is None:
        click.secho("Failed to retrieve RID or RTOE from BLAST response.", fg="red")
        sys.exit(1)
    else:
        console.print(f"RID: [bold]{rid}[/bold]", f"RTOE: [bold]{rtoe}[/bold]s")

    # Wait for search to complete with a dynamic progress bar
    console = Console()
    with Progress(
        "[progress.description]{task.description}",
        BarColumn(),
        SpinnerColumn(),
        TextColumn("{task.fields[status]}"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task(
            "[cyan]Waiting for BLAST results...",
            total=None,
            start=True,
            status="Preparing...",
        )
        # Initial waiting time split into smaller intervals
        initial_wait = rtoe // 2  # Half of estimated time to first check
        interval = 2  # Check every 2 seconds
        for _ in range(0, initial_wait, interval):
            time.sleep(interval)
            progress.update(task, advance=interval, status="Polling for results...")

        # Continuously check the status after initial wait
        while True:
            time.sleep(3)  # Shorter intervals for regular checks
            req = requests.get(f"{url}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}")
            if "Status=WAITING" in req.text:
                progress.advance(task)
                continue
            elif "Status=FAILED" in req.text:
                click.secho(
                    "Search failed; please report to blast-help@ncbi.nlm.nih.gov.",
                    fg="red",
                )
                sys.exit(4)
            elif "Status=UNKNOWN" in req.text:
                click.secho("Search RID expired.", fg="red")
                sys.exit(3)
            elif "Status=READY" in req.text:
                if "ThereAreHits=yes" in req.text:
                    progress.update(
                        task,
                        completed=progress.tasks[task].total,
                        status="Results ready!",
                    )
                    break
                else:
                    click.secho("No hits found.", fg="red")
                    sys.exit(2)
    # Retrieve and display results
    req = requests.get(f"{url}?CMD=Get&FORMAT_TYPE=Text&RID={rid}")
    result_text = req.text
    if cache:
        write_cache(input_hash, result_text)

    if not color:
        print(result_text)
    else:
        with console.pager(styles=True):
            console.print(result_text)


if __name__ == "__main__":
    main()
