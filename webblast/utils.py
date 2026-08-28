#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility helpers to read sequences from FASTA / FASTQ / BAM & build a query."""

from __future__ import annotations

import gzip
import os
from typing import Iterator, List, Optional, Tuple

try:  # optional: only needed for BAM/SAM input
    import pysam
except ImportError:  # pragma: no cover
    pysam = None


def open_text(path: str):
    """Open a (possibly gzipped) text file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def iter_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(name, sequence)`` records from a FASTA file (plain or gzipped)."""
    name, seq = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split(None, 1)[0].strip() if line[1:].strip() else "seq"
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def iter_fastq(path: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(name, sequence)`` records from a FASTQ file (plain or gzipped)."""
    name = None
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if i % 4 == 0:
                name = line[1:].split()[0] if line.startswith("@") and line[1:].split() else "seq"
            elif i % 4 == 1:
                yield name, line


def iter_bam(path: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(name, sequence)`` records from an alignment file (requires pysam)."""
    if pysam is None:
        raise ImportError("pysam is required to read BAM/SAM files: pip install pysam")
    with pysam.AlignmentFile(path, "r") as fh:
        for aln in fh:
            if aln.query_sequence:
                yield aln.query_name, aln.query_sequence


def iter_records(path: str) -> Iterator[Tuple[str, str]]:
    """Auto-detect format (by extension) and yield ``(name, sequence)`` records."""
    ext = os.path.splitext(path)[0].lower()
    if ext in (".fa", ".fasta", ".faa"):
        yield from iter_fasta(path)
    elif ext in (".fq", ".fastq"):
        yield from iter_fastq(path)
    elif ext in (".bam", ".sam"):
        yield from iter_bam(path)
    else:
        # default to FASTA
        yield from iter_fasta(path)


def build_query(paths: List[str], limit: Optional[int] = None) -> str:
    """Read sequence files and build a single multi-record FASTA query string.

    If no paths are given (or they are empty), read FASTA from stdin.
    """
    records: List[Tuple[str, str]] = []
    if not paths:
        import sys

        # parse a paste-in FASTA from stdin
        name, seq = None, []
        for raw in sys.stdin:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                    if limit and len(records) >= limit:
                        break
                name = line[1:].split(None, 1)[0].strip() if line[1:].strip() else "seq"
                seq = []
            else:
                seq.append(line)
        if name is not None:
            records.append((name, "".join(seq)))
    else:
        for path in paths:
            for name, seq in iter_records(path):
                records.append((name, seq))
                if limit and len(records) >= limit:
                    break

    if not records:
        raise ValueError("No sequences found in input.")

    return "\n".join(f">{name}\n{seq}" for name, seq in records)


# ---------------------------------------------------------------------- #
# Backwards-compatible helper (kept from the original package)
# ---------------------------------------------------------------------- #
def read_file(file_path, record_limit=None):
    """Return FASTA records as ``["name\\nseq", ...]`` strings.

    Retained for backward compatibility; prefer :func:`build_query`.
    """
    if file_path == "-":
        import sys
        paths = [sys.stdin]
        # not supported cleanly here; read as text
        return [f"stdin\n{''.join(l for l in sys.stdin if not l.startswith('>'))}"]

    out = []
    for name, seq in iter_records(file_path):
        out.append(f">{name}\n{seq}")
        if record_limit and len(out) >= record_limit:
            break
    return out
