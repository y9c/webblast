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


def read_records(paths: List[str], limit: Optional[int] = None) -> List[Tuple[str, str]]:
    """Read ``(name, sequence)`` records from files (or FASTA from stdin)."""
    records: List[Tuple[str, str]] = []
    if not paths:
        import sys

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
    return records


def join_records(records: List[Tuple[str, str]]) -> str:
    """Build a multi-record FASTA string from ``(name, sequence)`` pairs."""
    return "\n".join(f">{name}\n{seq}" for name, seq in records)


def chunk_records(records: List[Tuple[str, str]], n: int) -> List[List[Tuple[str, str]]]:
    """Split records into at most ``n`` roughly equal chunks."""
    n = max(1, n)
    if n == 1 or len(records) <= n:
        return [records]
    k, rem = divmod(len(records), n)
    out, idx = [], 0
    for i in range(n):
        size = k + (1 if i < rem else 0)
        out.append(records[idx : idx + size])
        idx += size
    return out



