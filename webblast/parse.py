#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Parse NCBI BLAST JSON2 / JSON2_S reports into clean structured records.

The NCBI BLAST JSON formats wrap everything in ``{"BlastOutput2": [...]}``.
This module walks that structure and produces lightweight dataclasses that the
CLI can render as a table, JSON, or tab-delimited output.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union


@dataclass
class HSP:
    """A single high-scoring pair (alignment)."""

    num: Optional[int] = None
    bit_score: Optional[float] = None
    score: Optional[float] = None
    evalue: Optional[float] = None
    identity: Optional[int] = None
    positive: Optional[int] = None
    align_len: Optional[int] = None
    gaps: Optional[int] = None
    query_from: Optional[int] = None
    query_to: Optional[int] = None
    hit_from: Optional[int] = None
    hit_to: Optional[int] = None
    query_strand: Optional[str] = None
    hit_strand: Optional[str] = None
    qseq: str = ""
    hseq: str = ""
    midline: str = ""

    @property
    def pident(self) -> Optional[float]:
        """Percent identity (identity / alignment length * 100)."""
        if self.identity is None or not self.align_len:
            return None
        return 100.0 * self.identity / self.align_len


@dataclass
class Hit:
    """One database hit for a query (a matched sequence)."""

    num: Optional[int] = None
    accession: Optional[str] = None
    id: Optional[str] = None
    title: Optional[str] = None
    taxid: Optional[int] = None
    sciname: Optional[str] = None
    len: Optional[int] = None
    hsps: List[HSP] = field(default_factory=list)


@dataclass
class QueryResult:
    """The search result for a single query sequence."""

    query_id: Optional[str] = None
    query_title: Optional[str] = None
    query_len: Optional[int] = None
    hits: List[Hit] = field(default_factory=list)
    message: Optional[str] = None


@dataclass
class BlastReport:
    """A whole BLAST run."""

    program: Optional[str] = None
    version: Optional[str] = None
    db: Optional[str] = None
    params: Dict = field(default_factory=dict)
    queries: List[QueryResult] = field(default_factory=list)

    @property
    def total_hits(self) -> int:
        return sum(len(q.hits) for q in self.queries)


def _first_hit_desc(hit: dict) -> dict:
    """Return the first description entry as a normalised dict."""
    descs = hit.get("description") or []
    if not descs:
        return {}
    d = descs[0]
    return {
        "id": d.get("id"),
        "accession": d.get("accession"),
        "title": d.get("title"),
        "taxid": d.get("taxid"),
        "sciname": d.get("sciname"),
    }


def _parse_hsp(h: dict) -> HSP:
    return HSP(
        num=h.get("num"),
        bit_score=h.get("bit_score"),
        score=h.get("score"),
        evalue=h.get("evalue"),
        identity=h.get("identity") if isinstance(h.get("identity"), (int, float)) else None,
        positive=h.get("positive") if isinstance(h.get("positive"), (int, float)) else None,
        align_len=h.get("align_len"),
        gaps=h.get("gaps"),
        query_from=h.get("query_from"),
        query_to=h.get("query_to"),
        hit_from=h.get("hit_from"),
        hit_to=h.get("hit_to"),
        query_strand=h.get("query_strand"),
        hit_strand=h.get("hit_strand"),
        qseq=h.get("qseq") or "",
        hseq=h.get("hseq") or "",
        midline=h.get("midline") or "",
    )


def _parse_hit(h: dict) -> Hit:
    d = _first_hit_desc(h)
    hsps = [_parse_hsp(x) for x in (h.get("hsps") or [])]
    return Hit(
        num=h.get("num"),
        accession=d.get("accession"),
        id=d.get("id"),
        title=d.get("title"),
        taxid=d.get("taxid"),
        sciname=d.get("sciname"),
        len=h.get("len"),
        hsps=hsps,
    )


def _parse_search(search: dict, query_title: str, query_id: str, query_len: int) -> QueryResult:
    hits = [_parse_hit(h) for h in (search.get("hits") or [])]
    return QueryResult(
        query_id=search.get("query_id", query_id),
        query_title=search.get("query_title", query_title),
        query_len=search.get("query_len", query_len),
        hits=hits,
        message=search.get("message"),
    )


def parse_blast_json(raw: Union[bytes, bytearray, str]) -> BlastReport:
    """Parse NCBI BLAST JSON (``JSON2`` / ``JSON2_S``) into a BlastReport."""
    if isinstance(raw, (bytes, bytearray)):
        # JSON2 (full) may be served as a zip; JSON2_S is plain JSON.
        if raw[:2] == b"PK":
            import io
            import zipfile

            with zipfile.ZipFile(io.BytesIO(raw)) as zf:
                name = zf.namelist()[0]
                text = zf.read(name).decode("utf-8")
        else:
            text = raw.decode("utf-8")
    else:
        text = raw

    data = json.loads(text)
    if not isinstance(data, dict) or "BlastOutput2" not in data:
        raise ValueError("Input is not a BLAST JSON2 report (missing BlastOutput2).")

    report = BlastReport()
    for entry in data["BlastOutput2"]:
        rep = entry.get("report") or {}
        if not report.program:
            report.program = rep.get("program")
            report.version = rep.get("version")
        search_target = rep.get("search_target") or {}
        report.db = search_target.get("db") or report.db
        report.params = rep.get("params") or report.params

        results = rep.get("results") or {}
        # results is a dict (keyed per search) or a list
        for search in _iter_results(results):
            qid = search.get("query_id")
            qtitle = search.get("query_title")
            qlen = search.get("query_len")
            report.queries.append(
                _parse_search(search, qtitle or "", qid or "", qlen or 0)
            )
    return report


def _iter_results(results):
    if isinstance(results, dict):
        for value in results.values():
            if isinstance(value, dict):
                yield value
    elif isinstance(results, list):
        for value in results:
            if isinstance(value, dict):
                yield value
