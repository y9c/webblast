#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""NCBI BLAST client using the fast *interactive/website* submission path.

Why this is faster than the plain URL API
------------------------------------------
The official BLAST URL API (``Blast.cgi?CMD=Put``) submits automated searches,
which NCBI places on a lower-priority queue. Interactive web users get
**priority**. The website browser talks to the *same* ``Blast.cgi`` endpoint but
goes through the interactive path:

  * it establishes a real NCBI session (``ncbi_sid`` cookie),
  * it sends a full web-form parameter set (``PAGE``, ``BLAST_PROGRAMS``,
    ``MAX_NUM_SEQ`` ...), and
  * it uses a genuine browser ``User-Agent``.

This module reproduces exactly that, so searches land on the fast interactive
queue. Results are then pulled back as ``JSON2_S`` (compact, machine-readable)
instead of the multi-megabyte HTML the ``Text`` format returns.
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from typing import Callable, Optional

import requests

# The endpoint is identical for submission and retrieval; only ``CMD`` differs.
BASE_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

# A realistic browser User-Agent so NCBI treats us as an interactive client.
BROWSER_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
)

# Program -> (PROGRAM, BLAST_PROGRAMS, PAGE, SERVICE).
# * PROGRAM is the underlying blast binary.
# * BLAST_PROGRAMS selects the web algorithm variant (megablast vs blastn ...).
# * PAGE selects the search form (drives how NCBI prioritizes/routes the job).
# * SERVICE is only needed for rpsblast.
PROGRAMS: dict[str, dict] = {
    # The web server needs MEGABLAST=on (in addition to BLAST_PROGRAMS=megaBlast)
    # to actually select the megablast algorithm.
    "megablast": dict(program="blastn", blast_programs="megaBlast", page="Nucleotides", service=None,
                      params={"MEGABLAST": "on"}),
    "discontiguous-megablast": dict(program="blastn", blast_programs="discoMegablast", page="Nucleotides",
                                    service=None, params={"MEGABLAST": "on"}),
    "blastn": dict(program="blastn", blast_programs="blastn", page="Nucleotides", service=None,
                   params=None),
    "blastp": dict(program="blastp", blast_programs="blastp", page="Proteins", service=None,
                   params=None),
    "quickblastp": dict(program="blastp", blast_programs="kmerBlastp", page="Proteins", service=None,
                        params=None),
    "blastx": dict(program="blastx", blast_programs="blastx", page="Nucleotides", service=None,
                   params=None),
    "tblastn": dict(program="tblastn", blast_programs="tblastn", page="Proteins", service=None,
                    params=None),
    "tblastx": dict(program="tblastx", blast_programs="tblastx", page="Nucleotides", service=None,
                    params=None),
    "rpsblast": dict(program="blastp", blast_programs=None, page="Proteins", service="rpsblast",
                     params=None),
}

# A curated set of common databases (the web form accepts many more).
DATABASES: list[str] = [
    "core_nt",
    "nt",
    "refseq_rna",
    "refseq_genomes",
    "refseq_representative_genomes",
    "nr",
    "refseq_protein",
    "swissprot",
    "pdbaa",
    "pdbnt",
    "est",
    "wgs",
    "16S_ribosomal_RNA",
]

# Output formats we can parse natively.
FORMAT_TYPES = {
    "json": "JSON2_S",        # simplified JSON (compact, what we parse)
    "json-full": "JSON2",     # full JSON (larger, richer)
    "xml": "XML2",
    "text": "Text",
    "tabular": "Text",        # Text + ALIGNMENT_VIEW=Tabular
}

# NCBI status strings from SearchInfo.
STATUS_WAITING = "WAITING"
STATUS_READY = "READY"
STATUS_FAILED = "FAILED"
STATUS_UNKNOWN = "UNKNOWN"

_DEFAULT_EMAIL = "yech1990@gmail.com"
_DEFAULT_TOOL = "webblast"


@dataclass
class Submission:
    """The result of submitting a search."""

    rid: str
    rtoe: int

    def __post_init__(self):
        if not self.rid:
            raise ValueError("empty RID")


@dataclass
class BlastClient:
    """Thin client that talks to NCBI BLAST through the interactive web path."""

    program: str = "megablast"
    database: str = "nt"
    email: str = _DEFAULT_EMAIL
    tool: str = _DEFAULT_TOOL
    api_key: str = ""
    timeout: float = 60.0
    # poll cadence (seconds) & backoff. NCBI's RTOE is a very conservative
    # over-estimate (fresh searches finish in ~3.4s), so we poll quickly from the
    # start to catch the result at the NCBI floor, then back off politely.
    poll_interval: float = 1.5
    max_poll_interval: float = 8.0
    backoff: float = 1.35
    last_response: str = ""

    def __post_init__(self):
        if self.program not in PROGRAMS:
            raise ValueError(f"Unknown program {self.program!r}. Choose from: {sorted(PROGRAMS)}")
        self.cfg = PROGRAMS[self.program]
        self._page = self.cfg["page"]
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": BROWSER_UA})
        self._init_session()

    # ------------------------------------------------------------------ #
    # session / transport
    # ------------------------------------------------------------------ #
    def _init_session(self) -> None:
        """Prime the session so NCBI assigns us an ``ncbi_sid`` cookie."""
        try:
            self.session.get(
                BASE_URL,
                params={"PAGE": self._page, "PROGRAM": self.cfg["program"]},
                timeout=self.timeout,
            )
        except requests.RequestException:
            # A session isn't strictly required; submission will still work.
            pass

    def _build_params(self, query: str, format_type: str, **extra) -> dict:
        params = {
            "CMD": "Put",
            "PROGRAM": self.cfg["program"],
            "DATABASE": self.database,
            "QUERY": query,
            "PAGE": self._page,
            "EMAIL": self.email,
            "TOOL": self.tool,
            "HITLIST_SIZE": str(extra.pop("hitlist_size", 100)),
            "MAX_NUM_SEQ": str(extra.pop("max_num_seq", 100)),
            "ALIGNMENTS": str(extra.pop("alignments", 100)),
            "DESCRIPTIONS": str(extra.pop("descriptions", 100)),
            "FORMAT_TYPE": format_type,
            "USER_TYPE": "2",
        }
        if self.cfg.get("blast_programs"):
            params["BLAST_PROGRAMS"] = self.cfg["blast_programs"]
        if self.cfg.get("service"):
            params["SERVICE"] = self.cfg["service"]
        if self.cfg.get("params"):
            params.update(self.cfg["params"])
        if self.api_key:
            params["api_key"] = self.api_key
        params.update(extra)  # allow arbitrary passthrough (EXPECT, WORD_SIZE...)
        # drop Nones / empties
        return {k: v for k, v in params.items() if v not in (None, "")}

    # ------------------------------------------------------------------ #
    # public API
    # ------------------------------------------------------------------ #
    def submit(self, query: str, format_type: str = "JSON2_S", **extra) -> Submission:
        """Submit one FASTA string. Returns (rid, rtoe)."""
        params = self._build_params(query, format_type, **extra)
        resp = self.session.post(BASE_URL, data=params, timeout=self.timeout)
        resp.raise_for_status()
        self.last_response = resp.text
        rid, rtoe = _parse_rid_rtoe(resp.text)
        return Submission(rid=rid, rtoe=rtoe)

    def get_status(self, rid: str) -> str:
        """Return WAITING / READY / FAILED / UNKNOWN for a RID."""
        resp = self.session.get(
            BASE_URL,
            params={"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"},
            timeout=self.timeout,
        )
        resp.raise_for_status()
        m = re.search(r"Status=(\w+)", resp.text)
        if not m:
            # a blank / non-standard page usually means still processing
            return STATUS_WAITING
        return m.group(1).upper()

    def wait(
        self,
        rid: str,
        rtoe: Optional[int] = None,
        timeout: float = 1200.0,
        on_status: Optional[Callable[[int, float, str], None]] = None,
    ) -> str:
        """Block until READY (raise on FAILED/UNKNOWN), polling early & adaptively.

        ``rtoe`` from NCBI is a conservative estimate that is usually much larger
        than the real completion time (fresh searches often finish in 3-10s), so
        we start polling after a short delay and back off geometrically.

        ``on_status(attempt, elapsed, status)`` is called after every poll so a
        caller can drive a live progress animation.
        """
        start = time.time()
        delay = self.poll_interval
        attempt = 0
        while time.time() - start < timeout:
            time.sleep(delay)
            attempt += 1
            elapsed = time.time() - start
            status = self.get_status(rid)
            if on_status is not None:
                on_status(attempt, elapsed, status)
            if status != STATUS_WAITING:
                return status
            delay = min(delay * self.backoff, self.max_poll_interval)
        raise TimeoutError(f"BLAST search {rid} not ready within {timeout:.0f}s")

    def fetch(self, rid: str, format_type: str = "JSON2_S", **extra) -> bytes:
        """Download the final report as bytes."""
        params = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": format_type}
        params.update(extra)
        resp = self.session.get(BASE_URL, params=params, timeout=self.timeout)
        resp.raise_for_status()
        return resp.content

    def run(self, query: str, format_type: str = "JSON2_S", timeout: float = 1200.0, **extra):
        """Submit, wait, fetch, and return the raw report bytes."""
        sub = self.submit(query, format_type=format_type, **extra)
        status = self.wait(sub.rid, sub.rtoe, timeout=timeout)
        if status != STATUS_READY:
            raise RuntimeError(f"BLAST {sub.rid} ended with status {status}")
        return self.fetch(sub.rid, format_type=format_type)


# ---------------------------------------------------------------------- #
# helpers
# ---------------------------------------------------------------------- #
def _parse_rid_rtoe(html: str):
    """Extract (rid, rtoe) from a submission response page."""
    rid_m = re.search(r"RID\s*=\s*(\S+)", html)
    rtoe_m = re.search(r"RTOE\s*=\s*(\d+)", html)
    if not rid_m:
        snippet = re.sub(r"<[^>]+>", "", html)
        raise RuntimeError(
            "Could not get a RID from NCBI BLAST. "
            f"Response (first 400 chars): {snippet[:400]!r}"
        )
    rid = rid_m.group(1)
    rtoe = int(rtoe_m.group(1)) if rtoe_m else 10
    return rid, rtoe
