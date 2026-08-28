#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Resolve a species' Chinese common name (中文名) from the bundled dictionary.

The species-name domain is finite and small, so the complete mapping is stored
losslessly as an lzma-compressed dict (~1.3 MB) and looked up in O(1) — far
smaller, faster and more accurate than any trained model.

You can extend coverage without touching the package by dropping more
``latin_name<TAB>chinese`` rows into a user mapping file (see ``_USER_FILE``);
user rows override the bundled ones.
"""

from __future__ import annotations

import lzma
import os
import pickle
from pathlib import Path
from typing import Dict, Optional, Tuple

_DATA_FILE = Path(__file__).parent / "data" / "species.zh.pkl.xz"

# Optional user overrides / additions: latin_name<TAB>chinese per line.
_USER_FILE = Path(
    os.environ.get(
        "WEBLAST_SPECIES_TSV",
        str(Path.home() / ".config" / "webblast" / "species.tsv"),
    )
)

_cache: Optional[Dict[str, str]] = None
_lower_index: Optional[Dict[str, str]] = None


def _load() -> Dict[str, str]:
    """Load the authoritative dict (cached), merged with any user overrides."""
    global _cache
    if _cache is None:
        with lzma.open(_DATA_FILE, "rb") as fh:
            _cache = pickle.load(fh)
        # user mapping wins; lets the table grow as you collect more names
        if _USER_FILE.exists():
            try:
                for line in _USER_FILE.read_text(encoding="utf-8").splitlines():
                    if "\t" in line:
                        key, val = line.split("\t", 1)
                        key, val = key.strip(), val.strip()
                        if key:
                            _cache.setdefault(key, val)
            except OSError:
                pass
    return _cache


def _lower() -> Dict[str, str]:
    global _lower_index
    if _lower_index is None:
        _lower_index = {k.casefold(): v for k, v in _load().items()}
    return _lower_index


def _binomial(name: str) -> str:
    parts = name.split()
    return " ".join(parts[:2]) if len(parts) >= 2 else name


def lookup(name: Optional[str]) -> Optional[str]:
    """Exact / robust lookup. Returns the Chinese name or None."""
    if not name:
        return None
    d = _load()
    key = name.strip()
    # 1) exact
    hit = d.get(key)
    if hit is not None:
        return hit
    # 2) case-insensitive
    hit = _lower().get(key.casefold())
    if hit is not None:
        return hit
    # 3) two-word binomial (handles author suffixes like "Gray")
    bino = _binomial(key)
    if bino != key:
        hit = d.get(bino)
        if hit is None:
            hit = _lower().get(bino.casefold())
        if hit is not None:
            return hit
    return None


class SpeciesResolver:
    """Dictionary-only resolver. ``resolve()`` never calls an external service.

    The mapping is baked into the package; resolve() is a pure in-memory dict
    lookup (fast, deterministic, offline).
    """

    def __bool__(self) -> bool:
        return True

    def resolve(self, sci_name: Optional[str]) -> Tuple[Optional[str], bool]:
        """Return ``(Chinese name or None, is_ai_guess=False)``."""
        return lookup(sci_name), False
