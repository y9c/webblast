#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Backward-compatible entry point (previously the whole CLI).

The old package exposed the CLI as ``webblast.blast:main``; it now re-exports
the new :func:`webblast.cli.main`. Importing this module still works for anyone
that did ``from webblast.blast import main``.
"""

from .cli import main  # noqa: F401

__all__ = ["main"]
