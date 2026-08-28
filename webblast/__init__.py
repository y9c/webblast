"""webblast — fast NCBI BLAST searches through the interactive web queue.

``webblast`` submits nucleotide/protein queries to NCBI BLAST using the exact
submission path the website uses (real session + browser User-Agent + web
parameter set), which lands on NCBI's fast interactive queue instead of the
slower automated API queue. Results are parsed from the compact ``JSON2_S``
format into clean records.
"""

from .client import BlastClient, Submission, BASE_URL, PROGRAMS, DATABASES
from .parse import BlastReport, Hit, HSP, QueryResult, parse_blast_json
from .utils import build_query, iter_records

__version__ = "0.1.0"

__all__ = [
    "BlastClient",
    "Submission",
    "BlastReport",
    "QueryResult",
    "Hit",
    "HSP",
    "parse_blast_json",
    "build_query",
    "iter_records",
    "BASE_URL",
    "PROGRAMS",
    "DATABASES",
    "__version__",
]
