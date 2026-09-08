"""Reusable annotation components, independent of PlasmiCord and PlasBench."""
__version__ = "0.1.0"
from .engine import annotate, parse_gff, parse_amrfinder, to_plasbench

__all__ = ["annotate", "parse_gff", "parse_amrfinder", "to_plasbench"]
