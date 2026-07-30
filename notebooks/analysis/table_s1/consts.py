"""Paths written by the Table S1 analysis scripts."""
from pathlib import Path

TABLE_S1_PATH = Path(__file__).resolve().parent
OUT_PATH = TABLE_S1_PATH / "out"

NAN_PERCENTAGE_CSV = OUT_PATH / "nan_percentage.csv"
