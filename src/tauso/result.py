import pandas as pd

from .consts import EXPERIMENT_RESULTS


def save_results_organism(df: pd.DataFrame, organism: str, experiment_name: str, algorithm: str):
    EXPERIMENT_PATH = EXPERIMENT_RESULTS / f'{experiment_name}'
    ORGANISM_RESULTS_PATH = EXPERIMENT_PATH / f'{organism}_results'
    ALGORITHM_PATH = ORGANISM_RESULTS_PATH / f'{algorithm}.csv'

    EXPERIMENT_PATH.mkdir(exist_ok=True)
    ORGANISM_RESULTS_PATH.mkdir(exist_ok=True)

    df.to_csv(str(ALGORITHM_PATH), index=False)


def save_results_on_target(df: pd.DataFrame, experiment_name: str, algorithm: str):
    csv_path = EXPERIMENT_RESULTS / f'{experiment_name}/antisense_results/{algorithm}.csv'
    df.to_csv(csv_path, index=False)
