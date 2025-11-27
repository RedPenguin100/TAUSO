import os

def find_project_root(target_folder='scripts'):
    """
    Walk up the directory tree to find the project root,
    which is defined as the folder that contains `scripts/`.
    """
    current_path = os.path.abspath(os.getcwd())
    while True:
        if os.path.isdir(os.path.join(current_path, target_folder)):
            return current_path
        parent = os.path.dirname(current_path)
        if parent == current_path:
            raise FileNotFoundError(f"Could not find folder '{target_folder}' in any parent directory.")
        current_path = parent

def save_feature(df, feature_name):
    project_root = find_project_root('scripts')
    feature_dir = os.path.join(project_root, 'scripts', 'features')
    os.makedirs(feature_dir, exist_ok=True)
    sub_df = df[['index', feature_name]]
    file_path = os.path.join(feature_dir, f'{feature_name}.csv')
    sub_df.to_csv(file_path, index=False)