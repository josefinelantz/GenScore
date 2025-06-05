import pandas as pd

def load_pathway_summary_csv_to_dict(csv_path):
    """
    Load a CSV file with pathway scores and return a dict: pathway → {total_score, genes}
    """
    df = pd.read_csv("data/Pathway_Impact_Summary.csv")
    df = df.set_index("pathway")
    return df.to_dict(orient="index")