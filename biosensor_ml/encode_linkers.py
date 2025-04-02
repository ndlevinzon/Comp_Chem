import os
import json
import pandas as pd
from tqdm import tqdm
from protlearn.features import aaindex1


def extract_linker_sequence(json_path):
    """Extracts the 4-residue linker from positions 70, 71, 313, and 314."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    sequence = data["sequences"][0]["protein"]["sequence"]
    linker = ''.join([sequence[i - 1] for i in [70, 71, 313, 314]])
    return linker


def encode_linker(linker):
    """Returns (feature names, feature values) for a 4-residue linker."""
    features, feature_names = aaindex1([linker])
    return feature_names, features[0]


def process_all_jsons(json_dir, output_csv):
    """Process all JSON files in a directory tree, extract linker, encode, and save."""
    encoded_data = []
    feature_names = None

    for root, _, files in os.walk(json_dir):
        for filename in tqdm(files, desc="Encoding linkers"):
            if filename.endswith(".json"):
                path = os.path.join(root, filename)
                try:
                    linker = extract_linker_sequence(path)
                    if feature_names is None:
                        feature_names, features = encode_linker(linker)
                    else:
                        _, features = encode_linker(linker)
                    row = {
                        "filename": filename,
                        "linker": linker,
                    }
                    row.update(dict(zip(feature_names, features)))
                    encoded_data.append(row)
                except Exception as e:
                    print(f"[WARNING] Failed on {filename}: {e}")

    df = pd.DataFrame(encoded_data)
    df.to_csv(output_csv, index=False)
    print(f"\nSaved AAindex-encoded linker features to: {output_csv}")


if __name__ == "__main__":
    # Customize your input and output paths here:
    INPUT_JSON_DIR = "C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake_test/json_test/"
    OUTPUT_CSV_PATH = "C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake_test/biosensor_ml/data/encoded_linkers.csv"

    process_all_jsons(INPUT_JSON_DIR, OUTPUT_CSV_PATH)
