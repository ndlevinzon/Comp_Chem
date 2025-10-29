#!/usr/bin/env python3
import argparse
import os
import sys
import json
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import List, Tuple, Dict

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.feature_selection import f_classif

# ----------------------------- utilities -----------------------------

META_COLS = {
    "replicate", "ligand", "mutant"
}

RES_PATTERNS = [
    re.compile(r"^res[_\-]?\d+$", re.IGNORECASE),   # res_372, res395
    re.compile(r"^resi[_\-]?\d+$", re.IGNORECASE),  # resi_372
    re.compile(r"^R\d+$", re.IGNORECASE),           # R372
]
ATOM_PAT = re.compile(r"^atom[_\-]?\d+$", re.IGNORECASE)


def _is_feature_col(col: str) -> bool:
    if col in META_COLS:
        return False
    # numeric-like RMSF columns (residue or atom)
    if any(p.match(col) for p in RES_PATTERNS):
        return True
    if ATOM_PAT.match(col):
        return True
    # generic numeric columns can also be treated as features
    # but we filter later by dtype
    return False


def collect_rmsf_features(df: pd.DataFrame) -> Tuple[np.ndarray, List[str]]:
    """
    Collect numeric RMSF feature columns. Prefer residue-style names if present;
    otherwise fall back to atom_* columns (or any numeric non-meta columns).
    """
    # First pass: candidate names by pattern
    cand = [c for c in df.columns if _is_feature_col(c)]
    # Enforce numeric dtype
    num_cols = [c for c in cand if pd.api.types.is_numeric_dtype(df[c])]
    if not num_cols:
        # fallback: any numeric column not in meta
        num_cols = [c for c in df.columns
                    if c not in META_COLS and pd.api.types.is_numeric_dtype(df[c])]
    if not num_cols:
        raise SystemExit("[rmsf-deeplda] No numeric RMSF feature columns found.")

    X = df[num_cols].astype("float32").values
    return X, num_cols


# ----------------------------- deep model -----------------------------
def _torch_available():
    try:
        import torch  # noqa
        return True
    except Exception:
        return False


def train_deep_classifier(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                          n_splits: int = 3, hidden: int = 64,
                          max_epochs: int = 25, lr: float = 1e-3,
                          weight_decay: float = 1e-4, seed: int = 123):
    """
    Small MLP with a linear last layer; returns:
      - folds (cv accuracies)
      - W_mean: class prototypes in the last latent (for importance/hist)
      - last_layer_W: (n_classes, hidden) for importance
      - scaler used for standardization
    """
    import torch
    import torch.nn as nn
    import torch.optim as optim

    gkf = GroupKFold(n_splits=n_splits)
    rng = np.random.RandomState(seed)
    n_classes = len(np.unique(y))

    folds = []
    W_list = []

    class MLP(nn.Module):
        def __init__(self, in_dim, hidden, out_dim):
            super().__init__()
            self.f = nn.Sequential(
                nn.Linear(in_dim, hidden),
                nn.ReLU(),
                nn.BatchNorm1d(hidden),
                nn.Linear(hidden, hidden),
                nn.ReLU(),
            )
            self.head = nn.Linear(hidden, out_dim)

        def forward(self, x):
            h = self.f(x)
            z = self.head(h)
            return z, h

    # We'll also pool class means in the latent space across folds
    class_means_accum = np.zeros((n_classes, hidden), dtype=np.float32)
    class_counts = np.zeros((n_classes,), dtype=np.int64)

    for tr, te in gkf.split(X, y, groups):
        scaler = StandardScaler().fit(X[tr])
        Xtr = scaler.transform(X[tr]).astype(np.float32)
        Xte = scaler.transform(X[te]).astype(np.float32)

        model = MLP(X.shape[1], hidden, n_classes)
        torch.manual_seed(seed)
        model.train()

        opt = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
        ce = nn.CrossEntropyLoss()

        Xt = torch.from_numpy(Xtr)
        yt = torch.from_numpy(y[tr]).long()
        Xe = torch.from_numpy(Xte)
        ye = torch.from_numpy(y[te]).long()

        for epoch in range(max_epochs):
            opt.zero_grad()
            logits, _ = model(Xt)
            loss = ce(logits, yt)
            loss.backward()
            opt.step()

        # Eval
        model.eval()
        with torch.no_grad():
            logits, h_tr = model(Xt)
            _, h_te = model(Xe)
            pred = logits.argmax(dim=1).cpu().numpy()
            acc = (pred == y[tr]).mean()  # training accuracy just for sanity
            # test predictions
            logits_te, _ = model(Xe)
            pred_te = logits_te.argmax(dim=1).cpu().numpy()
            acc_te = (pred_te == y[te]).mean()
            folds.append(float(acc_te))

            # accumulate class means in latent for importance/hist projection
            h_tr = h_tr.cpu().numpy()
            for c in range(n_classes):
                mask = (y[tr] == c)
                if np.any(mask):
                    class_means_accum[c] += h_tr[mask].mean(axis=0)
                    class_counts[c] += 1

        # Save last-layer weights
        W = model.head.weight.detach().cpu().numpy()  # (n_classes, hidden)
        W_list.append(W)

    W_mean = np.stack(W_list, axis=0).mean(axis=0)  # (n_classes, hidden)
    # If we computed latent class means across folds, average them
    with np.errstate(invalid="ignore"):
        class_means = (class_means_accum / np.maximum(class_counts[:, None], 1)).astype(np.float32)

    return {
        "folds": folds,
        "W_mean": W_mean,             # class weights in the last layer
        "latent_class_means": class_means,  # not strictly needed
        "hidden_dim": hidden,
    }


def train_linear_baseline(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                          n_splits: int = 3, C: float = 1.0, seed: int = 123):
    """
    Multinomial logistic regression with replica-aware GroupKFold CV.
    """
    gkf = GroupKFold(n_splits=n_splits)
    folds = []
    W_list = []

    for tr, te in gkf.split(X, y, groups):
        scaler = StandardScaler().fit(X[tr])
        Xtr = scaler.transform(X[tr])
        Xte = scaler.transform(X[te])

        clf = LogisticRegression(
            penalty="l2", C=C, solver="lbfgs",
            multi_class="multinomial", max_iter=4000,
            random_state=seed
        )
        clf.fit(Xtr, y[tr])
        acc = clf.score(Xte, y[te])
        folds.append(float(acc))
        W_list.append(clf.coef_)  # (n_classes, n_features)

    W_mean = np.stack(W_list, axis=0).mean(axis=0)
    return {"folds": folds, "W_mean": W_mean}


def feature_importance_from_W(W_mean: np.ndarray, feat_names: List[str]) -> pd.DataFrame:
    """
    Importance per feature from class weights:
      - (1, d): binary logistic returns one separating vector -> use |w|
      - (2, d): two-class multinomial style -> use |w1 - w0|
      - (k>=3, d): center across classes then L2 norm
    """
    if W_mean.ndim != 2:
        raise ValueError(f"Expected 2D W_mean, got {W_mean.shape}")
    k, d = W_mean.shape
    if k == 1:
        imp = np.abs(W_mean[0, :])
    elif k == 2:
        imp = np.abs(W_mean[1, :] - W_mean[0, :])
    else:
        Wc = W_mean - W_mean.mean(axis=0, keepdims=True)
        imp = np.sqrt((Wc ** 2).sum(axis=0))

    out = pd.DataFrame({"feature": feat_names, "importance": imp})
    out = out.sort_values("importance", ascending=False).reset_index(drop=True)
    return out



# ----------------------------- grouping (same syntax as your deeplda) -----------------------------
def parse_set_token(tok: str) -> Tuple[str, str, str, List[str] | None]:
    parts = tok.split("|")
    if len(parts) not in (3, 4):
        raise argparse.ArgumentTypeError(
            f"--set must be 'GROUP|LIGAND|MUTANT|rep1,rep2,...' or 'GROUP|LIGAND|MUTANT'; got: {tok}"
        )
    group, lig, mut = parts[0].strip(), parts[1].strip(), parts[2].strip()
    reps = None
    if len(parts) == 4:
        reps_str = parts[3].strip()
        if reps_str and reps_str != "*":
            reps = [r.strip() for r in reps_str.split(",") if r.strip()]
    return group, lig, mut, reps


def build_group_labels(df: pd.DataFrame, set_tokens: List[str]) -> pd.DataFrame:
    df = df.copy()
    need = {"ligand", "mutant", "replicate"}
    if not need.issubset(df.columns):
        raise SystemExit("CSV must contain columns: ligand, mutant, replicate")

    # Normalize columns for robust matching
    df["ligand_norm"]    = df["ligand"].astype(str).str.strip().str.casefold()
    df["mutant_norm"]    = df["mutant"].astype(str).str.strip().str.casefold()
    df["replicate_norm"] = df["replicate"].astype(str).str.strip().str.casefold()

    # Use a simple (object) dtype to avoid pandas NA string semantics
    df["group_label"] = ""

    total_matched = 0
    for tok in set_tokens:
        group, lig, mut, reps = parse_set_token(tok)
        lig_n = lig.strip().casefold()
        mut_n = mut.strip().casefold()
        mask = (df["ligand_norm"] == lig_n) & (df["mutant_norm"] == mut_n)
        if reps is not None:
            reps_n = [r.strip().casefold() for r in reps]
            mask &= df["replicate_norm"].isin(reps_n)

        n_match = int(mask.sum())
        print(f"[select] {tok} -> matched {n_match} rows")
        total_matched += n_match
        df.loc[mask, "group_label"] = group

    df_sel = df.loc[df["group_label"] != ""].copy()
    if df_sel.empty:
        # Helpful diagnostics
        ligs = sorted(df["ligand"].astype(str).unique().tolist())
        muts = sorted(df["mutant"].astype(str).unique().tolist())
        reps = sorted(df["replicate"].astype(str).unique().tolist())
        raise SystemExit(
            "[rmsf-deeplda] No rows matched your --set selections.\n"
            f"  unique ligand values: {ligs}\n"
            f"  unique mutant values: {muts}\n"
            f"  unique replicate values: {reps}\n"
            "  (Tip: matching is case-insensitive with whitespace stripped.)"
        )

    # Clean up helper columns
    df_sel.drop(columns=["ligand_norm", "mutant_norm", "replicate_norm"], inplace=True)
    return df_sel



# ----------------------------- plotting (same as your style) -----------------------------
def compute_cv1_scores(X: np.ndarray, y: np.ndarray,
                       model_name: str, W_mean: np.ndarray,
                       class_names: List[str]) -> Tuple[np.ndarray, str]:
    from sklearn.preprocessing import StandardScaler
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    import numpy as np

    # sanity: do we actually have 2 classes?
    if np.unique(y).size < 2:
        raise ValueError("compute_cv1_scores: only one class present after selection")

    scaler = StandardScaler().fit(X)
    Xs = scaler.transform(X)

    if len(class_names) == 2 and W_mean is not None:
        if W_mean.ndim != 2 or W_mean.shape[1] != X.shape[1]:
            raise ValueError(f"Unexpected W_mean shape {W_mean.shape} for X {X.shape}")

        if W_mean.shape[0] == 1:
            w = W_mean[0, :]
        elif W_mean.shape[0] == 2:
            w = W_mean[1, :] - W_mean[0, :]
        else:
            # very unusual; treat as multiclass
            w = None

        if w is not None:
            n = np.linalg.norm(w)
            if n > 0:
                scores = Xs @ (w / n)
                # If projection is degenerate (all same or NaN), fall back
                if np.allclose(np.nanstd(scores), 0.0) or not np.isfinite(scores).all():
                    w = None
            else:
                w = None

            if w is not None:
                return scores, "Canonical variable (CV1) — binary (logistic projection)"

    # Fallback: 1D LDA for visualization (works for 2+ classes)
    lda = LinearDiscriminantAnalysis(n_components=1)
    scores = lda.fit_transform(Xs, y).ravel()
    return scores, "Canonical variable (CV1) — LDA projection"




def plot_cv_hist(scores: np.ndarray, y: np.ndarray, class_names: List[str], out_png: str):
    import numpy as np
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 6))
    bins = 60

    # Map each encoded label actually present to its name
    present_labels = np.unique(y)               # e.g., array([0, 1])
    # guard: if length mismatch, fall back to indices
    if len(present_labels) == len(class_names):
        label_to_name = {lab: class_names[i] for i, lab in enumerate(present_labels)}
    else:
        # typical case should not hit this; still keep plotting usable
        label_to_name = {lab: f"class {int(lab)}" for lab in present_labels}

    plotted_any = False
    for lab in present_labels:
        mask = (y == lab)
        xs = scores[mask]
        xs = xs[~np.isnan(xs)]
        if xs.size == 0:
            continue
        plt.hist(xs, bins=bins, density=True, alpha=0.35,
                 label=f"{label_to_name[lab]} (n={mask.sum()})")
        plotted_any = True

    if not plotted_any:
        print("[warn] plot_cv_hist: no non-NaN scores to plot for any class")

    plt.xlabel("Canonical variable (CV1)")
    plt.ylabel("Density")
    plt.title("Separation along CV1 by selection")
    plt.legend(frameon=False, loc="best")
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()
    print(f"[plot] wrote: {out_png}")



def plot_topk_bar(imp_df: pd.DataFrame, out_png: str, k: int = 20, title: str = "Top Discriminating Residues/Atoms"):
    top = imp_df.head(k)
    plt.figure(figsize=(7.5, 6))
    plt.barh(top["feature"], top["importance"])
    plt.gca().invert_yaxis()
    plt.xlabel("Importance (L2 across classes)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()
    print(f"[plot] wrote: {out_png}")


# ------------------- greedy subset (same idea as your dihedral code) -------------------
def _cv_accuracy_for_subset(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                            cols: List[int], seed: int = 42, n_splits: int = 3) -> float:
    Xs = X[:, cols]
    gkf = GroupKFold(n_splits=n_splits)
    accs = []
    for tr, te in gkf.split(Xs, y, groups):
        scaler = StandardScaler().fit(Xs[tr])
        Xtr = scaler.transform(Xs[tr]); Xte = scaler.transform(Xs[te])
        clf = LogisticRegression(
            penalty="l2", C=1.0, multi_class="multinomial",
            solver="lbfgs", max_iter=2000, random_state=seed
        )
        clf.fit(Xtr, y[tr])
        accs.append(clf.score(Xte, y[te]))
    return float(np.mean(accs))


def greedy_select_best_features(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                                feat_names: List[str], imp_df: pd.DataFrame,
                                max_k: int = 12, min_improve: float = 1e-3,
                                seed: int = 42) -> Dict:
    # candidate order by importance
    name_to_idx = {n: i for i, n in enumerate(feat_names)}
    candidates = [n for n in imp_df["feature"].tolist() if n in name_to_idx]
    if not candidates:
        return {"best_subset": [], "best_acc": 0.0, "curve": []}

    best_subset: List[str] = []
    best_cols: List[int] = []
    curve = []

    # start with best single
    first = candidates[0]
    best_subset = [first]
    best_cols = [name_to_idx[first]]
    best_acc = _cv_accuracy_for_subset(X, y, groups, best_cols, seed=seed)
    curve.append({"k": 1, "acc": best_acc, "subset": best_subset.copy()})

    remaining = candidates[1:]
    while len(best_subset) < max_k and remaining:
        improvements = []
        for n in remaining:
            cols = best_cols + [name_to_idx[n]]
            acc = _cv_accuracy_for_subset(X, y, groups, cols, seed=seed)
            improvements.append((acc, n))
        improvements.sort(reverse=True, key=lambda x: x[0])
        top_acc, top_n = improvements[0]
        if top_acc >= best_acc + min_improve:
            best_subset.append(top_n)
            best_cols.append(name_to_idx[top_n])
            best_acc = top_acc
            curve.append({"k": len(best_subset), "acc": best_acc, "subset": best_subset.copy()})
            remaining = [n for n in remaining if n != top_n]
        else:
            break

    return {"best_subset": best_subset, "best_acc": best_acc, "curve": curve}


# ----------------------------- pipeline -----------------------------
def run_pipeline(
    in_csv: str,
    out_prefix: str,
    set_tokens: List[str],
    stride: int,
    model: str,
    seed: int,
    hidden: int,
    epochs: int,
    lr: float,
    weight_decay: float,
):
    df = pd.read_csv(in_csv)

    # Thin (if requested)
    if stride and stride > 1:
        df = df.iloc[::stride].copy()

    # Group labels (same syntax as your dihedral script)
    df = build_group_labels(df, set_tokens)

    # Collect features
    X_all, feat_names = collect_rmsf_features(df)

    # Labels and replica groups
    le = LabelEncoder()
    y = le.fit_transform(df["group_label"].astype(str).values)
    class_names = le.classes_.tolist()
    groups = df["replicate"].astype(str).values

    # Train
    if model.lower() == "deep":
        if not _torch_available():
            print("[warn] PyTorch not found; falling back to linear baseline.")
            res = train_linear_baseline(X_all, y, groups, n_splits=3, C=1.0, seed=seed)
            model_used = "linear"
        else:
            res = train_deep_classifier(
                X_all, y, groups,
                n_splits=3, hidden=hidden, max_epochs=epochs,
                lr=lr, weight_decay=weight_decay, seed=seed
            )
            # For importance on original features, fit a linear layer on X → classes
            # so we have W in feature space comparable to linear baseline
            # (quick refit with CV and average)
            lin = train_linear_baseline(X_all, y, groups, n_splits=3, C=1.0, seed=seed)
            res["W_mean"] = lin["W_mean"]
            model_used = "deep"
    elif model.lower() == "linear":
        res = train_linear_baseline(X_all, y, groups, n_splits=3, C=1.0, seed=seed)
        model_used = "linear"
    else:
        raise SystemExit("model must be 'deep' or 'linear'.")

    # Importance
    imp = feature_importance_from_W(res["W_mean"], feat_names)

    # ANOVA (how RMSF changes across groups): effect stats for each feature
    try:
        F, p = f_classif(X_all, y)
        anova_df = pd.DataFrame({"feature": feat_names, "F": F, "pval": p})
        anova_df = anova_df.sort_values("F", ascending=False).reset_index(drop=True)
    except Exception as e:
        print(f"[warn] f_classif failed: {e}")
        anova_df = pd.DataFrame({"feature": feat_names})

    # Outputs
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)

    imp_csv = f"{out_prefix}_rmsf_importance.csv"
    imp.to_csv(imp_csv, index=False)

    anova_csv = f"{out_prefix}_rmsf_anova.csv"
    anova_df.to_csv(anova_csv, index=False)

    # Plots (same style as your other code)
    plot_topk_bar(imp, f"{out_prefix}_top20.png", k=20, title="Top Discriminating Residues/Atoms (RMSF)")

    scores, axis_label = compute_cv1_scores(
        X_all, y, model_used, res.get("W_mean", None), class_names
    )
    hist_png = f"{out_prefix}_cv1_hist.png"
    plot_cv_hist(scores, y, class_names, hist_png)

    # Greedy subset of residues
    sel = greedy_select_best_features(
        X=X_all, y=y, groups=groups,
        feat_names=feat_names, imp_df=imp,
        max_k=12, min_improve=1e-3, seed=seed
    )
    with open(f"{out_prefix}_best_subset.txt", "w") as f:
        f.write("Greedy best RMSF feature subset (residues/atoms)\n")
        f.write(f"Mean CV accuracy: {sel['best_acc']:.4f}\n")
        f.write("Subset (order of addition):\n")
        for k, n in enumerate(sel["best_subset"], 1):
            f.write(f"  {k:2d}. {n}\n")
        f.write("\nLearning curve (k, acc):\n")
        for entry in sel["curve"]:
            f.write(f"  k={entry['k']:2d}, acc={entry['acc']:.4f}\n")
    with open(f"{out_prefix}_best_subset.json", "w") as f:
        json.dump(sel, f, indent=2)

    # Minimal CV summary
    with open(f"{out_prefix}_cv.json", "w") as f:
        json.dump({"fold_acc": res["folds"], "mean_acc": float(np.mean(res["folds"]))}, f, indent=2)

    # Frames used (for parity with your other script; if 'time' missing, omit)
    cols_keep = [c for c in ["time", "ligand", "mutant", "replicate", "group_label"] if c in df.columns]
    if cols_keep:
        df.loc[:, cols_keep].to_csv(f"{out_prefix}_frames_used.csv", index=False)

    print(f"[rmsf-deeplda] groups: {class_names}")
    print(f"[rmsf-deeplda] model: {model_used}")
    print(f"[rmsf-deeplda] CV folds: {res['folds']} (mean={np.mean(res['folds']):.3f})")
    print(f"[rmsf-deeplda] wrote: {imp_csv}")
    print(f"[rmsf-deeplda] wrote: {anova_csv}")
    print(f"[rmsf-deeplda] wrote: {out_prefix}_top20.png")
    print(f"[rmsf-deeplda] wrote: {hist_png}")
    print(f"[rmsf-deeplda] wrote: {out_prefix}_cv.json")
    print(f"[rmsf-deeplda] wrote: {out_prefix}_best_subset.txt")
    print(f"[rmsf-deeplda] wrote: {out_prefix}_best_subset.json")


# ----------------------------- CLI -----------------------------
def build_argparser():
    ap = argparse.ArgumentParser(
        description="DeepLDA-style analysis on residue/atom RMSF CSV with replica-aware CV, importance, and CV1 histogram."
    )
    ap.add_argument("--in-csv", required=True, help="Input CSV with RMSF features + ligand, mutant, replicate")
    ap.add_argument("--out-prefix", required=True, help="Output prefix (path/prefix for results)")
    ap.add_argument("--set", dest="sets", action="append", required=True,
                    help="Selection 'GROUP|LIGAND|MUTANT|rep1,rep2,...' (rep list optional or '*' for all). "
                         "Repeat to add more sets; GROUP names label the classes.")
    ap.add_argument("--stride", type=int, default=1, help="Keep every Nth row (thin data)")
    ap.add_argument("--model", choices=["deep", "linear"], default="deep", help="Deep MLP (if torch) or linear baseline")
    ap.add_argument("--seed", type=int, default=123)

    # Deep hyperparams
    ap.add_argument("--hidden", type=int, default=64, help="Hidden width (deep model)")
    ap.add_argument("--epochs", type=int, default=25, help="Training epochs (deep model)")
    ap.add_argument("--lr", type=float, default=1e-3, help="Learning rate (deep model)")
    ap.add_argument("--weight-decay", type=float, default=1e-4, help="Weight decay (deep model)")
    return ap


def main():
    ap = build_argparser()
    args = ap.parse_args()

    run_pipeline(
        in_csv=args.in_csv,
        out_prefix=args.out_prefix,
        set_tokens=args.sets,
        stride=args.stride,
        model=args.model,
        seed=args.seed,
        hidden=args.hidden,
        epochs=args.epochs,
        lr=args.lr,
        weight_decay=args.weight_decay,
    )


if __name__ == "__main__":
    main()
