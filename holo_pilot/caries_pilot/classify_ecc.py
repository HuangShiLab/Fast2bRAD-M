#!/usr/bin/env python3
"""Compare microbiome-only, host-only, and combined classifiers for ECC."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder, StandardScaler


def encode_labels(y: pd.Series) -> tuple[np.ndarray, LabelEncoder]:
    le = LabelEncoder()
    return le.fit_transform(y), le


def evaluate(X: pd.DataFrame, y: np.ndarray, name: str, n_splits: int = 5) -> dict:
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    aucs = []
    accs = []
    importances = np.zeros(X.shape[1])
    for train_idx, test_idx in cv.split(X, y):
        clf = RandomForestClassifier(
            n_estimators=500, max_depth=None, min_samples_leaf=2,
            class_weight="balanced", n_jobs=-1, random_state=42,
        )
        clf.fit(X.iloc[train_idx], y[train_idx])
        prob = clf.predict_proba(X.iloc[test_idx])[:, 1]
        pred = clf.predict(X.iloc[test_idx])
        aucs.append(roc_auc_score(y[test_idx], prob))
        accs.append(accuracy_score(y[test_idx], pred))
        importances += clf.feature_importances_
    importances /= n_splits
    return {
        "model": name,
        "AUC_mean": float(np.mean(aucs)),
        "AUC_std": float(np.std(aucs)),
        "Accuracy_mean": float(np.mean(accs)),
        "Accuracy_std": float(np.std(accs)),
        "importances": pd.Series(importances, index=X.columns).sort_values(ascending=False),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--feature-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    fdir = Path(args.feature_dir)
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    X_micro = pd.read_csv(fdir / "X_microbiome.tsv", sep="\t", index_col=0)
    X_host = pd.read_csv(fdir / "X_host.tsv", sep="\t", index_col=0)
    y = pd.read_csv(fdir / "y.tsv", sep="\t", index_col=0)["phenotype"]

    # Impute missing host SNPs (coded as -1) with mean dosage
    imp = SimpleImputer(strategy="mean")
    X_host_imp = pd.DataFrame(
        imp.fit_transform(X_host), index=X_host.index, columns=X_host.columns
    )

    # Optional: filter low-variance / sparse microbiome features
    keep = X_micro.columns[(X_micro > 0).sum() >= 3]
    X_micro = X_micro[keep]

    y_enc, le = encode_labels(y)

    results = []
    results.append(evaluate(X_micro, y_enc, "microbiome_only"))
    results.append(evaluate(X_host_imp, y_enc, "host_only"))

    combined = pd.concat([X_micro, X_host_imp], axis=1)
    results.append(evaluate(combined, y_enc, "microbiome_plus_host"))

    summary = pd.DataFrame([
        {
            "model": r["model"],
            "AUC_mean": r["AUC_mean"],
            "AUC_std": r["AUC_std"],
            "Accuracy_mean": r["Accuracy_mean"],
            "Accuracy_std": r["Accuracy_std"],
        }
        for r in results
    ])
    summary.to_csv(outdir / "model_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))

    for r in results:
        r["importances"].head(30).to_csv(
            outdir / f"feature_importance_{r['model']}.tsv", sep="\t"
        )


if __name__ == "__main__":
    main()
