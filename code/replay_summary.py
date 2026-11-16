from pathlib import Path
import io, zipfile
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

root = Path(__file__).resolve().parents[1]
with zipfile.ZipFile(root / "derived_data.zip") as z:
    def table(name):
        return pd.read_csv(io.BytesIO(z.read("results/" + name + ".tsv")), sep="\t")
    lineage = table("larry_bootstrap")
    print(lineage[(lineage["mode"] == "original") & (lineage["model"] == "ridge")][["dataset", "n_clones", "brier", "ci_low", "ci_high"]].to_string(index=False))
    print(table("axis_permutation_summary")[["target", "observed", "null_mean", "qvalue"]].to_string(index=False))
    clinical = table("clinical_primary_bootstrap")
    print(clinical[clinical.metric == "auroc"][["cohort", "model", "estimate", "ci_low", "ci_high"]].to_string(index=False))
    predictions = table("ispy2_arm_adjustment_predictions")
    scores = predictions.pivot(index="patient", columns="model", values="predicted")
    y = predictions.drop_duplicates("patient").set_index("patient").loc[scores.index, "outcome"]
    increment = roc_auc_score(y, scores["hr_her2_plus_arm_plus_five_gene"]) - roc_auc_score(y, scores["hr_her2_plus_arm"])
    paired = table("ispy2_arm_adjustment_paired_bootstrap_increments")
    row = paired[(paired.model == "hr_her2_plus_arm_plus_five_gene") & (paired.baseline == "hr_her2_plus_arm") & (paired.metric == "auroc")].iloc[0]
    assert np.isclose(increment, row.delta, atol=1e-12)
    print("I-SPY2 five-gene AUROC increment after regimen adjustment:", increment, (row.ci_low, row.ci_high))
    print(table("mixseq_molecular_compact").to_string(index=False))
    survival = table("clinical_survival")
    assert set(survival["n"]) == {317, 303}
