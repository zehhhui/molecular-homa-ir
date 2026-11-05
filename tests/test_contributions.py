"""Scientific invariants of the contribution calculation, using artificial data."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "treatment_contributions", ROOT / "scripts" / "treatment_contributions.py"
)
contribution = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(contribution)


def artificial_trial():
    rng = np.random.default_rng(7103)
    participants = pd.DataFrame({
        "arm": np.repeat(["CG32", "CG48", "PB"], 30),
        "sex": np.tile([1, 2], 45),
        "age": rng.integers(35, 72, 90),
        "baseline_mol": rng.normal(1.5, 0.5, 90),
    })
    design = contribution.build_design(participants, "PB")
    beta = np.array([
        [0.1, 0.2, -0.1, 0.0],
        [-0.3, -0.2, 0.1, 0.0],
        [1.0, 0.7, -0.4, 0.0],
        [0.2, -0.3, 0.1, 0.0],
        [0.1, 0.2, 0.2, 0.0],
        [0.1, -0.1, 0.1, 0.0],
    ])
    outcomes = design.to_numpy() @ beta
    outcomes[:, :3] += rng.normal(0, 0.08, (90, 3))
    return participants, design, outcomes, beta


class LoweringShareTests(unittest.TestCase):
    def test_denominator_is_gross_lowering_not_net_effect(self):
        # Net reduction is 2.2, but the lowering denominator must be 3.0.
        effects = np.array([[-2.0, -1.0, 0.8]])
        shares = contribution.lowering_shares(effects)
        np.testing.assert_allclose(shares, [[200 / 3, 100 / 3, 0]])
        np.testing.assert_allclose(shares.sum(axis=1), 100)

    def test_display_filter_does_not_renormalize_contributions(self):
        # A displayed subset does not explain 100% of the full denominator.
        effects = np.array([-2.0, -1.0, -1.0, 0.5])
        all_shares = contribution.lowering_shares(effects)
        np.testing.assert_allclose(all_shares[:2], [50, 25])
        self.assertEqual(all_shares[:2].sum(), 75)

    def test_no_lowering_and_nonfinite_effects_are_rejected(self):
        for effects in ([0, 1, 2], [0, 0], [-1, np.nan]):
            with self.subTest(effects=effects):
                with self.assertRaises(ValueError):
                    contribution.lowering_shares(effects)


class EligibilityTests(unittest.TestCase):
    @staticmethod
    def records():
        return pd.DataFrame([
            {"feature_id": feature, "repetition": repeat, "n_selected": count,
             "n_folds": 10}
            for feature, counts in [("stable", [1, 1, 1]),
                                    ("many_but_inconsistent", [10, 10, 0]),
                                    ("all_folds", [10, 10, 10])]
            for repeat, count in enumerate(counts, 1)
        ])

    def test_at_least_one_fold_in_every_repetition(self):
        actual = contribution.eligible_features(self.records()).set_index("feature_id")
        self.assertTrue(actual.loc["stable", "eligible"])
        self.assertFalse(actual.loc["many_but_inconsistent", "eligible"])
        self.assertTrue(actual.loc["all_folds", "eligible"])
        stringent = contribution.eligible_features(self.records(), min_folds=8)
        self.assertEqual(stringent.loc[stringent.eligible, "feature_id"].tolist(), ["all_folds"])

    def test_missing_repeat_is_not_treated_as_selected(self):
        missing = self.records().iloc[1:]
        with self.assertRaises(ValueError):
            contribution.eligible_features(missing)

    def test_duplicate_records_are_rejected(self):
        records = self.records()
        with self.assertRaises(ValueError):
            contribution.eligible_features(pd.concat([records, records.iloc[[0]]]))


class ContrastTests(unittest.TestCase):
    def test_molecular_effects_sum_to_total_score_contrast(self):
        _, design, outcomes, _ = artificial_trial()
        contrasts = contribution.contrast_matrix(design.columns, "PB")
        individual = contribution.fit_contrasts(design, outcomes, contrasts)
        total = contribution.fit_contrasts(design, outcomes.sum(axis=1), contrasts)
        np.testing.assert_allclose(individual["effect"].sum(axis=1),
                                   total["effect"].ravel(), atol=1e-12)

    def test_contrast_sign_and_order_against_known_effects(self):
        _, design, _, beta = artificial_trial()
        contrasts = contribution.contrast_matrix(design.columns, "PB")
        actual = contribution.fit_contrasts(design, design.to_numpy() @ beta, contrasts)
        np.testing.assert_allclose(actual["effect"],
            [[-1.0, -0.7, 0.4, 0.0], [-1.3, -0.9, 0.5, 0.0]], atol=1e-12)

    def test_seeded_arm_bootstrap_preserves_joint_decomposition(self):
        participants, design, outcomes, _ = artificial_trial()
        args = (design, outcomes, participants.arm, "PB", 75, 9107)
        first = contribution.bootstrap_contrasts(*args)
        second = contribution.bootstrap_contrasts(*args)
        self.assertEqual(first["effects"].shape, (75, 2, 4))
        np.testing.assert_array_equal(first["effects"], second["effects"])
        np.testing.assert_array_equal(first["shares"], second["shares"])
        np.testing.assert_allclose(first["shares"].sum(axis=2), 100, atol=1e-12)
        np.testing.assert_array_equal(first["effects"][:, :, 3], 0)
        # Same resampling of all molecules must reproduce the summed outcome.
        total = contribution.bootstrap_contrasts(
            design, outcomes.sum(axis=1, keepdims=True), participants.arm, "PB", 75, 9107)
        np.testing.assert_allclose(first["effects"].sum(axis=2),
                                   total["effects"][:, :, 0], atol=1e-12)
        np.testing.assert_allclose(first["se"], first["shares"].std(axis=0, ddof=1))


if __name__ == "__main__":
    unittest.main()
