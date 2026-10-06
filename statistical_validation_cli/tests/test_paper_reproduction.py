from pathlib import Path
import json

import pandas as pd
import pytest

from statistical_validation_cli import paper
from statistical_validation_cli.verify_paper import compare_frames, gene_digest, gene_table_digest, verify_outputs
from statistical_validation_cli.common import config_matches


def test_same_total_does_not_hide_changed_tissue_or_p_value():
    expected = pd.DataFrame({"tissue": ["Liver", "Ovary"], "count": [2254, 3300], "p": [.0026, .0072]})
    swapped = expected.copy()
    swapped["count"] = swapped["count"].iloc[::-1].to_numpy()
    with pytest.raises(AssertionError):
        compare_frames(swapped, expected, ["tissue"])
    changed = expected.copy()
    changed.loc[0, "p"] += .0001
    with pytest.raises(AssertionError):
        compare_frames(changed, expected, ["tissue"])


def test_missing_q_value_and_duplicate_tissue_are_rejected():
    expected = pd.DataFrame({"tissue": ["Liver"], "q": [.0216666666666667]})
    with pytest.raises(AssertionError):
        compare_frames(expected.assign(q=float("nan")), expected, ["tissue"])
    with pytest.raises(AssertionError, match="Duplicate"):
        compare_frames(pd.concat([expected, expected]), expected, ["tissue"])


def test_gene_identity_is_checked_even_when_stable_gene_total_matches():
    frame = pd.DataFrame({"gene_id": ["ENSG1", "ENSG2"], "observed_switching": [True, True],
                          "bootstrap_switch_frequency": [.9, .8],
                          "same_bracket_and_direction_frequency": [.9, .7]})
    assert gene_digest(frame) == gene_digest(frame.iloc[::-1])
    changed = frame.copy()
    changed["same_bracket_and_direction_frequency"] = [.7, .9]
    assert gene_digest(frame) != gene_digest(changed)


def test_observed_direction_and_smote_gene_frequencies_are_checked():
    calls = pd.DataFrame({"gene_id": ["ENSG1"], "switch_bracket": ["30-39"], "direction": ["up"]})
    assert gene_table_digest(calls) != gene_table_digest(calls.assign(direction="down"))
    stability = pd.DataFrame({"gene_id": ["ENSG1", "ENSG2"], "smote_switch_frequency": [.8, .9]})
    assert gene_table_digest(stability) == gene_table_digest(stability.iloc[::-1])
    changed = stability.copy()
    changed['smote_switch_frequency'] = [.9, .8]
    assert gene_table_digest(stability) != gene_table_digest(changed)


def test_modified_reference_cannot_be_accepted(tmp_path, monkeypatch):
    reference = tmp_path / "reference"
    reference.mkdir()
    table = reference / "primary.csv"
    table.write_text("tissue,count\nLiver,2254\n")
    manifest = {"files": {"primary.csv": paper.sha256(table)}}
    (reference / "manifest.json").write_text(json.dumps(manifest))
    table.write_text("tissue,count\nLiver,2255\n")
    monkeypatch.setattr(paper, "REFERENCE", reference)
    with pytest.raises(RuntimeError, match="modified"):
        paper.check_reference()


def test_incomplete_run_is_never_reported_as_a_full_replication(tmp_path):
    with pytest.raises(FileNotFoundError, match="Incomplete reproduction"):
        verify_outputs(tmp_path, tmp_path)


def test_input_hash_mismatch_is_rejected(tmp_path, monkeypatch):
    root = tmp_path / "project"
    (root / "reproducibility").mkdir(parents=True)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    path = inputs / "metadata.parquet"
    path.write_bytes(b"reference")
    manifest = {"full": {"files": {path.name: {"bytes": path.stat().st_size, "sha256": paper.sha256(path)}}}}
    (root / "reproducibility/inputs.json").write_text(json.dumps(manifest))
    path.write_bytes(b"different")
    monkeypatch.setattr(paper, "ROOT", root)
    with pytest.raises(RuntimeError, match="SHA-256"):
        paper.check_inputs(inputs)


@pytest.mark.parametrize("force", [False, True])
def test_every_full_run_command_is_accepted_by_its_actual_parser(tmp_path, force):
    plan = paper.full_command_plan(tmp_path, jobs=2, force=force)
    paper.validate_command_plan(plan)
    paper.validate_threshold_command(tmp_path)
    smote = [arguments for module, arguments in plan if module == "smote"]
    assert len(smote) == 2
    assert all("--allow-expensive-test" in arguments for arguments in smote)


def test_copied_checkpoints_without_provenance_cannot_be_certified(tmp_path):
    (tmp_path / "copied-checkpoint.json").write_text('{"count": 123456}')
    with pytest.raises(RuntimeError, match="no provenance"):
        paper.prepare_output(tmp_path, {"inputs": "current"}, False)
    assert not (tmp_path / "paper_provenance.json").exists()
    paper.prepare_output(tmp_path, {"inputs": "current"}, True)


def test_changed_inputs_or_code_reject_resume_and_old_success_is_removed(tmp_path):
    old = {"inputs": "original", "source": "original"}
    paper.prepare_output(tmp_path, old, False)
    marker = tmp_path / "paper_verified.json"
    marker.write_text("old success")
    with pytest.raises(RuntimeError, match="changed"):
        paper.prepare_output(tmp_path, {**old, "source": "changed"}, False)
    paper.prepare_output(tmp_path, old, False)
    assert not marker.exists()


def test_parameter_match_alone_does_not_reuse_a_checkpoint():
    config = {"seed": 42}
    assert not config_matches({"config": config}, config)
    config = {"seed": 42, "provenance": {"inputs": "original", "code": "original"}}
    assert config_matches({"config": config}, config)
    assert not config_matches({"config": config}, {**config, "provenance": {"inputs": "changed"}})


def test_execution_summary_distinguishes_resumed_and_fresh_tissues(tmp_path):
    for label, step in (("paper", "02_primary"), ("paper", "04_bootstrap"),
                        ("paper-smote-fixed", "03_smote"), ("paper-smote-adaptive", "03_smote")):
        folder = tmp_path / "v10/statistical_validation_cli/runs" / label / step / "execution"
        folder.mkdir(parents=True)
        for index in range(50):
            (folder / f"{index}.json").write_text(json.dumps({"tissue": str(index),
                           "status": "recomputed" if index == 0 else "reused"}))
    summary = paper.execution_summary(tmp_path)
    assert all(len(value["recomputed"]) == 1 and len(value["reused"]) == 49 for value in summary.values())


def test_reference_p_and_bh_checks_do_not_call_the_inference_helpers(monkeypatch):
    from statistical_validation_cli import common
    from statistical_validation_cli.verify_paper import validate_reference_statistics
    def fail(*args, **kwargs):
        raise AssertionError("Inference helper must not validate itself")
    monkeypatch.setattr(common, 'bh_qvalues', fail)
    monkeypatch.setattr(common, 'empirical_right_tail_p', fail)
    validate_reference_statistics(paper.REFERENCE)
