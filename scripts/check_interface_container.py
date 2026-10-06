"""Check packaged pages, in-memory exports and download isolation inside the Linux container."""
from pathlib import Path
from importlib.util import find_spec
import hashlib
import os
import sys
import zipfile
from io import BytesIO
from unittest.mock import patch
import streamlit as st

sys.path.insert(0, "/app/gui")
import matplotlib.pyplot as plt
import seaborn
from streamlit.testing.v1 import AppTest
from data_loader import GTEX_VERSION, COMPLETE_AGE_BINS

assert GTEX_VERSION == "v10" and COMPLETE_AGE_BINS
assert os.getuid() != 0
assert sorted(p.name for p in Path('/app/output').iterdir()) == ['v10_complete', 'v8_complete']
for version in ('v10','v8'):
    assert sorted(p.name for p in Path(f'/app/output/{version}_complete').iterdir()) == ['jaccard','normalized','sets']
for excluded in ("scripts", "statistical_validation_cli", "data", "reproducibility/reference", "gui/_local", "gui/_archive"):
    assert not (Path("/app") / excluded).exists(), excluded
assert find_spec("statistical_validation_cli") is None
assert (Path("/app/reproducibility/README.md")).is_file()
assert (Path("/app/requirements-reproduction.lock.txt")).is_file()
assert not Path('/app/gui/data/ad_panel_top1000.txt').exists()
atlas = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in Path('/app/output').rglob('*') if p.is_file()}
assert not os.access('/app/output/v10_complete/sets/Artery_Aorta_sets.txt', os.W_OK)
assert len(list(Path("/app/output/v10_complete/normalized").glob("*.parquet"))) == 50
assert len(list(Path("/app/output/v8_complete/normalized").glob("*.parquet"))) == 49

assert not Path('/app/gui/output').exists()
server_files = {p: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in Path('/app').rglob('*') if p.is_file()}

app = AppTest.from_file("/app/gui/main.py", default_timeout=60).run()
assert not app.exception, [e.message for e in app.exception]
buttons = [b.label for b in app.button]
assert not any("AD Threshold" in label for label in buttons)
labels = [b.label for b in app.get("download_button")]
assert "Download reproduction guide" in labels
assert "Download dependency versions" in labels
for page in ("stamp_generator", "upload_analysis", "tissue_comparison", "multi_tissue", "age_specific", "gene_sharing", "single_gene", "group_comparison", "panel_explorer"):
    app.session_state.current_page = page
    app.run()
    assert not app.exception, (page, [e.message for e in app.exception])
    assert not app.error, (page, [e.value for e in app.error])
    print("PASS page", page)

# Warm up library font caches before checking /tmp and the home directory
# for writes caused by generation and downloads.
writable_roots = (Path('/tmp'), Path('/home/stamp'))
def writable_snapshot():
    return {p: hashlib.sha256(p.read_bytes()).hexdigest()
            for root in writable_roots for p in root.rglob('*') if p.is_file()}
writable_files = writable_snapshot()

app.session_state.current_page = "stamp_generator"
app.run()
app.slider[0].set_value(0.7)
app.run()
for control in app.selectbox:
    if "Choose one tissue" in control.label:
        control.select("Artery_Aorta")
        break
app.run()
for button in app.button:
    if "Generate STAMP Files" in button.label:
        button.click()
        break
app.run()
assert not app.exception, [e.message for e in app.exception]
assert not app.error, [e.value for e in app.error]
assert app.session_state.processed
exports = dict(app.session_state.generator_exports)
assert len(exports) == 2
assert all(len(e['data'].decode().splitlines()) == 5 for e in exports.values())
assert all(not e['data'].split() for e in exports.values())
assert app.session_state.generator_threshold == 0.7
assert all('tau_0.7' in value for value in app.dataframe[-1].value['File Name'])
for label in ('Select All', 'Create ZIP Archive'):
    next(b for b in app.button if label in b.label).click()
    app.run()
    assert not app.exception and not app.error
with zipfile.ZipFile(BytesIO(app.session_state.generator_zip[1])) as archive:
    assert set(archive.namelist()) == set(exports)
    for name, export in exports.items():
        assert archive.read(name) == export['data']
assert any('Download ZIP Archive' in b.label for b in app.get('download_button'))
# A normal rerun keeps the prepared ZIP available.
app.run()
assert any('Download ZIP Archive' in b.label for b in app.get('download_button'))
# The second session must have its own exports.
second = AppTest.from_file('/app/gui/main.py', default_timeout=60)
second.session_state.current_page = 'stamp_generator'
second.run()
for control in second.selectbox:
    if 'Choose one tissue' in control.label:
        control.select('Artery_Aorta')
second.run()
for button in second.button:
    if 'Generate STAMP Files' in button.label:
        button.click()
second.run()
assert not second.exception and not second.error
assert any(e['data'].split() for e in second.session_state.generator_exports.values())
assert app.session_state.generator_exports == exports
# Generating again replaces exports and clears the old ZIP.
app.slider[0].set_value(0.5)
app.run()
next(b for b in app.button if 'Generate STAMP Files' in b.label).click()
app.run()
assert not app.exception and not app.error
assert any(e['data'].split() for e in app.session_state.generator_exports.values())
assert all('tau_0.5' in name for name in app.session_state.generator_exports)
assert 'generator_zip' not in app.session_state
for label in ('Select All', 'Create ZIP Archive'):
    next(b for b in app.button if label in b.label).click()
    app.run()
assert 'generator_zip' in app.session_state
# Reset one session while preserving the other session's exports.
second_exports = dict(second.session_state.generator_exports)
next(b for b in app.button if 'Start New Generation' in b.label).click()
with patch.object(st, 'download_button', wraps=st.download_button) as download:
    app.run()
    assert download.call_count == 0, 'Reset must not register the previous ZIP again'
assert not app.get('download_button')
assert not app.exception and not app.error
assert not app.session_state.processed
assert 'generator_exports' not in app.session_state
assert 'generator_zip' not in app.session_state
assert second.session_state.generator_exports == second_exports
assert atlas == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in Path('/app/output').rglob('*') if p.is_file()}
assert server_files == {p: hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in Path('/app').rglob('*') if p.is_file()}
assert not Path('/app/gui/output').exists()
assert writable_files == writable_snapshot(), 'Generation wrote into /tmp or /home/stamp'
print("PASS memory-only exports, ZIP contents, isolated sessions, replacement and reset; no files created or changed under /app, /tmp or /home/stamp")
plt.close("all")
print("PASS packaged interface, absent private pages and absent statistical CLI")
