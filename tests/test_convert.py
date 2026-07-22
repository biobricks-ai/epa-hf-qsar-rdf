import gzip,importlib.util,json
from pathlib import Path
R=Path(__file__).parents[1];s=importlib.util.spec_from_file_location('c',R/'stages/convert.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);SRC=Path('/mnt/raid2/biobricks/catalog-qsar-dataset-catalog-data-gov/brick')
def test_full_coverage(tmp_path):
 c=m.convert(SRC,tmp_path);assert c['source_rows']==6769 and c['represented_rows']==6769 and c['row_coverage']==c['mapped_cell_coverage']==1;assert len(c['tables'])==11 and all(t['field_disposition_coverage']==1 for t in c['tables'].values())
def test_prediction_reference_separation(tmp_path):
 c=m.convert(SRC,tmp_path);x=gzip.open(tmp_path/'epa-hf-qsar.nt.gz','rt').read();assert c['predicted_datums']>0 and c['reference_datums']>0;assert 'model-derived prediction or decision-analysis value' in x and 'experimental or authoritative reference value' in x
def test_n_hexane(tmp_path):
 m.convert(SRC,tmp_path);x=gzip.open(tmp_path/'epa-hf-qsar.nt.gz','rt').read();assert 'https://identifiers.org/cas/110-54-3' in x and '437.672' in x and 'Confidence Classification' not in x
def test_reports():assert json.loads((R/'health/license.json').read_text())['status']=='pass' and json.loads((R/'health/field-disposition.json').read_text())['all_fields_dispositioned']
