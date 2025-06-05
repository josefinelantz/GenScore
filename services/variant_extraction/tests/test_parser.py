import json
from services.variant_extraction.main import extract_gene_variants

def test_parser_creates_expected_output(path):
    vcf_path = "/Users/lantzan/python/scilifelab/GenScore/data/SNV.somatic.sweetelf.merged.clinical.ranked.vcf.gz"
    output_file = tmp_path / "genes.json"

    extract_gene_variants(vcf_path, output_file)

    with open(output_file) as f:
        data = json.load(f)
        assert "GNB1" in data or "AJAP1" in data
        assert isinstance(data, dict)