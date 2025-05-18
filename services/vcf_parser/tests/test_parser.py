import json
from app.parser import extract_gene_variants

def test_parser_creates_expected_output(path):
    vcf_path = ""
