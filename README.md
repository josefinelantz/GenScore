# Project Structure 

GenScore/
├── data/                # Sample VCFs, annotations, test sets
├── genscore/            # Python package with modules
│   ├── __init__.py
│   ├── io.py            # Input parsing (VCF reader)
│   ├── gene_mapper.py   # Map variants to genes
│   ├── scorer.py        # Score genes by variant effect
│   ├── pathway.py       # Map genes to pathways
│   └── viz.py           # Create pathway visualizations
├── tests/               # Unit tests (pytest)
├── cli.py               # CLI entry point
├── requirements.txt
├── Dockerfile
└── README.md