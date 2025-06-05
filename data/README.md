# Genmod Score Evaluatioin

## Overview

This project evaluates the correctnes off Genmod scores for somatic variants in cancer genomics VCF files. It parses VCF files, groups variants based on clinical significance, and visualizes the distribution and composition of Genmode scores. 

## Directory Structure 

genmod_evaluation/
├── data/
├── scripts/
├── tests/
├── docs/
├── requirements.txt
├── setup.py
└── main.py

## Requirements

- Python ^3.6
- Packages listed in `requirements.txt`

## Setup

1. **Clone the repository:**

```bash
git clone <repo_url>
cd GenScore
```

2. Create and activate a virtual environment 

```bash 
python3 -m venv venv 
source venc/bin/activate # windows: venv\Scripts\activate
```

3. Install dependencies

```bash
pip install -r requirements.txt
```

# Usage

```bash
python3 main.py path_to_vcf.vcf --output results/ 
```
This will generate boxplot.png and stackedplot.png in the specified output directory. 

# Testing 

```bash
pytest
``` 

# Recommendations for Further Development 

- Error Handling: Implement robust error handling for file parsin and data processing. 
- Modularity: Refactor scripts into smaller, reusable modules. 
- Configuration: Use configuration files (e.g YAML) for pararmeter settings. 
- Integration: Develop the scripts as a package to integrate seamlessly with other pipelines. 
- Documentation: Expand documentation with usage examples and detailed explanations. 

# Licence

MIT Licence # add link here