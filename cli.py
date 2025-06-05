import typer 

app = typer.Typer() 

@app.command() 
def extract_genes(vcf_path: str, output_path: str):
    """Extract genes from VCF and save to file."""