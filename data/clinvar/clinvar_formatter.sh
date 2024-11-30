#!/bin/bash
# Script that formats clinvar dataset to fit to MIP pipeline v12.0.0 prerequisites to allow variant annotation.
# 1. Strips and removes all clinvar annotations to have a clean dataset.
# 2. Re-adds and renames Clinvar INFO/CLNSIG as INFO/CLINVAR_GROUND_TRUTH as training label
# 3. Adds FORMAT/ID dummy data to enable pipeline processing
# 4. Adds dummy sample info column name, 'NOTASAMPLE'
set -e
set -x
export DATAFILE=`realpath $1`
stat $DATAFILE &>/dev/null
export OUTFILE=`echo $DATAFILE | sed 's/\.vcf//g'`-cleaned.vcf

# Setup header
cat <<EOF > $OUTFILE
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=2>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=MT>
##contig=<ID=X>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
EOF
echo -e "##INFO=<ID=CLINVAR_GROUND_TRUTH,Number=.,Type=String,Description=\"Clinical significance ground truth from CLINVAR database.\">" >> $OUTFILE

echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNOTASAMPLE" >> $OUTFILE

# Filter data
## NOTE to self: You can do multiple INFO fields addition, like so: \t%INFO/CLNSIG;%INFO/CLNVC
bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t.\t.\tCLINVAR_GROUND_TRUTH=%INFO/CLNSIG\tGT:DP:AD:GQ\t1/1:30:4,26:38\n' $DATAFILE >> $OUTFILE

# Compress
bgzip -c -l 9 -@ 16 $OUTFILE > $OUTFILE.gz

# Add index
tabix -f -p vcf $OUTFILE.gz > $OUTFILE.gz.tbi

echo Output file: $OUTFILE
