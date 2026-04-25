#!/bin/bash
#SBATCH --mem 16G
#SBATCH --partition uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/projecy_dark_genes/logs/outputs
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

# sbatch /uoa/home/r02hw22/sharedscratch/scripts/genome_inspect.sh

cd /uoa/home/r02hw22/sharedscratch/project_dark_genes/00_raw

# Look at headers
grep "^>" equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa | head
grep "^>" equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa | head

# Check sequence lengths
awk '
/^>/{
    if (seqlen){print seqlen}
    seqlen=0
    next
}
{seqlen+=length($0)}
END{if(seqlen) print seqlen}
' equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa | sort -n | awk '
BEGIN{count=0; sum=0}
{a[++count]=$1; sum+=$1}
END{
    print "count="count
    print "min="a[1]
    print "max="a[count]
    print "mean="sum/count
}'

# nuncleotide file
	# count=55607
	# min=94
	# max=31372
	# mean=1354.16


# amino acid file
awk '
/^>/{
    if (seqlen){print seqlen}
    seqlen=0
    next
}
{seqlen+=length($0)}
END{if(seqlen) print seqlen}
' equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa | sort -n | awk '
BEGIN{count=0; sum=0}
{a[++count]=$1; sum+=$1}
END{
    print "count="count
    print "min="a[1]
    print "max="a[count]
    print "mean="sum/count
}'

# count=55607
# min=32
# max=10458
# mean=452.05

# QC the protein set
awk '
/^>/{
    if (name!="") print name"\t"seqlen
    name=substr($0,2)
    seqlen=0
    next
}
{seqlen+=length($0)}
END{
    if (name!="") print name"\t"seqlen
}
' equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa > ../01_qc/protein_lengths.tsv

# mark short proteins
awk '$2 < 50' ../01_qc/protein_lengths.tsv  > ../01_qc/proteins_lt50aa.tsv
awk '$2 < 100' ../01_qc/protein_lengths.tsv > ../01_qc/proteins_lt100aa.tsv
awk '$2 < 150' ../01_qc/protein_lengths.tsv > ../01_qc/proteins_lt150aa.tsv

# cds lengths
awk '
/^>/{
    if (name!="") print name"\t"seqlen
    name=substr($0,2)
    seqlen=0
    next
}
{seqlen+=length($0)}
END{
    if (name!="") print name"\t"seqlen
}
' equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa > ../01_qc/cds_lengths.tsv

