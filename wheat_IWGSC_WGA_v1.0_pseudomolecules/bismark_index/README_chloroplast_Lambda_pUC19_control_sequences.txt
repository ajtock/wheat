Wheat chloroplast, Lamda and pUC19 sequences are appended to the wheat
reference genome to enable calculation of non-conversion rates and to
test whether methylation percentage is sufficiently high for the control
CpG methylated plasmid pUC19. Lambda and pUC19 sequences can be used for
these purposes where the NEBNext Enzymatic Methyl-seq Kit was used
( https://international.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit#Product%20Information ,
with "Control DNA CpG unmethylated Lambda" and "Control DNA CpG methylated pUC19").

Lambda and pUC19 sequences are obtainable from
https://international.neb.com/tools-and-resources/interactive-tools/dna-sequences-and-maps-tool

Retrieved with:
curl "https://international.neb.com/-/media/nebus/page-images/tools-and-resources/interactive-tools/dna-sequences-and-maps/text-documents/lambdafsa.txt?la=en&rev=c0c6669b9bd340ddb674ebfd9d55c691&hash=CB90FA368C445B65C29FE2F35BF47E22FD635809" > lambda_NEB.fa
curl "https://international.neb.com/-/media/nebus/page-images/tools-and-resources/interactive-tools/dna-sequences-and-maps/text-documents/puc19fsa.txt?la=en&rev=6e10f4c4a4234d638e401cd2f4578ef0&hash=F349B10F8358E1F12AC282122E17831C7EC23D67" > pUC19_NEB.fa


The wheat chloroplast sequence can be used for calculating non-coversion
rates where whole-genome bisulfite sequencing (WGBS) is applied.

The wheat chloroplast genome sequence (GenBank accession AB042240.3) is as described in Ogihara et al. (2002) Mol. Genet. Genomics 266: 740-746.
Retrieved with:
wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=AB042240.3&rettype=fasta" -O wheat_CS_chloroplast_genome.fa

The wheat mitochondrial genome sequence (GenBank accession AP008982.1) is as described in Ogihara et al. (2005) Nucleic Acids Res. 33: 6235-6250.
Retrieved with:
wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=AP008982.1&rettype=fasta" -O wheat_CS_mitochondrial_genome.fa

