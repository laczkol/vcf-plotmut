# vcf-plotmut
Quick and dirty python script to plot the different types of mutations by sample. The input is a `.vcf` file.

## Dependencies
### Python libraries
- `argparse` to parse CLI arguments
- `scikit-allel` to parse `.vcf` files
- `pandas` and `numpy` to create and manipulate data tables
- `matplotlib` to plot the results

## Terminology
- homref: homozygous position that matches the reference sequence
- het: heterozygous position
- homalt: homozygous position that does not match the reference sequence

## Usage
There are two positional arguments to be supplied. The first is the input `.vcf` file and the second is the stem name of the output files.
e.g. `python plotmut.py input.vcf outstem`

## Output

### `outstem_{count,ratio}.tsv`
Two tab-separated tables, one of which contains the raw counts of different mutations, and the other contains their ratios. The homref, het and homalt values are relative to all positions. The ratios in the rest of columns are relative to the number of mutations (i.e. homref is exluded). The number of total and variable positions are counted for each sample and is used to assess the ratios of different mutations.

### `outstem_{counts,hom_het_ratio,Tis_ratio,Tvs_ratio}.pdf`
`.pdf` files showing the count, the ratio of hom and het regions relative to all positions, the ratio of different transitions and the ratio of different transversions relative to variable positions.

## Limitations
- Not very fast
- Only uses diploid sites and SNPs (i.e. indels and MNPs are excluded)

