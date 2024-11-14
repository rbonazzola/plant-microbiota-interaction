# Plant-microbiota interaction

## Bacterial genomes
### Genomic sequence alignment
Use BWA to perform alignment.
- Input: FASTQ and FASTA files.
- Output: SAM files.

### Variant calling
Use BCFtools to perform variant calling from the SAM files.
- Input: SAM files.
- Intermediate files: BAM and BAM indices.
- Output: VCF files.

### Analysis of evolution


___
## Phenotype analysis
**Description of experiments:** 
- Two treatments: MS and K (three replicates each).
- Control: NB (no bacteria)

### Plant phenotypes
Examine principal root length, number of lateral roots and lateral root density (the ratio between the first two) as a function of the generation.

### Bacterial fitness
Examine optical density (OD) as a function of time, for each generation.
