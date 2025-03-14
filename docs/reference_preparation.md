# Reference Preparation Workflow

### 1. Download the FASTA and GFF files from the [*P. patens* Gene Model Lookup Database](https://peatmoss.plantcode.cup.uni-freiburg.de/ppatens_db/downloads.php) and the CORE collection of plant motiffs from [JASPAR2022](https://jaspar2022.genereg.net/)
[Accessed 2024/09/04]

FASTA: https://peatmoss.plantcode.cup.uni-freiburg.de/downloads/Genome%20sequences/P.patens_v6/P.patens_v6_genome.fna.gz

GFF: https://peatmoss.plantcode.cup.uni-freiburg.de/downloads/Annotations/P.patens_v6/P.patens_genome_annotation_v6.gff.gz

JASPAR Motifs: https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt

### 2. Remove Unneeded Attributes
Command:
```bash
agat_sp_manage_attributes.pl --gff P.patens_genome_annotation_v6.gff --tag all_attributes --outfile P.patens_genome_annotation_v6_essentialAttr.gff
```
### 3. Convert to GTF format
Command:
```bash
agat_convert_sp_gff2gtf.pl --gff P.patens_genome_annotation_v6_essentialAttr.gff --outfile P.patens_genome_annotation_v6.gtf
```
### 4. Generate Reference with 10x Makeref
Create ppatensv6.config file:
```json
{
   organism: "P. patens",
   genome : ["Ppatens_v6"],
   input_fasta : ["../Physcogenome/Version6/P.patens_v6_genome.fna"],
   input_gtf : ["../Physcogenome/Version6/P.patens_genome_annotation_v6.gtf"],
   input_motifs : "../Physcogenome/Motifs/JASPAR2022_CORE_plants_non-redundant_pfms_jaspar.txt"
}
```
Command:
```bash
10xGenomics/cellranger-arc-2.0.2/cellranger-arc mkref --config=../ppatensv6.config --memgb=40
```






