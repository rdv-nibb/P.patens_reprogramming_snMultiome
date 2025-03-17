# Reference Preparation and Cellranger-ARC workflow

### 1. Download the FASTA and GFF files from the [*P. patens* Gene Model Lookup Database](https://peatmoss.plantcode.cup.uni-freiburg.de/ppatens_db/downloads.php) and the CORE collection of plant motiffs from [JASPAR2022](https://jaspar2022.genereg.net/)
[Accessed 2024/09/04]

FASTA: [P.patens_v6_genome.fna.gz](https://peatmoss.plantcode.cup.uni-freiburg.de/downloads/Genome%20sequences/P.patens_v6/P.patens_v6_genome.fna.gz)

GFF: [P.patens_genome_annotation_v6.gff.gz](https://peatmoss.plantcode.cup.uni-freiburg.de/downloads/Annotations/P.patens_v6/P.patens_genome_annotation_v6.gff.gz)

JASPAR Motifs: [JASPAR2022_CORE_non-redundant_pfms_jaspar.txt](https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt)

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
### 4. Generate Reference with 10x mkref
Create ppatensv6.config file:
```
{
   organism: "P. patens"
   genome : ["Ppatens_v6"]
   input_fasta : ["/mnt/gpfsA/home/devil/Physcogenome/Version6/P.patens_v6_genome.fna"]
   input_gtf : ["/mnt/gpfsA/home/devil/Physcogenome/Version6/P.patens_genome_annotation_v6.gtf"]
   input_motifs : "/mnt/gpfsA/home/devil/Physcogenome/Motifs/JASPAR2022_CORE_plants_non-redundant_pfms_jaspar.txt"
}
```
Command:
```bash
../10xGenomics/cellranger-arc-2.0.2/cellranger-arc mkref --config=../ppatensv6.config --memgb=40
```
### 5.  Libraries CSV files

#### Data Access

The raw data used in this project is too large to be hosted on GitHub. You can access the raw data files from the following location:
[to be added, URL to ddbj]

Make a directory of CSV files specifying the type and location of sequencing data for executing Cellranger-ARC 2.0.2

e.g. sample_Wt_T0.csv
```csv
fastqs,sample,library_type
../home/MultiomeProject/Raw_Data/F22FTSAPJT0200_LIBrlvkR/Wt_T0_GEx,Wt_T0_GEx,Gene Expression
../home/MultiomeProject/Raw_Data/F22FTSAPJT0200_LIBrlvkR/Wt_T0_ATAC,Wt_T0_ATAC,Chromatin Accessibility
```
Sample File Names

+ sample_Wt_T00.csv
+ sample_Wt_T06.csv
+ sample_Wt_T12.csv
+ sample_Wt_T24.csv
+ sample_ste_T0.csv
+ sample_ste_T06.csv
+ sample_ste_T12.csv
+ sample_ste_T24.csv

### 6. Run Cellranger-ARC

Prepare PBS script `runCellRanger-arc.pbs`:
```txt
#!/bin/bash
#PBS -q large
#PBS -V
#PBS -l select=1:ncpus=1:mem=90gb
#PBS -J 1-9

SAMPLE=(Wt_T0 Wt_T6 Wt_T12 Wt_T24 Wt_T0_2021 ste_T0 ste_T6 ste_T12 ste_T24)
SAMPLE2=(Wt_T00 Wt_T06 Wt_T12 Wt_T24 Wt_T00_2021 ste_T00 ste_T06 ste_T12 ste_T24)
DATE=$(date +%Y-%m-%d)
echo $DATE
cd ../home/MultiomeProject/
mkdir -p runCellRanger-arc_Count_${DATE}
cd runCellRanger-arc_Count_${DATE}
../home/10xGenomics/cellranger-arc-2.0.2/cellranger-arc count \
			--id=${SAMPLE2[$PBS_ARRAY_INDEX - 1]} \
                        --libraries=../home/MultiomeProject/Sample_Files/sample_${SAMPLE[$PBS_ARRAY_INDEX - 1]}.txt \
                        --reference=../home/MultiomeProject/reference/Ppatens_v6 \
                        --localcores=3 \
                        --localmem=90 \
```
Submit PBS job.
```bash
runCellRanger-arc.pbs
```
Following above, run mainRScript.R, to generate normalised and batch corrected Seurat/Signac object.



