## Ensembl Regulatory Annotation

This folder contains scripts for annotating **Regulatory Features (RFs)** as part of the *Regulatory Annotation Pipeline* of Ensembl Regulation (EMBL-EBI).  

The scripts support:

-   Extraction of **TSSs and CDSs** annotation used for RFs identification.

-   **Regulatory Features identification**: *promoters*, *enhancers*, and *open chromatin regions*.

-   **Features liftover** between two assemblies.

---

### Table of Contents
- [Requirements](#requirements)  
- [Scripts overview](#scripts-overview)  
- [Preparing genome annotation](#preparing-genome-annotation)  
- [Annotating Regulatory Features](#annotating-regulatory-features)  
- [Lifting over Regulatory Features between assemblies](#lifting-over-regulatory-features-between-assemblies)  
- [Citation](#citation)  
- [License](#license)  

---

### Requirements
- [bedtools](https://bedtools.readthedocs.io/en/latest/) (>= 2.30) 
- [WiggleTools](https://github.com/Ensembl/WiggleTools) (>= 1.20) 
- [UCSC liftOver](https://hgdownload.soe.ucsc.edu/admin/exe/)  
- [R](https://cran.r-project.org/) (>= 4.0) with packages: 
  - [argparse](https://cloud.r-project.org/web/packages/argparse/index.html) 
  - [dplyr](https://cloud.r-project.org/web/packages/dplyr/index.html)
  - [stringr](https://cloud.r-project.org/web/packages/stringr/index.html)
  - [tidyr](https://cloud.r-project.org/web/packages/tidyr/index.html) 
  - [bedtoolsr](http://phanstiel-lab.med.unc.edu/bedtoolsr.html)

---

### Scripts overview

- [`AnnotationExtractor.R`](./scripts/AnnotationExtractor.R): Auxiliary script to extract **TSSs and CDSs** from a genome annotation (Ensembl/GENCODE).  

- [`RegulatoryBuild.R`](./scripts/RegulatoryBuild.R): Main script for **regulatory feature identification** across a genome assembly.  

- [`FeaturesLiftOver.R`](./scripts/FeaturesLiftOver.R): Script for **mapping features between assemblies** using UCSC liftOver.  

---

### Preparing genome annotation

#### Description

The [AnnotationExtractor.R](./scripts/AnnotationExtractor.R) script extracts *TSSs* and *CDSs* annotations that are useful for *promoters* and *enhancers* identification.

Input:

- `annotation.gff3`: Genome annotation produced by *Ensembl* or *GENCODE*.
- `chrsList.txt`: File listing the chromosomes of interest (i.e. canonical chromosomes) where features will be annotated. The order of the chromosomes in this file will determine the ordering of the **RFs** file as well. 

Output:

- `annotation_TSS.bed`: A bed-like file containing *TSS* coordinates, *chr*, *start*, *end*, in the first three columns (as any bed file), followed by the *strand*, *Ensembl gene id*, *gene biotype*, and *gene name*.

- `annotation_CDS.bed`: A bed file containing *CDS* coordinates, *chr*, *start*, *end*.

#### How to use

By default, all *TSSs* and *CDSs* from all transcripts annotated in the input file will be exported. 

The parameters `--transcriptsBiotypes`, `--flagsToKeep`, and `--flagsToRemove` allow to filter transcripts. For instance, to keep only *protein_coding* transcripts we can specify `--transcriptsBiotypes protein_coding`.

For more details, run the script with the `-h` option:

<details>

<summary> Code:</summary>

```         
Rscript AnnotationExtractor.R -h
usage: Rscript AnnotationExtractor.R [-h]
                                     [--transcriptsBiotypes TRANSCRIPTSBIOTYPES]
                                     [--flagsToKeep FLAGSTOKEEP]
                                     [--flagsToRemove FLAGSTOREMOVE]
                                     ChrFile annotationGFF3 isGENCODE

Extract TSSs and CDSs annotation from Ensembl GFF3 file

positional arguments:
  ChrFile               Full path to the file listing the canonical
                        chromosomes / chromosomes of interest (one per line).
  annotationGFF3        Full path to the Ensembl gff3 file containing gene
                        annotation.
  isGENCODE             Logical indicating if the genome annotation file has
                        been produced by GENCODE.

options:
  -h, --help            show this help message and exit
  --transcriptsBiotypes TRANSCRIPTSBIOTYPES
                        Comma-separated list of transcript biotypes to keep.
                        Default value is `NULL` indicating that all biotypes
                        should be kept.
  --flagsToKeep FLAGSTOKEEP
                        Comma-separated list of transcript flags to be keep.
                        Default value is `NULL` indicating that all
                        transcripts with any flag are kept.
  --flagsToRemove FLAGSTOREMOVE
                        Comma-separated list of flags indicating transcripts
                        having them should be removed. Default value is NULL,
                        all transcripts are kept.
```

</details>

#### Example: Extracting *GRCh38* *TSSs* and *CDSs* annotation

Ensembl 114 **RFs** are based on the *GENCODE Basic* gene annotation [release 48](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.basic.annotation.gff3.gz). 

We created a file, `GRCh38.chromList.txt` containing the canonical chromosomes (excluding *MT*).

Each chromosome is listed on a separate line, ordered numerically from *1* to *22*, followed by *X* and *Y*.  

For human **RFs**, we:

- Considered the following transcript biotypes: *protein_coding*, *lncRNA*, *retained_intron*, and *protein_coding_LoF*.

- Discarded transcripts flagged as *mRNA_start_NF*, since their *TSSs* are not defined.

<details>

<summary>Code:</summary>

```         
Rscript AnnotationExtractor.R \
'/path/to/chr/list/file/GRCh38.chromList.txt' '/path/to/GENCODE/v48/gencode.v48.basic.annotation.gff3' 'TRUE' \
--transcriptsBiotypes 'protein_coding,lncRNA,retained_intron,protein_coding_LoF' --flagsRemove 'mRNA_start_NF'
```

</details>

<details>

<summary>Output:</summary>

As result, we obtain the *TSSs* and *CDSs* annotation files exported to the working directory, looking like:

```
$ head -n 3 gencode.v48.annotation_TSS.bed
1	11425	11426	+	ENSG00000290825	lncRNA	DDX11L16
1	29366	29367	-	ENSG00000310526	lncRNA	WASH7P
1	29345	29346	-	ENSG00000310526	lncRNA	WASH7P

$ head -n 3 gencode.v48.annotation_CDS.bed
1	65564	65573
1	69036	70008
1	450739	451678
```
</details>

---

### Annotating Regulatory Features

#### Description

The [RegulatoryBuild.R](./scripts/RegulatoryBuild.R) identifies a set of **Regulatory Features** across
the whole genome of a species (`S`), based on a specific assembly (`A`).

Required input:

- A set of open chromatin peaks (obtained from `ATAC-seq` and/or `DNase-seq` experiments)

- *TSSs* and *CDSs* annotation files generated by the [AnnotationExtractor.R](./scripts/AnnotationExtractor.R) script

Whit these, the script can generates *promoter* and *open chromatin region* annotation.
 
Optional input: 

- `H3K4me1` and/or `H3K27ac` ChIP-seq peaks for adding *enhancers* identification.
- Histone ChIP-seq (`H3K4me1`, `H3K27ac`, or `H3K4me3`) peaks for annotating *EMARs* which are regions identified as accessible and containing epigenetic modifications.

The full strategy for identifying regulatory elements is described
[here](https://regulation.ensembl.org/help/regulatory_build).

#### How to use

Required and optional parameters can be checked by using the `-h` option:

<details>

<summary>Code:</summary>

```      
Rscript RegulatoryBuild.R -h
usage: Rscript RegulatoryBuild.R [-h] [--path2K4me1Files PATH2K4ME1FILES] [--path2K27acFiles PATH2K27ACFILES] 
[--path2K4me3Files PATH2K4ME3FILES] [--oldAnnotation OLDANNOTATION] [--ocPattern OCPATTERN] 
[--k4me1Pattern K4ME1PATTERN] [--k27acPattern K27ACPATTERN] [--k4me3Pattern K4ME3PATTERN] 
[--TSSwindowUp TSSWINDOWUP] [--TSSwindowDown TSSWINDOWDOWN] [--coreWUp COREWUP] [--coreWDown COREWDOWN]
[--featureOverlap FEATUREOVERLAP] [--histoneOverlap HISTONEOVERLAP] [--OCRsHistoneOverlap OCRSHISTONEOVERLAP]
[--cOCRsCDSOverlap COCRSCDSOVERLAP] [--maxDistOCR MAXDISTOCR] [--maxDistH MAXDISTH]
[--maxPromExtLen MAXPROMEXTLEN] [--minOCRLen MINOCRLEN] [--maxEnhancerLen MAXENHANCERLEN]
[--maxOCRLen MAXOCRLEN] [--exportActivity EXPORTACTIVITY] [--path2EpigenomesOCFiles EPIOCFILES]
[--epiPattern EPIPATTERN] S A EPrefix startID TSSsFile CDSFile output OCPFiles 

Regulatory features annotation

positional arguments: 

- S Species name in lower case and using underscore as a
separator (e.g Xxx_yyy).

- A Assembly name. EPrefix Ensembl prefix for the species. 

- startID Integer specifying the starting value for feature IDs. 

- TSSsFile Full path to the bed file containing TSSs annotation (chr, start, 
end, strand).

- CDSFile Full path to the bed file containing CDSs annotation (chr, start, end).

- output Full path to the output folder. 

- OCPFiles Full path to the folder where atac-seq/dnase-seq peak files are.

options: 
-h, --help show this help message and exit

--path2K4me1Files PATH2K4ME1FILES Full path to the folder where H3K4me1 ChiP-seq peaks files are.

--path2K27acFiles PATH2K27ACFILES Full path to the folder where H3K27ac ChiP-seq peaks files are. 

--path2K4me3Files PATH2K4ME3FILES Full path to the folder where H3K4me3 ChIP-seq peaks files are. 

--oldAnnotation OLDANNOTATION Full path to the gff3 file with previous regulatory annotation. 

--ocPattern OCPATTERN Pattern characterizing the subset of atac-seq/dnase- seq bed files to be used.

--k4me1Pattern K4ME1PATTERN Pattern characterizing the subset of H3K4me1 bed files to be used. 

--k27acPattern K27ACPATTERN Pattern characterizing the subset of H3K27ac bed files to be used. 

--k4me3Pattern K4ME3PATTERN Pattern characterizing the subset of H3K4me3 bed files to be used. 

--TSSwindowUp TSSWINDOWUP Window size (bp) upstream TSSs for promoters identification (default 90). 

--TSSwindowDown TSSWINDOWDOWN Window size (bp) downstream TSSs for promoters identification (default 10). 

--coreWUp COREWUP Upstream core size (bp) for potential promoters (default 490). 

--coreWDown COREWDOWN Downstream core size (bp) for potential promoters (default 10). 

--featureOverlap FEATUREOVERLAP Fraction of overlapping between RFs and open chromatin peaks for activity annotation (default 0.2). 

--histoneOverlap HISTONEOVERLAP Fraction of overlap between histone and open chromatin peaks, for features identification (default 0.5). 

--OCRsHistoneOverlap OCRSHISTONEOVERLAP Fraction of overlap between cOCRs and histone peaks for enhancers identification (default 0.5). 

--cOCRsCDSOverlap COCRSCDSOVERLAP Fraction of overlap between cOCRs and CDSs for enhancers identification (default 0.1). 

--maxDistOCR MAXDISTOCR Maximum distance between peaks to be merged when defining cOCRs (default 100). 

--maxDistH MAXDISTH Maximum distance between histone peaks to be merged (default 100).

--maxPromExtLen MAXPROMEXTLEN Maximum size for promoters extended regions (default 1000). 

--minOCRLen MINOCRLEN Minimum size for cOCRs (default 100). 

--maxEnhancerLen MAXENHANCERLEN Maximum length for enhancers (default 2500).

--maxOCRLen MAXOCRLEN Maximum length for open chromatin regions (default 2500).

--exportActivity EXPORTACTIVITY Should activity files be generated?(default FALSE) 

--path2EpigenomesOCFiles EPIOCFILES Full path to the folder containing the open chromatin peak files to be used for activity annotation. 

--epiPattern EPIPATTERN Pattern characterizing the subset of bed files to be used for activity annotation.
```

</details>

#### Example: Generating *GRCh38* RFs

For human, we have analyzed 61 `ATAC-seq`/`DNase-seq` experiments from different tissues and cell lines that passed our
quality control and were then used for RFs identification.

Because of the high heterogeneity observed in the open chromatin peaks set, we decided to retain only
those regions that were detected as open in at least 10% of the available experiments (six in this case).

For obtaining these regions, we

-   Combine all open chromatin peaks across experiments, and sort them by position using [bedtools
    merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html)

-   Compute coverage per base across experiment using [WiggleTools](https://github.com/Ensembl/WiggleTools).

-   Filter regions, keeping only those with coverage $\geq$ 10%.

<details>

<summary>Code:</summary>

```         
# /path/to/bed/files/ contain all the open chromatin and ChIP-seq peak files

cd /path/to/bed/files/

find . -type f \( -name '*atac*' -o -name '*dnase*' \) -exec cat {} + | bedtools sort -i - > all.bed

wiggletools write_bg allCov.bed coverage all.bed

awk 'BEGIN{FS="\t";OFS="\t"}{if($4 > 5){print $0}}' allCov.bed | bedtools merge -i - > OCPeaksGT10Perc_6.bed

rm all.bed allCov.bed
```

</details>
<br>
The `OCPeaksGT10Perc_6.bed` contains all the open chromatin regions that were used for the **RFs** definition.

Assuming we have:

-   The `OCPeaksGT10Perc_6.bed` file and all `H3K4me1`, `H3K27ac`, and `H3K4me3` peak files in the `/path/to/bed/files/`
    folder

-   The *TSSs* and *CDSs* annotation files in `/path/to/gencodeV48/`

-   The 113 version of the **GRCh38 Regulatory annotation** file in `/path/to/113/GRCh38/annotation/`

-   The output folder, `/path/to/output/folder` where the annotated features and the log file will be saved

The command line we used to run the [RegulatoryBuild.R](./scripts/RegulatoryBuild.R) script is:

<details>

<summary>Code:</summary>

```         
Rscript /path/to/the/RBscript/RegulatoryBuild.R "homo_sapiens" "GRCh38" "ENS" 1 \
"/path/to/gencodeV48/gencode.v48.annotation_TSS.bed" \
"/path/to/gencodeV48/gencode.v48.annotation_CDS.bed" \
"/path/to/output/folder/" "/path/to/bed/files/" \
--oldAnnotation "/path/to/113/GRCh38/annotation/homo_sapiens.GRCh38.regulatory_features.113.gff3" \
--ocPattern "OCPeaksGT10Perc_6" --path2H3K4me1Files /path/to/bed/files/" \
--path2H3K27acFiles "/path/to/bed/files/" \
--path2H3K4me3Files "/path/to/bed/files/" \
--h3k4me1pattern "k4me1" --h3k27acpattern "k27ac" --h3k4me3pattern "k4me3" \
--exportActivity FALSE > /path/to/output/folder/homo_sapiens.GRCh38.v114.log
```

</details>

---

### Lifting over Regulatory Features between assemblies

#### Description

The [FeaturesLiftOver.R](./scripts/FeaturesLiftOver.R) script implements the procedure we designed to translate features from one assembly to
another.

Starting from a set of **RFs** and Transcription Factor (TF) binding **motif features** annotated in assembly *A*, it uses [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to project these features onto assembly *B* using the appropriate chain file that describes the conversion between these assemblies.

<details>
<summary>The processing steps are:</summary>
<br> 

   a. Load current **motifs features** and **RFs** identified in assembly *A* and expand the list to include extended *promoters* when present.
   
   b. Map chromosome IDs from Ensembl style to the one used in the liftOver chain file. 
   
   c. Translate the **RFs** by running liftOver to assembly *B*.
   
   d. Load lifted **RFs** and map chromosome IDs back to *Ensembl* style 
   
   e. Use the extended *promoter* feature coordinates to fill the attribute columns of *promoters* when possible.
   
   f. Remove **RFs** mapping to non-canonical chromosomes and **RFs** mapped to different chromosomes (if specified).
   
   g. Update promoter-genes links matching `gene_name` or `gene_id` and keep only the links matching gene attributes between *A* and *B*.
   
   h. Solving overlapping features applying the following rules:
   
      - Keep *promoters* linked to *protein coding* over *lncRNA* genes and then remove *EMARs* overlapping the discarded *promoters*.
      
      - Keep *EMARs* with the largest number of *promoters* linked to *protein_coding* genes and remove *promoters* overlapping discarded *EMARs*.
      
      - Keep *EMARs* with the largest region overlapping *promoters*, or the largest number of *promoters* and remove *promoters* overlapping discarded *EMARs*.
      
      - For *EMARs* that do not overlap *promoters*, keep the first one and remove *enhancers* overlapping discarded *EMARs*.
      
      - Keep *open chromatin regions* that do not overlap *promoters* or *enhancers*, selecting the first in location when there is any overlapping.
      
      - Remove **motifs** and *CTCFs* overlapping any discarded *promoter*, *enhancer*, *open chromatin region*, or *EMAR*. 
      
      - Keep only *CTCFs* that overlap kept *CTCF* **motifs**. If there is more than one, then pick the first one.

</details>

#### How to use

In addition to the **RFs** and **motif *features** files from assembly *A*, and the corresponding chain file for `liftOver`, the script also requires the assembly report files for both assemblies *A* and *B* in order to correctly map chromosome names between the different files. 

The *TSSs* annotation for assembly *B* is also needed. It can be created by running the [AnnotationExtractor.R](./scripts/AnnotationExtractor.R) script.

A detailed list of required and optional parameters can be checked by using the `-h` option:

<details>

<summary>Code:</summary>

```         
Rscript FeaturesLiftOver.R -h
usage: Rscript FeaturesLiftOver.R [-h] [--keepRFsInDiffChr keepRFsInDiffChr] 
                                  [--updateGenePromLink updateGenePromLink]
                                  [--keepNonTSSsProms keepNonTSSsProms]
                                  [--promColor promColor]
                                  fromAnnotation fromMotifs toAnnotation
                                  toMotifs chainFile outPath fromAssemblyInfo
                                  toAssemblyInfo TSSsFile toAssemblyChr

Lifting over regulatory features between two assemblies

positional arguments:
  fromAnnotation        Full path to the gff3 file containing the regulatory
                        features to be lifted over.
  fromMotifs            Full path to the gff3 file containing the TF motifs
                        features to be lifted over.
  toAnnotation          Name of the gff3 file where lifted regulatory features
                        will be exported.
  toMotifs              Name of the gff3 file where lifted TF motifs features
                        will be exported.
  chainFile             Full path to the chain file.
  outPath               Full path to the output folder.
  fromAssemblyInfo      Full path to the report file of the original assembly
                        used for RFs identification, useful for chromosomes ID
                        mapping.
  toAssemblyInfo        Full path to the report file of the new assembly,
                        useful for chromosomes ID mapping.
  TSSsFile              Full path to the bed file containing TSSs (chr, start,
                        end, strand) and gene annotations of the new assembly.
  toAssemblyChr         Full path to the chromsomes file of the new assembly.

options:
  -h, --help            show this help message and exit
  --keepRFsInDiffChr keepRFsInDiffChr
                        Should features mapping to different chromosomes be
                        kept? Default is to remove them.
  --updateGenePromLink updateGenePromLink
                        Should the gene-promoter links be updated? Default is
                        yes, keep the links that are still valid in the new
                        assembly. If it is FALSE, then new links will be
                        created.
  --keepNonTSSsProms keepNonTSSsProms
                        Should promoters do not overlapping TSSs (without
                        genes linked to) be kept? Default is FALSE.
  --promColor promColor
                        Color of promoters.
```

</details>


#### Example: Mapping features from *GRCh38* to *GRCh37*

Required files:

- *GRCh38* [RFs](https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.regulatory_features.v114.gff3.gz), [EMARs](https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.EMARs.v114.gff.gz), and [motif features](https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.motif_features.v114.gff3.gz) file.

- Chain file [hg38ToHg19.over.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz)

- Assembly report files:
  
  - [GRCh38 assembly report](GCF_000001405.40_GRCh38.p14_assembly_report.txt)
  
  - [GRCh37 assembly report](GCF_000001405.25_GRCh37.p13_assembly_report.txt)
  
- *GENCODE Basic* gene annotation for [GRCh37 release 48](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh37_mapping/gencode.v48lift37.basic.annotation.gff3.gz). 

<details>

<summary>Pre-processing:</summary>

- Combine **RFs** and *EMARs* in a single file, `Homo_sapiens.GRCh38.RFs.v114.gff3`.

- Extract *TSSs* annotation in *GRCh37*:

    ```
    Rscript AnnotationExtractor.R '/path/to/chr/list/file/GRCh37.chromList.txt' \
    'gencode.v48lift37.basic.annotation.gff3' 'TRUE' \
    --transcriptsBiotypes 'protein_coding,lncRNA,lincRNA,retained_intron,protein_coding_LoF' \
    --flagsRemove 'mRNA_start_NF'
    ```
</details>


<details>

<summary>Code:</summary>

```
Rscript FeaturesLiftOver.R /path/to/GRCh38/RFs/file/Homo_sapiens.GRCh38.RFs.v114.gff3 \
/path/to/GRCh38/motifs/file/Homo_sapiens.GRCh38.motif_features.v114.gff3 \
Homo_sapiens.GRCh37.regulatory_features.v114.gff3 Homo_sapiens.GRCh37.motif_features.v114.gff3 \
/path/to/chain/file/hg38ToHg19.over.chain.gz /path/to/GRCh37/output/ \
/path/to/GRCh38/assembly/report/GCF_000001405.40_GRCh38.p14_assembly_report.txt \
/path/to/GRCh37/assembly/report/GCF_000001405.25_GRCh37.p13_assembly_report.txt \
/path/to/GRCh37/TSS/annotation/gencode.v48lift37.basic.annotation_TSS.bed \
'/path/to/chr/list/file/GRCh37.chromList.txt'
```

</details>

