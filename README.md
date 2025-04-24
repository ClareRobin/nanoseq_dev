

## For Demonstration Purposes only

This forked version of nanoseq is for demonstration purposes only. It contains a pared-down version of the main workflow (nanoseq.nf) called nanoseq_interview.nf, that only does: 
   - alignment of long read RNA-Seq data using minimap2 
   - transcript discovery & gene expression estimation using bambu
   - differential gene expression using DEseq2. 

The pipeline was run on downsampled data from the SG-NEx (Singapore Nanopore Expression) project found at http://sg-nex-data.s3-website-ap-southeast-1.amazonaws.com/#data/data_tutorial/fastq/, using downsampled reference Hg38 genome. 

The pipeline was run with a samplesheet formatted as follows: 
```
group,replicate,barcode,input_file,fasta,gtf
A549_directRNA,1,,data/A549_directRNA_sample1.fastq.gz,data/hg38_chr22.fa,data/hg38_chr22.gtf
A549_directRNA,2,,data/A549_directRNA_sample2.fastq.gz,data/hg38_chr22.fa,data/hg38_chr22.gtf
A549_directRNA,3,,data/A549_directRNA_sample3.fastq.gz,data/hg38_chr22.fa,data/hg38_chr22.gtf
HepG2_directRNA,1,,data/HepG2_directRNA_sample1.fastq.gz,data/hg38_chr22.fa,data/hg38_chr22.gtf
HepG2_directRNA,2,,data/HepG2_directRNA_sample2.fastq.gz,data/hg38_chr22.fa,data/hg38_chr22.gtf
HepG2_directRNA,3,,data/HepG2_directRNA_sample3.fastq.gz,data/hg38_chr22.fa,data/hg38_chr22.gtf
```

And using the following command (using docker): 
```
nextflow run main.nf --input data/samplesheet.csv -profile docker --protocol directRNA
```

Changes and deviations from nanoseq are documented in the changelog, bambu gene and transcript count and deseq2 results are uploaded for reference. 


## Begin Nanoseq READMe

## Introduction

**nfcore/nanoseq** is a bioinformatics analysis pipeline for Nanopore DNA/RNA sequencing data that can be used to perform basecalling, demultiplexing, QC, alignment, and downstream analysis.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a [full-sized dataset](https://github.com/nf-core/test-datasets/tree/nanoseq#full-sized-test-data) obtained from the [Singapore Nanopore Expression Consortium](https://github.com/GoekeLab/sg-nex-data) on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/nanoseq/results).

## Pipeline Summary

1. Demultiplexing ([`qcat`](https://github.com/nanoporetech/qcat); _optional_)
2. Raw read cleaning ([NanoLyse](https://github.com/wdecoster/nanolyse); _optional_)
3. Raw read QC ([`NanoPlot`](https://github.com/wdecoster/NanoPlot), [`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Alignment ([`GraphMap2`](https://github.com/lbcb-sci/graphmap2) or [`minimap2`](https://github.com/lh3/minimap2))
   - Both aligners are capable of performing unspliced and spliced alignment. Sensible defaults will be applied automatically based on a combination of the input data and user-specified parameters
   - Each sample can be mapped to its own reference genome if multiplexed in this way
   - Convert SAM to co-ordinate sorted BAM and obtain mapping metrics ([`samtools`](http://www.htslib.org/doc/samtools.html))
5. Create bigWig ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/)) and bigBed ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedToBigBed`](http://hgdownload.soe.ucsc.edu/admin/exe/)) coverage tracks for visualisation
6. DNA specific downstream analysis:
   - Short variant calling ([`medaka`](https://github.com/nanoporetech/medaka), [`deepvariant`](https://github.com/google/deepvariant), or [`pepper_margin_deepvariant`](https://github.com/kishwarshafin/pepper))
   - Structural variant calling ([`sniffles`](https://github.com/fritzsedlazeck/Sniffles) or [`cutesv`](https://github.com/tjiangHIT/cuteSV))
7. RNA specific downstream analysis:
   - Transcript reconstruction and quantification ([`bambu`](https://bioconductor.org/packages/release/bioc/html/bambu.html) or [`StringTie2`](https://ccb.jhu.edu/software/stringtie/))
     - bambu performs both transcript reconstruction and quantification
     - When StringTie2 is chosen, each sample can be processed individually and combined. After which, [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/) will be used for both gene and transcript quantification.
   - Differential expression analysis ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and/or [`DEXSeq`](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html))
   - RNA modification detection ([`xpore`](https://github.com/GoekeLab/xpore) and/or [`m6anet`](https://github.com/GoekeLab/m6anet))
   - RNA fusion detection ([`JAFFAL`](https://github.com/Oshlack/JAFFA))
8. Present QC for raw read and alignment results ([`MultiQC`](https://multiqc.info/docs/))

### Functionality Overview

A graphical overview of suggested routes through the pipeline depending on the desired output can be seen below.

<p align="center">
    <img src="docs/images/nanoseq_subwaymap_v3.1.png" alt="nf-core/nanoseq metro map" width="90%"
</p>

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/nanoseq -profile test,YOURPROFILE
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity` and are persistently observing issues downloading Singularity images directly due to timeout or network issues, then you can use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, you can use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

## Documentation

The nf-core/nanoseq pipeline comes with documentation about the pipeline [usage](https://nf-co.re/nanoseq/usage), [parameters](https://nf-co.re/nanoseq/parameters) and [output](https://nf-co.re/nanoseq/output).

```bash
nextflow run nf-core/nanoseq \
    --input samplesheet.csv \
    --protocol DNA \
    --barcode_kit SQK-PBK004 \
    -profile <docker/singularity/podman/institute>
```

See [usage docs](https://nf-co.re/nanoseq/usage) for all of the available options when running the pipeline.

An example input samplesheet for performing both basecalling and demultiplexing can be found [here](assets/samplesheet.csv).

## Credits

nf-core/nanoseq was originally written by [Chelsea Sawyer](https://github.com/csawye01) and [Harshil Patel](https://github.com/drpatelh) from [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) for use at [The Francis Crick Institute](https://www.crick.ac.uk/), London. Other primary contributors include [Laura Wratten](https://github.com/lwratten), [Ying Chen](https://github.com/cying111), [Yuk Kei Wan](https://github.com/yuukiiwa) and [Jonathan Goeke](https://github.com/jonathangoeke) from the [Genome Institute of Singapore](https://www.a-star.edu.sg/gis), [Christopher Hakkaart](https://github.com/christopher-hakkaart) from [Institute of Medical Genetics and Applied Genomics](https://www.medizin.uni-tuebingen.de/de/das-klinikum/einrichtungen/institute/medizinische-genetik-und-angewandte-genomik), Germany, and [Johannes Alneberg](https://github.com/alneberg) and [Franziska Bonath](https://github.com/FranBonath) from [SciLifeLab](https://www.scilifelab.se/), Sweden.

Many thanks to others who have helped out along the way too, including (but not limited to): [@crickbabs](https://github.com/crickbabs), [@AnnaSyme](https://github.com/AnnaSyme), [@ekushele](https://github.com/ekushele).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/channels/nanoseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
