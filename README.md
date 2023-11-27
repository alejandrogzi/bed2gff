![Crates.io](https://img.shields.io/crates/v/bed2gff?color=green)
![GitHub](https://img.shields.io/github/license/alejandrogzi/bed2gff?color=blue)

# **bed2gff**

A Rust BED-to-GFF3 parallel translator.


translates
```
chr7 56766360 56805692 ENST00000581852.25 1000 + 56766360 56805692 0,0,200 3 3,135,81, 0,496,39251,
```
into
```
chr7 bed2gff gene 56399404 56805892 . + . ID=ENSG00000166960;gene_id=ENSG00000166960

chr7 bed2gff transcript 56766361 56805692 . + . ID=ENST00000581852.25;Parent=ENSG00000166960;gene_id=ENSG00000166960;transcript_id=ENST00000581852.25

chr7 bed2gff exon 56766361 56766363 . + . ID=exon:ENST00000581852.25.1;Parent=ENST00000581852.25;gene_id=ENSG00000166960;transcript_id=ENST00000581852.25,exon_number=1

chr7 bed2gff CDS 56766361 56766363 . + 0 ID=CDS:ENST00000581852.25.1;Parent=ENST00000581852.25;gene_id=ENSG00000166960;transcript_id=ENST00000581852.25,exon_number=1

...

chr7 bed2gff start_codon 56766361 56766363 . + 0 ID=start_codon:ENST00000581852.25.1;Parent=ENST00000581852.25;gene_id=ENSG00000166960;transcript_id=ENST00000581852.25,exon_number=1

chr7 bed2gff stop_codon 56805690 56805692 . + 0 ID=stop_codon:ENST00000581852.25.3;Parent=ENST00000581852.25;gene_id=ENSG00000166960;transcript_id=ENST00000581852.25,exon_number=3

...
```

in a few seconds.

Converts
- *Homo sapiens* GRCh38 GENCODE 44 (252,835 transcripts) in 4.16 seconds.
- *Mus musculus* GRCm39 GENCODE 44 (149,547 transcritps) in 2.15 seconds.
- *Canis lupus* familiaris ROS_Cfam_1.0 Ensembl 110 (55,335 transcripts) in 1.30 seconds.
- *Gallus gallus* bGalGal1 Ensembl 110 (72,689 transcripts) in 1.51 seconds.

>**What's new on v.0.1.2**
>
> - Now bed2gtf works over a parallel algorithm that reduces computation time x3 (compared to the previous implementation).
> - Fixes a recently noted bug on gene line coordinates (wrong ends).
> - Due to the breaking changes,  `*.bed.gz` and `*.bed.zlib` support have been disable (will come back in future releases).
> - The library feature is temporarily disable and now bed2gtf only works a CLI tool
> - Disables the lexicograph-based algorithm implemented in the previous version and tries to outputs a somewhat sorted .gtf file (chromosome + start). Note that features will not be in order, if user needs that it is recommended to use [gtfsort](https://github.com/alejandrogzi/gtfsort)


## Usage
``` rust
Usage: bed2gff[EXE] --bed <BED> --isoforms <ISOFORMS> --output <OUTPUT>

Arguments:
    --bed <BED>: a .bed file
    --isoforms <ISOFORMS>: a tab-delimited file
    --output <OUTPUT>: path to output file

Options:
    --help: print help
    --version: print version
    --threads/-t: number of threads (default: max cpus)
```

>[!WARNING] 
>
>All the transcripts in .bed file should appear in the isoforms file.
#### crate: [https://crates.io/crates/bed2gff](https://crates.io/crates/bed2gff)

<details>
<summary>click for detailed formats</summary>
<p>
bed2gff just needs two files:

1. a .bed file

    tab-delimited files with 3 required and 9 optional fields:

    ```
    chrom   chromStart  chromEnd      name    ...
      |         |           |           |
    chr20   50222035    50222038    ENST00000595977    ...
    ```

    see [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) for more information

2. a tab-delimited .txt/.tsv/.csv/... file with genes/isoforms (all the transcripts in .bed file should appear in the isoforms file):

    ```
    > cat isoforms.txt

    ENSG00000198888 ENST00000361390
    ENSG00000198763 ENST00000361453
    ENSG00000198804 ENST00000361624
    ENSG00000188868 ENST00000595977
    ```

    you can build a custom file for your preferred species using [Ensembl BioMart](https://www.ensembl.org/biomart/martview). 

</p>
</details>

## Installation
to install bed2gff on your system follow this steps:
1. get rust: `curl https://sh.rustup.rs -sSf | sh` on unix, or go [here](https://www.rust-lang.org/tools/install) for other options
2. run `cargo install bed2gff` (make sure `~/.cargo/bin` is in your `$PATH` before running it)
4. use `bed2gff` with the required arguments
5. enjoy!

## Build
to build bed2gff from this repo, do:

1. get rust (as described above)
2. run `git clone https://github.com/alejandrogzi/bed2gff.git && cd bed2gff`
3. run `cargo run --release -- -b <BED> -i <ISOFORMS> -o <OUTPUT>`


## Output

bed2gff will send the output directly to the same .bed file path if you specify so

```
bed2gff annotation.bed isoforms.txt output.gff

.
├── ...
├── isoforms.txt
├── annotation.bed
└── output.gff3
```
where `output.gff3` is the result.

## FAQ
### Why?

Converting formats is a daily practice in bioinformatics. This is way more common while working with gene annotations as tools differ in input/output layouts. GTF/GFF/BED are the most used structures to store gene-related annotations and the conversion needs are not well covered by available software. 

A considerable portion of genomic tools reduce the software space by accepting GTF/GFF3 files only, directing BED users to translate their files into different formats. While some of this issues have already been covered (e.g. [bed2gtf](https://github.com/alejandrogzi/bed2gtf)) with GTF files, the GFF3 layout lacks stable converting tools (1, 2).

bed2gff is presented as a straightforward option to convert BED files into ready-to-use GFF3 files, closing that gap.  


### How?
bed2gff, takes the base code of [bed2gtf](https://github.com/alejandrogzi/bed2gtf), that basically is the reimplementation of UCSC's C binaries merged in 1 step (bedToGenePred + genePredToGtf). This tool evaluates the position of exons and other features (CDS, stop/start, UTRs), preserving reading frames and adjusting the indexing count. The main approach now is a parallel algorithm that significantly reduces computation times. 

Following the rationale of [bed2gtf](https://github.com/alejandrogzi/bed2gtf), bed2gff is able to produce a ready-to-use gff3 file by using an isoforms file, that works as the refTable in C binaries to map each transcript to their respective gene. 


## References

1. https://bioinformatics.stackexchange.com/questions/2242/how-to-convert-bed-to-gff3
2. https://www.biostars.org/p/2/
