# Translational Profiler - a tool for bacterial ribosome footprint analysis

This tool is designed to aid in the analysis of ribosomal footprint data (i.e., riboseq) for non-model archaea and bacteria and maximize input flexibility. As a starting point, this data handles raw, sorted `bam` files, `fasta` inputs for reference genome(s) that can be one or many contigs, and gene call positions in `gff3` (prodigal) output.

Note that this tool **requires** perfect match between the contig names in the `bam` files and the `gff` / prodigal file with gene annotations. Different tools have different sensitivities to certain characters like whitespace or colons which may result in altered contig / reference names, so best case use letters, numbers, and the underscore for contig names but if other characters must be used, ensure they are appropriate for the upstream tools.


## Installation

The simplest way is to through conda, like 

```
conda create -y -n translational_profiler python pysam=0.22.1 biopython pandas matplotlib numpy
```

Older versions of pysam should work but have not been tested.

Currently, all functions are in the `translational_profiler.py` file which can be made available to any python script as 

```
from translational_profiler.py import *
```


## Usage

Examples and the code used to make the figures from the manuscript are found in `figures.ipynb` with more full description of the usage to follow.

Much appreciation for the developers of other ribosomal footprinting tools [mRibo](https://github.com/dgelsin/mRibo) and [HRIBO](https://github.com/RickGelhausen/HRIBO) whose work influenced this project. 


### Quickstart reference
Arguments in square brackets optional. Arguments in double-square brackets optional but likely relevant (e.g. footprint-specific options)

#### Common arguments
* `bamfile="PATH"` - path to a bam file
* `annfile="PATH"` - path to a annotation file (GFF3 or prodigal output), or similar pandas dataframe
* `window=INT` - number of bases around region / gene
* `codon_resolved=BOOL` - Assign coverage to a single reference base, not all bases spanned by the mapped read
* `offset=INT` - Default 0, offset from 5' end for assigning coverage when `codon_resolved=True`
* `fiveprime=BOOL` - Default True, count from 5' end of read (False would be 3' end, automatically flips `offset`)
* `rpm=BOOL` - apply reads per million normalization, usually True/False
* `fasta_path=PATH` - path to a fasta file. Deflines should match bam contig names exactly

#### Main functions
* `metagene(bamfile, annfile, window, [[codon_resolved=False, footprint_size=24, offset=0, fiveprime=True,]] [num_missing=0, meannorm=False, maxnorm=False, rpm=False, verbose=False])`
  - Returns two 1d np arrays: first is start, second is stop.
  
* `plot_metagene(covs_start, covs_stop, bg_start=None, bg_stop=None, window=None, codonlines=True, axis_units=None, colors=['black','#dddddddd'], bg_fill=True)`
  - `bg_fill` is bool for whether background is a trace (good for foots / single-nt) or filled (default)

* `plot_genes(query, bamfile, annfile, footprint_size, txfile=None, window=50, offset=0, rpm=False,
codon_resolved=True, cols=['red','#919100','green','blue','black'], axis_units=False, 
colors=['black','#dddddddd'], bg_fill=True, norenorm=False, y1lab='Footprints', y2lab='Coverage', genelabs=False)`
  - `colors` control foreground/background colors
  - `cols` control gene colors
  - `genelabs` defaults false to not label genes, but can be a list like `['mcrA', 'mcrD']` etc from left to right.

* `codonhalves(bamfile, annfile, fasta_path, defline_trim=True, codon='M', query='', minrpm=0, rpm=True, fiveprime=True, 
    codon_resolved=False, offset=0, footprint_size=range(18,30), normalgenes=True)`
  - Returns two arrays, first being coverages before codon and second being coverages after codon
  - `codon` is the codon searched for. i.e. `M` would find coverage in windows around all `M` codons
  - `normalgenes` subsets to genes with an `M` stop and a `*` stop (UAG is not treated as a `*`)

* `codonslice(x, y, window=250)`
  - Normalize output of `codonhalves()` such that the cumulative sum of each `x` and `y` pair add up to 1
  - `window` is the range around, works best at relatively short windows
    
*  `plotPrePostAA(residue_dict, window=100, 
                  colors={'K':'green','W':'brown','P':'purple','O':'red', 'L':'gray', 'S':'orange', 'F':'pink','R':'#7a1a00','*':'#303030'},
                 n_best=10000, traces=False)`
   - Takes a dict in the format `{'M': codonslice(), 'K': ...}`and creates a plot
   - `n_best` is the number of 'best' codon windows to use if there are more than `n_best`. 'Best' windows are those that have a high read density, i.e. reads mapping aligning to different positions within the window, i.e. cumulative coverage is not flat.


#### Helpers
* `readDist(bams, [verbose=False])` - histogram of read lengths
* `read_prodigal(annfile)` - pandas dataframe from GFF or prodigal file, corrects to 0-index
* `featureCounts(bamfile, annfile, [paired=False, minlength=0, maxlength=np.inf])` - mimic featureCounts for rpf
* `checkGeneSpacing(annfile, window, [inverse=False])` - subset genes based on distance (`window`) to neighbor
* `geneFromBam(bamfile, contig, start, stop, d, window, [[offset=0, footprint_size=24, fiveprime=True, codon_resolved=False]])` - slice gene coverage from bam by coordinates
* `transformCoverages(cov, d, window, [rpm=False, meannorm=False, maxnorm=False])` - transform coverage reported by `geneFromBam`
* `read_fasta(fasta, trim=False)` - read fasta into a dict, `trim=True` cuts contig deflines at whitespace to emulate bbmap bams
* `add_stop_start_to_gff(annfile, fasta_path, trim=False)` - adds columns 'stop_cdn', etc for each gene
* `whereAreX(annfile, fasta_path, codon='O', trim=False)` - adds list of CODON position(s) where `codon` is found