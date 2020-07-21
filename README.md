
-----------------------------------------------------------
!!! WARNING !!!
-----------------------------------------------------------

This pipeline is under active development and NOT stable!

-----------------------------------------------------------
Overview
-----------------------------------------------------------

'linCNV' is a pipeline that can be called as a single
executable program. This executable can be incorporated
into the workflow management system of your preference.

The pipeline takes the output of 10x CellRanger DNA and
uses the bam file to execute a different method of genome
binning, based on bins with a fixed statistical weight.

The objective of the pipeline is to find high confidence
copy number variant (CNV) calls in individual cells, and
common CNVs in groups of cells, and to use that CNV
information to assemble a cell lineage tree.

The pipeline calls a series of scripts in Perl, Bash and R,
which sometimes call other programs on your system.

-----------------------------------------------------------
Installation
-----------------------------------------------------------

No installation is necessary beyond simply cloning or copying
this folder of scripts to your computer.

Program dependencies can be examined using the executable
(see help, below), which must be installed prior to use if
not already present.

-----------------------------------------------------------
Usage
-----------------------------------------------------------

Use './linCNV' for command line help.

The pipeline supports parallel processing via options
'--n-cpu' and '--ram-per-cpu'.

Option '--dry-run' allows a test of the command and
option configuration prior to actual execution.

Commands 'options' and 'dependencies' provide further
configuration assistance, while 'status' and 'rollback'
provide support for monitoring and manipulating pipeline
execution status.

-----------------------------------------------------------
Required Data Inputs
-----------------------------------------------------------

First, the pipeline expects that you have already run the
10x Genomics Cell Ranger pipeline to complete at least
read alignment, duplicate purging, and cell assignments.
Point the pipeline to that Cell Ranger output using option
'cell-ranger-dir', where the code expects to find files
named 'possorted_bam.bam' and 'per_cell_summary_metrics.csv'.

Additional required  files define information about the
genome in use. Option 'genome' might typically be 'hg38'.
Option 'genome-dir' tells the pipeline where to find the
genome reference files including '\<genome\>.fa' and
its index '\<genome\>.fa.fai'.

The pipeline ignores low quality regions of the genome, which
you must provide using option 'bad-regions-file', a BED format
file that lists those offending regions. See
https://pubmed.ncbi.nlm.nih.gov/31249361/ for details and
where to obtain a suitable file for human and mouse.

Finally, 'gc-file' and 'mappability-file' must be set to
two BED files that define genome bins that were assigned
fraction GC and fraction mappability scores, respectively.
These are used during bin data normalization. The size of
bins in these files is not important, e.g. 1 kb bins.

-----------------------------------------------------------
Output Naming Conventions
-----------------------------------------------------------

Data files written by the pipeline are always placed into
a folder you specify using option 'output-dir'. The file
names are always prefixed with the value of option 'data-name'.

'output-dir' and 'data-name' are thus universally required
options, and their values should be the same between 
sequential commands applied to the same data. 

-----------------------------------------------------------
Workflow
-----------------------------------------------------------

The general workflow occurs in the following conceptual steps:

1) The 'bin' command is run at the command line to create an
initial set of genome bins.

2) The user individually examines and manually marks every cell
for quality using the 'mark_cells' R Shiny web tool.

3) The 'analyze' command is run at the command line to refine
the bins based on the user's cell marks and to complete data
normalization and CNV calling.

4) The user examines the results using the 'heat_map' and other
R Shiny web tools.

-----------------------------------------------------------
R Shiny Web Tools
-----------------------------------------------------------

Instructions pending.

