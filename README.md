ens-sandbox
===========

A sandbox of development for the Ensembl project. See each section for more information on what's available.

## ensembl_chains

This is a single script which connects to an Ensembl database and writes out the mapping held in that DB into a UCSC formatted chain file. This allows you to use Ensembl generated mappings with tools such as liftover (http://genome-euro.ucsc.edu/cgi-bin/hgLiftOver) or crossmap (http://crossmap.sourceforge.net/). More information on the chain format can be found http://genome.ucsc.edu/goldenPath/help/chain.html.

Ensembl holds only 1 mapping and internally can provide bi-directional mapping (internally held as asm & cmp); these tools require a chain file per mapping direction. The script writes out two files:

- asmTocmp.chain
- cmpToasm.chain

What these are requires you to know what is available in the Ensembl database in question. By convention Ensembl holds the current live assembly as `asm` and the older assembly as `cmp`. So you could invoke the script like so:

```
./ensembl_assembly_to_chain.pl -host ensembldb.ensembl.org -port 3306 -user anonymous -version 75 -asm GRCh37 -cmp NCBI36
```

This will produce 2 files in your current working directory:

- NCBI36ToGRCh37.chain
- GRCh37ToNCBI36.chain

To use this with liftover first download the tool locally, format a bed file (with Ensembl chromsome names not UCSC names) and run it like so (this will work with this repo):

```
cd ensembl_chains
liftOver example.ensembl.36.bed example_chains/NCBI36ToGRCh37.chain test.new test.unmapped
```

You can compare your results to the `example.ensembl.mapped37.bed` file which is a pre-computed mapping identical to one generated from UCSC's live lift-over tool.