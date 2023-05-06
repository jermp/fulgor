Fulgor
======

Fulgor is a *colored compacted de Bruijn graph* index for large-scale matching and color queries, powered by [SSHash](https://github.com/jermp/sshash).


#### Table of contents
* [Compiling the code](#compiling-the-code)
* [Tools](#tools)
* [Demo](#Demo)
* [Indexing an example Salmonella pan-genome](#indexing-an-example-pan-genome)


Compiling the code
------------------

    git clone --recursive https://github.com/jermp/fulgor.git
    mkdir build; cd build
    cmake ..
    make

If you forgot `--recursive`, do

    git submodule update --init --recursive

before compiling.

Tools
-----


Demo
----

First, download [Cuttlefish](https://github.com/COMBINE-lab/cuttlefish).
Then, from withing `cuttlefish`, do

	git checkout inverted-colors 
	
then compile the Cufflefish code. After compilation,
from within `cuttlefish/build`, do

    ulimit -n 2048
    ./src/cuttlefish build -d [RELATIVE-FULGOR-PATH]/fulgor/test_data/salmonella_10 -k 31 -t 2 -o [RELATIVE-FULGOR-PATH]/fulgor/test_data/salmonella_10 --extract-inverted-colors

on our example data in `test_data/salmonella_10`, where `[RELATIVE-FULGOR-PATH]` is the path to `fulgor` relative to your machine.

Cuttlefish will then generate the files
`salmonella_10.fa` and `salmonella_10.cf_inv_col`.

Then, from the parent directory `fulgor`, do

	python3 scripts/build_index.py --bin-dir build -k 31 -m 17 --tmp-dir build -g 1 test_data/salmonella_10 

which will run, in order, the sub-tools `invert`, `sort_unique`, `permute_unitigs`, and `build`.

<!--Then, from within `fulgor/build`, invert the reference to unitig mapping with:

    ./invert -i ../test_data/salmonella_10 -g 1 -d tmp_dir --verbose

Deduplicate the color classes and build the map from unitig ids to color classes:

    ./sort_unique -i ../test_data/salmonella_10 -g 1 -d tmp_dir --verbose

Then permute the unitigs by color class:

    ./permute_unitigs -i ../test_data/salmonella_10 -g 1 -d tmp_dir --verbose

And finally build the index with:

    ./build_index -i ../test_data/salmonella_10 -k 31 -m 17 -d tmp_dir --verbose

Check correctness of colors:

    ./check_colors -i ../test_data/salmonella_10-->


Indexing an example Salmonella pan-genome
-----------------------------------------

The genomes can be downloaded from [here](https://zenodo.org/record/1323684). We will build a Fulgor index for the 4,546 Salmonella genomes.

Then create a list of all `.fasta` filenames with

	find $(pwd)/Salmonella_enterica/Genomes/*.fasta > salmonella_4546_filenames.txt
	
and run Cuttlefish with

    ./src/cuttlefish build -l salmonella_4546_filenames.txt -k 31 -t 8 -w tmp_dir/ -o salmonella_4546 --extract-inverted-colors
    
Then from `fulgor`, do

	python3 scripts/build_index.py --bin-dir build -k 31 -m 20 --tmp-dir build/tmp_dir -g 8 build/salmonella_4546
