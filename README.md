Fulgor
======

Fulgor is a *colored compacted de Bruijn graph* index for large-scale matching and color queries, powered by [SSHash](https://github.com/jermp/sshash).

#### Table of contents
* [Compiling the code](#compiling-the-code)
* [Tools](#tools)
* [Demo](#Demo)
* [Indexing an example Salmonella pan-genome](#indexing-an-example-salmonella-pan-genome)


Compiling the code
------------------

    git clone --recursive https://github.com/jermp/fulgor.git
    mkdir build; cd build
    cmake ..
    make

If you forgot `--recursive`, do

    git submodule update --init --recursive

before compiling.

To compile the code for a release environment (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following:

    mkdir build
    cd build
    cmake ..
    make -j

For a testing environment, use the following instead:

    mkdir debug_build
    cd debug_build
    cmake .. -D CMAKE_BUILD_TYPE=Debug -D FULGOR_USE_SANITIZERS=On
    make -j
    
Tools
-----

There is one executable called `fulgor` after the compilation, which can be used to run a tool.
Run `./fulgor` to see a list of available tools.

	== Fulgor: a colored compacted de Bruijn graph index ==========================
	
	Usage: ./fulgor <tool> ...
	
	Available tools:
	  invert          	 invert the reference-to-unitig Cuttlefish's file 
	  sort-unique     	 deduplicate the color sets 
	  permute-unitigs 	 permute unitigs according to the color they map to 
	  build           	 build a fulgor index 
	  pseudoalign     	 pseudoalign reads to references using a fulgor index 
	  stats           	 print index statistics 
	  print-filenames 	 print all reference filenames 


Demo
----

First, download [Cuttlefish](https://github.com/COMBINE-lab/cuttlefish).
Then, from within `cuttlefish`, do

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

which will run, in order, the tools `invert`, `sort-unique`, `permute-unitigs`, and `build`.

We now have an index serialized to the file `test_data/salmonella_10.hybrid.index`.



<!--Then, from within `fulgor/build`, invert the reference to unitig mapping with:

    ./fulgor invert -i ../test_data/salmonella_10 -g 1 -d tmp_dir --verbose

Deduplicate the color classes and build the map from unitig ids to color classes:

    ./fulgor sort_unique -i ../test_data/salmonella_10 -g 1 -d tmp_dir --verbose

Then permute the unitigs by color class:

    ./fulgor permute_unitigs -i ../test_data/salmonella_10 -g 1 -d tmp_dir --verbose

And finally build the index with:

    ./fulgor build -i ../test_data/salmonella_10 -k 31 -m 17 -d tmp_dir --verbose

Check correctness of colors:

    ./check_colors -i ../test_data/salmonella_10-->


Indexing an example Salmonella pan-genome
-----------------------------------------

In this example, we will build a Fulgor index for the 4,546 Salmonella genomes that can be downloaded from [here](https://zenodo.org/record/1323684).

We assume all commands are issue from within the home (`~/`) directory.

After download,
create a list of all `.fasta` filenames with

	find $(pwd)/Salmonella_enterica/Genomes/*.fasta > salmonella_4546_filenames.txt
	
and run Cuttlefish with

    ./src/cuttlefish build -l ~/salmonella_4546_filenames.txt -k 31 -t 8 -w tmp_dir/ -o ~/Salmonella_enterica/salmonella_4546 --extract-inverted-colors
    
Then from the `fulgor` parent directory, do

	python3 scripts/build_index.py --bin-dir build -k 31 -m 20 --tmp-dir build/tmp_dir -g 8 ~/Salmonella_enterica/salmonella_4546

which will create an index with the following stats:

	total index size: 0.266677 GB
	SPACE BREAKDOWN:
	  K2U: 66154784 bytes / 0.0661548 GB (24.8071%)
	  CCs: 199938374 bytes / 0.199938 GB (74.9741%)
	  Other: 583396 bytes / 0.000583396 GB (0.218765%)
	    U2C: 298192 bytes / 0.000298192 GB (0.111818%)
	    filenames: 285204 bytes / 0.000285204 GB (0.106948%)
	Color id range 0..4545
	Number of distinct color classes: 972178
	Number of ints in distinct color classes: 2139057102 (0.747763 bits/int)
	k: 31
	m: 20 (minimizer length used in K2U)
	Number of kmers in dBG: 43788757 (12.0862 bits/kmer)
	Number of unitigs in dBG: 1908149

Lastly, we can now pseudoalign the reads from [SRR801268](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR801/SRR801268/SRR801268_1.fastq.gz) with:

	./build/fulgor pseudoalign -i ~/Salmonella_enterica/salmonella_4546.hybrid.index -q ~/SRR801268_1.fastq.gz -t 8 -o /dev/null

using 8 parallel threads and writing the output to `/dev/null`.
