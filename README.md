Fulgor
======

Fulgor is a *colored compacted de Bruijn graph* index for large-scale matching and color queries, powered by [SSHash](https://github.com/jermp/sshash) and [GGCAT](https://github.com/algbio/GGCAT).

The index is described in the following paper.

[**Fulgor: A Fast and Compact k-mer Index for Large-Scale Matching and Color Queries**](https://drops.dagstuhl.de/opus/volltexte/2023/18644/)
(WABI 2023)

Please, cite this paper if you use Fulgor.

### Table of contents
* [Dependencies](#dependencies)
* [Compiling the code](#compiling-the-code)
* [Tools](#tools)
* [Demo](#Demo)
* [Indexing an example Salmonella pan-genome](#indexing-an-example-salmonella-pan-genome)


Dependencies
------------

#### GGCAT

The code uses the [GGCAT](https://github.com/algbio/GGCAT) Rust library,
so make sure you have Rust installed. If not, Rust can be installed as recommended [here](https://www.rust-lang.org/tools/install), with

	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

It is also recommended to use the [Nightly](https://doc.rust-lang.org/book/appendix-07-nightly-rust.html#rustup-and-the-role-of-rust-nightly) release channel.

	rustup toolchain install nightly

#### zlib

If you do not have `zlib` installed, you can do

    sudo apt-get install zlib1g

if you are on Linux/Ubuntu, or

    brew install zlib

if you are using MacOS.


Compiling the code
------------------

The code is tested on Linux with `gcc` and on MacOS with `clang`.
To build the code, [`CMake`](https://cmake.org/) is required.

First clone the repository with

    git clone --recursive https://github.com/jermp/fulgor.git

If you forgot `--recursive` when cloning, do

    git submodule update --init --recursive

before compiling.

To compile the code for a release environment (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following, within the parent `fulgor` directory:

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
	  build           	 build a fulgor index
	  pseudoalign     	 pseudoalign reads to references using a fulgor index
	  stats           	 print index statistics
	  print-filenames 	 print all reference filenames


Demo
----

This short demo shows how to index the 10-genome collection
in the folder `test_data/salmonella_10` with Fulgor.
We will use the standard value k = 31.

First create a list of filenames (with absolute paths) for the files in `test_data/salmonella_10`.
From `fulgor/test_data`, do

	find $(pwd)/salmonella_10/* > salmonella_10_filenames.txt

Then, from `fulgor/build`, run

	./fulgor build -l ../test_data/salmonella_10_filenames.txt -o ../test_data/salmonella_10 -k 31 -m 19 -d tmp_dir -g 1 -t 1 --verbose --check

to build an index that will be serialized to the file `test_data/salmonella_10.hybrid.index`.


Indexing an example Salmonella pan-genome
-----------------------------------------

In this example, we will build a Fulgor index, with k = 31, for the 4,546 Salmonella genomes that can be downloaded from [here](https://zenodo.org/record/1323684).

We assume all commands are issue from within the home (`~/`) directory.

After download,
create a list of all `.fasta` filenames with

	find $(pwd)/Salmonella_enterica/Genomes/*.fasta > salmonella_4546_filenames.txt

and, from `fulgor/build`, run

	./fulgor build -l ~/salmonella_4546_filenames.txt -o ~/Salmonella_enterica/salmonella_4546 -k 31 -m 20 -d tmp_dir -g 8 -t 8 --verbose --check

which will create an index with the following stats:

	total index size: 0.266428 GB
	SPACE BREAKDOWN:
	  K2U: 65891238 bytes / 0.0658912 GB (24.7314%)
	  CCs: 199938374 bytes / 0.199938 GB (75.0442%)
	  Other: 597956 bytes / 0.000597956 GB (0.224435%)
	    U2C: 294568 bytes / 0.000294568 GB (0.110562%)
	    filenames: 303388 bytes / 0.000303388 GB (0.113873%)
	Color id range 0..4545
	Number of distinct color classes: 972178
	Number of ints in distinct color classes: 2139057102 (0.747763 bits/int)
	k: 31
	m: 20 (minimizer length used in K2U)
	Number of kmers in dBG: 43788757 (12.038 bits/kmer)
	Number of unitigs in dBG: 1884865

We can now pseudoalign the reads from [SRR801268](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR801/SRR801268/SRR801268_1.fastq.gz) with:

	./fulgor pseudoalign -i ~/Salmonella_enterica/salmonella_4546.hybrid.index -q ~/SRR801268_1.fastq.gz -t 8 -o /dev/null

	mapped 6584304 reads
	elapsed = 130133 millisec / 130.133 sec / 2.16888 min / 19.7641 musec/read
	num_mapped_reads 5797119/6584304 (88.0445%)

using 8 parallel threads and writing the mapping output to `/dev/null`.
