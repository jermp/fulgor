Fulgor
======

Fulgor is a *(meta-) colored compacted de Bruijn graph* index for large-scale matching and color queries, powered by [SSHash](https://github.com/jermp/sshash) and [GGCAT](https://github.com/algbio/GGCAT).

The Fulgor index is described in the following paper.

[**Fulgor: A Fast and Compact k-mer Index for Large-Scale Matching and Color Queries**](https://drops.dagstuhl.de/opus/volltexte/2023/18644/)
(WABI, 2023)

And the meta-colored version in this pre-print: [**Meta-colored compacted de Bruijn graphs: overview and challenges**](https://www.biorxiv.org/content/10.1101/2023.07.21.550101v1) (bioRxiv, 2023).

Please, cite these papers if you use Fulgor.

### Table of contents
* [Dependencies](#dependencies)
* [Compiling the code](#compiling-the-code)
* [Tools](#tools)
* [Demo](#Demo)
* [Indexing an example Salmonella pangenome](#indexing-an-example-salmonella-pan-genome)

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

    git clone https://github.com/jermp/fulgor.git

and checkout the `mac-dbg` branch with

	git checkout mac-dbg

Then do

    git submodule update --init --recursive

to pull all necessary submodules before compilation.

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

	== Fulgor: a (meta-) colored compacted de Bruijn graph index =============================

	Usage: ./fulgor <tool> ...

	Available tools:
	  build              build a Fulgor index
	  pseudoalign        pseudoalign reads to references
	  stats              print index statistics
	  print-filenames    print all reference filenames
	  partition          partition a Fulgor index and build a meta-colored Fulgor index
	  permute            permute the reference names of a Fulgor index

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

to build an index that will be serialized to the file `test_data/salmonella_10.fur`.


Indexing an example Salmonella pangenome
----------------------------------------

In this example, we will build a Fulgor index, with k = 31, for the 4,546 Salmonella genomes that can be downloaded from [here](https://zenodo.org/record/1323684).

We assume all commands are issue from within the home (`~/`) directory.

After download,
create a list of all `.fasta` filenames with

	find $(pwd)/Salmonella_enterica/Genomes/*.fasta > salmonella_4546_filenames.txt

and, from `fulgor/build`, run

	./fulgor build -l ~/salmonella_4546_filenames.txt -o ~/Salmonella_enterica/salmonella_4546 -k 31 -m 20 -d tmp_dir -g 8 -t 8 --verbose --check

which will create an index named `~/Salmonella_enterica/salmonella_4546.fur` of 0.266 GB.

We can now pseudoalign the reads from SRR801268, as follows.

First, download the reads in `~/` with (assuming you have `wget` installed):

	cd
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR801/SRR801268/SRR801268_1.fastq.gz

and then process them with:

	./fulgor pseudoalign -i ~/Salmonella_enterica/salmonella_4546.fur -q ~/SRR801268_1.fastq.gz -t 8 -o /dev/null

	mapped 6584304 reads
	elapsed = 130133 millisec / 130.133 sec / 2.16888 min / 19.7641 musec/read
	num_mapped_reads 5797119/6584304 (88.0445%)

using 8 parallel threads and writing the mapping output to `/dev/null`.

To partition the index to obtain a meta-colored Fulgor index, then do:

	./fulgor partition -i ~/Salmonella_enterica/salmonella_4546.fur -d tmp_dir --check

The meta-colored index will be serialized to the file `~/Salmonella_enterica/salmonella_4546.mfur`
and will take 0.104 GB (2.55X smaller than the `.fur` index).
