
# The program HLAbaits generates a target sequence tiling to design probes for targeted enrichment. 
If you want to have an existing design for HLA enrichment look into the design folder.


## Installing the tool

To install this pipeline, simply clone the repository to a location:

`git clone https://github.com/MiWitt/HLAbaits.git`

Move to the project directory:

`cd HLAbaits`

To compile the tool run make:

`make`

To update the code, run git update inside of the local clone:

`git update`

## Usage
run like: 
    HLAbaits probe_length mis_matches tiling_step ref_seq_fasta scratch [existing_probes_fasta]

Parameters
    probe_length:  define the length of the probes here (recommended 120)
    mis_matches:   for each new fasta sequence all existing probes are aligned to it using blat. Define the number of allowed mismatches here (recommended 5)
    tiling_step:   for nun overlapping probes set tiling step = probe_length. Set it to probe_length/2 for half overlapping (recommended 60)
    ref_seq_fasta: Provide your collection of fasta sequences here. Set the path to a multiple fasta file
    scratch:       Set the path to a scratch directory where HLAbaits has read write access
    existing_probes_fasta]: Provide a file of already existing probes here. Be aware of the important note below


All output goes to stdout, so redirect the output to your destination file. (e.g. > /foo/bar/my_baits.fa)

IMPORTANT: if you provide an existing probes (baits) file, provide it in following format. First line unique probe name, second line probe sequence, third line unique probe name, etc. NO BLANK LINES, NO REDUNDANT PROBE NAMES.

IMPORTANT: This software depends on a preinstalled blat (https://genome.ucsc.edu/goldenPath/help/blatSpec.html) that must be executable by simply typing 'blat' in a terminal!


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details



