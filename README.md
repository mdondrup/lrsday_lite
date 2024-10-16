# lrsday_lite
A snakemakified version of LRSDAY. The annotation section only, not the assembly part.
The shell scripts have been implemented as snakemake rules, and the dependencies have been ported to
conda environments. Several tools are downloaded, compiled, installed, and configured outside of conda.

Installs all dependencies:

snakemake -c 8 --use-conda setup

RepeatMasker libraries will be downloaded during setup and RepeatMaster configuration will be run interactively. When asked to select the search engine, select RPMBlast and continue.

To run the annotation, place the assembly in the directory 07.Supervised_Final_Assembly, by default, the workflow will look for CPG_1a.assembly.final.fa (same as LRSDAY). The prefix can be configured in config/config.yaml or on the commandline. Mitochondrial chromosomes are expected to be present with the sequence id prefixed with "MT" (configuration parameter chrMT).




