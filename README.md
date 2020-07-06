#Phylogenetic-pipeline

In order to run the script you will need following programs installed and added to your computers PATH:
-blastp
-mcl
-clustalo
-clann
also you may want to install pandas, biopython, ete and all their dependencies if you don't have them already.

Phylogenetic pipeline used in the python script is fully automatic, which means that all you need to do is provide a txt file with standard NCBI ids and script will do the rest. "The rest" means download files, blast them, cluster them, perform MSA, build NJ trees and at the end build supertree. After each bigger or smaller step, user will get suitable communicate just to know on which stage of work pipeline is at.
Despite this, that the pipeline is meant to be run in the "full" mode, you can stil choose if you don't want to run some clue parts of it (and have right files in suitable folders). More about it below, in argument No.6 description.

Default usage:

python3 phylogenetic_pipeline.py test test.txt group 3 nj pdbcmts

Arguments:
1) Name of your analysis used to create a suitable directory.
2) Txt file with standard NCBI ids that can by separated by any of this: '!@/#$|:'.
3) Whatever you want to name your files e.g. escherichia_coli, coronaviruses etc.
4) Cutoff - clusters with less genes than it are dropped.
5) Supertree reconstruction method, you can use any from the available in clann manual.
6) Mode - which part of pipeline do you want to use?
   - p - precheck of needed directories and creating them if necessary
   - d - download the files corresponding to ids you provided
   - b - make blast database and blast sequences form "data" directory
   - c - cluster using mcl, create unique fasta file for every cluster
   - m - perform MSA with clustalo on the files from clusters directory
   - t - build NJ trees using MSAs from MSA directory
   - s - build supertree using chosen method using NJ trees
