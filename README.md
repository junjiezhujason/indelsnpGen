# indelsnpGen Source Code

This Code is used to simulate SNPs and INDELs in a reference genome to create a target sequence. The input is the fasta file of a chormosome and the output is a fasta file of the target sequence, as well as a BCF/VCF file that records the simulated SNPs and INDELs. 

The input parameters into this program include snp rate, indel rate and the bases (AGTC).

This code was built by Ki Sing (Ben) Chan, and modified by Junjie (Jason) Zhu. 

Modfications include:

- documentation in the source code 

- Keeping the ambigous bases (N) in the target sequence, and handling the lower-case/ upper-case representation of the bases.

- Recalcuation of the position random variant generation



