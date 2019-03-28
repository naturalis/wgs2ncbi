---
title: 'Toolkit for preparing genomes for submission to NCBI'
tags:
  - bioinformatics
  - computational biology
authors:
 - name: Rutger Aldo Vos
   orcid: 0000-0001-9254-7318
   affiliation: Research Group 'Understanding Evolution', Naturalis Biodiversity Center, Leiden, The Netherlands
 - name: Nnadi Nnaemeka Emmanuel
   affiliation: Department of Microbiology, Faculty of Natural and Applied Science, Plateau State University, Bokkos, Plateau State, Nigeria.
date: 28 March 2019
bibliography: paper.bib
---

# Summary

The process of submitting annotated genomes to NCBI GenBank [@GenBank2015] - and having them pass review - is a labour intensive, iterative process due to the stringent quality requirements that NCBI imposes [@Pirovano2015]. These requirements cannot typically be met on the first iteration because they involve checks based on the entirety of NCBI's presently held data (e.g. contamination checks), checks for the continuously moving target of vendor-specific adaptor sequences, and checks for the validity of gene product names. In many genome annotation projects, these product names are copied over from homologous sequences in related model organisms. This may introduce terminology that is not, or no longer, permitted by NCBI, such as molecular weights and protein structure, organism names, database identifiers, and so on. 

__wgs2ncbi__ is a standalone Perl package for preparing the results of whole genome sequencing (WGS) and annotation projects to NCBI GenBank. The purpose of the package is to automate responding to NCBI's reviews by allowing batch corrections to detected problems. The functionality consists of a command line program that takes sub-commands for the various steps of the process: 

1. Preparing the input data for rapid processing downstream (sub-command `prepare`).
2. Generating valid scaffolds and feature tables from the sequence data and genome annotation (`process`).
3. Converting the processing results to ASN.1/seqin files for upload to NCBI (`convert`).
4. Compressing the converted files to archives acceptable to NCBI's upload service (`compress`). 

In step 2., the toolkit allows for easy masking of detected contaminations and adaptors, a generalized mapping between invalid product names (as detected by NCBI) and valid alternatives, and automatic conversion of putative-but-invalid genes (e.g. those with introns that are 'too short') to pseudogenes. This functionality has helped make the submission of the genome of the King Cobra [@Vonk2013] and ([AZIM00000000.1](https://www.ncbi.nlm.nih.gov/nuccore/AZIM00000000.1)), and that of the Velvet Bean ([QJKJ00000000.1](https://www.ncbi.nlm.nih.gov/nuccore/QJKJ00000000.1)).

# References
