---
title: 'Toolkit for preparing genomes for submission to NCBI'
tags:
  - bioinformatics
  - computational biology
authors:
 - name: Rutger A. Vos
   orcid: 0000-0001-9254-7318
   affiliation: Naturalis Biodiversity Center, P.O. Box 9517, 2300 RA, Leiden, The Netherlands
date: 28 March 2019
bibliography: paper.bib
---

# Summary

__wgs2ncbi__ is a standalone Perl package for preparing the results of whole genome sequencing (WGS) and annotation projects to NCBI. The functionality consists of a command line program that takes sub-commands for the various steps of the process: i) preparing the input data for rapid processing downstream (sub-command `prepare`); ii) generating valid scaffolds and feature tables from the sequence data and genome annotation (`process`); iii) converting the processing results to ASN.1/seqin files for upload to NCBI (`convert`); iv) compressing the converted files to archives acceptable to NCBI's upload service (`compress`). The process of submitting the results - and having them pass review - is a labour intensive, iterative process due to the stringent quality requirements that NCBI imposes [@Pirovano2015]. These requirements cannot typically be met on the first submission because they involve checks based on the entirety of NCBI's presently held data (e.g. contamination checks), checks for the continuously moving target of vendor-specific adaptor sequences, and checks for the validity of gene product names. In many genome annotation projects, these product names are copied over from homologous sequences in related model organisms. This may introduce terminology that is not, or no longer, permitted by NCBI, such as molecular weights and protein structure, organism names, database identifiers, and so on. The toolkit allows for easy masking of detected contaminations and adaptors, a generalized mapping between invalid product names (as detected by NCBI) and valid alternatives, and automatic conversion of putative-but-invalid genes (e.g. those with introns that are 'too short') to pseudogenes.

# References
