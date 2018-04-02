The scripts and file in this repository pertain to comparative chromatin (or gene expression) reprogramming profiling analysis during cellular differentiation process. The concept of chromatin reprogramming level (CRL) - or analogous GRL, for gene expression - measures are discussed below. These scripts and results have been used to identify reprogrammed and non-reprogrammed sites of chromatin accessibility and their relation to gene expression during myogenic conversion of primary human fibroblasts - see this paper (https://academic.oup.com/nar/article/45/20/11684/4107215). 

## Chromatin Reprogramming Level (CRL)
CRL, short for "Chromatin Reprogramming Level", is a measure of the degree that chromatin accessibility is "reprogrammed" during some cellular differentiation (or protein knockout) process. Specifically, given a starting, target and current chromatin accessibility measures (for example, as indicated by normalized DNAase readcounts), CRL specifies the degree that the chromatin accessibility has changed in relation to the starting and target conditions. In relation to chromatin accessibility, CRL is computed only for the DHS sites that are significantly differentially accessible between the starting and target conditions (eg. cell lines).

CRL is computed as follows: 

$\frac{signal_b - signal_a}{signal_c - signalb}$, where
signal_a, signal_b and signal_c represent the normalized reads in starting state "a", current state "b" and target state "c".

CRL is a continuous real number, and its density is typically highest in the range (0, 1). A CRL score of 1 means that the chromatin accessibility matches exactly that in the target condition, and a score of 0 means that chromatin accessibility is not remodeled compared to the starting condition. 


## Gene expression Reprogramming level (GRL)
Similar to CRL in concept, GRL pertains to quantifying the degree of reprogramming in the expression level of a gene. For comparison, transcripts per million (TPM) expression signal should be used instead of RPKM or FPKM values, that are not directly comparable across cell types or conditions.

## Application in time-course studies

The CRL, GRL or other similar measures for histone modification data, for example, can also be used in time-course data. If the initial (t=0) and final condition are fixed, the CRL (or GRL) measures for intermediate signals would correspond to **the level that chromatin (or gene expression) is reprogrammed** at those time points. However, as mentioned above, only those data points that are significantly differential between these initial and final conditions can be included given the definition of CRL and GRL measures.

## Contact

Please do not hesitate to contact dinesh(dot)mdh01(at)gmail.com to discuss about these measures and their applications.
