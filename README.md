# chromatin_remodeling

## CRL
CRL stands for "Chromatin Reprogramming Level". Given a starting, target and current chromatin accessibility states (for example, DNAase readcounts) across three conditions (for example, cell lines), the degree that the chromatin is remodeled in current condition is computed, after normalizing the count matrix with DESeq. CRL is computed only for the DHS sites that are significantly differentially accessible between the starting and target cell lines.

CRL is a continuous real number, and its density is typically highest in the range (0, 1). A CRL score of 1 means that the chromatin accessibility matches exactly that in the target condition, and a score of 0 means that chromatin accessibility is not remodeled compared to the starting condition. 
