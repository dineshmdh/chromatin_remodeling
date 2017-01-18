# Created on 8/3/2015

# Idea on the protocol:
#   Estimate library size for all four cell lines together, and use only the sizes for Fu and Mb replicates
#   when nbinomTest is performed (between Fu and Mb). FibCell (Fd or Fc) replicates are normalized using the sizeFactor 
#   estimated. Then filtering based on the CRL values, fold change and p-value of significance gives the desired output.
# Load the result image for quick access to the results.
#
# chrM and chrX/Y are removed right after the count matrix is loaded.
# Just checked the countMatrix. There are no chrs other than the standard chromosomes. There is also no chrM.

# Naming of the objects.
# In the data below, mbOpen_well and mbOpen_badly (and corresponding mbClosed_* objects) are only those that are filtered
# for their CRL_wFd values. So, any pk that has CRL_wFd of > 1.1 or < -0.1 are not include in the list. 
#
# What is the best number of mbOpen/mbClosd "well" and "badly" reprogrammed sites to take for downstream analysis?
# Using FDR analysis, for an FDR of 0.01 it seems that 1k well reprogrammed sites is reasonable for mbOpen set, 
# and ~1.5k for mbClosed. for an FDR of 0.05, the number is ~2.5k for both sets. 

require(DESeq)
require(ggplot2)

# set some parameters
pVal_cutOff = 0.01 # using pVal_cutOff as the threshold padj value to call differentially expressed peaks
pkSize = "top100k" # either "all" or "top100k" # used in the directory specification..
sizeFactorByModel = TRUE
getTopPksOnly = TRUE

if (getTopPksOnly){
  getTopXPks = 1000 # get top X peaks in each DHS category (eg. MbOpen_well etc)
  wellReprog_upperLim = 1.1
  badlyReprog_lowerLim = -0.1
}else{ # use the following constraints as before.
  wellReprog_upperLim = 1.20
  wellReprog_lowerLim = 0.75
  badlyReprog_upperLim = 0.20
  badlyReprog_lowerLim = -0.05
}


# Fix the input and output directories.
inputDir = "/Users/Dinesh/Dropbox/Lab_work/GordanLab/myoDNew/DNase_new/DESeq/DESeq_input"
outputDir = "/Users/Dinesh/Dropbox/Lab_work/GordanLab/myoDNew/DNase_new/DESeq/DESeq_output/wDESeq/using_FuMb_top100kPks/DESeq_outFiles_orig_afterDeconvolution"

# if to just load the result_dataset
#load(paste(outputDir,"/allResults_generated.RData", sep=""))

setwd(inputDir)

# ============================================
#
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Run DESeq algorithm. 
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================

counts = read.table(paste("FuMbFdFc.allreps.dnaseReadCounts.usingWholePk.wFlankBp0.forFu_", pkSize, "_Mb_", pkSize, "_dnasepks.rep0.q0.01.mergedwD0.txt", sep=""))
names(counts) = c("Chr", "ss", "es", "Fu_r1", "Fu_r2", "Fu_r3", "Mb_r1", "Mb_r2", "Mb_r3", "Mb_r4", "Fd_r1", "Fd_r2", "Fd_r3", "Fc_r1", "Fc_r2", "Fc_r3")

# remove chrM, chrX/Y; there are no non-standard chromosomes. 
autosomalChrIndices = (counts$Chr != "chrM") & (counts$Chr != "chrX") & (counts$Chr != "chrY")
counts = counts[autosomalChrIndices,]


# Change the working directory to output directory now (that we have the count matrix loaded).
#=====================================================
setwd(outputDir)
#=====================================================

# making dhsLocs the rowNames
rownames(counts) = paste(counts$Chr, paste(counts$ss, counts$es, sep="-"), sep=":")
countsOnly = subset(counts, select=4:16)
FuVsMb_countsOnly = subset(countsOnly, select=1:7)

#=====================================================
#
#Run the nBinomTest between Fu and Mb dataset
#
#=====================================================
conds = factor(c("Fu", "Fu", "Fu", "Mb", "Mb", "Mb", "Mb", "Fd", "Fd", "Fd", "Fc", "Fc", "Fc"))
cds = newCountDataSet(countsOnly, conds)
cds = estimateSizeFactors(cds)
sizeFactors = cds$sizeFactor # has to be saved now before the nbinom test here below. 
df_normedFdFc_only = counts(cds, normalized=TRUE)[,8:13] # saving now to normalize the Fd and Fc data below.

conds = factor(c("Fu", "Fu", "Fu", "Mb", "Mb", "Mb", "Mb"))
cds = newCountDataSet(FuVsMb_countsOnly, conds) # instantiate a new countDataSet object.
cds$sizeFactor = sizeFactors[1:7]
cds = estimateDispersions(cds) 
pdf("Dispersion.pdf")
plotDispEsts( cds )
dev.off()
result = nbinomTest(cds, "Fu", "Mb")

# Some visualization of the results :: MA plot requires 'x' to be named "baseMean".
pdf("MA_plot.pdf")
plotMA(result, ylim=c(-7, 7))
dev.off()
pdf("pVal_adj.pdf")
hist(result$padj, breaks=200, col="skyblue", border="slateblue", main="", xlab="Adjusted pValues")
dev.off()

colnames(result) = c("dhs_locs", "baseMean_FuMb", "baseMean_Fu", "baseMean_Mb", "FC_MbByFu", "log2FC_MbByFu", "pval_FuMb", "padj_FuMb")
rownames(result) = result$dhs_locs # set dhs_locs as the rownames (useful in df merging below)
result = result[-1] # now we don't need the dhs_locs column

# ==========================================
#
# Now merge all: the normalized data and the result
# Also add the baseMean_Fd and baseMean_Fc columns
#
# ==========================================

normed_fumb_df = counts(cds, normalized=TRUE)
stopifnot(rownames(normed_fumb_df) == rownames(result))
stopifnot(rownames(normed_fumb_df) == rownames(df_normedFdFc_only))
dhs_df = cbind(normed_fumb_df, df_normedFdFc_only, result)
dhs_df$baseMean_Fd = rowMeans(cbind(dhs_df$Fd_r1, dhs_df$Fd_r2, dhs_df$Fd_r3))
dhs_df$baseMean_Fc = rowMeans(cbind(dhs_df$Fc_r1, dhs_df$Fc_r2, dhs_df$Fc_r3))

# ============================================
#
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Get all significant results along w/ CRL_Fd and CRL_Fc info.
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================

# Now get only the significantly differential dhs peaks (USE FILTERS)
dhs_df_noNA = is.na(dhs_df$padj_FuMb) == FALSE
dhs_df_sigPVal = dhs_df$padj_FuMb < pVal_cutOff # For Fd; ~300 when pVal is 0.00001 vs. 2200 when it is 10E-3.
dhs_df_condFC = abs(dhs_df$log2FC_MbByFu) >= 1 # condFC == condition on Fold Change;
dhs_df_sig = dhs_df[dhs_df_noNA & dhs_df_sigPVal & dhs_df_condFC,]

# ==========================================
#
# Get the CRL info for all DHSs in dhs_df_sig
#
# ==========================================
dhs_df_sig$crl_wFd = (dhs_df_sig$baseMean_Fd - dhs_df_sig$baseMean_Fu) / (dhs_df_sig$baseMean_Mb - dhs_df_sig$baseMean_Fu)
dhs_df_sig$crl_wFc = (dhs_df_sig$baseMean_Fc - dhs_df_sig$baseMean_Fu) / (dhs_df_sig$baseMean_Mb - dhs_df_sig$baseMean_Fu)

# ============================================
#
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ===== Get MbOpen and MbClosed dataframes ===
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================
dhs_df_sig_backup = dhs_df_sig

dhs_df_sig$mbStatus = rep("unss", dim(dhs_df_sig)[1])
dhs_df_sig$mbStatus[dhs_df_sig$log2FC_MbByFu >= 1] = "mbOpen"
dhs_df_sig$mbStatus[dhs_df_sig$log2FC_MbByFu <= -1] = "mbClosed"
stopifnot(sum(dhs_df_sig$mbStatus == "mbOpen") + sum(dhs_df_sig$mbStatus == "mbClosed") == dim(dhs_df_sig)[1])

# bound the dhs_df_sig sites by well and badly reprog upper and lower lim conditions
# this removes about 340  dhs sites genomewide (belonging to both mbOpen and mbClosed combined)
dhs_df_bounded = dhs_df_sig$crl_wFd <= wellReprog_upperLim & dhs_df_sig$crl_wFd >= badlyReprog_lowerLim # bounded by upper and lower lim
mbOpen_df_bounded = dhs_df_sig[dhs_df_bounded & dhs_df_sig$mbStatus == "mbOpen",] # 7489 by 25
mbOpen_df_bounded = mbOpen_df_bounded[order(mbOpen_df_bounded$crl_wFd, decreasing = TRUE),]

mbClosed_df_bounded = dhs_df_sig[dhs_df_bounded & dhs_df_sig$mbStatus == "mbClosed",] # 4782 by 25
mbClosed_df_bounded = mbClosed_df_bounded[order(mbClosed_df_bounded$crl_wFd, decreasing = TRUE),]

# ============================================
#
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# === CALCULATE THE FDR VALUES ====
# assuming we are taking the top 1k or bottom 1k sites
# as well or badly reprog sites respectively
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================

# ========= For MbOpen_well ==========
fdrs_mbOpen_well = rep(-10, 10)
getTheseMany = seq(1:20) * 100 + 1000 # used this to narrow down for fdr 0.05; first started w/ seq(1:10) * 100 + 1000 
getTheseMany = c(getTopXPks, getTheseMany)

for (i in seq(1:21)){
  mbOpen_well = mbOpen_df_bounded[1:getTheseMany[i],] # the crl_fd order is decreasing from the top
  fp_mbOpen_well = sum(mbOpen_well$crl_wFc >= min(mbOpen_well$crl_wFd)) # if we use mbOpen_allInfo$CRL_wFc instead, fdr=0.018
  fdr_mbOpen_well = fp_mbOpen_well/getTheseMany[i]
  fdrs_mbOpen_well[i] = fdr_mbOpen_well
}
pdf("how_many_mbOpen_well_dhss_to_take.pdf")
plot(x = getTheseMany, y = fdrs_mbOpen_well, type="l", main="mbOpen", xlab="Number of DHS sites below CRL 1.1", ylab="FDR")
abline(h=0.05, col="red"); abline(h=0.01, col = "red")
dev.off()
# From the plot it looks like for an fdr of 0.05, we could take ~2500 top DHS sites. For fdr 0.01, ~1k.


# ======== For MbClosed_well ==========
fdrs_mbClosed_well = rep(-10, 10)
getTheseMany = seq(1:20) * 100 + 1000
getTheseMany = c(getTopXPks, getTheseMany)

for (i in seq(1:21)){
  mbClosed_well = mbClosed_df_bounded[1:getTheseMany[i],] # the crl_fd order is decreaseing from the top
  fp_mbClosed_well = sum(mbClosed_well$crl_wFc >= min(mbClosed_well$crl_wFd)) # if we use mbOpen_allInfo$CRL_wFc instead, fdr=0.005 still.
  fdr_mbClosed_well = fp_mbClosed_well/getTheseMany[i]
  fdrs_mbClosed_well[i] = fdr_mbClosed_well
}
pdf("how_many_mbClosed_well_dhss_to_take.pdf")
plot(x = getTheseMany, y = fdrs_mbClosed_well, type="l", main="mbClosed", xlab="Number of DHS sites below CRL 1.1", ylab="FDR")
abline(h=0.05, col="red");  abline(h=0.01, col = "red")
dev.off()
# From the plot it looks like for an fdr of 0.05, we could take ~2800 top DHS sites. 



# ============================================
#
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# === NOW SAVING DATA AND WRITING FILES ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================
save.image("allResults_generated.RData")
write.table(dhs_df_sig, file="dhs_df_sig.txt", quote=F, sep="\t")
write.table(mbOpen_df_bounded, file="mbOpen_df_bounded.txt", quote=F, sep="\t")
write.table(mbClosed_df_bounded, file="mbClosed_df_bounded.txt", quote=F, sep="\t")

# ============================================
#
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ========= NOW SAVING THE FIGURES ===========
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================

pdf("crls_mbOpen.pdf")
plot(density(mbOpen_df_bounded$crl_wFc), col="pink", lwd=2, main="mbOpen", xlim=c(-1,1))
lines(density(mbOpen_df_bounded$crl_wFd), col="darkred", lwd=2)
dev.off()

pdf("crls_mbClosed.pdf")
plot(density(mbClosed_df_bounded$crl_wFc), col="lightblue", lwd=2, main="mbClosed", xlim=c(-1,1))
lines(density(mbClosed_df_bounded$crl_wFd), col="darkblue", lwd=2)
dev.off()

# ============================================
#
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# === Are non-reprogrammed sites wider? =====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================

# working on mbOpen
sses = unlist(lapply(strsplit(rownames(mbOpen_df_bounded),":"), function(x) x[2]))
sses_v = as.numeric(unlist(strsplit(sses, "-"))) # as vector
sses_m = matrix(sses_v, ncol=2, byrow=T) # as matrix
dhs_len_mbOpen = sses_m[,2] - sses_m[,1]
pdf("mbOpen_crlFd_vs_dhsLen.pdf")
scatter.smooth(x = mbOpen_df_bounded$crl_wFd, y= dhs_len_mbOpen, ylim=c(0,1400), col="#CCCCCC")
#plot(x = mbOpen_df_bounded$crl_wFd, y= dhs_len_mbOpen, pch=16, cex=0.2)
#qplot(x = mbOpen_df_bounded$crl_wFd, y= dhs_len_mbOpen, geom='smooth', span =5)
dev.off()
# the correlation seems to ~0..

# working on mbClosed
sses = unlist(lapply(strsplit(rownames(mbClosed_df_bounded),":"), function(x) x[2]))
sses_v = as.numeric(unlist(strsplit(sses, "-"))) # as vector
sses_m = matrix(sses_v, ncol=2, byrow=T) # as matrix
dhs_len_mbClosed = sses_m[,2] - sses_m[,1]
pdf("mbClosed_crlFd_vs_dhsLen.pdf")
scatter.smooth(x = mbClosed_df_bounded$crl_wFd, y= dhs_len_mbClosed, ylim=c(0,1400), col="#CCCCCC")
dev.off()


# ============================================
#
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# === GENERATE PLOTS USING MEDIAN VALUES===
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================
# fu_med = apply(cds_new[, 1:3], 1, median)





# ============================================
#
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# === Assuming Fd is x% fib ===
# removing the Fib "contribution" to recompute the CRLs
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================

percent_fib = c(80, 80, 80) # for fd1, fd2 and fd3 respectively
percent_fib_info = paste(percent_fib[1], percent_fib[2], percent_fib[3], sep="_")

# first work on mbOpen
mbOpen_df_bounded$Fd_r1_new = ((100*mbOpen_df_bounded$Fd_r1) - (percent_fib[1]*mbOpen_df_bounded$baseMean_Fu)) / (100 - percent_fib[1])
mbOpen_df_bounded$Fd_r2_new = ((100*mbOpen_df_bounded$Fd_r2) - (percent_fib[2]*mbOpen_df_bounded$baseMean_Fu)) / (100 - percent_fib[2])
mbOpen_df_bounded$Fd_r3_new = ((100*mbOpen_df_bounded$Fd_r3) - (percent_fib[3]*mbOpen_df_bounded$baseMean_Fu)) / (100 - percent_fib[3])
mbOpen_df_bounded$baseMean_Fd_new = rowMeans(cbind(mbOpen_df_bounded$Fd_r1_new, mbOpen_df_bounded$Fd_r2_new, mbOpen_df_bounded$Fd_r3_new))
mbOpen_df_bounded$crl_wFd_new = (mbOpen_df_bounded$baseMean_Fd_new - mbOpen_df_bounded$baseMean_Fu) / (mbOpen_df_bounded$baseMean_Mb - mbOpen_df_bounded$baseMean_Fu)

pdf(paste("crls_mbOpen_",percent_fib_info,".pdf", sep=""))
plot(density(mbOpen_df_bounded$crl_wFc), col="pink", lwd=2, 
     main=paste("mbOpen: %Fib is ", percent_fib_info, sep=""), xlim=c(-0.5,4))
lines(density(mbOpen_df_bounded$crl_wFd), col="darkred", lwd=2)
lines(density(mbOpen_df_bounded$crl_wFd_new), col="black", lwd=2)
dev.off()


# now working with mbClosed
mbClosed_df_bounded$Fd_r1_new = ((100*mbClosed_df_bounded$Fd_r1) - (percent_fib[1]*mbClosed_df_bounded$baseMean_Fu)) / (100 - percent_fib[1])
mbClosed_df_bounded$Fd_r2_new = ((100*mbClosed_df_bounded$Fd_r2) - (percent_fib[2]*mbClosed_df_bounded$baseMean_Fu)) / (100 - percent_fib[2])
mbClosed_df_bounded$Fd_r3_new = ((100*mbClosed_df_bounded$Fd_r3) - (percent_fib[3]*mbClosed_df_bounded$baseMean_Fu)) / (100 - percent_fib[3])
mbClosed_df_bounded$baseMean_Fd_new = rowMeans(cbind(mbClosed_df_bounded$Fd_r1_new, mbClosed_df_bounded$Fd_r2_new, mbClosed_df_bounded$Fd_r3_new))
mbClosed_df_bounded$crl_wFd_new = (mbClosed_df_bounded$baseMean_Fd_new - mbClosed_df_bounded$baseMean_Fu) / (mbClosed_df_bounded$baseMean_Mb - mbClosed_df_bounded$baseMean_Fu)

pdf(paste("crls_mbClosed_",percent_fib_info,".pdf", sep=""))
plot(density(mbClosed_df_bounded$crl_wFc), col="lightblue", lwd=2, 
     main=paste("mbClosed: %Fib is ", percent_fib_info, sep=""), xlim=c(-0.5,4))
lines(density(mbClosed_df_bounded$crl_wFd), col="darkblue", lwd=2)
lines(density(mbClosed_df_bounded$crl_wFd_new), col="black", lwd=2)
dev.off()

