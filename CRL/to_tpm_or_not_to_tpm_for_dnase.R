# Created on Jan 23, 2016
# Inspiration for this script: http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

# Some Interesting notes
# 1. Since "tpm" conversion is relative (to 1 million total reads), the absolute measure of 
#    accessibility is lost with this conversion.
#
# 2. Reordering levels in cds_all does not change the sizeFactor plot, but the x-axis labels is not asserted;
#    But since cds_all$condition is the same vector of all samples (despite different levels), the x-axis labels
#    should also be the same. 
#
# 3. The p-value plot is "uniform", but there are also a few "peaks"?

library(dplyr)

wellReprog_upperLim = 2
badlyReprog_lowerLim = -1


# set some parameters
pVal_cutOff = 0.01 # using pVal_cutOff as the threshold padj value to call differentially expressed peaks
pkSize = "top100k" # either "all" or "top100k" # used in the directory specification..

# Fix the input and output directories.
inputDir = "/Users/Dinesh/Dropbox/Lab_work/GordanLab/myoDNew/DNase_new/DESeq/DESeq_input"
outputDir = "/Users/Dinesh/Dropbox/Lab_work/GordanLab/myoDNew/DNase_new/DESeq/DESeq_output/wDESeq/using_FuMb_top100kPks/noDESeq_justTPM"
setwd(inputDir)

counts = read.table(paste("FuMbFdFc.allreps.dnaseReadCounts.usingWholePk.wFlankBp0.forFu_", pkSize, "_Mb_", pkSize, "_dnasepks.rep0.q0.01.mergedwD0.txt", sep=""))
names(counts) = c("Chr", "ss", "es", "Fu_r1", "Fu_r2", "Fu_r3", "Mb_r1", "Mb_r2", "Mb_r3", "Mb_r4", "Fd_r1", "Fd_r2", "Fd_r3", "Fc_r1", "Fc_r2", "Fc_r3")

# Change the working directory to output directory (now that we have the count matrix loaded).
#=====================================================
setwd(outputDir)
#=====================================================

# remove chrM, chrX/Y; there are no non-standard chromosomes. 
autosomalChrIndices = (counts$Chr != "chrM") & (counts$Chr != "chrX") & (counts$Chr != "chrY")
counts = counts[autosomalChrIndices,]

# making dhsLocs the rowNames
counts$dhsLen = (counts$es - counts$ss)/100 # in 100bp scale
rownames(counts) = paste(counts$Chr, paste(counts$ss, counts$es, sep="-"), sep=":")
counts = subset(counts, select=4:ncol(counts))

# plot the boxplot before TPM normalization
pdf("boxPlot_counts_beforeNorm.pdf")
boxplot(log2(counts), ylab=c("log2 counts before norm"))
dev.off()

#=====================================================
# ==== Normalize as in TPM ===== 
#=====================================================

# Step1: Normalize by DHS length (in the scale of 100bp)
counts_scaledByLen = apply( select(counts, Fu_r1:Fc_r3), 2, function(i){i/counts$dhsLen})

# Step2: Normalize by the sum of the columns
scaleBy = as.numeric(colSums(counts_scaledByLen))
counts_tpm = data.frame(t(apply(counts_scaledByLen, 1, function(i){i/scaleBy})))
counts_tpm = counts_tpm * 1000000
colnames(counts_tpm) = colnames(counts)[1:13] # upto Fc_r3

pdf("counts_tpm_colSums.pdf")
plot(colSums(counts_tpm))
dev.off()

pdf("counts_ceilingOfTPM_colSums.pdf")
plot(colSums(ceiling(counts_tpm)))
dev.off()

pdf("boxPlot_countsTPM.pdf")
boxplot(log2(counts_tpm), ylab=c("log2 counts_tpm"))
dev.off()

counts_orig = counts # for backup
counts = counts_tpm

#=====================================================
#
#Run the nBinomTest between Fu and Mb dataset
#
#=====================================================

#=====================================================
# see what the new sizeFactors look like (# see note # 2 above)
#=====================================================
conds_all = factor(c("Fu", "Fu", "Fu", "Mb", "Mb", "Mb", "Mb", "Fd", "Fd", "Fd", "Fc", "Fc", "Fc"), levels = c("Fu", "Mb", "Fd", "Fc"))
cds_all = newCountDataSet(ceiling(counts_tpm), conds_all)
cds_all = estimateSizeFactors(cds_all)
sizeFactors = cds_all$sizeFactor 
pdf("sizeFactors_afterTPM_withDESeq.pdf") # see note # 2 above (saved a file *_levelsReordered.pdf)
plot(sizeFactors)
dev.off()

# now make the sizeFactors all 1 (b/c we don't want to normalize by the size factor)
sizeFactors = rep(1, ncol(counts_tpm))

#=====================================================
# run DESeq against Fu and Mb
#=====================================================

conds = factor(c("Fu", "Fu", "Fu", "Mb", "Mb", "Mb", "Mb"), levels = c("Fu", "Mb"))
cds = newCountDataSet(ceiling(counts_tpm[,1:7]), conds) # instantiate a new countDataSet object.
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

stopifnot(rownames(counts_tpm) == rownames(result))
dhs_df = cbind(counts_tpm, result)
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

# plot the CRLs
pdf("CRLs_using_all_sigDHSs.pdf")
plot(density(dhs_df_sig$crl_wFc), col="black", xlim=c(-0.7, 1.5), main="CRL using all sig DHS")
lines(density(dhs_df_sig$crl_wFd), col="red", lwd=2)
legend (0.5,4, c("CRL_wFc", "CRL_wFd"), lty=c(1,1), lwd=c(2,2), col=c("black", "red")) 
dev.off()

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
dhs_df_bounded = dhs_df_sig$crl_wFd <= wellReprog_upperLim & dhs_df_sig$crl_wFd >= badlyReprog_lowerLim # bounded by upper and lower lim
mbOpen_df_bounded = dhs_df_sig[dhs_df_bounded & dhs_df_sig$mbStatus == "mbOpen",] # 7489 by 25 if w/o Fdh and Fdl
mbOpen_df_bounded = mbOpen_df_bounded[order(mbOpen_df_bounded$crl_wFd, decreasing = TRUE),] # 7123 by 37 if w/ Fdh and Fdl

mbClosed_df_bounded = dhs_df_sig[dhs_df_bounded & dhs_df_sig$mbStatus == "mbClosed",] # 4782 by 25 if w/o Fdh and Fdl
mbClosed_df_bounded = mbClosed_df_bounded[order(mbClosed_df_bounded$crl_wFd, decreasing = TRUE),] # 6836 by 37 if w/ Fdh and Fdl

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
# === NOW SAVING DATA AND WRITING FILES ====
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#
# ============================================
save.image("allResults_generated.RData")
write.table(dhs_df_sig, file="dhs_df_sig.txt", quote=F, sep="\t")
write.table(mbOpen_df_bounded, file="mbOpen_df_bounded.txt", quote=F, sep="\t")
write.table(mbClosed_df_bounded, file="mbClosed_df_bounded.txt", quote=F, sep="\t")
