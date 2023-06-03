
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

############################################
#
# Write data tables from DADA2 results
#
############################################

#---------------------------------------------
# Description:
# This script prepares DADA2 results to work with R
#  - Separates ASVs sequences and count reads data.
#  - Creates ASVs names.
#  - Remove data of not Prokaryotic species: Eukaryota, Chloroplast and 
#    Mitochondria.
#---------------------------------------------


############################################
# Prepare DADA2 results 
############################################

# load DADA2 results
load("raw_data/DADA2_CR30.No_primers.RData")

#  count data matrix
asvs_raw <- as.data.frame(seqtab.nochim)

# separate sequences from count data
asvs_seqs <- data.frame(seq = names(asvs_raw))

# create ASVs names
names(asvs_raw) <- paste("asv",sprintf("%04d", 1:ncol(asvs_raw)), sep="")
asvs_seqs$asv <-names(asvs_raw)
asvs_seqs <- asvs_seqs[,c("asv","seq")]

# sample names to a column
asvs_raw$sample <- rownames(asvs_raw)
asvs_raw <- asvs_raw[, c(ncol(asvs_raw), 2:ncol(asvs_raw)-1)]

# Load taxonomy
taxonomy <- as.data.frame(taxa.sp.print)
taxonomy$asv <- names(asvs_raw)[2:ncol(asvs_raw)]


##########################################################
# Clean data of Eukaryotes, Chloroplasts and Mitochondria
##########################################################

# Check where they are in the taxonomy
lapply(taxonomy, function(x){grep("Eukaryota", x)})
# 
# $Kingdom
# [1] 5210 6520 6707 7670 7829
# 
# $Phylum
# integer(0)
# 
# $Class
# integer(0)
# 
# $Order
# integer(0)
# 
# $Family
# integer(0)
# 
# $Genus
# integer(0)
# 
# $Species
# integer(0)
# 
# $asv
# integer(0)


lapply(taxonomy, function(x){grep("Chloroplast", x)})
# $Kingdom
# integer(0)
# 
# $Phylum
# integer(0)
# 
# $Class
# integer(0)
# 
# $Order
# [1]    4    7   20   43   44   59   80   84   85  137  138  203  206  207  216  238  282  316  331  338  354  381  384  392  401  412  417  421  424  443  450  468  469  485  491  525  547  557  565  567
# [41]  580  591  599  603  634  642  643  663  668  689  717  723  727  730  735  748  754  811  849  867  882  884  896  930  936  940  951  957  979  996  998 1021 1022 1049 1062 1078 1082 1124 1154 1164
# [81] 1203 1211 1237 1266 1277 1312 1317 1318 1320 1347 1350 1354 1362 1375 1377 1396 1420 1424 1428 1447 1455 1457 1467 1479 1498 1513 1514 1519 1534 1575 1595 1629 1630 1647 1660 1673 1700 1708 1770 1810
# [121] 1814 1815 1816 1854 1875 1876 1883 1894 1900 1914 1941 1942 1948 1949 2044 2064 2070 2077 2089 2106 2116 2118 2123 2147 2208 2213 2225 2238 2239 2333 2351 2371 2376 2377 2383 2393 2395 2413 2446 2453
# [161] 2504 2505 2523 2534 2537 2568 2584 2618 2628 2647 2657 2682 2712 2728 2731 2765 2771 2773 2779 2814 2859 2864 2922 2991 3002 3016 3018 3049 3050 3066 3076 3079 3100 3130 3134 3138 3144 3145 3182 3183
# [201] 3184 3188 3189 3199 3200 3211 3238 3241 3254 3294 3301 3303 3305 3310 3353 3358 3364 3365 3382 3402 3406 3414 3439 3459 3480 3484 3490 3523 3524 3544 3545 3595 3612 3630 3648 3651 3652 3657 3670 3709
# [241] 3713 3719 3720 3736 3755 3791 3860 3889 3902 3903 3917 3918 3932 3952 3957 4021 4022 4023 4029 4055 4069 4100 4106 4115 4119 4120 4142 4170 4174 4187 4216 4218 4219 4220 4242 4246 4250 4259 4271 4272
# [281] 4309 4313 4315 4324 4326 4338 4339 4352 4370 4375 4380 4382 4393 4443 4463 4465 4468 4482 4490 4509 4519 4534 4535 4549 4550 4555 4583 4616 4622 4644 4646 4649 4650 4661 4690 4711 4721 4722 4732 4735
# [321] 4737 4752 4753 4761 4767 4788 4797 4798 4802 4805 4814 4815 4839 4841 4842 4845 4848 4876 4878 4881 4886 4887 4897 4934 4939 4950 4998 5007 5009 5013 5016 5077 5079 5082 5083 5102 5107 5132 5140 5154
# [361] 5155 5176 5178 5202 5205 5208 5212 5219 5220 5234 5235 5286 5289 5293 5302 5312 5325 5354 5356 5364 5384 5385 5396 5430 5431 5436 5448 5460 5461 5463 5467 5526 5531 5532 5533 5537 5580 5629 5663 5681
# [401] 5689 5694 5696 5700 5711 5740 5765 5767 5769 5771 5776 5834 5877 5879 5919 5920 5922 5928 5963 5965 6010 6011 6014 6016 6057 6069 6082 6093 6138 6167 6168 6170 6174 6182 6223 6235 6241 6245 6254 6260
# [441] 6274 6275 6276 6283 6329 6357 6367 6381 6387 6432 6440 6460 6463 6526 6545 6611 6614 6620 6636 6640 6656 6711 6712 6809 6820 6840 6842 6845 6849 6850 6852 6854 6863 6915 6924 6927 6942 6945 6947 6948
# [481] 6961 6968 6983 7025 7039 7049 7054 7072 7100 7175 7176 7184 7185 7209 7244 7286 7303 7304 7309 7312 7313 7329 7340 7343 7371 7424 7431 7446 7447 7450 7474 7497 7507 7509 7510 7534 7621 7633 7634 7640
# [521] 7649 7662 7730 7749 7770 7780 7815 7831 7843 7856 7866 7867 7892 7956 7973 7979 8170
# 
# $Family
# integer(0)
# 
# $Genus
# integer(0)
# 
# $Species
# integer(0)
# 
# $asv
# integer(0)



lapply(taxonomy, function(x){grep("Mitochondria", x)})
# $Kingdom
# integer(0)
# 
# $Phylum
# integer(0)
# 
# $Class
# integer(0)
# 
# $Order
# integer(0)
# 
# $Family
# [1]  700  742  860  926  992 2215 2469 2493 2548 2601 2671 3077 3455 3542 3636 3951 4114 4307 4931 4966 5152 5173 5487 5545 5658 5697 5759 5831 5866 6005 6056 6385 6513 6551 6552 6605 6606 6634 6670 6730
# [41] 6748 6909 6922 6953 7016 7028 7051 7087 7149 7230 7239 7302 7318 7346 7360 7452 7453 7496 7627 7630 7651 7667 7706 7787 7811 7818 7819 7833 7835 7886 7902 7949 7957
# 
# $Genus
# integer(0)
# 
# $Species
# integer(0)
# 
# $asv
# integer(0)


# Count how many not prokariotic ASVs there are:
length(grep("Eukaryota", taxonomy$Kingdom))
#[1] 5
length(grep("Chloroplast", taxonomy$Order))
#[1] 537
length(grep("Mitochondria", taxonomy$Family))
#[1] 73

# Get ASVs to remove
toremove <- subset(taxonomy, Kingdom == "Eukaryota" | Order == "Chloroplast" | Family == "Mitochondria")$asv

# filter all tables
asvs_raw <- subset(asvs_raw, select = names(asvs_raw)[!(names(asvs_raw) %in% toremove)])
asvs_seqs <- subset(asvs_seqs, !(asv %in% toremove))
taxonomy <- subset(taxonomy, !(asv %in% toremove))

# save to files
# counts
write.table(asvs_raw, file = "raw_data/asvs_raw_counts.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) 
# sequences
write.table(asvs_seqs, file = "raw_data/asvs_seqs.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# taxonomy
write.table(taxonomy, file = "raw_data/asvs_taxonomy.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
