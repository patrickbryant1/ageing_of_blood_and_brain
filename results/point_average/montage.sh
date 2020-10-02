#Figure 1
#Add labels
convert  ag_distribution.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  ag_distribution.png
convert  ag_point_cutoffs.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  ag_point_cutoffs.png
montage ag_distribution.png ag_point_cutoffs.png -tile 2x1 -geometry +2+2 Figure1.png


#Figure 2
#Add labels
convert  grad_diff_distribution.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" grad_diff_distribution.png
convert  grad_diff_vs_FC.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  grad_diff_vs_FC.png
convert  gradients.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C"  gradients.png
montage grad_diff_distribution.png grad_diff_vs_FC.png gradients.png -tile 3x1 -geometry +2+2 Figure2.png


#Figure 3
convert  pos_ra.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" pos_ra.png
convert  neg_ra.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" neg_ra.png
montage pos_ra.png neg_ra.png -tile 2x1 -geometry +2+2 Figure3.png

#Figure 4
montage AFAP1.png     FBXL16.png   KLF14.png         NPY.png      PCDHA3.png   PEX5L.png   TRIM58.png ARHGAP22.png  FBXL7.png    LOC100271715.png  NRIP3.png    PCDHA4.png   PODXL2.png  TUBGCP5.png ASXL3.png     FILIP1L.png  LOC375196.png     PCDHA10.png  PCDHA5.png   PRKAA2.png  WDR17.png C3orf26.png  FOXG1.png    MARCH4.png        PCDHA11.png  PCDHA6.png  RNF32.png   ZNF518B.png -tile 4x7 -geometry +2+2 all1.png
montage    C7orf13.png   GSX1.png     MIR2277.png       PCDHA12.png  PCDHA7.png   SDHAP3.png CAT.png       INA.png   MIR548G.png       PCDHA13.png  PCDHA8.png   SGPP2.png COBL.png      IRF8.png     MOSC2.png         PCDHA1.png   PCDHA9.png   TBX20.png FAM172A.png   ISOC2.png    MYRIP.png         PCDHA2.png   PCDHAC1.png  TOX2.png -tile 4x6 -geometry all2.png
