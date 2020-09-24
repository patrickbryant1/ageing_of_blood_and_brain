# #Montage in a 7x7 grid having the pvalue distributions on the uppertriangular and volcano on lower
# montage whitebox.png pval_0_1.png pval_0_2.png pval_0_3.png pval_0_4.png pval_0_5.png pval_0_6.png -tile 7x1 -geometry +2+2 row1.png
# montage volcano_0_1.png whitebox.png pval_1_2.png pval_1_3.png pval_1_4.png pval_1_5.png pval_1_6.png -tile 7x1 -geometry +2+2 row2.png
# montage volcano_0_2.png volcano_1_2.png whitebox.png pval_2_3.png pval_2_4.png pval_2_5.png pval_2_6.png -tile 7x1 -geometry +2+2 row3.png
# montage volcano_0_3.png volcano_1_3.png volcano_2_3.png whitebox.png pval_3_4.png pval_3_5.png pval_3_6.png -tile 7x1 -geometry +2+2 row4.png
# montage volcano_0_4.png volcano_1_4.png volcano_2_4.png volcano_3_4.png whitebox.png pval_4_5.png pval_4_6.png -tile 7x1 -geometry +2+2 row5.png
# montage volcano_0_5.png volcano_1_5.png volcano_2_5.png volcano_3_5.png volcano_4_5.png whitebox.png pval_5_6.png -tile 7x1 -geometry +2+2 row6.png
# montage volcano_0_6.png volcano_1_6.png volcano_2_6.png volcano_3_6.png volcano_4_6.png volcano_5_6.png whitebox.png -tile 7x1 -geometry +2+2 row7.png
# #All rows together
# montage row1.png row2.png row3.png row4.png row5.png row6.png row7.png -tile 1x7 -geometry +2+2 matrix.png
#
#
# #qval volcano
# #Montage in a 7x7 grid having the qvalue distributions on the uppertriangular and qval_volcano on lower
# montage whitebox_45.png qval_0_1.png qval_0_2.png qval_0_3.png qval_0_4.png qval_0_5.png qval_0_6.png -tile 7x1 -geometry +2+2 row1.png
# montage qval_volcano_0_1.png whitebox_45.png qval_1_2.png qval_1_3.png qval_1_4.png qval_1_5.png qval_1_6.png -tile 7x1 -geometry +2+2 row2.png
# montage qval_volcano_0_2.png qval_volcano_1_2.png whitebox_45.png qval_2_3.png qval_2_4.png qval_2_5.png qval_2_6.png -tile 7x1 -geometry +2+2 row3.png
# montage qval_volcano_0_3.png qval_volcano_1_3.png qval_volcano_2_3.png whitebox_45.png qval_3_4.png qval_3_5.png qval_3_6.png -tile 7x1 -geometry +2+2 row4.png
# montage qval_volcano_0_4.png qval_volcano_1_4.png qval_volcano_2_4.png qval_volcano_3_4.png whitebox_45.png qval_4_5.png qval_4_6.png -tile 7x1 -geometry +2+2 row5.png
# montage qval_volcano_0_5.png qval_volcano_1_5.png qval_volcano_2_5.png qval_volcano_3_5.png qval_volcano_4_5.png whitebox_45.png qval_5_6.png -tile 7x1 -geometry +2+2 row6.png
# montage qval_volcano_0_6.png qval_volcano_1_6.png qval_volcano_2_6.png qval_volcano_3_6.png qval_volcano_4_6.png qval_volcano_5_6.png whitebox_45.png -tile 7x1 -geometry +2+2 row7.png
# #All rows together
# montage row1.png row2.png row3.png row4.png row5.png row6.png row7.png -tile 1x7 -geometry +2+2 qval_matrix.png


# #Montage fold changes for genes
montage fold_changes/genes/ARHGAP22.png  fold_changes/genes/LOH12CR2.png fold_changes/genes/ARHGEF4.png   fold_changes/genes/LRRK1.png fold_changes/genes/ARSG.png      fold_changes/genes/MAGI2.png fold_changes/genes/B3GALT6.png   fold_changes/genes/MOSC2.png fold_changes/genes/BMP8A.png     fold_changes/genes/MYO10.png fold_changes/genes/C12orf65.png  fold_changes/genes/PRDXDD1P.png fold_changes/genes/CACNA1A.png   fold_changes/genes/PRKG1.png fold_changes/genes/DOCK1.png     fold_changes/genes/RABGGTB.png fold_changes/genes/DUSP5.png     fold_changes/genes/SDF4.png fold_changes/genes/EPB41L4A.png  fold_changes/genes/SFRS6.png  -tile 4x5 -geometry +2+2 fold_changes/genes/all_genes1.png

montage fold_changes/genes/FLJ32810.png  fold_changes/genes/SHANK1.png fold_changes/genes/GGNBP2.png    fold_changes/genes/SHANK2.png fold_changes/genes/HLA-B.png     fold_changes/genes/SLC16A6.png fold_changes/genes/JAKMIP3.png   fold_changes/genes/SLC7A5.png fold_changes/genes/KIAA1409.png  fold_changes/genes/SNORD45C.png fold_changes/genes/LEPR.png      fold_changes/genes/TMEM181.png fold_changes/genes/LHX3.png      fold_changes/genes/TP73.png fold_changes/genes/LIMCH1.png    fold_changes/genes/TYMP.png fold_changes/genes/LMTK3.png     fold_changes/genes/VWC2.png fold_changes/genes/LOH12CR1.png  fold_changes/genes/ZMAT4.png -tile 4x5 -geometry +2+2 fold_changes/genes/all_genes2.png
#
# #Montage beta value changes for markers
#montage fold_changes/gap10_even.png fold_changes/gap10_uneven.png -tile 2x1 -geometry +2+2 fold_changes/all_gap10.png

montage fold_changes/markers/cg02725290.png fold_changes/markers/cg14074924.png fold_changes/markers/cg16651196.png fold_changes/markers/cg02725290_vs_age.png fold_changes/markers/cg14074924_vs_age.png fold_changes/markers/cg16651196_vs_age.png -tile 3x2 -geometry +2+2 fold_changes/markers/all_makers.png
