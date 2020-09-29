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
BASE=fold_changes/genes
montage $BASE/ARHGAP22.png  $BASE/LOH12CR2.png $BASE/ARHGEF4.png  $BASE/ARHGAP22_vs_age.png  $BASE/LOH12CR2_vs_age.png $BASE/ARHGEF4_vs_age.png  $BASE/LRRK1.png $BASE/ARSG.png  $BASE/MAGI2.png $BASE/LRRK1_vs_age.png $BASE/ARSG_vs_age.png  $BASE/MAGI2_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes1.png
montage $BASE/B3GALT6.png   $BASE/MOSC2.png $BASE/BMP8A.png  $BASE/B3GALT6_vs_age.png   $BASE/MOSC2_vs_age.png $BASE/BMP8A_vs_age.png  $BASE/MYO10.png $BASE/C12orf65.png  $BASE/PRDXDD1P.png $BASE/MYO10_vs_age.png $BASE/C12orf65_vs_age.png  $BASE/PRDXDD1P_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes2.png
montage $BASE/CACNA1A.png $BASE/PRKG1.png $BASE/DOCK1.png $BASE/CACNA1A_vs_age.png $BASE/PRKG1_vs_age.png $BASE/DOCK1_vs_age.png $BASE/RABGGTB.png $BASE/DUSP5.png $BASE/SDF4.png $BASE/RABGGTB_vs_age.png $BASE/DUSP5_vs_age.png $BASE/SDF4_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes3.png
montage $BASE/EPB41L4A.png  $BASE/SFRS6.png  $BASE/FLJ32810.png $BASE/EPB41L4A_vs_age.png  $BASE/SFRS6_vs_age.png  $BASE/FLJ32810_vs_age.png $BASE/SHANK1.png $BASE/GGNBP2.png    $BASE/SHANK2.png $BASE/SHANK1_vs_age.png $BASE/GGNBP2_vs_age.png    $BASE/SHANK2_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes4.png
montage $BASE/HLA-B.png $BASE/SLC16A6.png $BASE/JAKMIP3.png $BASE/HLA-B_vs_age.png $BASE/SLC16A6_vs_age.png $BASE/JAKMIP3_vs_age.png $BASE/SLC7A5.png $BASE/KIAA1409.png  $BASE/SNORD45C.png $BASE/SLC7A5_vs_age.png $BASE/KIAA1409_vs_age.png  $BASE/SNORD45C_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes5.png
montage $BASE/LEPR.png $BASE/TMEM181.png $BASE/LHX3.png $BASE/LEPR_vs_age.png $BASE/TMEM181_vs_age.png $BASE/LHX3_vs_age.png $BASE/TP73.png $BASE/LIMCH1.png $BASE/TYMP.png  $BASE/TP73_vs_age.png $BASE/LIMCH1_vs_age.png $BASE/TYMP_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes6.png
montage $BASE/LMTK3.png $BASE/VWC2.png $BASE/LOH12CR1.png $BASE/LMTK3_vs_age.png $BASE/VWC2_vs_age.png $BASE/LOH12CR1_vs_age.png  $BASE/ZMAT4.png $BASE/ZMAT4_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes7.png
#

# #Montage beta value changes for markers
#montage fold_changes/gap10_even.png fold_changes/gap10_uneven.png -tile 2x1 -geometry +2+2 fold_changes/all_gap10.png

#montage fold_changes/markers/cg02725290.png fold_changes/markers/cg14074924.png fold_changes/markers/cg16651196.png fold_changes/markers/cg02725290_vs_age.png fold_changes/markers/cg14074924_vs_age.png fold_changes/markers/cg16651196_vs_age.png -tile 3x2 -geometry +2+2 fold_changes/markers/all_makers.png
