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
# BASE=fold_changes/genes
# montage $BASE/ARHGAP22.png  $BASE/LOH12CR2.png $BASE/ARHGEF4.png  $BASE/ARHGAP22_vs_age.png  $BASE/LOH12CR2_vs_age.png $BASE/ARHGEF4_vs_age.png  $BASE/LRRK1.png $BASE/ARSG.png  $BASE/MAGI2.png $BASE/LRRK1_vs_age.png $BASE/ARSG_vs_age.png  $BASE/MAGI2_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes1.png
# montage $BASE/B3GALT6.png   $BASE/MOSC2.png $BASE/BMP8A.png  $BASE/B3GALT6_vs_age.png   $BASE/MOSC2_vs_age.png $BASE/BMP8A_vs_age.png  $BASE/MYO10.png $BASE/C12orf65.png  $BASE/PRDXDD1P.png $BASE/MYO10_vs_age.png $BASE/C12orf65_vs_age.png  $BASE/PRDXDD1P_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes2.png
# montage $BASE/CACNA1A.png $BASE/PRKG1.png $BASE/DOCK1.png $BASE/CACNA1A_vs_age.png $BASE/PRKG1_vs_age.png $BASE/DOCK1_vs_age.png $BASE/RABGGTB.png $BASE/DUSP5.png $BASE/SDF4.png $BASE/RABGGTB_vs_age.png $BASE/DUSP5_vs_age.png $BASE/SDF4_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes3.png
# montage $BASE/EPB41L4A.png  $BASE/SFRS6.png  $BASE/FLJ32810.png $BASE/EPB41L4A_vs_age.png  $BASE/SFRS6_vs_age.png  $BASE/FLJ32810_vs_age.png $BASE/SHANK1.png $BASE/GGNBP2.png    $BASE/SHANK2.png $BASE/SHANK1_vs_age.png $BASE/GGNBP2_vs_age.png    $BASE/SHANK2_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes4.png
# montage $BASE/HLA-B.png $BASE/SLC16A6.png $BASE/JAKMIP3.png $BASE/HLA-B_vs_age.png $BASE/SLC16A6_vs_age.png $BASE/JAKMIP3_vs_age.png $BASE/SLC7A5.png $BASE/KIAA1409.png  $BASE/SNORD45C.png $BASE/SLC7A5_vs_age.png $BASE/KIAA1409_vs_age.png  $BASE/SNORD45C_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes5.png
# montage $BASE/LEPR.png $BASE/TMEM181.png $BASE/LHX3.png $BASE/LEPR_vs_age.png $BASE/TMEM181_vs_age.png $BASE/LHX3_vs_age.png $BASE/TP73.png $BASE/LIMCH1.png $BASE/TYMP.png  $BASE/TP73_vs_age.png $BASE/LIMCH1_vs_age.png $BASE/TYMP_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes6.png
# montage $BASE/LMTK3.png $BASE/VWC2.png $BASE/LOH12CR1.png $BASE/LMTK3_vs_age.png $BASE/VWC2_vs_age.png $BASE/LOH12CR1_vs_age.png  $BASE/ZMAT4.png $BASE/ZMAT4_vs_age.png -tile 3x4 -geometry +2+2 $BASE/genes7.png
#

# #Montage beta value changes for markers
#montage fold_changes/gap10_even.png fold_changes/gap10_uneven.png -tile 2x1 -geometry +2+2 fold_changes/all_gap10.png

#montage fold_changes/markers/cg02725290.png fold_changes/markers/cg14074924.png fold_changes/markers/cg16651196.png fold_changes/markers/cg02725290_vs_age.png fold_changes/markers/cg14074924_vs_age.png fold_changes/markers/cg16651196_vs_age.png -tile 3x2 -geometry +2+2 fold_changes/markers/all_makers.png

#Montage probes not picked up by correlation
# BASE=/home/pbryant/ageing_of_blood/results/fold_changes/diff
# montage $BASE/cg07527631_vs_age.png  $BASE/cg15553911_vs_age.png $BASE/cg00225623_vs_age.png  $BASE/cg07624948_vs_age.png  $BASE/cg16595484_vs_age.png $BASE/cg00287122_vs_age.png  $BASE/cg08346664_vs_age.png  $BASE/cg17105886_vs_age.png $BASE/cg00532901_vs_age.png  $BASE/cg08657228_vs_age.png  $BASE/cg17728503_vs_age.png $BASE/cg00625131_vs_age.png  $BASE/cg08858875_vs_age.png  $BASE/cg17741712_vs_age.png $BASE/cg00786635_vs_age.png  $BASE/cg09620491_vs_age.png  $BASE/cg17765304_vs_age.png $BASE/cg01251835_vs_age.png  $BASE/cg09821905_vs_age.png  $BASE/cg18314882_vs_age.png $BASE/cg01270001_vs_age.png  $BASE/cg10016610_vs_age.png  $BASE/cg18455249_vs_age.png $BASE/cg01637175_vs_age.png  $BASE/cg10684905_vs_age.png  $BASE/cg18606723_vs_age.png $BASE/cg02725290_vs_age.png  $BASE/cg11701583_vs_age.png  $BASE/cg18857216_vs_age.png $BASE/cg02943497_vs_age.png  $BASE/cg11808100_vs_age.png  $BASE/cg19627864_vs_age.png $BASE/cg03031660_vs_age.png  $BASE/cg11905617_vs_age.png  $BASE/cg19726630_vs_age.png $BASE/cg03249590_vs_age.png  $BASE/cg12014753_vs_age.png  $BASE/cg20540608_vs_age.png $BASE/cg03919694_vs_age.png  $BASE/cg12030941_vs_age.png  -tile 5x8 -geometry +2+2 $BASE/all1.png
#
# montage $BASE/cg20894495_vs_age.png $BASE/cg04315227_vs_age.png  $BASE/cg13021192_vs_age.png  $BASE/cg21012737_vs_age.png $BASE/cg04633600_vs_age.png  $BASE/cg13275176_vs_age.png  $BASE/cg21245975_vs_age.png $BASE/cg05268155_vs_age.png  $BASE/cg13332925_vs_age.png  $BASE/cg22133704_vs_age.png $BASE/cg05273049_vs_age.png  $BASE/cg13353717_vs_age.png  $BASE/cg23109721_vs_age.png $BASE/cg05979320_vs_age.png  $BASE/cg13372488_vs_age.png  $BASE/cg23243080_vs_age.png $BASE/cg06069187_vs_age.png  $BASE/cg13529619_vs_age.png  $BASE/cg23484755_vs_age.png $BASE/cg06103394_vs_age.png  $BASE/cg13682223_vs_age.png  $BASE/cg24613083_vs_age.png $BASE/cg06343869_vs_age.png  $BASE/cg14074924_vs_age.png  $BASE/cg25307371_vs_age.png $BASE/cg06960514_vs_age.png  $BASE/cg14300823_vs_age.png  $BASE/cg26195366_vs_age.png $BASE/cg07089783_vs_age.png  $BASE/cg14390683_vs_age.png  $BASE/cg26466970_vs_age.png $BASE/cg07193998_vs_age.png  $BASE/cg14464244_vs_age.png  $BASE/cg26657240_vs_age.png $BASE/cg07518837_vs_age.png  $BASE/cg14795253_vs_age.png  $BASE/cg27281836_vs_age.png -tile 5x8 -geometry +2+2 $BASE/all2.png

#Montage top10 correlations
#montage /home/pbryant/ageing_of_blood/results/fold_changes/top10/*.png -tile 5x2 -geometry +2+2 /home/pbryant/ageing_of_blood/results/fold_changes/top10/all.png

#Montage overlapping probes with 10 year gaps
montage /home/pbryant/ageing_of_blood/results/fold_changes/overlap_10_year_gaps/*.png -tile 4x2 -geometry +2+2 /home/pbryant/ageing_of_blood/results/fold_changes/overlap_10_year_gaps/all.png
