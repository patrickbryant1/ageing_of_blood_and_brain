#Montage in a 7x7 grid having the pvalue distributions on the uppertriangular and volcano on lower
montage whitebox.png pval_0_1.png pval_0_2.png pval_0_3.png pval_0_4.png pval_0_5.png pval_0_6.png -tile 7x1 -geometry +2+2 row1.png
montage volcano_0_1.png whitebox.png pval_1_2.png pval_1_3.png pval_1_4.png pval_1_5.png pval_1_6.png -tile 7x1 -geometry +2+2 row2.png
montage volcano_0_2.png volcano_1_2.png whitebox.png pval_2_3.png pval_2_4.png pval_2_5.png pval_2_6.png -tile 7x1 -geometry +2+2 row3.png
montage volcano_0_3.png volcano_1_3.png volcano_2_3.png whitebox.png pval_3_4.png pval_3_5.png pval_3_6.png -tile 7x1 -geometry +2+2 row4.png
montage volcano_0_4.png volcano_1_4.png volcano_2_4.png volcano_3_4.png whitebox.png pval_4_5.png pval_4_6.png -tile 7x1 -geometry +2+2 row5.png
montage volcano_0_5.png volcano_1_5.png volcano_2_5.png volcano_3_5.png volcano_4_5.png whitebox.png pval_5_6.png -tile 7x1 -geometry +2+2 row6.png
montage volcano_0_6.png volcano_1_6.png volcano_2_6.png volcano_3_6.png volcano_4_6.png volcano_5_6.png whitebox.png -tile 7x1 -geometry +2+2 row7.png
#All rows together
montage row1.png row2.png row3.png row4.png row5.png row6.png row7.png -tile 1x7 -geometry +2+2 matrix.png


#qval volcano
#Montage in a 7x7 grid having the qvalue distributions on the uppertriangular and qval_volcano on lower
montage whitebox_45.png qval_0_1.png qval_0_2.png qval_0_3.png qval_0_4.png qval_0_5.png qval_0_6.png -tile 7x1 -geometry +2+2 row1.png
montage qval_volcano_0_1.png whitebox_45.png qval_1_2.png qval_1_3.png qval_1_4.png qval_1_5.png qval_1_6.png -tile 7x1 -geometry +2+2 row2.png
montage qval_volcano_0_2.png qval_volcano_1_2.png whitebox_45.png qval_2_3.png qval_2_4.png qval_2_5.png qval_2_6.png -tile 7x1 -geometry +2+2 row3.png
montage qval_volcano_0_3.png qval_volcano_1_3.png qval_volcano_2_3.png whitebox_45.png qval_3_4.png qval_3_5.png qval_3_6.png -tile 7x1 -geometry +2+2 row4.png
montage qval_volcano_0_4.png qval_volcano_1_4.png qval_volcano_2_4.png qval_volcano_3_4.png whitebox_45.png qval_4_5.png qval_4_6.png -tile 7x1 -geometry +2+2 row5.png
montage qval_volcano_0_5.png qval_volcano_1_5.png qval_volcano_2_5.png qval_volcano_3_5.png qval_volcano_4_5.png whitebox_45.png qval_5_6.png -tile 7x1 -geometry +2+2 row6.png
montage qval_volcano_0_6.png qval_volcano_1_6.png qval_volcano_2_6.png qval_volcano_3_6.png qval_volcano_4_6.png qval_volcano_5_6.png whitebox_45.png -tile 7x1 -geometry +2+2 row7.png
#All rows together
montage row1.png row2.png row3.png row4.png row5.png row6.png row7.png -tile 1x7 -geometry +2+2 qval_matrix.png
