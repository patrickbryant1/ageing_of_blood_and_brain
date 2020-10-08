#Figure 1
#Add labels
convert  ag_distribution.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  ag_distribution.png
convert  ag_point_cutoffs.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  ag_point_cutoffs.png
montage ag_distribution.png ag_point_cutoffs.png -tile 2x1 -geometry +2+2 -title 'Blood' -pointsize 30 Figure1.png


#Figure 2
#Add labels
convert  pval_distribution.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" pval_distribution.png
convert  grad_diff_distribution.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" grad_diff_distribution.png
convert  grad_diff_vs_FC.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C"  grad_diff_vs_FC.png
montage pval_distribution.png grad_diff_distribution.png grad_diff_vs_FC.png -tile 3x1 -geometry +2+2 -title 'Blood' -pointsize 30  Figure2.png


#Figure 3
convert  pos_ra.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" pos_ra.png
convert  neg_ra.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" neg_ra.png
convert  gradients.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C"  gradients.png
montage pos_ra.png neg_ra.png gradients.png -tile 3x1 -geometry +2+2 -title 'Blood' -pointsize 30 Figure3.png

#Figure 4
convert  hannum_markers.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  hannum_markers.png
convert  correlations.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" correlations.png
montage hannum_markers.png correlations.png -tile 2x1 -geometry +2+2 Figure4.png

#Figure 5
convert  ./genes/go_ra.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  ./genes/go_ra.png
convert  ./genes/Hannum/go_hannum.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./genes/Hannum/go_hannum.png
montage ./genes/go_ra.png ./genes/Hannum/go_hannum.png -tile 1x2 -geometry +2+2 Figure5.png
