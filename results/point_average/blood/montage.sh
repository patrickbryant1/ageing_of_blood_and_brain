#Figure 1
#Add labels
convert  ag_distribution.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  ag_distribution.png
# convert  ag_point_cutoffs.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  ag_point_cutoffs.png
montage ag_distribution.png ag_point_cutoffs.png -tile 1x2 -geometry +2+2 -title 'Blood' -pointsize 30 Figure1.png
#
#
# #Figure 2
# #Add labels
# convert  pval_distribution.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" Figure2.png
#
#
# #Figure 3
# convert  pos_ra.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" pos_ra.png
# convert  neg_ra.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" neg_ra.png
# convert  gradients.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C"  gradients.png
# montage pos_ra.png neg_ra.png gradients.png -tile 3x1 -geometry +2+2 -title 'Blood' -pointsize 30 Figure3.png


#Figure 7
convert  FC_correlations.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  FC_correlations.png
convert  abs_correlations.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  abs_correlations.png
montage FC_correlations.png abs_correlations.png -tile 1x2 -geometry +2+2 -title 'Blood' -pointsize 30  Figure7.png

# #Figure 8
convert  FC_hannum_markers.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  FC_hannum_markers.png
convert  abs_hannum_markers.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  abs_hannum_markers.png
montage FC_hannum_markers.png abs_hannum_markers.png -tile 1x2 -geometry +2+2 -title 'Blood' -pointsize 30  Figure8.png
#Figure 9
montage ./models/FC_cv_results.png ./models/abs_cv_results.png ./models/overlap_cv_results.png -tile 3x1 -geometry +2+2 -title 'Blood' -pointsize 30 CV_blood.png
# #Figure 5
# convert  ./genes/go_ra.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  ./genes/go_ra.png
# convert  ./genes/Hannum/go_hannum.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./genes/Hannum/go_hannum.png
# montage ./genes/go_ra.png ./genes/Hannum/go_hannum.png -tile 1x2 -geometry +2+2 Figure5.png
