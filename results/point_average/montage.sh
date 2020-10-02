
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
