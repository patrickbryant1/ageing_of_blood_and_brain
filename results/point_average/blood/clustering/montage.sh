#!/usr/bin/env
#Add labels
convert FC1.png -pointsize 50 -gravity NorthWest -annotate +0+0 "A" FC1.png
convert FC2.png -pointsize 50 -gravity NorthWest -annotate +0+0 "B" FC2.png
convert FC_tsne.png -pointsize 50 -gravity NorthWest -annotate +0+0 "C" FC_tsne.png
convert FC_gradients.png -pointsize 50 -gravity NorthWest -annotate +0+0 "D" FC_gradients.png
montage FC1.png FC2.png FC_tsne.png FC_gradients.png -title 'FC' -pointsize 30 -tile 3x2 -geometry +2+2 FC_montage.png
#Add labels
convert abs1.png -pointsize 50 -gravity NorthWest -annotate +0+0 "E" abs1.png
convert abs_tsne.png -pointsize 50 -gravity NorthWest -annotate +0+0 "F" abs_tsne.png
convert abs_gradients.png -pointsize 50 -gravity NorthWest -annotate +0+0 "G" abs_gradients.png
montage abs1.png abs_tsne.png abs_gradients.png -title 'Absolute value' -pointsize 30 -tile 3x1 -geometry +2+2 abs_montage.png
#Add labels
convert overlap1.png -pointsize 50 -gravity NorthWest -annotate +0+0 "H" overlap1.png
convert overlap_tsne.png -pointsize 50 -gravity NorthWest -annotate +0+0 "I" overlap_tsne.png
convert overlap_gradients.png -pointsize 50 -gravity NorthWest -annotate +0+0 "J" overlap_gradients.png
montage overlap1.png overlap_tsne.png overlap_gradients.png -title 'Overlap' -pointsize 30 -tile 3x1 -geometry +2+2 overlap_montage.png
#All
montage FC_montage.png abs_montage.png overlap_montage.png -tile 1x3 -geometry +2+2 all_montage.png
