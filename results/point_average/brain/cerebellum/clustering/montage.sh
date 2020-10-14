#!/usr/bin/env
convert  1.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  1.png
convert  2.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  2.png
convert  3.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C"  3.png
convert  4.png -pointsize 60 -gravity NorthWest -annotate +0+0 "D"  4.png
convert  tsne.png -pointsize 60 -gravity NorthWest -annotate +0+0 "E" tsne.png
convert gradients.png -pointsize 60 -gravity NorthWest -annotate +0+0 "F" gradients.png
montage 1.png 2.png 3.png 4.png tsne.png gradients.png -tile 3x2 -geometry +2+2  Figure5.png
#-title 'Cerebellum' -pointsize 30


convert  1_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  1_unnormalized.png
convert  2_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  2_unnormalized.png
convert  3_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C"  3_unnormalized.png
convert  4_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "D"  4_unnormalized.png
montage 1_unnormalized.png 2_unnormalized.png 3_unnormalized.png 4_unnormalized.png -tile 2x2 -geometry +2+2  FigureS4.png
