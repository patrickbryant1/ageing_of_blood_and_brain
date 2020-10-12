#!/usr/bin/env
#Add labels
convert  1.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  1.png
convert  2.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  2.png
convert  tsne.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" tsne.png
convert gradients.png -pointsize 60 -gravity NorthWest -annotate +0+0 "D" gradients.png
montage 1.png 2.png tsne.png gradients.png -tile 2x2 -geometry +2+2 -title 'Blood' -pointsize 30 Figure3.png
