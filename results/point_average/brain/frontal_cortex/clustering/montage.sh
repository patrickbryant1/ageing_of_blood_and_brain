#!/usr/bin/env
convert  1.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  1.png
convert  tsne.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" tsne.png
convert gradients.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" gradients.png
montage 1.png tsne.png gradients.png -tile 3x1 -geometry +2+2 Figure4.png
#-title 'Frontal cortex' -pointsize 30
