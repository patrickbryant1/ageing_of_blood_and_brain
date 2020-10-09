#Figure 1
#Add labels
montage ./blood/Figure1.png ./brain/frontal_cortex/Figure1.png ./brain/cerebellum/Figure1.png -tile 3x1 -geometry +2+2 Figure1.png

#Figure 2
montage ./blood/Figure2.png ./brain/frontal_cortex/Figure2.png ./brain/cerebellum/Figure2.png -tile 1x3 -geometry +2+2 Figure2.png

#Figure 3
montage ./blood/Figure3.png ./brain/frontal_cortex/Figure3.png ./brain/cerebellum/Figure3.png -tile 1x3 -geometry +2+2 Figure3.png

#Figure 4
montage ./blood/Figure4.png ./brain/frontal_cortex/Figure4.png ./brain/cerebellum/Figure4.png -tile 3x1 -geometry +2+2 Figure4.png


#Figure S1
#Add labels
convert  ./blood/blood_betaplot.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/blood_betaplot.png
convert  './brain/frontal cortex_betaplot.png' -pointsize 60 -gravity NorthWest -annotate +0+0 "B" './brain/frontal cortex_betaplot.png'
convert  ./brain/cerebellum_betaplot.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/cerebellum_betaplot.png
montage ./blood/blood_betaplot.png './brain/frontal cortex_betaplot.png' ./brain/cerebellum_betaplot.png  -tile 2x2 -geometry +2+2 FigureS1.png

#Figure S6
convert  ./blood/genes/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/genes/go.png
convert  ./blood/genes/Hannum/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./blood/genes/Hannum/go.png
montage ./blood/genes/go.png ./blood/genes/Hannum/go.png  -tile 1x2 -geometry +2+2 FigureS6.png

#Figure S7
convert  ./brain/frontal_cortex/genes/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./brain/frontal_cortex/genes/go.png
convert  ./brain/cerebellum/genes/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./brain/cerebellum/genes/go.png
convert  ./brain/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/go.png
montage ./brain/frontal_cortex/genes/go.png ./brain/cerebellum/genes/go.png ./brain/go.png -tile 1x3 -geometry +2+2 FigureS7.png
