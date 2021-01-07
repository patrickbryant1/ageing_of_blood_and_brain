#Figure 1
#Add labels
montage ./blood/Figure1.png ./brain/frontal_cortex/Figure1.png ./brain/cerebellum/Figure1.png -tile 3x1 -geometry +2+2 Figure1.png

#Figure 2
montage ./blood/Figure2.png ./brain/frontal_cortex/Figure2.png ./brain/cerebellum/Figure2.png -tile 3x1 -geometry +2+2 Figure2.png

#Figure 6
montage ./blood/Figure4.png ./brain/frontal_cortex/Figure4.png ./brain/cerebellum/Figure4.png -tile 3x1 -geometry +2+2 Figure6.png

#Figure 7
convert ./blood/models/cv_results.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/models/cv_results.png
convert ./brain/frontal_cortex/models/cv_results.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./brain/frontal_cortex/models/cv_results.png
convert ./brain/cerebellum/models/cv_results.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/cerebellum/models/cv_results.png
montage ./blood/models/cv_results.png ./brain/frontal_cortex/models/cv_results.png ./brain/cerebellum/models/cv_results.png -tile 3x1 -geometry +2+2 Figure7.png
# #Figure S1
# #Add labels
# convert  ./blood/blood_betaplot.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/blood_betaplot.png
# convert  './brain/frontal cortex_betaplot.png' -pointsize 60 -gravity NorthWest -annotate +0+0 "B" './brain/frontal cortex_betaplot.png'
# convert  ./brain/cerebellum_betaplot.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/cerebellum_betaplot.png
# montage ./blood/blood_betaplot.png './brain/frontal cortex_betaplot.png' ./brain/cerebellum_betaplot.png  -tile 2x2 -geometry +2+2 FigureS1.png
#
# #Figure S2
# convert ./blood/clustering/1_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/clustering/1_unnormalized.png
# convert ./blood/clustering/2_unnormalized.png  -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./blood/clustering/2_unnormalized.png
# convert ./brain/frontal_cortex/clustering/1_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/frontal_cortex/clustering/1_unnormalized.png
# convert ./brain/cerebellum/clustering/1_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "D" ./brain/cerebellum/clustering/1_unnormalized.png
# convert ./brain/cerebellum/clustering/2_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "E" ./brain/cerebellum/clustering/2_unnormalized.png
# convert ./brain/cerebellum/clustering/3_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "F" ./brain/cerebellum/clustering/3_unnormalized.png
# convert ./brain/cerebellum/clustering/4_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "G" ./brain/cerebellum/clustering/4_unnormalized.png
# montage ./blood/clustering/1_unnormalized.png ./blood/clustering/2_unnormalized.png ./brain/frontal_cortex/clustering/1_unnormalized.png ./brain/cerebellum/clustering/1_unnormalized.png  ./brain/cerebellum/clustering/2_unnormalized.png  ./brain/cerebellum/clustering/3_unnormalized.png  ./brain/cerebellum/clustering/4_unnormalized.png -tile 3x3 -geometry +2+2 FigureS2.png
# #Figure S6
# convert  ./blood/genes/go1.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/genes/go1.png
# convert  ./blood/genes/go2.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./blood/genes/go2.png
# convert  ./blood/genes/Hannum/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./blood/genes/Hannum/go.png
# montage ./blood/genes/go1.png ./blood/genes/go2.png ./blood/genes/Hannum/go.png  -tile 2x2 -geometry +2+2 FigureS6.png
#
# #Figure S7
# convert  ./brain/frontal_cortex/genes/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./brain/frontal_cortex/genes/go.png
# convert  ./brain/cerebellum/genes/go1.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./brain/cerebellum/genes/go1.png
# convert  ./brain/cerebellum/genes/go2.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/cerebellum/genes/go2.png
# convert  ./brain/cerebellum/genes/go3.png -pointsize 60 -gravity NorthWest -annotate +0+0 "D" ./brain/cerebellum/genes/go3.png
# montage ./brain/frontal_cortex/genes/go.png ./brain/cerebellum/genes/go1.png ./brain/cerebellum/genes/go2.png ./brain/cerebellum/genes/go3.png -tile 2x2 -geometry +2+2 FigureS7.png
#
# #Figure S8
# convert  ./brain/cerebellum/genes/go4.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./brain/cerebellum/genes/go4.png
# convert  ./brain/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./brain/go.png
# montage ./brain/cerebellum/genes/go4.png ./brain/go.png -geometry +2+2  -tile 2x1 FigureS8.png
