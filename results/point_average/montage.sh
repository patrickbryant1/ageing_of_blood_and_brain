#Figure 1
#Add labels
#montage ./blood/Figure1.png ./brain/frontal_cortex/Figure1.png ./brain/cerebellum/Figure1.png -tile 3x1 -geometry +2+2 Figure1.png

#Figure 2
#montage ./blood/Figure2.png ./brain/frontal_cortex/Figure2.png ./brain/cerebellum/Figure2.png -tile 3x1 -geometry +2+2 Figure2.png

# #Figure 7
montage ./blood/Figure7.png ./brain/frontal_cortex/Figure7.png ./brain/cerebellum/Figure7.png -tile 3x1 -geometry +2+2 Figure7.png

# #Figure 8
montage ./blood/Figure8.png ./brain/frontal_cortex/Figure8.png ./brain/cerebellum/Figure8.png -tile 3x1 -geometry +2+2 Figure8.png

#Figure 9
convert ./blood/CV_blood.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/CV_blood.png
convert ./brain/frontal_cortex/CV_frontal_cortex.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./brain/frontal_cortex/CV_frontal_cortex.png
convert ./brain/cerebellum/CV_cerebellum.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/cerebellum/CV_cerebellum.png
montage ./blood/CV_blood.png ./brain/frontal_cortex/CV_frontal_cortex.png ./brain/cerebellum/CV_cerebellum.png -tile 1x3 -geometry +2+2 Figure9.png

#Figure 8
# convert ./blood/models/cv_results.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/models/cv_results.png
# convert ./brain/frontal_cortex/models/cv_results.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./brain/frontal_cortex/models/cv_results.png
# convert ./brain/cerebellum/models/cv_results.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/cerebellum/models/cv_results.png
# montage ./blood/models/cv_results.png ./brain/frontal_cortex/models/cv_results.png ./brain/cerebellum/models/cv_results.png -tile 3x1 -geometry +2+2 Figure7.png
# #Figure S1
# #Add labels
# convert  ./blood/blood_betaplot.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/blood_betaplot.png
# convert  './brain/frontal cortex_betaplot.png' -pointsize 60 -gravity NorthWest -annotate +0+0 "B" './brain/frontal cortex_betaplot.png'
# convert  ./brain/cerebellum_betaplot.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/cerebellum_betaplot.png
# montage ./blood/blood_betaplot.png './brain/frontal cortex_betaplot.png' ./brain/cerebellum_betaplot.png  -tile 2x2 -geometry +2+2 FigureS1.png
#
# #Figure S2
convert ./blood/clustering/FC1_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/clustering/FC1_unnormalized.png
convert ./blood/clustering/FC2_unnormalized.png  -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./blood/clustering/FC2_unnormalized.png
convert ./brain/frontal_cortex/clustering/FC1_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/frontal_cortex/clustering/FC1_unnormalized.png
convert ./brain/cerebellum/clustering/FC1_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "D" ./brain/cerebellum/clustering/FC1_unnormalized.png
convert ./brain/cerebellum/clustering/FC2_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "E" ./brain/cerebellum/clustering/FC2_unnormalized.png
montage ./blood/clustering/FC1_unnormalized.png ./blood/clustering/FC2_unnormalized.png ./brain/frontal_cortex/clustering/FC1_unnormalized.png ./brain/cerebellum/clustering/FC1_unnormalized.png  ./brain/cerebellum/clustering/FC2_unnormalized.png  -tile 3x2 -geometry +2+2 FigureS2.png

# #Figure S6
convert  ./blood/genes/go1.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./blood/genes/go1.png
convert  ./blood/genes/go2.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./blood/genes/go2.png
convert  ./blood/genes/Hannum/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./blood/genes/Hannum/go.png
montage ./blood/genes/go1.png ./blood/genes/go2.png ./blood/genes/Hannum/go.png  -tile 2x2 -geometry +2+2 FigureS6.png
#
# #Figure S7
convert  ./brain/frontal_cortex/genes/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A" ./brain/frontal_cortex/genes/go.png
convert  ./brain/cerebellum/genes/go1.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B" ./brain/cerebellum/genes/go1.png
convert  ./brain/cerebellum/genes/go2.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C" ./brain/cerebellum/genes/go2.png
convert  ./brain/go.png -pointsize 60 -gravity NorthWest -annotate +0+0 "D" ./brain/go.png
montage ./brain/frontal_cortex/genes/go.png ./brain/cerebellum/genes/go1.png ./brain/cerebellum/genes/go2.png ./brain/go.png -tile 2x2 -geometry +2+2 FigureS7.png
