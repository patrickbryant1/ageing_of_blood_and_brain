#!/usr/bin/env

#Montage
for mode in 'abs' 'FC' 'overlap'
do
convert  $mode'1.png' -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  $mode'1.png'
convert  $mode'2.png' -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  $mode'2.png'
convert  $mode'3.png' -pointsize 60 -gravity NorthWest -annotate +0+0 "C"  $mode'3.png'
convert  $mode'4.png' -pointsize 60 -gravity NorthWest -annotate +0+0 "D"  $mode'4.png'
convert  $mode'_tsne.png' -pointsize 60 -gravity NorthWest -annotate +0+0 "E" $mode'_tsne.png'
convert $mode'_gradients.png' -pointsize 60 -gravity NorthWest -annotate +0+0 "F" $mode'_gradients.png'
montage $mode'1.png' $mode'2.png' $mode'3.png' $mode'4.png' $mode'_tsne.png' $mode'_gradients.png' -tile 3x2 -geometry +2+2  $mode'_montage.png'
#-title 'Cerebellum' -pointsize 30


# convert  1_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "A"  1_unnormalized.png
# convert  2_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "B"  2_unnormalized.png
# convert  3_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "C"  3_unnormalized.png
# convert  4_unnormalized.png -pointsize 60 -gravity NorthWest -annotate +0+0 "D"  4_unnormalized.png
# montage 1_unnormalized.png 2_unnormalized.png 3_unnormalized.png 4_unnormalized.png -tile 2x2 -geometry +2+2  FigureS4.png
done
