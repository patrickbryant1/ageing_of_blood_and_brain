###Cerebellum
SEL_MARKERS=../../results/point_average/brain/cerebellum/ra_sig_markers.csv
MV=../../results/point_average/brain/cerebellum_marker_values.npy
AGES=../../results/point_average/brain/cerebellum_ages.csv
OUTDIR=../../results/point_average/brain/cerebellum/models/
./age_rf.py --selected_markers $SEL_MARKERS --marker_values $MV --ages $AGES --outdir $OUTDIR
