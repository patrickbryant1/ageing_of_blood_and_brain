###Blood
#FC
SEL_MARKERS=../../results/point_average/blood/FC_ra_sig_markers.csv
MV=../../results/point_average/blood/marker_values.npy
AGES=../../results/point_average/blood/ages.csv
HANNUM_MARKERS=../../point_average/blood/hannum_markers.csv
FCDF=../../results/point_average/blood/marker_max_FC_pval.csv
MODE=FC
OUTDIR=../../results/point_average/blood/models/
./age_rf.py --selected_markers $SEL_MARKERS --marker_values $MV --ages $AGES --hannum_markers $HANNUM_MARKERS --max_fold_change_df $FCDF --mode $MODE --outdir $OUTDIR

#AV
SEL_MARKERS=../../results/point_average/blood/abs_ra_sig_markers.csv
MV=../../results/point_average/blood/marker_values.npy
AGES=../../results/point_average/blood/ages.csv
HANNUM_MARKERS=../../point_average/blood/hannum_markers.csv
FCDF=../../results/point_average/blood/marker_max_FC_pval.csv
MODE=abs
OUTDIR=../../results/point_average/blood/models/
./age_rf.py --selected_markers $SEL_MARKERS --marker_values $MV --ages $AGES --hannum_markers $HANNUM_MARKERS --max_fold_change_df $FCDF --mode $MODE --outdir $OUTDIR

#Overlap
SEL_MARKERS=../../results/point_average/blood/overlap_ra_sig_markers.csv
MV=../../results/point_average/blood/marker_values.npy
AGES=../../results/point_average/blood/ages.csv
HANNUM_MARKERS=../../point_average/blood/hannum_markers.csv
FCDF=../../results/point_average/blood/marker_max_FC_pval.csv
MODE=overlap
OUTDIR=../../results/point_average/blood/models/
./age_rf.py --selected_markers $SEL_MARKERS --marker_values $MV --ages $AGES --hannum_markers $HANNUM_MARKERS --max_fold_change_df $FCDF --mode $MODE --outdir $OUTDIR
