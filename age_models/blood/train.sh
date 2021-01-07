###Blood
SEL_MARKERS=../../results/point_average/blood/ra_sig_markers.csv
MV=/home/pbryant/results/methylation/aging_of_blood/10/median/marker_values.npy
AGES=../../results/point_average/blood/ages.csv
HANNUM_MARKERS=../../point_average/blood/hannum_markers.csv
FCDF=../../results/point_average/blood/marker_max_FC_pval.csv
OUTDIR=../../results/point_average/blood/models/
./age_rf.py --selected_markers $SEL_MARKERS --marker_values $MV --ages $AGES --hannum_markers $HANNUM_MARKERS --max_fold_change_df $FCDF --outdir $OUTDIR
