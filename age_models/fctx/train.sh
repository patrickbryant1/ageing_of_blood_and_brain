###frontal_cortex
#FC
SEL_MARKERS=../../results/point_average/brain/frontal_cortex/FC_ra_sig_markers.csv
MV=../../results/point_average/brain/frontal_cortex_marker_values.npy
AGES=../../results/point_average/brain/frontal_cortex_ages.csv
HORVATH_MARKERS=../../point_average/brain/datMiniAnnotation3.csv
FCDF=../../results/point_average/brain/frontal_cortex_marker_max_FC_pval.csv
MODE=FC
OUTDIR=../../results/point_average/brain/frontal_cortex/models/
./age_rf.py --selected_markers $SEL_MARKERS --marker_values $MV --ages $AGES --horvath_markers $HORVATH_MARKERS --max_fold_change_df $FCDF --mode $MODE --outdir $OUTDIR

#AV
SEL_MARKERS=../../results/point_average/brain/frontal_cortex/abs_ra_sig_markers.csv
MV=../../results/point_average/brain/frontal_cortex_marker_values.npy
AGES=../../results/point_average/brain/frontal_cortex_ages.csv
HORVATH_MARKERS=../../point_average/brain/datMiniAnnotation3.csv
FCDF=../../results/point_average/brain/frontal_cortex_marker_max_FC_pval.csv
MODE=abs
OUTDIR=../../results/point_average/brain/frontal_cortex/models/
./age_rf.py --selected_markers $SEL_MARKERS --marker_values $MV --ages $AGES --horvath_markers $HORVATH_MARKERS --max_fold_change_df $FCDF --mode $MODE --outdir $OUTDIR

#Overlap
SEL_MARKERS=../../results/point_average/brain/frontal_cortex/overlap_ra_sig_markers.csv
MV=../../results/point_average/brain/frontal_cortex_marker_values.npy
AGES=../../results/point_average/brain/frontal_cortex_ages.csv
HORVATH_MARKERS=../../point_average/brain/datMiniAnnotation3.csv
FCDF=../../results/point_average/brain/frontal_cortex_marker_max_FC_pval.csv
MODE=overlap
OUTDIR=../../results/point_average/brain/frontal_cortex/models/
./age_rf.py --selected_markers $SEL_MARKERS --marker_values $MV --ages $AGES --horvath_markers $HORVATH_MARKERS --max_fold_change_df $FCDF --mode $MODE --outdir $OUTDIR
