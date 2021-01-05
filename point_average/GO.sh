#Plot the GO
BLOODGO1=../results/point_average/blood/genes/pantherChart_cluster1.txt
BLOODGO2=../results/point_average/blood/genes/pantherChart_cluster2.txt
FCTXGO=../results/point_average/brain/frontal_cortex/genes/pantherChart.txt
CEREBELLUMGO=../results/point_average/brain/cerebellum/genes/all_pantherChart.txt
HANNUMGO=../results/point_average/blood/genes/Hannum/pantherChart.txt
HORVATHGO=../results/point_average/brain/cerebellum/genes/horvath/pantherChart.txt
OUTDIR=../results/point_average/
./gene_ontology.py --blood1 $BLOODGO1 --blood2 $BLOODGO2 --frontal_cortex $FCTXGO --cerebellum $CEREBELLUMGO --horvath $HORVATHGO --hannum $HANNUMGO --outdir $OUTDIR
