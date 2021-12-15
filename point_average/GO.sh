#Plot the GO
BLOODGO1=../results/point_average/blood/genes/pantherChart_cluster1.txt
BLOODGO2=../results/point_average/blood/genes/pantherChart_cluster2.txt
FCTXGO=../results/point_average/brain/frontal_cortex/genes/pantherChart_cluster1.txt
CEREBELLUMGO1=../results/point_average/brain/cerebellum/genes/pantherChart_cluster1.txt
CEREBELLUMGO2=../results/point_average/brain/cerebellum/genes/pantherChart_cluster2.txt
HANNUMGO=../results/point_average/blood/genes/Hannum/pantherChart.txt
HORVATHGO=../results/point_average/brain/cerebellum/genes/horvath/pantherChart.txt
OUTDIR=../results/point_average/
./gene_ontology.py --blood1 $BLOODGO1 --blood2 $BLOODGO2 --frontal_cortex $FCTXGO \
--cerebellum1 $CEREBELLUMGO1 --cerebellum2 $CEREBELLUMGO2 \
--horvath $HORVATHGO --hannum $HANNUMGO --outdir $OUTDIR
