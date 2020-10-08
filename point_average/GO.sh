#Plot the GO
BLOODGO=../results/point_average/blood/genes/pantherChart.txt
FCTXGO=../results/point_average/brain/frontal_cortex/genes/pantherChart.txt
CEREBELLUMGO=../results/point_average/brain/cerebellum/genes/pantherChart.txt
HANNUMGO=../results/point_average/blood/genes/Hannum/pantherChart.txt
HORVATHGO=../results/point_average/brain/cerebellum/genes/horvath/pantherChart.txt
OUTDIR=../results/point_average/
./gene_ontology.py --blood $BLOODGO --frontal_cortex $FTCXGO --cerebellum $CEREBELLUMGO --horvath $HORVATHGO --hannum $HANNUMGO --outdir $OUTDIR
