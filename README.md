# ageing_of_blood_and_brain
Scripts for assessing the relationship between ageing and changes in the human blood and brain methylomes.

1. Download the necessary data E-GEOD-40279 (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-40279/), E-GEOD-36194 (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-36194/) and E-GEOD-15745 (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-15745/) from https://www.ebi.ac.uk/arrayexpress/.
2. Run the scripts 
./point_average/blood/analyze.sh
./point_average/brain/analyze.sh
3. Fetch the GO info by pasting the result files
./results/point_average/blood/genes/unique_genes.txt 
./results/point_average/brain/cerebellum/genes/unique_genes.txt
./results/point_average/brain/frontal_cortex/genes/unique_genes.txt
into http://pantherdb.org/
4. The montage.sh scripts in the results folders will montage all figures as found in the publication.
