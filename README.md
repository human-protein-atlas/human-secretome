# The Human Secretome
The proteins secreted by human cells (the secretome) are important for the basic understanding of human biology, but also for the identification of potential targets for future diagnostics and therapy. Here, we present a genome-wide analysis of proteins predicted to be secreted from human cells providing information about their final localization in the human body, including the proteins actively secreted to peripheral blood. The proteins detected in human blood by mass spectrometry-based proteomics and antibody-based immune-assays are also presented to provide an open access Blood Atlas resource to facilitate exploration of individual proteins actively secreted by human cells.

## Code explanation
1.	The perl script proteinclasses.pl predicts transmembrane regions and signal peptides based on several different prediction methods, and use the resulting majority-based predictions to classify all transcripts belonging to ensemble genes as being secreted (signal peptide and no transmembrane region), membrane-bound (has at least one trans-membrane region) or intracellular (no signal peptide and no trans-membrane region).

2.	The php script classify_gene.php is used to classify all ensemble genes into secreted and non-secreted genes, based on the result of the proteinclasses.pl script above as well as annotation data regarding secretion from the UniProt database. All genes with at least one predicted or annotated secreted transcript is classified as secreted.

3.	The php script classify_function.php is used to assign one functional annotation to each of the secreted genes, mainly based on selected UniProt keyword terms for molecular function and biological process.

4.	The R script main.R is used to create the Figures 3 B-D in the paper, based on expression data from  The Human Protein Atlas and the functional annotations included in the paper.

5.	The R script function.R includes functions used in main.R. Please note only the functions “Aasa_facet_heatmap” and “chord_classification” were used from the script main.R.
