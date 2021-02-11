Please download the CyProduct.jar to run the tool and contact stian2@ualberta.ca if you have any questions or concerns regarding the repository.

To run the CyProduct tool 

the user should use command in the terminal as:
Java -jar CyProduct.jar QueryMolecule enzymeList OutputFolderPath

The QueryMolecule can be either a SMILES string or the path to a sdf file. Note that if the QueryMolecule is a SMILES string, it should use format:SMILES=yourSmilesString.
The enzymeList is the list of CYP450 enzymes you want to run CyProduct on. If the enzyme list contains more than one enzyme, CyProduct will predict the metabolites for all of them. For example, if the input enzyme list is 1A2,3A4, then metabolites predicted are catalyzed by either 1A2 or 3A4 or both of them.
The OutputFoldPath is the path to the folder where you want to store the predicted results.
One example of running the tool in the commandline is, you can run it as:
java -jar CyProduct.jar SMILES=CC(=O)Nc1ccc(O)cc1 1A2,2A6,2B6,2C8,2C9,2C19,2D6,2E1,3A4 E:\Users\cyproduct\

Note that there is a class called BioTransformerAPI in the jar file. It provides two static functions that take one molecule (or molecules), a list of cyp450 enzymes and a boolean variable useCypReact as input, and predict the corresponding metabolites. Note that if the useCypReact is set as false, then the CypReact filter module will be disabled. The user can call those two static functions directly when they want to use CyProduct in their own tool/software.

Please note that when you use the CyProduct tool or its API, please input enzyme without “CYP”. For example, please use 1A2 other than CYP1A2 in both cases.


-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Description of folders/files in the repository

1. The Datasets folder:

	1. EBoMD.sdf contains the compounds derived from the Zaretzki’s CYP450 substrate datasets of 679 compounds. It is used to train the CypBoM module.
	2. EBoMD2.sdf is the hold-out dataset used to test our model and generate hold-out test results
	3. HDFAME folder contains a 60Reactants.sdf which is a subset of EBoMD2 and a 30NonReactants.sdf file. These two datasets are used when we compare with FAME2 and GLORY.
	4.13321_2018_324_MOESM6_ESM.sdf is the biotransformer’s hold-out dataset mentioned in the CyProduct paper. Please find it in the biotransformer paper (DOI: 10.1186/s13321-018-0324-5).

2. PredictedResults folder:

	1. Reactants_Result folder:
		1. CyProductResult_ForCYP2C9_2D6_3A4 folder contains the predicted results for the union of CYP2C9, 2D6 and 3A4 for the 60Reactants.sdf.
		2. CyProductResults_For_9CYPs folder contains the predicted results for each of the 9 CYP450 enzymes for EBoMD2.
	2. NonReactants_Result folder:
		1. CyProducts_3CYPS  folder contains the predicted results for the union of CYP2C9, 2D6 and 3A4 for the 30NonReactants.sdf.
		2. CyProduct_AllCyps folder contains the predicted results for the union of the nine CYP450 enzymes for the 30NonReactants.sdf.
	3. CyProduct_Results_For59BiotransformerMolecules folder contains the predicted results for the 59 compounds of the biotransformer’s hold-out dataset for the union of the 9 CYP450 enzymes
	Please note that the results reported in the manuscript are produced by using different combinations of the predicted metabolites mentioned above (or their subsets).

3. Reaction_rules.csv contains the SMIRKS Strings used in the MetaboGen module.
