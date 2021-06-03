package cyProduct;

import java.io.File;
import java.util.ArrayList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import bioTransformerAPI.BioTransformerAPI;
import nu.xom.jaxen.function.SubstringAfterFunction;
import reactantpredictor.BioTransformerAPIs;
import reactantpredictor.ReactantPred;
import utils.ReadMolecules;
import utils.Utilities;

public class cyProductMain {
	/**
	 * Define the static variables that will be used multiple times later
	 */
	public static final String[] cyps = {"1A2","2A6","2B6","2C8","2C9","2C19","2D6","2E1","3A4"};
	public static final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);	
	public static void main(String[] args) throws Exception{
		if(args.length != 3) {
			System.out.println("Input arguemtns don't match: queriedCompund, cyp450enzyme and outputPath");
			return;
		}
		cyProductMain cyProduct = new cyProductMain();
		String inputPath = args[0];
		String enzyme = args[1];
		String outputPathFolder = args[2];
		File outputFoler = new File(outputPathFolder);
		if(!outputFoler.exists()) outputFoler.mkdirs();
		//useCypReact is set as true by default, because cyProduct alone is a tool that predicts metabolites for reactants. 
		boolean useCypReact = true;//true;
		ArrayList<String> enzymeList = new ArrayList<>();
		if(enzyme.split(",").length > 1){
			String[] enzymeParse = enzyme.split(",");
			for(int i = 0; i < enzymeParse.length; i++){
				enzymeList.add(enzymeParse[i]);
			}
		}
		InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
		IAtomContainerSet inputMoles = cyProduct.readSMILES(inputPath);
		if(inputMoles == null){
			inputMoles = cyProduct.readMolecule(inputPath);
		}
		for(int i = 0; i < inputMoles.getAtomContainerCount(); i++){
			IAtomContainer oneMole = inputMoles.getAtomContainer(i);
//			if(!oneMole.getProperty("cdk:Title").equals("(R)-mianserin_ 8-hydroxy")) continue;
//			if(!oneMole.getProperty("MolName").equals("Monodemethylated_mifepristone")) continue;
//			if(!oneMole.getProperty("NAME").equals("Benzo[a]pyrene")) continue;
//			boolean skip = true;
//			if(oneMole.getProperty("NAME").equals("Ondansetron")){
//				skip = false;
//			}
//			if(skip) continue;
			String outputPath;
			if(oneMole.getProperty("cdk:Title") != null){
				outputPath = outputPathFolder + "\\" + oneMole.getProperty("cdk:Title") + ".sdf";
		
			}
			else if(oneMole.getProperty("NAME") != null){
				outputPath = outputPathFolder + "\\" + oneMole.getProperty("NAME") + ".sdf";
			}
			else{
				String inChiKey = inchiFactory.getInChIGenerator(oneMole).getInchiKey();
				outputPath = outputPathFolder + "\\" + inChiKey + ".sdf";
			}
			if(!enzymeList.isEmpty() && enzymeList != null) makePredictionForEnzymeList(oneMole, enzymeList, outputPath, useCypReact);
			else{
				makePrediction(oneMole, enzyme, outputPath, useCypReact);
			}
		}
		
	}
	
	public IAtomContainerSet readMolecule(String inputPath) throws Exception{
		return ReadMolecules.extractMoleculesFromFile(inputPath);
	}
	/**
	 * This function will read a smiles string as input.
	 * If the input is a SMILES string, it returns a IAtomContainerSet that contains the molecule
	 * Otherwise return null;
	 * @param smiles
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet readSMILES(String smiles) throws Exception{
		try{
			SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
			IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			IAtomContainer oneMole = sp.parseSmiles(smiles);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(oneMole.getBuilder());
			Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
			aromaticity.apply(oneMole);
			adder.addImplicitHydrogens(oneMole);
			IChemObjectBuilder builder_1 = SilentChemObjectBuilder.getInstance();
			ModelBuilder3D mb3d = ModelBuilder3D.getInstance(builder_1);
			IAtomContainer oneMolecule = mb3d.generate3DCoordinates(oneMole, false);
//			StructureDiagramGenerator sdg = new StructureDiagramGenerator();
//			sdg.setMolecule(mol);
//			sdg.generateCoordinates();
//			IAtomContainer oneMolecule = sdg.getMolecule();
			AtomContainerManipulator.suppressHydrogens(oneMolecule);
			molecules.addAtomContainer(oneMolecule);
			return molecules;
		}catch(Exception e){
			return null;
		}
	}
	public static IAtomContainer readSMILES_oneMole(String smiles) throws Exception{
		try{
			SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
			IAtomContainer oneMole = sp.parseSmiles(smiles);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(oneMole.getBuilder());
			Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
			aromaticity.apply(oneMole);
			adder.addImplicitHydrogens(oneMole);
			IChemObjectBuilder builder_1 = SilentChemObjectBuilder.getInstance();
			ModelBuilder3D mb3d = ModelBuilder3D.getInstance(builder_1);
			IAtomContainer oneMolecule = mb3d.generate3DCoordinates(oneMole, false);
//			StructureDiagramGenerator sdg = new StructureDiagramGenerator();
//			sdg.setMolecule(mol);
//			sdg.generateCoordinates();
//			IAtomContainer oneMolecule = sdg.getMolecule();
			AtomContainerManipulator.suppressHydrogens(oneMolecule);
			return oneMolecule;
		}catch(Exception e){
			return null;
		}
	}
	public static IAtomContainerSet makePredictionForEnzymeList(IAtomContainer oneMole, ArrayList<String> cypList, String outputPath, boolean useCypReact) throws Exception{
		//We don't use score in this tool
		IAtomContainerSet results = BioTransformerAPI.runOnePrediction(oneMole, cypList, useCypReact,0.0);
		if(outputPath!=null) PredictorForAllThreeBoMs.outputResultIAtomContainerSet(results, outputPath);
		return results;
		
	}
	
	public static IAtomContainerSet makePrediction(IAtomContainer oneMole, String cyp, String outputPath, boolean useCypReact) throws Exception{
		//System.out.println(cyp);
		/**
		 * Setting up support files including feature file and model files
		 */
		if(useCypReact){
			if(!isReactant(oneMole, cyp)){
				IAtomContainerSet results = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
				if(outputPath!=null){
					PredictorForAllThreeBoMs.outputResultIAtomContainerSet(results, outputPath);
				}
				return results;
			}
		}
		PredictorForAllThreeBoMs bomPred_typeOne = new PredictorForAllThreeBoMs();
		PredictorForAllThreeBoMs bomPred_typeTwo = new PredictorForAllThreeBoMs();
		PredictorForAllThreeBoMs bomPred_typeThree = new PredictorForAllThreeBoMs();
		bomPred_typeOne.setup(cyp, 1);
		bomPred_typeTwo.setup(cyp, 2);
		bomPred_typeThree.setup(cyp, 3);
		/**
		 * Make prediction using the clone of the original molecule, so all the original information and properties will be kept			
		 */
		IAtomContainer tempMole = oneMole.clone();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tempMole);
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(tempMole);
		sdg.generateCoordinates();
		tempMole = sdg.getMolecule();		
		UniqueIDFunctionSet.assignUniqueID(tempMole);
		//System.out.println("Smiles of the original structure: " + sg.create(tempMole.clone()));
		/**
		 * The makePrediction function in the PredictorForAllThreeBoMs class was primarily designed for multiple molecules. The function can be modified to handle one molecule later
		 */

		//long start=System.currentTimeMillis();
		//TypeOne
		IAtomContainerSet moleculeSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		moleculeSet.addAtomContainer(tempMole.clone());
		ArrayList<MoleResultPair> predResultList_typeOne = bomPred_typeOne.makePrediction(moleculeSet, 1);	
		//TypeTwo
		moleculeSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		moleculeSet.addAtomContainer(tempMole.clone());
		ArrayList<MoleResultPair> predResultList_typeTwo = bomPred_typeTwo.makePrediction(moleculeSet, 2);	
		//TypeThree
		moleculeSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		moleculeSet.addAtomContainer(tempMole.clone());
		ArrayList<MoleResultPair> predResultList_Three = bomPred_typeThree.makePrediction(moleculeSet, 3);
		IAtomContainer subtrate = null;
		if(predResultList_typeOne!=null && !predResultList_typeOne.isEmpty()){
			subtrate = predResultList_typeOne.get(0).getMolecule();
		}
		ArrayList<ArrayList<IAtom>> typeOne_BoMs = PredictorForAllThreeBoMs.getTypeOneBoMList(predResultList_typeOne, subtrate);		
		if(predResultList_typeTwo!=null && !predResultList_typeTwo.isEmpty()) subtrate = predResultList_typeTwo.get(0).getMolecule();
		ArrayList<IAtom> typeTwo_BoMs = PredictorForAllThreeBoMs.getTypeTwoBoMList(predResultList_typeTwo, subtrate);
		
		if(predResultList_Three!=null && !predResultList_Three.isEmpty()) subtrate = predResultList_Three.get(0).getMolecule();
		ArrayList<IAtom> typeThree_BoMs = PredictorForAllThreeBoMs.getTypeThreeBoMList(predResultList_Three, subtrate);
		//long predict_time=System.currentTimeMillis();
		//System.out.println("Predicted TypeOne BoM");
		for(int t = 0; t < typeOne_BoMs.size(); t++){
			ArrayList<IAtom> oneBoM = typeOne_BoMs.get(t);
			int idx_one = predResultList_typeOne.get(0).getMolecule().indexOf(oneBoM.get(0));
			int idx_two = predResultList_typeOne.get(0).getMolecule().indexOf(oneBoM.get(1));
			//System.out.println(idx_one + "," + idx_two);
		}
		//System.out.println("Predicted TypeTwo BoM");
		for(int t = 0; t < typeTwo_BoMs.size(); t++){
			IAtom oneBoM = typeTwo_BoMs.get(t);
			int idx_one = predResultList_typeTwo.get(0).getMolecule().indexOf(oneBoM);
			//int idx_two = predResultList_typeOne.get(0).getMolecule().indexOf(oneBoM.get(1));
			//System.out.println(idx_one);
		}
		//System.out.println("Prediceted TypeThree BoM");
		for(int t = 0; t < typeThree_BoMs.size(); t++){
			IAtom oneBoM = typeThree_BoMs.get(t);
			int idx_one = predResultList_Three.get(0).getMolecule().indexOf(oneBoM);
			//int idx_two = predResultList_typeOne.get(0).getMolecule().indexOf(oneBoM.get(1));
			//System.out.println(idx_one);
		}
		//System.out.println("BoM prediction done");
		IAtomContainerSet results = ClarifyReactionType.arrangeReactionTypesAndPredictMetabolites(typeOne_BoMs, typeTwo_BoMs, typeThree_BoMs, moleculeSet.getAtomContainer(0));
		//long metabolite_time=System.currentTimeMillis();
		results = Utilities.removeDuplicates(results);
		//String outputPath = "C:/Users/Tian/Desktop/BioData/SOM-React/BioTransformerDB/WaitToMerge/Merged/" + "Molecule_" + i + "_metabolites.sdf";
		//String outputPath = "E:/CyProduct_results/" + "Molecule_" + i + "_metabolites.sdf";
		results = getValidMetabolites(results);
		for(int k = 0; k < results.getAtomContainerCount(); k++){
			results.getAtomContainer(k).setProperty("Enzyme", cyp);
		}
		if(outputPath!=null) PredictorForAllThreeBoMs.outputResultIAtomContainerSet(results, outputPath);
		//System.out.println("Prediction time: " + (predict_time - start));
		//System.out.println("Prediction time: " + (metabolite_time - predict_time));
		return results;
	}
	/**
	 * This function will use CypReact to check if the molecule is a reactant for the input enzyme
	 * @param oneMole
	 * @param cyp
	 * @return
	 * @throws Exception
	 */
	public static boolean isReactant(IAtomContainer oneMole, String cyp) throws Exception{
		//ReactantPred reactantPredictor = new ReactantPred();
		BioTransformerAPIs cypReactAPI = new BioTransformerAPIs();
		boolean isReactant = cypReactAPI.predictReactant(oneMole, cyp);
		//boolean isReactant = true;
		return isReactant;
	}
	/**
	 * This function will check every carbon atom within the molecule, see if its connection is valid, say valence.
	 * If all carbon atoms are valid, then the molecule is valid and true is returned. Otherwise it's not.
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static boolean validateMetabolites(IAtomContainer oneMole) throws Exception{
		for(int i = 0; i < oneMole.getAtomCount(); i++){
			IAtom oneAtom = oneMole.getAtom(i);
			if(oneAtom.getSymbol().equalsIgnoreCase("C") && (oneAtom.getImplicitHydrogenCount() + oneMole.getBondOrderSum(oneAtom)) > 4){
				return false;
			}
		}
		return true;
	}
	
	public static IAtomContainerSet getValidMetabolites(IAtomContainerSet molecules) throws Exception{
		IAtomContainerSet results = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			if(validateMetabolites(molecules.getAtomContainer(i))){
				results.addAtomContainer(molecules.getAtomContainer(i));
			}
		}
		return results;
	}
	
	
}
