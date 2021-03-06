package cyProduct;

import java.io.File;
import java.util.ArrayList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformerapicyproduct.BioTransformerAPI;
import utils.ReadMolecules;

public class cyProductMain {
	/**
	 * Define the static variables that will be used multiple times later
	 */
	public final String[] cyps = {"1A2","2A6","2B6","2C8","2C9","2C19","2D6","2E1","3A4"};
	public final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);	
	//public biotransformerapicyproduct.BioTransformerAPI biotransformeAPI = new biotransformerapicyproduct.BioTransformerAPI();
	
	public static void main(String[] args) throws Exception{
		if(args.length != 3) {
			throw new Exception("CyProduct: Input arguemtns don't match: queriedCompund, cyp450enzyme and outputPath");
			//return;
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
		if(enzyme.split(",").length > 0){
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
			if(!enzymeList.isEmpty() && enzymeList != null) cyProduct.makePredictionForEnzymeList(oneMole, enzymeList, outputPath, useCypReact);
			else{
				PredictionFunctions pf = new PredictionFunctions();
				pf.makePrediction(oneMole, enzyme, outputPath, useCypReact);
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
	public IAtomContainerSet readSMILES(String smiles) throws Exception{
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
	
	public IAtomContainerSet makePredictionForEnzymeList(IAtomContainer oneMole, ArrayList<String> cypList, String outputPath, boolean useCypReact) throws Exception{
		//We don't use score in this tool
		IAtomContainerSet results = BioTransformerAPI.runOnePrediction(oneMole, cypList, useCypReact,0.0);
		if(outputPath!=null) PredictorForAllThreeBoMs.outputResultIAtomContainerSet(results, outputPath);
		return results;
		
	}
	
	
}
