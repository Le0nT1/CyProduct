package cyProduct.predictionhelpers;

import java.util.ArrayList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformerapicypreact.BioTransformerAPI_cypreact;
import instances.GenerateInstance;
import utils.Utilities;

public class PredictionFunctions {
	SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
	SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Stereo);
	ArrayList<String> rawFeatureList_TypeOne;
	ArrayList<String> rawFeatureList_TypeTwo;
	ArrayList<String> rawFeatureList_TypeThree;
	PredictorForAllThreeBoMs common_usage_PredictorForAllThreeBoMs = new PredictorForAllThreeBoMs();
	public PredictionFunctions(IAtomContainer molecule) throws Exception {
		UniqueIDFunctionSet.assignUniqueID(molecule);
		GenerateInstance gi = new GenerateInstance();
		this.rawFeatureList_TypeOne = gi.generateRawInstances(molecule, 1);
		//System.out.println("TypeOne raw feature done");
		this.rawFeatureList_TypeTwo = gi.generateRawInstances(molecule, 2);
		//System.out.println("TypeTwo raw feature done");
		this.rawFeatureList_TypeThree = gi.generateRawInstances(molecule, 3);
		//System.out.println("TypeThree raw feature done");
	}
	public IAtomContainerSet makePrediction(IAtomContainer oneMole, String cyp, String outputPath, boolean useCypReact) throws Exception{
		/**
		 * Setting up support files including feature file and model files
		 */
		if(useCypReact){
			if(!isReactant(oneMole, cyp)){
				IAtomContainerSet results = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
				if(outputPath!=null){
					common_usage_PredictorForAllThreeBoMs.outputResultIAtomContainerSet(results, outputPath);
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
		
		//System.out.println("Smiles of the original structure: " + sg.create(tempMole.clone()));
		/**
		 * The makePrediction function in the PredictorForAllThreeBoMs class was primarily designed for multiple molecules. The function can be modified to handle one molecule later
		 */

		//long start=System.currentTimeMillis();
		//TypeOne
		IAtomContainerSet moleculeSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		moleculeSet.addAtomContainer(tempMole.clone());
		ArrayList<MoleResultPair> predResultList_typeOne = bomPred_typeOne.makePrediction(moleculeSet, this.rawFeatureList_TypeOne, 1);	
		//TypeTwo
		moleculeSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		moleculeSet.addAtomContainer(tempMole.clone());
		ArrayList<MoleResultPair> predResultList_typeTwo = bomPred_typeTwo.makePrediction(moleculeSet, this.rawFeatureList_TypeTwo, 2);	
		//TypeThree
		moleculeSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		moleculeSet.addAtomContainer(tempMole.clone());
		ArrayList<MoleResultPair> predResultList_Three = bomPred_typeThree.makePrediction(moleculeSet, this.rawFeatureList_TypeThree, 3);
		IAtomContainer subtrate = null;
		if(predResultList_typeOne!=null && !predResultList_typeOne.isEmpty()){
			subtrate = predResultList_typeOne.get(0).getMolecule();
		}
		ArrayList<ArrayList<IAtom>> typeOne_BoMs = this.common_usage_PredictorForAllThreeBoMs.getTypeOneBoMList(predResultList_typeOne, subtrate);		
		if(predResultList_typeTwo!=null && !predResultList_typeTwo.isEmpty()) subtrate = predResultList_typeTwo.get(0).getMolecule();
		ArrayList<IAtom> typeTwo_BoMs = this.common_usage_PredictorForAllThreeBoMs.getTypeTwoBoMList(predResultList_typeTwo, subtrate);
		
		if(predResultList_Three!=null && !predResultList_Three.isEmpty()) subtrate = predResultList_Three.get(0).getMolecule();
		ArrayList<IAtom> typeThree_BoMs = this.common_usage_PredictorForAllThreeBoMs.getTypeThreeBoMList(predResultList_Three, subtrate);
		//long predict_time=System.currentTimeMillis();
		IAtomContainerSet results = ClarifyReactionType.arrangeReactionTypesAndPredictMetabolites(typeOne_BoMs, typeTwo_BoMs, typeThree_BoMs, moleculeSet.getAtomContainer(0));
		results = Utilities.removeDuplicates(results);
		results = getValidMetabolites(moleculeSet.getAtomContainer(0),results);
//		for(int i = 0; i < results.getAtomContainerCount(); i++) {
//			System.out.println(sg.create(results.getAtomContainer(i)));
//		}
		for(int k = 0; k < results.getAtomContainerCount(); k++){
			results.getAtomContainer(k).setProperty("Enzyme", ("CYP" + cyp));
		}
		if(outputPath!=null) this.common_usage_PredictorForAllThreeBoMs.outputResultIAtomContainerSet(results, outputPath);
		//System.out.println("Prediction time: " + (predict_time - start));
		return results;
	}
	/**
	 * This function will use CypReact to check if the molecule is a reactant for the input enzyme
	 * @param oneMole
	 * @param cyp
	 * @return
	 * @throws Exception
	 */
	public boolean isReactant(IAtomContainer oneMole, String cyp) throws Exception{
		//ReactantPred reactantPredictor = new ReactantPred();
		BioTransformerAPI_cypreact cypReactAPI = new BioTransformerAPI_cypreact();
		boolean isReactant = cypReactAPI.predictReactant(oneMole, cyp);
		//System.out.println("is Reaction for " + cyp + " : " + isReactant); 
		return isReactant;
	}
	/**
	 * This function will check every carbon atom within the molecule, see if its connection is valid, say valence.
	 * If all carbon atoms are valid, then the molecule is valid and true is returned. Otherwise it's not.
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public boolean validateMetabolites(IAtomContainer substrate, IAtomContainer oneMole) throws Exception{
		//A valid metabolite must have a score assigned in CyProduct		
		if(oneMole.getProperty("Score")==null) return false;
		for(int i = 0; i < oneMole.getAtomCount(); i++){
			IAtom oneAtom = oneMole.getAtom(i);
			if(oneAtom.getSymbol().equals("R")) {
				if(oneAtom.getProperty("UniqueID")==null) {
					return false;
				}
				else {
					IAtom origin_atom = UniqueIDFunctionSet.getAtomByUniqueID(oneAtom.getProperty("UniqueID"), substrate);
					oneAtom.setAtomicNumber(origin_atom.getAtomicNumber());
					AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
				}
			}
			if(oneAtom.getSymbol().equalsIgnoreCase("C") && (oneAtom.getImplicitHydrogenCount() + oneMole.getBondOrderSum(oneAtom)) > 4){
				return false;
			}
			//System.out.println(oneAtom.getSymbol());
	
			
		}
		return true;
	}
	
	public IAtomContainerSet getValidMetabolites(IAtomContainer substrate, IAtomContainerSet molecules) throws Exception{
		IAtomContainerSet results = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			if(validateMetabolites(substrate, molecules.getAtomContainer(i))){				
				IAtomContainer temp = this.sp.parseSmiles(this.sg.create(molecules.getAtomContainer(i)));
				temp.addProperties(molecules.getAtomContainer(i).getProperties());
				results.addAtomContainer(temp);
			}
		}
		return results;
	}
}
