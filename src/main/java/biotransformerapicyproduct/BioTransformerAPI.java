/**
 *
 * @author Leon Pu , Anisha Jauhari
 *
 */
package biotransformerapicyproduct;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

import biotransformerapis.predictors.BioTransformerAPIs;
import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import cyProduct.cyProductMain;
//import reactantpredictor.utils.FileUtils;
import utils.Utilities;
import utils.Validator;

public class BioTransformerAPI implements BioTransformerAPIs {
	public static SmilesParser	smiParser	= new SmilesParser(SilentChemObjectBuilder.getInstance());
	public static SmilesGenerator smiGen 	= new SmilesGenerator().isomeric();

	@Override
	public Object predict(LinkedHashMap<String, Object> parameters) throws Exception {
		Validator validator = new Validator();
		Object results  = null;
		IAtomContainerSet inputMolecules = null;
		Validator.ioFormats inputFormat = Validator.ioFormats.valueOf( StringUtils.upperCase( (String) parameters.get("inputFormat")));
		Validator.ioFormats outputFormat = Validator.ioFormats.valueOf( StringUtils.upperCase( (String) parameters.get("outputFormat")));
		ArrayList<String> enzymes = (ArrayList<String>) parameters.get("properties");
		LinkedHashMap<String, Object> arguments = (LinkedHashMap<String, Object>) parameters.get("arguments");
		if(validator.validateParameters(parameters)) {
				if(inputFormat == Validator.ioFormats.IATOMCONTAINERSET) {
				inputMolecules = (IAtomContainerSet) parameters.get("input");
			}
			else if(inputFormat == Validator.ioFormats.SDFILE) {
				inputMolecules = Utilities.parseSdf(String.valueOf(parameters.get("input")) );
			}

			IAtomContainerSet molsWithPredictions = null;
			try {
				boolean useCypReact = (boolean) arguments.get("useCypReact");
				Double scoreThreshold = (Double) arguments.get("scoreThreshold");
				molsWithPredictions = runPredictions(inputMolecules, enzymes, useCypReact, scoreThreshold);
			} catch (Exception e) {
				e.printStackTrace();
			}

			if(outputFormat == Validator.ioFormats.IATOMCONTAINERSET) {
				results = molsWithPredictions;
			}
			else if(outputFormat == Validator.ioFormats.SDFILE) {
				String outputFileName = (String) parameters.get("output");

				if(molsWithPredictions != null) {
					Utilities.saveAtomContainerSetToSDF(molsWithPredictions, outputFileName);
					results = outputFileName;
				}
			}
		}
		return results;
	}

	/**
	 * This function This function will predict metabolites for each molecule in the input molecule set for the input list of CYP450 enzymes
	 * When useCypReact = true, the prediction will run CypReact first and predicts metabolites for reactants only
	 * When useCypReact = false, the input molecule will be treated as reactant.
	 * @param molecules
	 * @param enzymeNames
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet runPredictions(IAtomContainerSet molecules, ArrayList<String> enzymeNames, boolean useCypReact, Double scoreThreshold) throws Exception{
		IAtomContainerSet results = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			IAtomContainerSet oneResult = runOnePrediction(molecules.getAtomContainer(i), enzymeNames, useCypReact, scoreThreshold);
			results.add(oneResult);
		}
		return results;
	}
	/**
	 * This function will predict metabolites for the input molecule for the input list of CYP450 enzymes
	 * When useCypReact = true, the prediction will run CypReact first and predicts metabolites for reactants only
	 * When useCypReact = false, the input molecule will be treated as reactant.
	 * @param molecule
	 * @param enzymeNames
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet runOnePrediction(IAtomContainer molecule, ArrayList<String> enzymeNames, boolean useCypReact, Double scoreThreshold) throws Exception{
		System.out.println("CyProduct Working");
		InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
		IAtomContainerSet results = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		HashMap<String, IAtomContainer> existed = new HashMap();
		//ArrayList<String> checkExist = new ArrayList<>();
		for(int i = 0; i < enzymeNames.size(); i++){
			//ArrayList<IAtomContainer> results_enzyme = new ArrayList<>(); 
			IAtomContainerSet metabolites = cyProductMain.makePrediction(molecule, enzymeNames.get(i), null, useCypReact);	
			//System.out.println("CyProduct Prediction Done");
			for(int j = 0; j < metabolites.getAtomContainerCount(); j++){
				IAtomContainer oneMetabolite = metabolites.getAtomContainer(j);
				Double score_current = oneMetabolite.getProperty("Score");
				
				InChIGenerator inchiGen = inchiFactory.getInChIGenerator(oneMetabolite);
				//We only check the first 14 characters within the inChiKey
				String inChiKey = inchiGen.getInchiKey().split("-")[0];
				//If the score is lower than the given threshold, then 
				if(score_current <= scoreThreshold) continue;
				if(!existed.keySet().contains(inChiKey)){
					Double score_round = Math.round(score_current * 100.0) / 100.0;
					oneMetabolite.setProperty("Score", score_round);
					existed.put(inChiKey, oneMetabolite);
					
				}
				else{
					IAtomContainer storedMolecule = existed.get(inChiKey);
					Double score_previous = storedMolecule.getProperty("Score");
					//If the metabolite is produced by more than one enzymes, assign the higher score to the metabolite.
					if(score_current > score_previous){
						Double score_round = Math.round(score_current * 100.0) / 100.0;
						storedMolecule.setProperty("Score", score_round);
					}
					String enzymeList = storedMolecule.getProperty("Enzyme");
					enzymeList = enzymeList + " " + enzymeNames.get(i);
					storedMolecule.setProperty("Enzyme", enzymeList);
				}
			}
		}
		for(IAtomContainer metabolite : existed.values()){
			IAtomContainer mole = smiParser.parseSmiles(smiGen.create(metabolite));
			mole.setProperties(metabolite.getProperties());
			//results.addAtomContainer(metabolite);
			results.addAtomContainer(mole);
		}
		System.out.println("CyProduct Done");
		return results;
	}
}
