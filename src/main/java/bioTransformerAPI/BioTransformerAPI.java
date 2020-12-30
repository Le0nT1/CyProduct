package bioTransformerAPI;

import java.util.ArrayList;
import java.util.HashMap;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import cyProduct.cyProductMain;

public class BioTransformerAPI {
	/**
	 * This function This function will predict metabolites for each molecule in the input molecule set for the input list of CYP450 enzymes
	 * When useCypReact = true, the prediction will run CypReact first and predicts metabolites for reactants only
	 * When useCypReact = false, the input molecule will be treated as reactant.
	 * @param molecules
	 * @param enzymeNames
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet runPredictions(IAtomContainerSet molecules, ArrayList<String> enzymeNames, boolean useCypReact) throws Exception{
		IAtomContainerSet results = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			IAtomContainerSet oneResult = runOnePrediction(molecules.getAtomContainer(i), enzymeNames, useCypReact);
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
	public static IAtomContainerSet runOnePrediction(IAtomContainer molecule, ArrayList<String> enzymeNames, boolean useCypReact) throws Exception{
		System.out.println("CyProduct Working");
		InChIGeneratorFactory inchiFactory = InChIGeneratorFactory.getInstance();
		IAtomContainerSet results = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		HashMap<String, IAtomContainer> existed = new HashMap<>();
		//ArrayList<String> checkExist = new ArrayList<>();
		for(int i = 0; i < enzymeNames.size(); i++){
			//ArrayList<IAtomContainer> results_enzyme = new ArrayList<>(); 
			IAtomContainerSet metabolites = cyProductMain.makePrediction(molecule, enzymeNames.get(i), null, useCypReact);	
			System.out.println("CyProduct Prediction Done");
			for(int j = 0; j < metabolites.getAtomContainerCount(); j++){
				IAtomContainer oneMoetabolite = metabolites.getAtomContainer(j);
				InChIGenerator inchiGen = inchiFactory.getInChIGenerator(oneMoetabolite);
				//We only check the first 14 characters within the inChiKey
				String inChiKey = inchiGen.getInchiKey().split("-")[0];
				if(!existed.keySet().contains(inChiKey)){
					existed.put(inChiKey, oneMoetabolite);
				}
				else{
					IAtomContainer storedMolecule = existed.get(inChiKey);
					String enzymeList = storedMolecule.getProperty("Enzyme");
					enzymeList = enzymeList + " " + enzymeNames.get(i);
					storedMolecule.setProperty("Enzyme", enzymeList);
				}
			}
		}
		for(IAtomContainer metabolite : existed.values()){
			results.addAtomContainer(metabolite);
		}
		System.out.println("CyProduct Done");
		return results;
	}
}
