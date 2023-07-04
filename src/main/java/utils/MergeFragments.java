package utils;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import cyProduct.predictionhelpers.UniqueIDFunctionSet;

public class MergeFragments {
	/**
	 * This function is used to merge the modified structure into the original structure of the compound to generate the ring-rearrangement metabolite
	 * The final structure (fragment) will first take atoms and bonds within the original "reactant", then update the stored structure using "modified"
	 *  
	 * @param reactant
	 * @param modified
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainer mergeFragmentsForRearrangementReaction(IAtomContainer reactant, IAtomContainer modified, ArrayList<IAtom> removedAtoms) throws Exception{
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		AtomContainerManipulator.suppressHydrogens(reactant);
		AtomContainerManipulator.suppressHydrogens(modified);
		//Check the presence of atoms using its uniqueID. Then create new IAtom and add it into the fragmentStructure
		//Add all atoms and bonds within the reactant into the fragmentStructure. We create atoms and bonds manually instead of using the clone() method.
		for(int i = 0; i < reactant.getAtomCount(); i++){
			IAtom atom_one = new Atom(reactant.getAtom(i).getSymbol());
			if(reactant.getAtom(i).getProperty("ReactiveSite")==null){
				atom_one.setFormalCharge(reactant.getAtom(i).getFormalCharge());
			}
			atom_one.setProperties(reactant.getAtom(i).getProperties());
			fragmentStructure.addAtom(atom_one);
		}
		//Remove all the atoms that are in the removedAtomsList from the fragmentStructure
		for(int i = 0; i < removedAtoms.size(); i++){
			IAtom atomToRemove_temp = removedAtoms.get(i);
			if(UniqueIDFunctionSet.containAtom(atomToRemove_temp, fragmentStructure)){
				IAtom atomToRemove = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atomToRemove_temp), fragmentStructure);
				fragmentStructure.removeAtom(atomToRemove);
			}
		}
		//Add atoms that are newly added into the fragment 
		for(int i = 0; i < modified.getAtomCount(); i++){
			Integer uid =UniqueIDFunctionSet.getUniqID(modified.getAtom(i));
			//If that atom doesn't have a UID, then it is newly added and was not in the original molecule. We assign a uid to it.
			if(uid == null){
				IAtom atom_one = new Atom(modified.getAtom(i).getSymbol());
				int maxUid = UniqueIDFunctionSet.getLargestUID(fragmentStructure);
				UniqueIDFunctionSet.assingUID_toAtom(modified.getAtom(i), (maxUid + 1));
				UniqueIDFunctionSet.assingUID_toAtom(atom_one, (maxUid + 1));
				fragmentStructure.addAtom(atom_one);
			}
			uid =UniqueIDFunctionSet.getUniqID(modified.getAtom(i));
			IAtom updateAtom = UniqueIDFunctionSet.getAtomByUniqueID(uid, fragmentStructure);
			updateAtom.setFormalCharge(modified.getAtom(i).getFormalCharge());

		}
		//Explore the structures 
		for(int i = 0; i < fragmentStructure.getAtomCount(); i++){
			IAtom atom_one = fragmentStructure.getAtom(i);					
			for(int j = i+1; j < fragmentStructure.getAtomCount(); j++){
				IAtom atom_two = fragmentStructure.getAtom(j);
				//If the both atoms are within the modified structure, then use the bonds within the modified one.
				if(UniqueIDFunctionSet.containAtom(atom_one, modified) && UniqueIDFunctionSet.containAtom(atom_two, modified)){
					IAtom atom_one_temp = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_one), modified);
					IAtom atom_two_temp = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_two), modified); 
					IBond oneBond = modified.getBond(atom_one_temp, atom_two_temp);		
					if(oneBond !=null){
						IAtom[] atoms_for_bond = {atom_one, atom_two};
						IBond newBond = new Bond(atoms_for_bond, oneBond.getOrder());
						fragmentStructure.addBond(newBond);
					}					
				}
				//Otherwise, use the bonds within the original reactant molecule structure
				else{
					//If one of the atom is newly added, then it can not be found in the original reactant, so skip it.
					if(!UniqueIDFunctionSet.containAtom(atom_one, reactant) || !UniqueIDFunctionSet.containAtom(atom_two, reactant)) continue;
					IAtom atom_one_temp = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_one), reactant);
					IAtom atom_two_temp = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_two), reactant); 
					IBond oneBond = reactant.getBond(atom_one_temp, atom_two_temp);		
					if(oneBond !=null){
						IAtom[] atoms_for_bond = {atom_one, atom_two};
						IBond newBond = new Bond(atoms_for_bond, oneBond.getOrder());
						fragmentStructure.addBond(newBond);
					}				
				}
			}								
		}
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
		adder.addImplicitHydrogens(fragmentStructure);
		AtomContainerManipulator.suppressHydrogens(fragmentStructure);
		//System.out.println(sg.create(fragmentStructure));
		return fragmentStructure;
		
	}
	/**
	 * This function will check the reactive structure and the changed part (result part), and then
	 * put all the atoms that are in the reactive structure but removed in the changed part
	 * into an arrayList<IAtom>
	 * This function will return that arrayList
	 * @param reactiveSites
	 * @param changedPart
	 * @return
	 * @throws Excpetion
	 */
	public static ArrayList<IAtom> getRemovedAtomList(IAtomContainer reactiveSites, IAtomContainer changedPart) throws Exception{
		ArrayList<IAtom> removedAtomList = new ArrayList<>();
		for(int i = 0; i < reactiveSites.getAtomCount(); i++){
			IAtom atom_one = reactiveSites.getAtom(i);
			//If that atom is within the reactiveSites, it's not hydrogen and not included in the changedPart, then it's the atom that is removed
			if(!atom_one.getSymbol().equals("H") && !UniqueIDFunctionSet.containAtom(atom_one, changedPart)){
				IAtom atom_one_copy = atom_one.clone();
				removedAtomList.add(atom_one_copy);
			}
		}
		return removedAtomList;
	}
	/**
	 * This function will check if the fragment matches SMIRKS string
	 * return true if matches
	 * return false otherwise
	 * @param smirks
	 * @param fragment
	 * @return
	 * @throws Exception
	 */
	public static boolean matchSMIRKS(String smirks, IAtomContainer fragment) throws Exception{
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		Pattern smp = SmartsPattern.create(smirks, bldr);	
		boolean match = smp.matches(fragment);
		return match;
	}
	/**
	 * This function will find all atoms that match the SMIRKS query and return them as a ArrayList<IAtom>
	 * @param smirks
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<IAtom> findMatchedAtoms(String smirks, IAtomContainer oneMole) throws Exception{
		ArrayList<IAtom> matchedAtoms = new ArrayList<>();
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		Pattern smp = SmartsPattern.create(smirks, bldr);	
		int[] matchedAtomsIdx = smp.match(oneMole);
		IAtom oneAtom;
		for(int i = 0; i < matchedAtomsIdx.length; i++){
			oneAtom = oneMole.getAtom(matchedAtomsIdx[i]);
			matchedAtoms.add(oneAtom);
		}
		return matchedAtoms;
	}
	/**
	 * This function will find each list of atoms that matches the SMIRKS within the molecule
	 * Those atoms will the stored into a ArrayList<IAtom>
	 * It returns the ArrayList<ArrayList<IAtom>> of all such lists  
	 * @param smirks
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<ArrayList<IAtom>> findMatchedList(String smirks, IAtomContainer oneMole)  throws Exception{
		ArrayList<ArrayList<IAtom>> result = new ArrayList<>();
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		SMARTSQueryTool smartsPattern = new SMARTSQueryTool(smirks,builder); //Set the SMART Pattern
		//Pattern smp = SmartsPattern.create(smirks, builder);	
		boolean matched = smartsPattern.matches(oneMole);
		List<List<Integer>> allMatches = smartsPattern.getUniqueMatchingAtoms();
		IAtom oneAtom;
		ArrayList<IAtom> oneCandidateList = new ArrayList<>();
	
		for(List<Integer> mapping : allMatches){
			for(int i = 0; i < mapping.size(); i++){
				oneAtom = oneMole.getAtom(mapping.get(i));
				oneCandidateList.add(oneAtom);
			}
			result.add(oneCandidateList);
			oneCandidateList = null;
			oneAtom = null;
			return result;
		}
		return result;
	}
	
	
}
