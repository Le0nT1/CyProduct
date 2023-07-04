package typeOneFeatures;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.SybylAtomTypeMatcher;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class GenerateBondFeatures_TypeOne {	
	public static String generateCurrentBondAtomDescriptor(IAtomContainer molecule, IBond oneBond) throws Exception{
		List<IAtom> atoms = reorderAtomsOfBond(molecule, oneBond);
		String[] firstAtomDes = GenerateAtomFeatures_ForBond_TypeOne.generateAtomDescriptorFeatures(molecule,atoms.get(0));
		String[] secondAtomDes = GenerateAtomFeatures_ForBond_TypeOne.generateAtomDescriptorFeatures(molecule,atoms.get(1));
		String[] temp = new String[firstAtomDes.length + secondAtomDes.length];
		for(int i = 0; i < firstAtomDes.length; i++){
			temp[i] = firstAtomDes[i];
		}
		for(int i = 0; i< secondAtomDes.length; i++){
			temp[i+firstAtomDes.length] = secondAtomDes[i];
		}
		StringBuilder builder = new StringBuilder();
		for(int i = 0; i < temp.length; i++){
			if(i == temp.length-1){
				builder.append(temp[i]);
				break;
			}
			builder.append(temp[i] + ",");
		}
		return builder.toString();
		
	}
	/**
	 * This function is used to generate the a String consists of "0" and "1" that indicates the atom types of atoms of a <i,j> in the molecule.
	 * 
	 * @param molecule
	 * @return a 24-tuple binary string. Each tuple is for a AtomType in the AtomTypeLookupTable. The value of every tuple is either "1" or "0".
	 * @throws Exception 
	 */
	public static String generateCurrentBondAtomType(IAtomContainer molecule, IBond oneBond) throws Exception {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(bldr);
		
		IAtom leftAtom = oneBond.getAtom(0);
		IAtom rightAtom = oneBond.getAtom(1);
		String[] atomTypeLookupTable = GenerateNeighborFeatures_TypeOne.atomTypeLookupTable;
		String[] resultStrArray = new String[atomTypeLookupTable.length];
		//Initialize all 24 bits as "0"
		for(int j = 0; j < resultStrArray.length; j++){
			resultStrArray[j] = "0";
		}
		IAtomType atomType = typeMatcher.findMatchingAtomType(molecule, leftAtom); 
		String type = leftAtom.getSymbol();
		try{
			type = atomType.getAtomTypeName();	
		}catch (Exception e) {
			
		}
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(type.equalsIgnoreCase(atomTypeLookupTable[i])){
				resultStrArray[i] = "1";
				break;
			}
		}
		type = rightAtom.getSymbol();
		try {
			atomType = typeMatcher.findMatchingAtomType(molecule, rightAtom); 
			type = atomType.getAtomTypeName();	
		}catch(Exception e) {
			
		}
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(type.equalsIgnoreCase(atomTypeLookupTable[i])){
				resultStrArray[i] = "1";
				break;
			}
		}
					
		StringBuffer resultStrBuffer = new StringBuffer();
		for(int k = 0; k <resultStrArray.length; k++){
			if(k==0){
				resultStrBuffer.append(resultStrArray[0]);
			}
			else{
				resultStrBuffer = resultStrBuffer.append(",").append(resultStrArray[k]);
			}
		}
		return resultStrBuffer.toString();

	}
	/**
	 * This function generates fingerprint for any given bond in the given molecule.
	 * For example, bond <i, j> matches a substructure if both atom i and j are within the substructure
	 * @param mole
	 * @param bond
	 * @return
	 * @throws Exception
	 */
	public static String generateBondFingeprint(IAtomContainer molecule, IBond bond) throws Exception {
		IAtomContainer mole = molecule;
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mole);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mole.getBuilder());		
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(mole);
		adder.addImplicitHydrogens(mole);
		String[] matched = new String[FingerPrintQueries_TypeOne.queriesList.size()];
		for(int i = 0; i < matched.length; i++){
			matched[i] = "0";
		}
		/*
		 * Get the two atoms of the bond
		 */
		IAtom leftAtom = bond.getAtom(0);
		IAtom rightAtom = bond.getAtom(1);
		
		Integer leftAtomIdx = mole.indexOf(leftAtom);
		Integer rightAtomIdx = mole.indexOf(rightAtom);
		//Iterate through all queries in order
		Iterator<Entry<String, String>> allQueries = FingerPrintQueries_TypeOne.queriesList.entrySet().iterator();
		int counter = 0; //Track the index of a query
		while(allQueries.hasNext()) { // Exploit all querries
			Entry<String,String> onePattern = allQueries.next();
			IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(onePattern.getValue(),builder); //Set the SMART Pattern
			boolean occurrences = smartsPattern.matches(mole);//match the pattern with the molecule
			List<List<Integer>> matchingAtoms = smartsPattern.getUniqueMatchingAtoms();//Get the atoms that matches this pattern
			for(int i = 0; i < matchingAtoms.size(); i++){
				List<Integer> oneMatchedSet = matchingAtoms.get(i);
				if(oneMatchedSet.contains(leftAtomIdx)&&oneMatchedSet.contains(rightAtomIdx)){//If the two atoms of the bond are both in the set of this pattern
					matched[counter] = "1";		
					break;
				}				
			}
			//System.out.println(onePattern.getKey());
			counter++;
			
		}
		StringBuilder tempMatched = new StringBuilder();
		for(int i = 0; i < matched.length; i++){
			if(i==matched.length-1){
				tempMatched.append(matched[i]);
				break;
			}
			tempMatched.append(matched[i] + ",");
		}		
		return tempMatched.toString();
	}
	/**
	 * This function returns a List of IAtoms that follows a specific order. The one with smaller atomic number comes first.
	 * @param oneBond
	 * @return
	 * @throws Exception
	 */
	public static List<IAtom> reorderAtomsOfBond(IAtomContainer molecule, IBond oneBond) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		List<IAtom> reordered = new ArrayList<>();
		IAtom leftAtom = oneBond.getAtom(0);
		IAtom rightAtom = oneBond.getAtom(1);
		Integer leftNumber = leftAtom.getAtomicNumber();
		Integer rightNumber = rightAtom.getAtomicNumber();
		ArrayList<IAtom> checkVisitedAtoms = new ArrayList<>();
		checkVisitedAtoms.add(rightAtom);
		int length_left = findTotalConnected(molecule, leftAtom, checkVisitedAtoms)-1;
		checkVisitedAtoms = new ArrayList<>();
		checkVisitedAtoms.add(leftAtom);
		int length_right = findTotalConnected(molecule, rightAtom, checkVisitedAtoms)-1;
		if(leftNumber == rightNumber){
			if(length_left == length_right){
				String[] firstAtomDes = GenerateAtomFeatures_ForBond_TypeOne.generateAtomDescriptorFeatures(molecule,leftAtom);
				String[] secondAtomDes = GenerateAtomFeatures_ForBond_TypeOne.generateAtomDescriptorFeatures(molecule,rightAtom);
				boolean leftFirst = leftFirst(firstAtomDes,secondAtomDes);
				if(leftFirst){
					reordered.add(leftAtom);
					reordered.add(rightAtom);
				}
				else{
					reordered.add(rightAtom);
					reordered.add(leftAtom);
				}
			}
			else if(length_left < length_right){
				reordered.add(leftAtom);
				reordered.add(rightAtom);
			}
			else{
				reordered.add(rightAtom);
				reordered.add(leftAtom);
			}
		}
		else if(leftNumber < rightNumber){
			reordered.add(leftAtom);
			reordered.add(rightAtom);
		}
		else{
			reordered.add(rightAtom);
			reordered.add(leftAtom);
		}
		return reordered;
	}
	/**
	 * This function is the wrapper of findLongestPath. 
	 * checkVisitedAtoms need to be reinitialize here
	 * This function returns the length longest path
	 */
	public static int[] getLongestPath(IAtomContainer molecule, IBond oneBond) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		List<IAtom> orderedAtoms = reorderAtomsOfBond(molecule, oneBond);
		IAtom leftAtom = orderedAtoms.get(0);
		IAtom rightAtom = orderedAtoms.get(1);
		/**
		 * Get the longest path from the left node
		 */
		//checkVisitedAtoms.add(leftAtom);
		ArrayList<IAtom> checkVisitedAtoms = new ArrayList<>();
		checkVisitedAtoms.add(rightAtom);
		Integer length_left = findTotalConnected(molecule, leftAtom,checkVisitedAtoms) -1;
		/**
		 * Get the longest path from the right node
		 */
		checkVisitedAtoms = new ArrayList<>();
		checkVisitedAtoms.add(leftAtom);
		//checkVisitedAtoms.add(rightAtom);
		Integer length_right = findTotalConnected(molecule, rightAtom,checkVisitedAtoms) -1;		
		int[] result = {length_left, length_right}; //Check Point

		return result;
	}
	/**
	 * Use recursive call to find the longest path. Circles are only counted once.
	 * Need to initialize(empty) the checking ArrayList every time this function is called.
	 * Better to wrap this function by another one
	 * @param molecule
	 * @param oneBond
	 * @return
	 * @throws Exception
	 */
	public static int findTotalConnected(IAtomContainer molecule, IAtom center,ArrayList<IAtom> checkVisitedAtoms) throws Exception{
		int count = 0;
		List<IAtom> neighbors = molecule.getConnectedAtomsList(center);
		for(int i = 0; i < neighbors.size(); i++){
			if(checkVisitedAtoms.contains(neighbors.get(i))){
				neighbors.remove(i);
			}
		}
		if(!checkVisitedAtoms.contains(center)){//if it's a neighbor that has not been visited before, add length by 1
			count++;
			checkVisitedAtoms.add(center);
		}
		else{
			return count;
		}
		if(neighbors.isEmpty()){//If all atoms in the neighbors are checked, terminate
			return count;
		}
		for(int i = 0; i < neighbors.size(); i++){

			int stackCount = findTotalConnected(molecule, neighbors.get(i),checkVisitedAtoms);
			count = stackCount + count;
		}
		return count;
		
	}
	/**
	 * Order the atom connected by bond using descriptors' value
	 */
	public static boolean leftFirst(String[] firstAtom, String[] secondAtom) throws Exception{
		boolean leftFirst = true;
		for(int i = 0; i < firstAtom.length; i++){
			if(Double.parseDouble(firstAtom[i]) == Double.parseDouble(secondAtom[i])){
				continue;
			}
			if(Double.parseDouble(firstAtom[i]) < Double.parseDouble(secondAtom[i])){
				leftFirst = true;
				break;
			}
			else{
				leftFirst = false;
				break;
			}
		}
		return leftFirst;
	}
}
