package typeOneFeatures;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.SybylAtomTypeMatcher;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class GenerateNeighborFeatures_TypeOne {
	
	/**
	 * definition for all the adopted SYBYL atom types chosen by us
	 * 23 in total
	 */
	public static final String[] atomTypeLookupTable = { "C.1", "C.2", "C.3", "C.ar", "C.cat", "N.1",
			"N.2", "N.3", "N.4", "N.ar", "N.am", "N.pl3", "O.2", "O.3", "O.co2", "S.2",
			"S.3", "S.O", "S.O2", "P.3", "F", "Cl", "Br", "I" };
	/*
	 * 14 self defined atomTypes
	 */
	public static final String[] roughAtomTypeTable = {"C", "C.ar", "N", "N.ar", "O", "O.co2", "S", "S.O", "S.O2", "P", "F", "Cl", "Br", "I"};
	/**
	 * This function convert the neighborAtomDescriptorMap to a String.
	 * @param molecule
	 * @param oneBond
	 * @param upToDepth
	 * @return
	 * @throws Exception
	 */
	public static String generateNeighborAtomDescriptor(IAtomContainer molecule, IBond oneBond, int upToDepth) throws Exception{
		StringBuilder resultBuilder = new StringBuilder();
		for(int depth = 1; depth < (upToDepth+1); depth++){
			HashMap<String, String[]> neighborDescriptors = generateNeighborAtomDescriptorMap(molecule, oneBond, depth);			
			for(int i = 0; i < roughAtomTypeTable.length; i++){
				String atomType = roughAtomTypeTable[i];
				String[] tempList = neighborDescriptors.get(atomType);
				for(int j = 0; j < tempList.length; j++){
					if((depth == upToDepth)&&(i == (roughAtomTypeTable.length-1)&& (j == (tempList.length-1)))){//Last depth, last atomType, last descriptor
						resultBuilder.append(tempList[j]);
						break;
					}
					resultBuilder.append(tempList[j] + ",");
				}
			}
		}
		return resultBuilder.toString();
	}
	
	/**
	 * Return the HashMap<AtomType,DescriptorValueList> storing the neighAtomEnvironment of a bond within a molecule, at a certain depth
	 * @param molecule
	 * @param oneBond
	 * @param depth
	 * @return
	 * @throws Exception
	 */
	public static HashMap<String, String[]> generateNeighborAtomDescriptorMap(IAtomContainer molecule, IBond oneBond, int depth) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(bldr);
		
		ArrayList<IAtom> neighborAtoms = getAtomsAtRadius(oneBond,molecule, depth);
		HashMap<String,String[]> atomTypeValueMap = new HashMap<>();
		HashMap<String,Integer> atomCounterMap = new HashMap<>();
		String[] initStringList = new String[10];//10 molecular descriptors
		//int[] atomsCounterForEachType = new int[roughAtomTypeTable.length];
		for(int i = 0; i < initStringList.length; i++){
			initStringList[i] = "0.0";
			//atomsCounterForEachType[i] = 0;
		}
		for(int i = 0; i < roughAtomTypeTable.length; i++){//Initialize the hashmap
			atomTypeValueMap.put(roughAtomTypeTable[i], initStringList);
			atomCounterMap.put(roughAtomTypeTable[i], 0);
		}
		/**
		 * Get the total value of each descriptor for each atomType.
		 */
		for(int i = 0; i < neighborAtoms.size(); i++){
			IAtom oneAtom = neighborAtoms.get(i);
			IAtomType atomType = typeMatcher.findMatchingAtomType(molecule, oneAtom); 
			String type = oneAtom.getSymbol();
			try {
				type = atomType.getAtomTypeName();
			}catch (Exception e) {
				
			}
			String roughType;
			if(type.equalsIgnoreCase("Any")){
				type = oneAtom.getSymbol();
				roughType = type;
			}
			else roughType = getConvertType(type);
			String[] newValue = GenerateAtomFeatures_ForBond_TypeOne.generateAtomDescriptorFeatures(molecule, oneAtom);			
			String[] currentValue = atomTypeValueMap.get(roughType);
			String[] updatedValue = updateTypeValueList(currentValue, newValue);
			atomTypeValueMap.put(roughType, updatedValue);
			
			Integer counter = atomCounterMap.get(roughType);
			atomCounterMap.put(roughType, (counter+1));
			
			
		}
		/**
		 * Get the average value of each descriptor for each atomType
		 */
		for(int i = 0; i < roughAtomTypeTable.length; i++){
			String roughAtomType = roughAtomTypeTable[i];
			Integer counter = atomCounterMap.get(roughAtomType);
			String[] descriptorValueList = atomTypeValueMap.get(roughAtomType);
			for(int j = 0; j < descriptorValueList.length; j++){
				Double currentValue = Double.parseDouble(descriptorValueList[j]);
				Double avgValue = currentValue/counter;
				if(avgValue < 0.0001||counter==0){
					avgValue = 0.0;
				}
				descriptorValueList[j] = Double.toString(avgValue);
				
			}
			atomTypeValueMap.put(roughAtomType, descriptorValueList);
		}		
		return atomTypeValueMap;
	}
	/**
	 * This function generate the string of atom types for the neighbor atoms up to a certain depth
	 * @param molecule
	 * @param oneBond
	 * @param upToDepth
	 * @return
	 * @throws Exception
	 */
	public static String generateNeighborAtomType_upToDepth(IAtomContainer molecule, IBond oneBond, int upToDepth) throws Exception{
		StringBuilder tempBuilder = new StringBuilder();
		for(int d = 1; d < upToDepth+1; d++){
			String oneResult = generateNeighborAtomType(molecule, oneBond, d);
			if(d==upToDepth){
				tempBuilder.append(oneResult);
				break;
			}			
			tempBuilder.append(oneResult + ",");
		}
		return tempBuilder.toString();
	}
	/**
	 * This function is used to generate the a String consists of "0" and "1" that reflects the atom type for a atom in the molecule.
	 * 
	 * @param molecule
	 * @return a 24-tuple binary string. Each tuple is for a AtomType in the AtomTypeLookupTable. The value of every tuple is either "1" or "0".
	 * @throws Exception 
	 */
	public static String generateNeighborAtomType(IAtomContainer molecule, IBond oneBond, int depth) throws Exception {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(bldr);
		
		ArrayList<IAtom> neighborAtoms = getAtomsAtRadius(oneBond,molecule, depth);

		String[] resultStrArray = new String[atomTypeLookupTable.length];
		//Initialize all 24 bits as "0"
		for(int j = 0; j < resultStrArray.length; j++){
			resultStrArray[j] = "0";
		}
		for(int atomIdx = 0; atomIdx < neighborAtoms.size(); atomIdx++){
			IAtom oneAtom = neighborAtoms.get(atomIdx);
			IAtomType atomType = typeMatcher.findMatchingAtomType(molecule, oneAtom); 
			String type = oneAtom.getSymbol();
			try{
				type = atomType.getAtomTypeName();	
			}catch (Exception e) {
				
			}
			for(int i = 0; i < atomTypeLookupTable.length; i++){
				if(type.equalsIgnoreCase(atomTypeLookupTable[i])){
					Integer current = Integer.parseInt(resultStrArray[i]);
					resultStrArray[i] = Integer.toString(current+1);
				}
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
	 * This function returns the list contains all neighbor atoms of a <i,j> bond by searching to the left/right from i/j by a certain depth
	 * @param oneBond
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<IAtom> getAtomsAtRadius(IBond oneBond, IAtomContainer molecule, int depth) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		ArrayList<IAtom> resultList = new ArrayList<>();
		IAtom leftAtom = oneBond.getAtom(0);
		IAtom rightAtom = oneBond.getAtom(1);
		List<List<IAtom>> atomOnPath_left = PathTools.getPathsOfLength(molecule, leftAtom, depth);
		List<List<IAtom>> atomOnPath_right = PathTools.getPathsOfLength(molecule, rightAtom, depth);
		/**
		 * check Left then Right
		 */
		for(int pathIdx = 0; pathIdx < atomOnPath_left.size(); pathIdx++){
			List<IAtom> path = atomOnPath_left.get(pathIdx);
			if(path.contains(rightAtom)){
				continue;
			}
			IAtom targetAtom = path.get(path.size()-1);//path.size = depth - 1
			if(!resultList.contains(targetAtom)){
				resultList.add(targetAtom);
			}
		}
		for(int pathIdx = 0; pathIdx < atomOnPath_right.size(); pathIdx++){
			List<IAtom> path = atomOnPath_right.get(pathIdx);
			if(path.contains(leftAtom)){
				continue;
			}
			IAtom targetAtom = path.get(path.size()-1);//path.size = depth - 1
			if(!resultList.contains(targetAtom)){
				resultList.add(targetAtom);
			}
		}
		//resultList.remove(leftAtom);
		//resultList.remove(rightAtom);
		return resultList;
	}
	/**
	 * Add up the corresponding currentList and newList and return the updated List
	 * @param currentList
	 * @param newList
	 * @return
	 */
	public static String[] updateTypeValueList(String[] currentList, String[] newList) {
		String[] updatedList = new String[currentList.length];
		for(int i = 0; i < updatedList.length; i++){
			Double oneValue = Double.parseDouble(currentList[i]) + Double.parseDouble(newList[i]);
			updatedList[i] = Double.toString(oneValue);
		}
		return updatedList;
	}
	/**
	 * Convert the exact atomType to the rough one. The rough atomType List is:
	 * {"C", "C.ar", "N", "N.ar", "O", "O.co2", "S", "S.O", "S.O2", "P", "F", "Cl", "Br", "I"};
	 * @param atomType
	 * @return
	 */
	public static String getConvertType(String atomType) throws Exception{
		String result;
		if(atomType.contains("C.") || atomType.equalsIgnoreCase("C")){
			if(atomType.equals("C.ar")){
				result = "C.ar";
			}
			else{
				result = "C";
			}
		}
		else if(atomType.contains("N.") || atomType.equalsIgnoreCase("N")){
			if(atomType.equals("N.ar")){
				result = "N.ar";
			}
			else{
				result = "N";
			}
		}
		else if(atomType.contains("O.") || atomType.equalsIgnoreCase("O")){
			if(atomType.equals("O.co2")){
				result = "O.co2";
			}
			else{
				result = "O";
			}
			
		}
		else if(atomType.contains("S.") || atomType.equalsIgnoreCase("S")){
			if(atomType.equals("S.O")){
				result = "S.O";
			}
			else if(atomType.equals("S.O2")){
				result = "S.O2";
			}
			else{
				result = "S";
			}
		}
		else if(atomType.contains("P.") || atomType.equalsIgnoreCase("P")){
			result = "P";
		}
		else if(atomType.contains("F")){
			result = "F";
		}
		else if(atomType.contains("Cl")){
			result = "Cl";
		}
		else if(atomType.contains("Br")){
			result = "Br";
		}
		else if(atomType.contains("I")){
			result = "I";
		}
		else{
			result = "Unknown";
			throw new Exception("The atomType is not included in the rough atomType List.");
		}
		
		return result;
	}
	
	
	
}
