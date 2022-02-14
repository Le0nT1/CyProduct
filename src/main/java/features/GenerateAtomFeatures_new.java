package features;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.SybylAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.NoSuchAtomTypeException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.atomic.AtomDegreeDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.AtomHybridizationDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.AtomValenceDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.EffectiveAtomPolarizabilityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PartialSigmaChargeDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PartialTChargeMMFF94Descriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PiElectronegativityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.SigmaElectronegativityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.StabilizationPlusChargeDescriptor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import utils.DistanceRelatedFunciton;



public class GenerateAtomFeatures_new {
	DistanceRelatedFunciton drf = new DistanceRelatedFunciton();
	/**
	 * definition for all the adopted SYBYL atom types chosen by us
	 * 23 in total
	 */
	public final String[] atomTypeLookupTable = { "C.1", "C.2", "C.3", "C.ar", "C.cat", "N.1",
			"N.2", "N.3", "N.4", "N.ar", "N.am", "N.pl3", "O.2", "O.3", "O.co2", "S.2",
			"S.3", "S.O", "S.O2", "P.3", "F", "Cl", "Br", "I" };
	
	/**
	 * This function will generate the "merged" fingerprint for the "numClosestAtoms many closest atoms" for the IAtom oneAtom
	 * The numClosestAtoms is determined by the 3D-coordinates
	 * @param molecule
	 * @param oneAtom
	 * @param numClosestAtoms
	 * @return
	 * @throws Exception
	 */
	public String generateClosestAtomsFPs(IAtomContainer molecule, IAtom oneAtom, Integer numClosestAtoms) throws Exception{
		ArrayList<IAtom> closestAtoms;
		if(oneAtom.getPoint3d() != null){
			closestAtoms = this.drf.getClosestAtomListInMolecule_3D(molecule, oneAtom, numClosestAtoms);
		}
		else closestAtoms = this.drf.getClosestAtomListInMolecule_NoCoordinates(molecule, oneAtom, numClosestAtoms);
		IAtomContainer mole = molecule;
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mole);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mole.getBuilder());		
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(mole);
		adder.addImplicitHydrogens(mole);
		String[] matched = new String[FingerPrintQueries_new.newQueriesList.size()];
		for(int i = 0; i < matched.length; i++){
			matched[i] = "0";
		}
		for(int j = 0; j < closestAtoms.size(); j++){
			Integer atomIdx = molecule.indexOf(closestAtoms.get(j));
			//Iterate through all queries in order
			Iterator<Entry<String, String>> allQueries = FingerPrintQueries_new.newQueriesList.entrySet().iterator();
			int counter = 0; //Track the index of a query
			while(allQueries.hasNext()) { // Exploit all querries
				Entry<String,String> onePattern = allQueries.next();
				IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
				SMARTSQueryTool smartsPattern = new SMARTSQueryTool(onePattern.getValue(),builder); //Set the SMART Pattern
				boolean occurrences = smartsPattern.matches(mole);//match the pattern with the molecule
				List<List<Integer>> matchingAtoms = smartsPattern.getUniqueMatchingAtoms();//Get the atoms that matches this pattern
				for(int i = 0; i < matchingAtoms.size(); i++){
					List<Integer> oneMatchedSet = matchingAtoms.get(i);
					if(oneMatchedSet.contains(atomIdx)){
						matched[counter] = "1";//String.valueOf(Integer.parseInt(matched[counter]) + 1);//"1";		
						break;
					}				
				}
				counter++;			
			}
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
	 * This function will generate the atom types features the numClosestAtoms for the oneAtom within the molecule
	 * @param oneMole
	 * @param oneAtom
	 * @param numClosestAtoms
	 * @return
	 * @throws Exception
	 */
	public String generateClosestAtomsAtomTypes(IAtomContainer molecule, IAtom oneAtom, Integer numClosestAtoms) throws Exception{
		ArrayList<IAtom> closestAtoms;
		if(oneAtom.getPoint3d() != null){
			closestAtoms = this.drf.getClosestAtomListInMolecule_3D(molecule, oneAtom, numClosestAtoms);
		}
		else closestAtoms = this.drf.getClosestAtomListInMolecule_NoCoordinates(molecule, oneAtom, numClosestAtoms);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(bldr);
		String[] resultStrArray = new String[atomTypeLookupTable.length];	
		//Initialize all 24 bits as "0"
		for(int j = 0; j < resultStrArray.length; j++){
			resultStrArray[j] = "0";
		}
		
		for(int k = 0; k < closestAtoms.size(); k++){
			Integer atomIdx = molecule.indexOf(closestAtoms.get(k));
			IAtomType[] atomTypeList = typeMatcher.findMatchingAtomTypes(molecule); 
			IAtomType atomType = atomTypeList[atomIdx];
			String type = atomType.getAtomTypeName();
			//Find the matched atomType and update the value to current + "1"
			for(int i = 0; i < atomTypeLookupTable.length; i++){
				if(type.equalsIgnoreCase(atomTypeLookupTable[i])){
					resultStrArray[i] = String.valueOf(Integer.parseInt(resultStrArray[i]) + 1);// "1";
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
	 * This function generates fingerprint for any given bond in the given molecule.
	 * For example, bond <i, j> matches a substructure if both atom i and j are within the substructure
	 * @param mole
	 * @param bond
	 * @return
	 * @throws Exception
	 */
	public String generateAtomFingeprint(IAtomContainer molecule, IAtom oneAtom) throws Exception {
		IAtomContainer mole = molecule;
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mole);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mole.getBuilder());		
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(mole);
		adder.addImplicitHydrogens(mole);
		String[] matched = new String[FingerPrintQueries_new.newQueriesList.size()];
		for(int i = 0; i < matched.length; i++){
			matched[i] = "0";
		}
		/*
		 * Get the two atoms of the bond
		 */
		Integer atomIdx = molecule.indexOf(oneAtom);
		//Iterate through all queries in order
		Iterator<Entry<String, String>> allQueries = FingerPrintQueries_new.newQueriesList.entrySet().iterator();
		int counter = 0; //Track the index of a query
		while(allQueries.hasNext()) { // Exploit all querries
			Entry<String,String> onePattern = allQueries.next();
			IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(onePattern.getValue(),builder); //Set the SMART Pattern
			boolean occurrences = smartsPattern.matches(mole);//match the pattern with the molecule
			List<List<Integer>> matchingAtoms = smartsPattern.getUniqueMatchingAtoms();//Get the atoms that matches this pattern
			for(int i = 0; i < matchingAtoms.size(); i++){
				List<Integer> oneMatchedSet = matchingAtoms.get(i);
				if(oneMatchedSet.contains(atomIdx)){
					matched[counter] = "1";		
					break;
				}				
			}
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
	 * This function is used to generate the a String consists of "0" and "1" that reflects the atom type for a atom in the molecule.
	 * 
	 * @param molecule
	 * @return a 24-tuple binary string. Each tuple is for a AtomType in the AtomTypeLookupTable. The value of every tuple is either "1" or "0".
	 * 
	 * @throws CDKException
	 */
	public String generateAtomType(IAtomContainer molecule, IAtom atom) throws CDKException {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(bldr);
		//IAtomContainer mol = molecule;
		//Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		//aromaticity.apply(mol);
		//IAtomType atomType = typeMatcher.findMatchingAtomType(molecule, atom);
		Integer atomIdx = molecule.indexOf(atom);
		IAtomType[] atomTypeList = typeMatcher.findMatchingAtomTypes(molecule); 
		IAtomType atomType = atomTypeList[atomIdx];
		String type = atomType.getAtomTypeName();
		String[] resultStrArray = new String[atomTypeLookupTable.length];
		//Initialize all 24 bits as "0"
		for(int j = 0; j < resultStrArray.length; j++){
			resultStrArray[j] = "0";
		}
		//Find the matched atomType and update the value to "1"
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(type.equalsIgnoreCase(atomTypeLookupTable[i])){
				resultStrArray[i] = "1";
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
	 * This function generates atomic features for a atom in the molecule
	 * @param molecule
	 * @param oneAtom
	 * @return String[] whose each tuple is a real value
	 * @throws CDKException
	 * @throws NoSuchAtomTypeException
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 */
	public String[] generateAtomDescriptorFeatures(IAtomContainer molecule, IAtom oneAtom) throws CDKException, NoSuchAtomTypeException, ClassNotFoundException, IOException {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);	
		adder.addImplicitHydrogens(molecule);				
		AtomDegreeDescriptor degree = new AtomDegreeDescriptor();
		AtomHybridizationDescriptor hy = new AtomHybridizationDescriptor();
		AtomValenceDescriptor va = new AtomValenceDescriptor();
		EffectiveAtomPolarizabilityDescriptor ep = new EffectiveAtomPolarizabilityDescriptor();	
		PartialSigmaChargeDescriptor psc = new PartialSigmaChargeDescriptor();
		//This PartialTChargeMMFF94Descriptor() may generate NaN for some Atoms in some Molecules 
		PartialTChargeMMFF94Descriptor ptc = new PartialTChargeMMFF94Descriptor();
		PiElectronegativityDescriptor pen = new PiElectronegativityDescriptor();
		SigmaElectronegativityDescriptor sen = new SigmaElectronegativityDescriptor();
		StabilizationPlusChargeDescriptor spc = new StabilizationPlusChargeDescriptor();	
//		CovalentRadiusDescriptor covR = new CovalentRadiusDescriptor();
//		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
//		covR.initialise(builder);
		
		String[] atomicFeatures =new String[9];		
		DescriptorValue d = degree.calculate(oneAtom, molecule);
		atomicFeatures[0] = d.getValue().toString();
		DescriptorValue h = hy.calculate(oneAtom, molecule);
		atomicFeatures[1] = h.getValue().toString();
		DescriptorValue v = va.calculate(oneAtom, molecule);
		atomicFeatures[2] = v.getValue().toString();
		DescriptorValue e = ep.calculate(oneAtom, molecule);
		atomicFeatures[3] = e.getValue().toString();
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		DescriptorValue ps = psc.calculate(oneAtom, molecule);
		atomicFeatures[4] = ps.getValue().toString();
		DescriptorValue pt = ptc.calculate(oneAtom, molecule);
		atomicFeatures[5] = pt.getValue().toString();
		DescriptorValue pe = pen.calculate(oneAtom, molecule);
		atomicFeatures[6]= pe.getValue().toString();
		DescriptorValue se = sen.calculate(oneAtom, molecule);
		atomicFeatures[7] = se.getValue().toString();
		//System.out.println(res[7]);
		DescriptorValue sp = spc.calculate(oneAtom, molecule);
		atomicFeatures[8] = sp.getValue().toString();
		
//		DescriptorValue cov = covR.calculate(oneAtom, molecule);
//		atomicFeatures[9] = cov.getValue().toString();
		AtomContainerManipulator.suppressHydrogens(molecule);		
		return atomicFeatures;
	}
	
	/**
	 * This function generates bond features for "CH"-ish and <i,S> bonds 
	 * @param molecule
	 * @param oneAtom
	 * @return
	 * @throws Exception
	 */
	public String generateHydrogenBondDescriptor(IAtomContainer molecule, IAtom oneAtom) throws Exception{
		IAtomContainer mole = molecule;
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mole);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mole.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(mole);
		adder.addImplicitHydrogens(mole);		
		//Iterate through all queries in order
		String specialFp;
		Integer numberHydrogen = oneAtom.getImplicitHydrogenCount();
		String atomType = oneAtom.getSymbol();//Here, the symbol is "S","N","C"; Also include "O"
		switch(numberHydrogen){
			case 0:{
				specialFp = atomType;
				break;
			}
			case 1:{
				specialFp = atomType + "H"; //eg CH
				break;
			}
			case 2:{
				specialFp = atomType + "H2"; // eg. can be changed to CH2
				break;
			}
			case 3:{
				specialFp = atomType + "H3"; // eg. can be changed to CH3
				break;
			}
			case 4:{
				specialFp = atomType + "H4"; //eg. Can be change to CH4
				break;
			}
			default:{
				specialFp = "UnkownError";
				System.out.println(specialFp);
				break;
			}							
		}	
		return specialFp;
	}

	
}
