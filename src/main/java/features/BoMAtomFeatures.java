package features;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.SybylAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.NoSuchAtomTypeException;
import org.openscience.cdk.geometry.surface.NumericalSurface;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
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
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class BoMAtomFeatures {
	/**
	 * definition for all the adopted SYBYL atom types chosen by us
	 * 23 in total
	 */
	public static final String[] atomTypeLookupTable = { "C.1", "C.2", "C.3", "C.ar", "C.cat", "N.1",
			"N.2", "N.3", "N.4", "N.ar", "N.am", "N.pl3", "O.2", "O.3", "O.co2", "S.2",
			"S.3", "S.O", "S.O2", "P.3", "F", "Cl", "Br", "I" };
	/**
	 * This function generates atomic features for Type Two and Type Three bond with the relative atomType as special character
	 * @param hyrogen
	 * @return
	 */
	public static String[] generateTypeTwoAndThreeHydrogenCountFeature(String atomType, Integer length){
		String[] result = new String[length];
		String specialChara;
		if(atomType.equals("H")){
			specialChara = "H"; // For "CH"-ish bond
		}
		else{
			specialChara = "None"; // For "Multi-Valence" bond
		}
		for(int i = 0; i < result.length; i++){
			result[i] = specialChara;
		}
		return result;
	}
	
	/**
	 * This function generates bond features for "CH"-ish and <i,S> bonds 
	 * @param molecule
	 * @param oneAtom
	 * @return
	 * @throws Exception
	 */
	public static String generateHydrogenBondDescriptor(IAtomContainer molecule, IAtom oneAtom) throws Exception{
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
	
	/**
	 * This function is used to generate the a String consists of "0" and "1" that indicates the atom type for a atom in the molecule.
	 * 
	 * @param molecule
	 * @return a 24-tuple binary string. Each tuple is for a AtomType in the AtomTypeLookupTable. The value of every tuple is either "1" or "0".
	 * 
	 * @throws CDKException
	 */
	public static String generateVariantAtomType(IAtomContainer molecule, IBond bond) throws CDKException {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		//IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(SilentChemObjectBuilder.getInstance());
		IAtom atomLeft = bond.getAtom(0);
		IAtom atomRight = bond.getAtom(1);
		String atomTypeNameLeft = typeMatcher.findMatchingAtomType(molecule, atomLeft).getAtomTypeName(); 
		String atomTypeNameRight = typeMatcher.findMatchingAtomType(molecule, atomRight).getAtomTypeName();
		
		//Integer atomIdx = molecule.getAtomNumber(atom);
		//IAtomType[] atomTypeList = typeMatcher.findMatchingAtomTypes(molecule); 
		//IAtomType atomType = atomTypeList[atomIdx];
		//String type = atomType.getAtomTypeName();
		String[] resultStrArray = new String[atomTypeLookupTable.length];
		//Initialize all 24 bits as "0"
		for(int j = 0; j < resultStrArray.length; j++){
			resultStrArray[j] = "0";
		}
		//Find the matched atomType and update the value to "1"
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(atomTypeNameLeft.equalsIgnoreCase(atomTypeLookupTable[i])){
				resultStrArray[i] = "1";
			}
		}
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(atomTypeNameRight.equalsIgnoreCase(atomTypeLookupTable[i])){
				if(resultStrArray[i].equals("1")){
					resultStrArray[i] = "2";
				}
				else{
					resultStrArray[i] = "1";
				}
			}
		}
		/**
		 * Convert the resultStrArray to the output string.
		 */
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
	 * Generate the AtomType String for the neighbors of the bond by searching towards left and right by depth
	 * @param molecule
	 * @param neighAtoms
	 * @return
	 * @throws Exception
	 */
	public static String generateNeighborAtomType(IAtomContainer molecule, List<IAtom> neighAtoms) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(bldr);
		String[] resultStrArray = new String[atomTypeLookupTable.length];
		for(int j = 0; j < resultStrArray.length; j++){
			resultStrArray[j] = "0";
		}
		for(int i = 0; i < neighAtoms.size(); i++){
			IAtom atom = neighAtoms.get(i);
			Integer atomIdx = molecule.indexOf(atom);
			IAtomType atomType = typeMatcher.findMatchingAtomType(molecule, atom); 
			String type = atomType.getAtomTypeName();			
			//Initialize all 24 bits as "0"
			//Find the matched atomType and update the value to "1"
			for(int j = 0; j < atomTypeLookupTable.length; j++){
				if(type.equalsIgnoreCase(atomTypeLookupTable[i])){
					resultStrArray[i] = "1";
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
	}	/**
	 * This function is used to generate the a String consists of "0" and "1" that reflects the atom type for a atom in the molecule.
	 * 
	 * @param molecule
	 * @return a 24-tuple binary string. Each tuple is for a AtomType in the AtomTypeLookupTable. The value of every tuple is either "1" or "0".
	 * 
	 * @throws CDKException
	 */
	public static String generateAtomTypeCurrentBond(IAtomContainer molecule, IBond oneBond) throws CDKException {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(bldr);

		IAtom leftAtom = oneBond.getAtom(0);
		IAtom rightAtom = oneBond.getAtom(1);
		
		IAtomType leftAtomType = typeMatcher.findMatchingAtomType(molecule, leftAtom);
		IAtomType rightAtomType = typeMatcher.findMatchingAtomType(molecule, rightAtom);

		String leftType = leftAtomType.getAtomTypeName();
		String rightType = rightAtomType.getAtomTypeName();
		
		String[] resultStrArray = new String[atomTypeLookupTable.length];
		//Initialize all 24 bits as "0"
		for(int j = 0; j < resultStrArray.length; j++){
			resultStrArray[j] = "0";
		}
		//Find the matched atomType and update the value to "1"
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(leftType.equalsIgnoreCase(atomTypeLookupTable[i])){
				resultStrArray[i] = "1";
			}
		}
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(rightType.equalsIgnoreCase(atomTypeLookupTable[i])){
				if(resultStrArray[i].equals("1")){
					resultStrArray[i] = "2";
				}
				else{
					resultStrArray[i] = "1";
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
	 * This function is used to generate the a String consists of "0" and "1" that reflects the atom type for a atom in the molecule.
	 * 
	 * @param molecule
	 * @return a 24-tuple binary string. Each tuple is for a AtomType in the AtomTypeLookupTable. The value of every tuple is either "1" or "0".
	 * 
	 * @throws CDKException
	 */
	public static String generateAtomType(IAtomContainer molecule, IAtom atom) throws CDKException {
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
	 * This function is used to generate the a String consists of "0" and "1" that reflects the atom type for a atom in the molecule.
	 * 
	 * @param molecule
	 * @return a 24-tuple binary string. Each tuple is for a AtomType in the AtomTypeLookupTable. The value of every tuple is either "1" or "0".
	 * 
	 * @throws CDKException
	 */
	public static String generateBondAtomType(IAtomContainer molecule, IBond oneBond) throws CDKException {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher typeMatcher = SybylAtomTypeMatcher.getInstance(bldr);
		
		Integer leftAtomIdx = molecule.indexOf(oneBond.getAtom(0));
		Integer rightAtomIdx = molecule.indexOf(oneBond.getAtom(0));
		IAtomType leftAtomType = typeMatcher.findMatchingAtomType(molecule, oneBond.getAtom(0)); 
		IAtomType rightAtomType = typeMatcher.findMatchingAtomType(molecule, oneBond.getAtom(1));
		String leftType = leftAtomType.getAtomTypeName();
		String rightType = rightAtomType.getAtomTypeName();
		String[] resultStrArray = new String[atomTypeLookupTable.length];
		//Initialize all 24 bits as "0"
		for(int j = 0; j < resultStrArray.length; j++){
			resultStrArray[j] = "0";
		}
		//Find the matched atomType and update the value to "1"
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(leftType.equalsIgnoreCase(atomTypeLookupTable[i])){
				resultStrArray[i] = "1";
			}
		}
		for(int i = 0; i < atomTypeLookupTable.length; i++){
			if(rightType.equalsIgnoreCase(atomTypeLookupTable[i])){
				if(resultStrArray[i].equals("1")){
					resultStrArray[i] = "2";
				}
				else{
					resultStrArray[i] = "1";
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
	 * This function generates atomic features for a atom in the molecule
	 * @param molecule
	 * @param oneAtom
	 * @return String[] whose each tuple is a real value
	 * @throws CDKException
	 * @throws NoSuchAtomTypeException
	 */
	public static String[] generateAtomDescriptorFeatures(IAtomContainer molecule, IAtom oneAtom) throws CDKException, NoSuchAtomTypeException {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
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
		NumericalSurface ns = new NumericalSurface(molecule);
		
		String[] atomicFeatures =new String[10];		
		DescriptorValue d = degree.calculate(oneAtom, molecule);
		atomicFeatures[0] = d.getValue().toString();
		DescriptorValue h = hy.calculate(oneAtom, molecule);
		atomicFeatures[1] = h.getValue().toString();
		DescriptorValue v = va.calculate(oneAtom, molecule);
		atomicFeatures[2] = v.getValue().toString();
		DescriptorValue e = ep.calculate(oneAtom, molecule);
		atomicFeatures[3] = e.getValue().toString();
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
		Integer atomIdx = molecule.indexOf(oneAtom);
		ns.calculateSurface();
		Double asa = ns.getSurfaceArea(atomIdx);
		atomicFeatures[9] = Double.toString(asa);

		return atomicFeatures;
	}

	
	
}
