package instances;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.SybylAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import otherFeatures.BoMAtomFeatures;
import otherFeatures.FingerPrintQueries_new;
import otherFeatures.GenerateAtomFeatures_new;
import otherFeatures.GenerateMolecularFeatures_new;
import utils.Utilities;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

public class CreateTypeTwoThreeInstances {
	public static final String[] cypList = {"1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4"};
	public static final String[]	str	= { "C.1", "C.2", "C.3", "C.ar", "C.cat", "N.1",
			"N.2", "N.3", "N.4", "N.ar", "N.am", "N.pl3", "O.2", "O.3", "O.co2", "S.2",
			"S.3", "S.O", "S.O2", "P.3", "F", "Cl", "Br", "I" };
	public int cypIdx, depth;
	String titleLine;
	ArrayList<String> attributes;
	GenerateAtomFeatures_new gaf_new = new GenerateAtomFeatures_new();
	public CreateTypeTwoThreeInstances(String cyp, ArrayList<String> attributes,int depth) throws Exception{
		loadCypIdx(cyp);
		this.depth = depth;
		this.attributes = attributes;
	}
	public CreateTypeTwoThreeInstances() {
		
	}
	
	public Instances generateTypeTwoInstances_OneMole(ArrayList<String> tempData) throws Exception{
		String titleLine = generateTitle(this.depth);		
		Instances data = createInstances(tempData,  titleLine);
		return data;
	}
	/**
	 * This function takes as input an ArrayList<String> for rawData and a String for titleLine
	 * @param rawData: the values for each eta-eta bond
	 * @param titleLine: the names of all attributes, including names and features, separated by comma 
	 * @return Instances: featured instances for each eta-eta bond within the the molecule. 
	 * @throws Exception
	 */
	public Instances createInstances(ArrayList<String> rawData, String titleLine) throws Exception{
		//Convert nameLsitInstances from "Instances" to "ArrayList" to make it more efficient in both computation and storage
		ArrayList<Attribute> attList = createAttributes(titleLine);		
		Instances bomInstances = new Instances("Rel", attList, 100000); 
		String oneLine;
		for(int i = 0; i < rawData.size(); i++){
			oneLine = rawData.get(i);
			Instance oneInstance = new DenseInstance(attList.size()); 
			oneInstance = createOneInstance(oneLine, oneInstance, attList);			
			bomInstances.add(oneInstance);			
		}
		if(bomInstances.classIndex()==-1){
			bomInstances.setClassIndex(bomInstances.numAttributes()-1);
		}
		/*
		 * Match attributes
		 */
		Instances matchedInstance = matchAttributes(bomInstances,this.attributes);
		//System.out.println("Instances of dataset have been created successfully");
		return matchedInstance;
	}
	
	/**
	 * This function assign values in the dataset file to a Instance.
	 * @param oneLine
	 * @param oneInstance
	 * @return
	 */
	public Instance createOneInstance(String oneLine, Instance oneInstance, ArrayList<Attribute> attList){
		String[] parseLine = oneLine.split(",");
		String label = parseLine[2+this.cypIdx];//need to consider cyp index
		Attribute att_Class = attList.get(attList.size()-1);
		oneInstance.setValue(att_Class, label);//The last attribute is the class attribute. Note that the label is "T" or "F"
		for(int i = 0; i < 11; i++){
			Attribute oneAtt = attList.get(i);
			oneInstance.setValue(oneAtt,parseLine[i]);
		}
		
		for(int i = 11; i < oneInstance.numAttributes()-1; i++){ //Set values for the other n-2 attributes. Assume n is the number of attributes 
			Attribute oneAtt = attList.get(i);//i = 10, i+1
			//System.out.println(oneAtt.toString() + "," + parseLine[i+2]);
			//Att[0] = par[0], Att[1] = par[1], Att[2] = par[3], ... Att[last] = par[2]
			//System.out.println(oneAtt.name());			
			if((!parseLine[i].contains("H")&&!parseLine[i].contains("S")&&!parseLine[i].contains("P")&&!parseLine[i].contains("N"))&&(Double.parseDouble(parseLine[i])< 1.0e-4)){
					parseLine[i] = "0";
					oneInstance.setValue(oneAtt, Double.parseDouble(parseLine[i]));
			}
			else if(oneAtt.isNumeric()){
				//String tt = parseLine[i];
				//Double db = Double.parseDouble(tt);
				//oneInstance.setValue(oneAtt, db);//parseLine[i]
				oneInstance.setValue(oneAtt, Double.parseDouble(parseLine[i]));
			}
			else oneInstance.setValue(oneAtt, parseLine[i]);
			if(oneInstance.isMissing(oneAtt)){
				String tt = parseLine[i];
				Double db = Double.parseDouble(tt);
				Integer it = Integer.parseInt(tt);
				oneInstance.setValue(oneAtt, db);
				oneInstance.setValue(oneAtt, it);
				//System.out.println("gg");
			}
		}		
		return oneInstance;
		
	}
	/**
	 * This function generates the attribute list for the dataset. The first attribute is the name of the molecule. The last attribute is the class attribute.
	 * @param titleLine: The string contains all the attribute names
	 * @param nameList: The names of molecules. Each name string is a nominal string in the name attribute.
	 * Note that the number of attributes is 486 = 26 + 17*26 + 18
	 * @return
	 */
	public ArrayList<Attribute> createAttributes(String titleLine){
		ArrayList<Attribute> attList = new ArrayList<>();
		String[] parseLine = titleLine.split(",");
		ArrayList<String> numHydrogenList = getNumHydrogenList();
		for(int i = 0; i<11; i++){
			Attribute oneAtt = new Attribute(parseLine[i],  (ArrayList<String>) null);
			attList.add(oneAtt);
		}
		//Feature attributes starts at index = 3.
		for(int i = 11; i < parseLine.length; i++){
			Attribute oneAtt = new Attribute(parseLine[i]);
			if(parseLine[i].equals("atomHydrogenCount")){
				oneAtt = new Attribute(parseLine[i],numHydrogenList);
				attList.add(oneAtt);
				continue;
			}
			attList.add(oneAtt);
		}
		List<String> class_nominal_values = new ArrayList<String>();
		class_nominal_values.add("F");// F = not BoM. 0
		class_nominal_values.add("T");// T = BoM. 1 
		Attribute classAttribute = new Attribute("BoM",class_nominal_values);
		attList.add(classAttribute);
		return attList;
	}
	
	
	public static String generateTitle(int depth) throws Exception{
		StringBuffer titleLineBuffer = new StringBuffer();
		titleLineBuffer.append("Name" + "," + "Bond");
		/**
		 * Real Class Label
		 */
		for(int i = 0; i < cypList.length; i++){
			if(cypList[i].equals("2E1")){
				titleLineBuffer.append("," + "2|E1");
			}
			else{
				titleLineBuffer.append("," + cypList[i]);
			}
		}
		String[] molecularDescriptorTitle = {"ALOGP", "APol", "HBondAcceptorCount", "HBondDonorCount", "MomentOfInertia", "RotatableBondsCount",
				"TPSA", "Weight", "XLogP"};
		for(int i = 0; i < molecularDescriptorTitle.length; i++){
			titleLineBuffer.append("," + molecularDescriptorTitle[i]);
		}
		String atomHydrogenCount = "atomHydrogenCount";
		titleLineBuffer.append("," + atomHydrogenCount);
		String[] atomTypeTitle = BoMAtomFeatures.atomTypeLookupTable;
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + "Current_" + atomTypeTitle[i]);
		}
		String[] atomDescriptorTitle = {"AtomDegree", "AtomHybridization", "AtomValence", "EffectiveAtomPolarizability", "PartialSigmaCharge",
										"PartialTChargeMMFF94", "PiElectronegativity", "SigmaElectronegativity", "StabilizationPlusCharge"};
		for(int i = 0; i < atomDescriptorTitle.length; i++){
			titleLineBuffer.append("," + atomDescriptorTitle[i]);
		}
		for(int i = 0; i < 28; i++){
			titleLineBuffer.append("," + "d=1_environmentFeature_" + i);			
		}
		for(int i = 0; i < 28; i++){
			titleLineBuffer.append("," + "d=2_environmentFeature_" + i);			
		}
		for(int i = 0; i < 28; i++){
			titleLineBuffer.append("," + "d=3_environmentFeature_" + i);			
		}
		for(int i = 0; i < 28; i++){
			titleLineBuffer.append("," + "d=4_environmentFeature_" + i);			
		}
		for(int i = 0; i < FingerPrintQueries_new.newQueriesList.size(); i++){
			titleLineBuffer.append("," + "closestAtomFP_" + i);		
		}
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + "closestAtomAtomType_" + atomTypeTitle[i]);
		}
		for(int d = 0; d < depth+1; d++){
			//System.out.println(ChemSearcher.getFingerprintPatterns().size());
			if(d== 0){
				for(int i = 0; i < FingerPrintQueries_new.newQueriesList.size(); i++){//GeneratePatterns.queriesList.size() + 
					titleLineBuffer.append("," + "currentFP_" + i);
				}
			}
			else{
				for(int i = 0; i < FingerPrintQueries_new.newQueriesList.size(); i++){
					titleLineBuffer.append("," + "d=" + d + "FP_" + i);
				}
			}
		}
		return titleLineBuffer.toString();
	}
	
	public ArrayList<String> generateTypeTwoAtomBasedFeatures(IAtomContainer molecule, int depth, boolean mergedFP) throws Exception {		
		HashMap<Integer, String> atomFPMap = generateAtomFPMap(molecule);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);		
		String molecularFeatures = GenerateMolecularFeatures_new.generateMolecularFeatures_new(molecule);
		ArrayList<String> allBonds = new ArrayList<String>();
		String moleName = (String) molecule.getProperties().get("cdk:Title");
		if(moleName==null){
			 SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Absolute);
			 moleName = sg.create(molecule); //			 
		}
		/**
		 * Generate atom features
		 */
		//go through all atoms in the molecule
		for(int i = 0; i < molecule.getAtomCount(); i++){
			/*
			 * Need to process the molecule every time an atom is processed because some operation changed the informatin of the molecule .
			 */
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
			adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
			aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
			aromaticity.apply(molecule);
			adder.addImplicitHydrogens(molecule);
			/**
			 * Generate atom features.
			 */
			IAtom oneAtom = molecule.getAtom(i);
			String atomType = oneAtom.getSymbol(); //Check whether it's a C, B, N, P, S or O atom
			int numberHydrogens = oneAtom.getImplicitHydrogenCount();
			//System.out.println(molecule.getAtomNumber(oneAtom) + ", " + molecule.getAtom(i).getImplicitHydrogenCount());
			if((numberHydrogens<1)|| (atomType.equals("O") && numberHydrogens > 0)){//Skip OH, and only process C-H -ish bons, if it has a Hydrogen attached, it will have hydroxylation instead of type Three
				continue;
			}
			ArrayList<String[]> boms = getBoMs(molecule);
			String bomLabels= realBoMLabels(boms, molecule, oneAtom);
			//If it's a not an Oxygen atom and has more than one Hydrogen atoms attached.
			String atomType_Feature = this.gaf_new.generateAtomType(molecule, oneAtom);
			String atomHydrogenCount = this.gaf_new.generateHydrogenBondDescriptor(molecule,oneAtom);//CH,CH2,CH3, etc
			String[] atomicDescriptorFeatures = this.gaf_new.generateAtomDescriptorFeatures(molecule,oneAtom); // Atomic Features for left Atom	
			String atomicDescriptorFeaturesString = Utilities.convertListToString(atomicDescriptorFeatures);
			String atomEnvironmentDepth_1 = generateEnviromentFeatures(molecule, oneAtom, 1);
			String atomEnvironmentDepth_2 = generateEnviromentFeatures(molecule, oneAtom, 2);
			String atomEnvironmentDepth_3 = generateEnviromentFeatures(molecule, oneAtom, 3);
			String atomEnvironmentDepth_4 = generateEnviromentFeatures(molecule, oneAtom, 4);

			//Here we use depth as the numClosestAtoms.
			String atomClosestAtomFP = this.gaf_new.generateClosestAtomsFPs(molecule, oneAtom, depth);
			String atomClosestAtomType = this.gaf_new.generateClosestAtomsAtomTypes(molecule, oneAtom, depth);

			/**
			 * generateAtomFingerPrints
			 */
			//long start=System.currentTimeMillis();
			String atomFP = generateAtomFingerPrints(molecule, oneAtom, depth, mergedFP, atomFPMap);
			
			//long end=System.currentTimeMillis();
			//System.out.println("Time cost for atom features: " + (end - start));
			
			String atomFeatures = atomHydrogenCount + "," + atomType_Feature + "," + atomicDescriptorFeaturesString + "," + 
							      atomEnvironmentDepth_1 + "," + atomEnvironmentDepth_2 + "," + atomEnvironmentDepth_3 + "," + atomEnvironmentDepth_4 + "," + 
							      //atomFP;
							      atomClosestAtomFP + "," + atomClosestAtomType + "," + atomFP;//Merged atom feature String
			String instanceFinal = moleName + "," + "<" + atomType + "." + (molecule.indexOf(oneAtom) + 1) + ";H" + ">" + ","
								+ bomLabels + "," + molecularFeatures + "," + atomFeatures;// + "," + bondFeatures + "," + mergedFeature_Left;
			allBonds.add(instanceFinal.replace("NaN", "-99.0")); //The bonds of the molecule is added into the array in the order of the bond
		}								
		return allBonds;
	}
	
	/**
	 * Generate bond features for a type Two bond <R,H>. R must has at leaset one Hydrogen atom attached.
	 * Note: 
	 * 1. Each molecule should have been processed by AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mole); 
	 * 2. Each molecule should have been processed by applying the aromatic model as below: 		
	 * 		2.1 Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
	 *		2.2 aromaticity.apply(mol);
	 * 3. Each bond must have two atoms, left and right, and CH4(methane) has 0 bond in this cdk system.
	 * 4. For every bond in the molecule, we will generate One instances for it.
	 * 		4.1	Instance left = currentBondFeature + leftBondFeature + rightBondFeature(leftBondFeatures)
	 * We only generate one instance for each bond in this case
	 * The patterns considered are "C-H", "B-H", "N-H", "S-H" and "P-H"
	 * Note that "O-H" is not considered here.
	 * @param molecule
	 * @param upToDepth is omitted. Could be used in the generateEnviromentFeatures() function later
	 * @return A HashMap<String, ArrayList<String>> ---> HashMap<"Left", Bond {1,2,3,4....}. And as well for "Right"
	 * @throws CDKException
	 * @throws IOException
	 */
	public ArrayList<String> generateTypeThreeAtomBasedFeatures(IAtomContainer molecule, int depth, boolean mergedFP) throws Exception {		
		HashMap<Integer, String> atomFPMap = generateAtomFPMap(molecule);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		String molecularFeatures = GenerateMolecularFeatures_new.generateMolecularFeatures_new(molecule);
		//HashMap<String, ArrayList<String>> allAtomsMap = new HashMap<String, ArrayList<String>>();		
		ArrayList<String> allBonds = new ArrayList<String>();
		String moleName = (String) molecule.getProperties().get("cdk:Title");
		if(moleName==null){
			 SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
			 moleName = sg.create(molecule); //			 
		}
		/**
		 * Generate atom features
		 */
		//go through all atoms in the molecule
		for(int i = 0; i < molecule.getAtomCount(); i++){
			/*
			 * Need to process the molecule every time an atom is processed because some operation changed the informatin of the molecule .
			 */
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
			adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
			aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
			aromaticity.apply(molecule);
			adder.addImplicitHydrogens(molecule);
			/**
			 * Generate atom features.
			 */
			IAtom oneAtom = molecule.getAtom(i);
			String atomType = oneAtom.getSymbol(); //Check whether it's a C, B, N, P, S or O atom
			int numberHydrogens = oneAtom.getImplicitHydrogenCount();
			//System.out.println(molecule.indexOf(oneAtom) + ", " + molecule.getAtom(i).getImplicitHydrogenCount());
//			if(molecule.indexOf(oneAtom) == 17) {
//				System.out.println("Check");
//			}
			if(!atomType.equals("S") && !atomType.equals("P") && !atomType.equals("N")){//Skip OH, and only process C-H -ish bons
				continue;
			}
			ArrayList<String[]> boms = getBoMs(molecule);
			String bomLabels= realBoMLabels(boms, molecule, oneAtom);
			//If it's a not an Oxygen atom and has more than one Hydrogen atoms attached.
			String atomType_Feature = this.gaf_new.generateAtomType(molecule, oneAtom);
			String atomHydrogenCount = this.gaf_new.generateHydrogenBondDescriptor(molecule,oneAtom);//CH,CH2,CH3, etc
			String[] atomicDescriptorFeatures = this.gaf_new.generateAtomDescriptorFeatures(molecule,oneAtom); // Atomic Features for left Atom	
			String atomicDescriptorFeaturesString = Utilities.convertListToString(atomicDescriptorFeatures);
			String atomEnvironmentDepth_1 = generateEnviromentFeatures(molecule, oneAtom, 1);
			String atomEnvironmentDepth_2 = generateEnviromentFeatures(molecule, oneAtom, 2);
			String atomEnvironmentDepth_3 = generateEnviromentFeatures(molecule, oneAtom, 3);
			String atomEnvironmentDepth_4 = generateEnviromentFeatures(molecule, oneAtom, 4);
			//Here we use depth as the numClosestAtoms.
			String atomClosestAtomFP = this.gaf_new.generateClosestAtomsFPs(molecule, oneAtom, depth);
			String atomClosestAtomType = this.gaf_new.generateClosestAtomsAtomTypes(molecule, oneAtom, depth);
			/**
			 * generateAtomFingerPrints
			 */
			String atomFP = generateAtomFingerPrints(molecule, oneAtom, depth, mergedFP, atomFPMap);
			
			String atomFeatures = atomHydrogenCount + "," + atomType_Feature + "," + atomicDescriptorFeaturesString + "," + 
							      atomEnvironmentDepth_1 + "," + atomEnvironmentDepth_2 + "," + atomEnvironmentDepth_3 + "," + atomEnvironmentDepth_4 + "," + 
							      //atomFP;
							      atomClosestAtomFP + "," + atomClosestAtomType + "," + atomFP;//Merged atom feature String
			String instanceFinal = moleName + "," + "<" + (molecule.indexOf(oneAtom) + 1) + ";" + atomType + ">" + ","
								+ bomLabels + "," + molecularFeatures + "," + atomFeatures;// + "," + bondFeatures + "," + mergedFeature_Left;
			allBonds.add(instanceFinal.replace("NaN", "-99.0")); //The bonds of the molecule is added into the array in the order of the bond
		}
									
		return allBonds;
	}
	/**
	 * This function check whether a <C,H>-ish bond or a <i,S>-ish bond is a BoM or not
	 * @param boms
	 * @param molecule
	 * @param leftAtom
	 * @return
	 * @throws Exception
	 */
	public String realBoMLabels(ArrayList<String[]> boms, IAtomContainer molecule, IAtom leftAtom) throws Exception{
		StringBuffer result = new StringBuffer();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		String leftIdx = String.valueOf((molecule.indexOf(leftAtom)+1));//atom index starting with 0
		for(int i = 0; i < boms.size(); i++){
			boolean isBoM = false;
			String[] bomOneCyp = boms.get(i);
			for(int j = 0; j < bomOneCyp.length; j++){
				String[] oneBoM = bomOneCyp[j].split(",");
				if(oneBoM[0].equals(leftIdx)&&oneBoM[1].equals("H")){
					//result.append("1");
					isBoM = true;
					break;
				}
				if(oneBoM[0].equals(leftIdx)&&(oneBoM[1].equals("S")||oneBoM[1].equals("N")||oneBoM[1].equals("P"))){
					//result.append("1");
					isBoM = true;
					break;
				}
			}
			if(i==0){
				if(isBoM){
					result.append("T");
				}
				else{
					result.append("F");
				}
			}
			else{
				if(isBoM){
					result.append(",T");
				}
				else{
					result.append(",F");
				}
			}
		}
		
		return result.toString();
	}
	/**
	 * This function uses a ArrayList contains the BoMs of a molecule for 9 CYPs
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public ArrayList<String[]> getBoMs(IAtomContainer molecule) throws Exception{
		ArrayList<String[]> boms9CYPs = new ArrayList<String[]>();
		for(int i = 0; i < cypList.length; i++){
			String cyp = "BOM_" + cypList[i];
			String bom = (String) molecule.getProperties().get(cyp);
			if(bom==null||bom.equals("")){
				String[] bomTemp = {"None"};
				boms9CYPs.add(bomTemp);
			}
			else{
				String[] bomList = bom.split("\n");
				String[] resultList = new String[bomList.length];
				//StringBuffer bomBuffer = new StringBuffer();
				for(int j = 0; j < bomList.length; j++){
						resultList[j] = bomList[j].split(";")[0].replace("<", "");
						resultList[j] = bomList[j].split(";")[0].replace("<", "");// 7,H   1,2    5,S											
				}
				boms9CYPs.add(resultList);
		
			}
			
		}
		return boms9CYPs;		
	}
	/**
	 * This function is used to generate the environmental features for a atom within the molecule.
	 * @param molecule
	 * @param depth
	 *            : the bond length of an atom to the centered atom.
	 * @return String representation of the environmental features for a single
	 *         molecule where each feature is separated by a comma.
	 *         A (number of atoms) * (24 + 4 tuple vector). Each of the first 24 tuple represents a atom type. Each of the last tuple is either 0 or 1. Tuple(Bond_) 1,2,3,4 means single, double, triple and aromatic 
	 * @throws CDKException
	 */
	public String generateEnviromentFeatures(IAtomContainer molecule, IAtom oneAtom, int depth) throws CDKException {
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		IAtomContainer mol = molecule;
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(mol);
		//int dep = depth;
		SybylAtomTypeMatcher satm = SybylAtomTypeMatcher.getInstance(bldr);
		int length = mol.getAtomCount();
		int o;
		int[] atomType = new int[24];
		String enviromentFeatures;
		for (int typeIdx = 0; typeIdx < 24; typeIdx++){
			atomType[typeIdx] = 0;
		}
		int singleBond = 0;
		int doubleBond = 0;
		int tripleBond = 0;
		int aromaticBond = 0;
		List<List<IAtom>> pathCollection = PathTools.getPathsOfLength(mol, oneAtom, depth);
		List<IAtom> atomsOnOnePath;
		IAtom lastSecondAtom, lastAtom;
		for (int j = 0; j < pathCollection.size(); j++) {
			int position = 0;
			atomsOnOnePath = pathCollection.get(j);
			lastSecondAtom = atomsOnOnePath.get(atomsOnOnePath.size() - 1);
			lastAtom = atomsOnOnePath.get(atomsOnOnePath.size() - 2);
			String type_1 = lastSecondAtom.getSymbol();
			try {
				IAtomType s_1 = satm.findMatchingAtomType(mol, lastSecondAtom);
				type_1 = s_1.getAtomTypeName();
			}catch(Exception e) {
				
			}
			for (int k = 0; k < 24; k++) {
				if (type_1.equalsIgnoreCase(str[k]))
					position = k;
			}
			atomType[position] = 1;
			IBond bond = mol.getBond(lastSecondAtom, lastAtom);
			IBond.Order order = bond.getOrder();
			o = order.numeric();
			boolean isAromatic = bond.getFlag(CDKConstants.ISAROMATIC);
			if (isAromatic)
				aromaticBond = 1;
			else if (o == 1)
				singleBond = 1;
			else if (o == 2)
				doubleBond = 1;
			else
				tripleBond = 1;
		}
		StringBuffer strbuf = new StringBuffer();
		for (int n = 0; n < 24; n++) {
			strbuf.append(atomType[n]).append(',');
		}
		strbuf.append(singleBond).append(",");
		strbuf.append(doubleBond).append(",");
		strbuf.append(tripleBond).append(",");
		strbuf.append(aromaticBond);
		enviromentFeatures = strbuf.toString();
			//System.out.println("enviromentFeatures's length" + enviromentFeatures.split(",").length);
		return enviromentFeatures;
	}
	/**
	 * Generate an HashMap<String, String> based on the atom's uniqueID within the molecule
	 */
	public HashMap<Integer, String> generateAtomFPMap(IAtomContainer molecule) throws Exception{
		HashMap<Integer, String> resultMap = new HashMap<>();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		
		for(int i = 0; i < molecule.getAtomCount(); i++){
			IAtom oneAtom = molecule.getAtom(i);
			Integer uniqueID = oneAtom.getProperty("UniqueID");
			if(!resultMap.containsKey(uniqueID)){
				String fp = this.gaf_new.generateAtomFingeprint(molecule,oneAtom);
				resultMap.put(uniqueID, fp);
			}
		
		}
		return resultMap;
		
	}
	/**
	 * Generate atom fingerprint that satisfies some specific patterns
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public String generateAtomFingerPrints(IAtomContainer molecule, IAtom oneAtom, int depth, boolean mergeFP, HashMap<Integer, String> resultMap) throws Exception{
		/**
		 * Always reset the the molecule before processsing it
		 */
		String[] resultFP = new String[(depth+1) * FingerPrintQueries_new.newQueriesList.size()];
		//Initialize resultFP
		for(int i = 0; i < resultFP.length; i++){
			resultFP[i] = "0";
		}
		for(int i = 0; i < depth+1; i++){			
			ArrayList<IAtom> neighAtoms = getNeighAtomsDepth(molecule, oneAtom, depth);
			for(int j = 0; j < neighAtoms.size(); j++){
				IAtom neighAtom = neighAtoms.get(j);
				String[] fp = resultMap.get(neighAtom.getProperty("UniqueID")).split(",");
				//String[] fp = GenerateAtomFeatures_new.generateAtomFingeprint(molecule,neighAtom).split(",");

				for(int k = 0; k < fp.length; k++){
					if(fp[k].equals("1")){
						//resultFP[k + FingerPrintQueries_new.newQueriesList.size()*i + GeneratePatterns.queriesList.size()] = "1";
						resultFP[k + FingerPrintQueries_new.newQueriesList.size()*i] = "1";
					}
				}
			}			
		}
		StringBuffer outPut = new StringBuffer();
		for(int i = 0; i < resultFP.length; i++){
			if(i == resultFP.length-1){//If reaches the end
				outPut.append(resultFP[i]);
			}
			else outPut.append(resultFP[i] + ",");
		}
		return outPut.toString();
	}
	
	public ArrayList<IAtom> getNeighAtomsDepth(IAtomContainer molecule, IAtom center, int depth){
		ArrayList<IAtom> resultList = new ArrayList<>();
		List<List<IAtom>> pathCollection = PathTools.getPathsOfLength(molecule, center, depth);
		//System.out.println("center is " + center.getSymbol() + (molecule.indexOf(center)+1));
		for(int i = 0; i < pathCollection.size(); i++){
			IAtom endpoint = pathCollection.get(i).get(pathCollection.get(i).size()-1);
			resultList.add(endpoint);//Get the atom at the given depth
			//System.out.println("End is " + endpoint.getSymbol() + (molecule.indexOf(endpoint)+1));
		}
		return resultList;
		
	}
	public void loadCypIdx(String cyp) throws Exception{
		int index = -1;
		for(int i = 0; i < cypList.length; i++){
			if(StringUtils.containsIgnoreCase(cyp, cypList[i])){
				index = i;
			}
		}
		if(index == -1){
			throw new Exception("The input CYP can not be identified.");
		}
		else{
			this.cypIdx = index;
		}
		
	}
	/**
	 * Reformat the instances of the input molecules so they can be input to the learned models by applying feature selection and normalization
	 * @param testSet: The instances of the user's input molecules
	 * @param attList: The attributes that are used in the learned model
	 * @param meanList: The mean values of each attribute
	 * @param maxList: The max values of each attribute
	 * @param minList: The min values of each attribute
	 * @return
	 * @throws Exception
	 */
	public Instances matchAttributes(Instances testSet, ArrayList<String> attListToMatch) throws Exception{
		/*
		 * Create the class atrribute
		 */
		List my_nominal_values = new ArrayList(3); 
		my_nominal_values.add("F"); 
		my_nominal_values.add("T"); 
		Attribute lastAttribute  = new Attribute(attListToMatch.get(attListToMatch.size()-1), my_nominal_values);
		/*
		 * Convert the attribute name strings to attribute
		 */
		ArrayList<Attribute> attsToMatch = new ArrayList<Attribute>();
		for(int i = 0; i<attListToMatch.size();i++){
			if(i<attListToMatch.size()-1){
				Attribute attribute = new Attribute(attListToMatch.get(i));
				attsToMatch.add(attribute);
			}
			else{
				
				attsToMatch.add(lastAttribute);
			}
		}

		int numAttributes = attListToMatch.size();
		Instances matched = new Instances("Rel", attsToMatch, testSet.size());
		Instance sample = new DenseInstance(numAttributes);
		/*
		 * Assign values for each testing instance
		 */
		for(int idxInstance = 0; idxInstance< testSet.size(); idxInstance++){
			//ArrayList misFeatures = new ArrayList();
			Instance temp = testSet.get(idxInstance);
			String label;
			//System.out.println("Check Predictor temp.classValue() == " + temp.classValue());
			if(temp.classValue()==0.0){
				label = "F";
			}
			else{
				label = "T";
			}
			
			for(int idxAttribute = 0; idxAttribute < numAttributes; idxAttribute++){
				  Attribute oneAttribute = attsToMatch.get(idxAttribute);
				  if(idxAttribute<numAttributes-1){					  
					  String nameOfAtt = oneAttribute.toString();
					  String[] formatAtt = nameOfAtt.split(" ");
					  
					  if(formatAtt.length!=3){
						  String combineSpaces = "";
						  for(int kk = 1; kk<formatAtt.length-1; kk++){
							  if(kk<formatAtt.length - 2){
								  combineSpaces = combineSpaces + formatAtt[kk]+ " ";
							  }
							  else combineSpaces = combineSpaces + formatAtt[kk];
						  }
						  formatAtt[1] = combineSpaces;
					  }
					 
					  if(formatAtt[1].contains("'")){
						  formatAtt[1] = formatAtt[1].replace("'", "");
					  }
					  //System.out.println(formatAtt[1]);
					  Attribute convAtt = testSet.attribute(formatAtt[1]);
					  double vle = temp.value(convAtt);
					  if(!Double.isNaN(vle)){					
						  sample.setValue(oneAttribute, vle);
					  }

				  }
				  else{
					  //Set class with random value. Does not matter
					  sample.setValue(lastAttribute, label);
					  
				  }
			}						
			matched.add(sample);
		}
		if(matched.classIndex() == -1){
			matched.setClassIndex(matched.numAttributes()-1);
		}
		
		return matched;
		
	}
	/**
	 * Create getNumHydrogenList nominal attribute.
	 * @return
	 */
	public ArrayList<String> getNumHydrogenList(){
		ArrayList<String> numHydrogenList = new ArrayList<>();
		numHydrogenList.add("C");//new for type 3
		numHydrogenList.add("S");//new for type 3
		numHydrogenList.add("N");//new for type 3 
		numHydrogenList.add("P");//new for type 3
		numHydrogenList.add("CH");
		numHydrogenList.add("CH2");
		numHydrogenList.add("CH3");
		numHydrogenList.add("NH");
		numHydrogenList.add("NH2");
		numHydrogenList.add("NH3");
		numHydrogenList.add("PH");
		numHydrogenList.add("PH2");
		numHydrogenList.add("PH3");
		numHydrogenList.add("SH");
		numHydrogenList.add("SH2");
		numHydrogenList.add("SH3");
		return numHydrogenList;
	}
}
