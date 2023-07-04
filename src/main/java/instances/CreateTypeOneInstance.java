
package instances;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import typeOneFeatures.*;
//import reactantpredictor.utils.FeatureGeneration;
//import reactantpredictor.utils.MoleculeExplorer;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

public class CreateTypeOneInstance {
	public final String[] cypList = {"1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4"};
	public int cypIdx, depth_neighborAtomType, depth_neighborAtomDescriptor;
	String titleLine;
	ArrayList<String> attributes;
	public CreateTypeOneInstance(String cyp, ArrayList<String> attributes,int depth_neighborAtomType, int depth_neighborAtomDescriptor)
			throws Exception{
		loadCypIdx(cyp);
		this.depth_neighborAtomType = depth_neighborAtomType;
		this.depth_neighborAtomDescriptor = depth_neighborAtomDescriptor;
		this.attributes = attributes;
	}
	
	
	public CreateTypeOneInstance() {
		// TODO Auto-generated constructor stub
	}


	public Instances generateTypeOneBondInstances_OneMole(ArrayList<String> tempData) throws Exception{
		GenerateTitleLine_TypeOne tl = new GenerateTitleLine_TypeOne(depth_neighborAtomType, depth_neighborAtomDescriptor);
		this.titleLine = tl.generateTitleLine();
		Instances data = createInstances(tempData,titleLine);
		return data;
	}
	
	/**
	 * Generate bond features for a Type ONE bond <R1,R2>
	 * Note:
	 * 1. Each molecule should have been processed by AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mole); 
	 * 2. Each molecule should have been processed by applying the aromatic model as below: 		
	 * 		2.1 Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
	 *		2.2 aromaticity.apply(mol);
	 * 3. Each bond must have two atoms, left and right, and CH4(methane) has 0 bond in this cdk system.
	 * 4. For every bond in the molecule, we will generate two instances for it.
	 * 		4.1	Instance left = currentBondFeature + leftBondFeature + rightBondFeature
	 * 		4.2	Instance right = currentBondFeature + rightBondFeature + leftBondFeature
	 * We do this because bond <1,2> and <2,1> should be the same bond and the order should not matter. Hence, we create two instances for the same bond to avoid this order issue.
	 * @param molecule
	 * @return A HashMap<String, ArrayList<String>> ---> HashMap<"Left", Bond {1,2,3,4....}. And as well for "Right"
	 * @throws CDKException
	 * @throws IOException
	 */
	public ArrayList<String> generateTypeOneBondFeatures_OneMole(IAtomContainer molecule, int depth_neighborAtomType, int depth_neighborAtomDescriptor) throws Exception {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		String molecularFeatures = GenerateMolecularFeatures_TypeOne.generateMolecularFeatures_new(molecule);
		String molecularFP = GenerateMolecularFeatures_TypeOne.generateMoleculeFP(molecule);
		ArrayList<String> allInstances = new ArrayList<String>();
		String moleName = (String) molecule.getProperties().get("cdk:Title");
		ArrayList<String[]> boms = getBoMs(molecule);
		//If the molecule doesn't have a name, generate a SMILES string for it.
		if(moleName==null){
			 SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Absolute);
			 moleName = sg.create(molecule); //			 
		}
		/**
		 * Type one bond is <i,j> bond and thus must be recorded in the sdf file
		 */
		for(int i = 0; i < molecule.getBondCount(); i++){
			//long start =System.currentTimeMillis();	
			/*
			 * Need to process the molecule every time an atom is processed because some operation changed the informatin of the molecule .
			 */
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
			adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
			aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
			aromaticity.apply(molecule);
			adder.addImplicitHydrogens(molecule);			
			/*
			 * Start to generate features.
			 */
			IBond oneBond = molecule.getBond(i);
			//System.out.println("Processing bond: " + oneBond.getAtom(0).getIndex() + ", " + oneBond.getAtom(1).getIndex());
			String bomLabels= realBoMLabels(boms, molecule, oneBond.getAtom(0), oneBond.getAtom(1));
			String currentBondAtomType = GenerateBondFeatures_TypeOne.generateCurrentBondAtomType(molecule, oneBond);
				
			String currentBondFP = GenerateBondFeatures_TypeOne.generateBondFingeprint(molecule, oneBond);
			String currentBondDescriptor = GenerateBondFeatures_TypeOne.generateCurrentBondAtomDescriptor(molecule, oneBond);

		
			String neighborAtomType = GenerateNeighborFeatures_TypeOne.generateNeighborAtomType_upToDepth(molecule, oneBond, depth_neighborAtomType);
			String neighborAtomDescriptor = GenerateNeighborFeatures_TypeOne.generateNeighborAtomDescriptor(molecule, oneBond, depth_neighborAtomDescriptor);
			int[] longestPath = GenerateBondFeatures_TypeOne.getLongestPath(molecule, oneBond);
			List<IAtom> orderedAtoms = GenerateBondFeatures_TypeOne.reorderAtomsOfBond(molecule, oneBond);
			
			//long end =System.currentTimeMillis();
			//System.out.println("Neighbour features time cost: " + (end - start));
			
			IAtom leftAtom = orderedAtoms.get(0);
			IAtom rightAtom = orderedAtoms.get(1);
			//IAtom leftAtom = oneBond.getAtom(0);
			//IAtom rightAtom = oneBond.getAtom(1);
			String oneInstance = moleName + "," +"<" + leftAtom.getSymbol() + "." + (molecule.indexOf(leftAtom) + 1) + ";" 
					+ rightAtom.getSymbol() + "." + (molecule.indexOf(rightAtom) + 1) + ">" + ","
					+ bomLabels + ","
					+ longestPath[0] + "," + longestPath[1] + "," 
					+ molecularFeatures + "," + molecularFP + "," + currentBondFP + "," + currentBondAtomType + "," + currentBondDescriptor + "," + neighborAtomType + "," + neighborAtomDescriptor;
							
			allInstances.add(oneInstance); //The bonds of the molecule is added into the array in the order of the bond			
		}
		return allInstances;
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
			if((!parseLine[i].contains("H"))&&(Double.parseDouble(parseLine[i])< 1.0e-4)){
					parseLine[i] = "0";
			}
			oneInstance.setValue(oneAtt, Double.parseDouble(parseLine[i]));
		}		
		return oneInstance;
		
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
	 * This function checks whether the <i,j> bond connected with leftAtom and rightAtom is a BoM
	 * @param molecule
	 * @param leftAtom
	 * @param rightAtom
	 * @return A 9-tuple string vector consists of "0"(No) or "1"(Yes). Each bit for one of the nine CYP450 enzymes.
	 * @throws Exception
	 */
	public String realBoMLabels(ArrayList<String[]> boms, IAtomContainer molecule, IAtom leftAtom, IAtom rightAtom) throws Exception{
		StringBuffer result = new StringBuffer();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		String leftIdx = String.valueOf((molecule.indexOf(leftAtom)+1));
		String rightIdx = String.valueOf((molecule.indexOf(rightAtom)+1));
		for(int i = 0; i < boms.size(); i++){
			boolean isBoM = false;
			String[] bomOneCyp = boms.get(i);
			for(int j = 0; j < bomOneCyp.length; j++){
				String[] oneBoM = bomOneCyp[j].split(",");
				//System.out.println(Arrays.toString(bomOneCyp));
				if(oneBoM[0].equals(leftIdx)&&oneBoM[1].equals(rightIdx)){
					//result.append("1");
					isBoM = true;
					break;
				}
				if(oneBoM[0].equals(rightIdx)&&oneBoM[1].equals(leftIdx)){
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
	 * This function generates the attribute list for the dataset. The first attribute is the name of the molecule. The last attribute is the class attribute.
	 * @param titleLine: The string contains all the attribute names
	 * @param nameList: The names of molecules. Each name string is a nominal string in the name attribute.
	 * Note that the number of attributes is 486 = 26 + 17*26 + 18
	 * @return
	 */
	public ArrayList<Attribute> createAttributes(String titleLine){
		ArrayList<Attribute> attList = new ArrayList<>();
		String[] parseLine = titleLine.split(",");
		for(int i = 0; i<11; i++){
			Attribute oneAtt = new Attribute(parseLine[i],  (ArrayList<String>) null);
			attList.add(oneAtt);
		}
		//Feature attributes starts at index = 3.
		for(int i = 11; i < parseLine.length; i++){
			Attribute oneAtt = new Attribute(parseLine[i]);
			attList.add(oneAtt);
		}
		List<String> class_nominal_values = new ArrayList<String>();
		class_nominal_values.add("F");// F = not BoM. 0
		class_nominal_values.add("T");// T = BoM. 1 
		Attribute classAttribute = new Attribute("BoM",class_nominal_values);
		attList.add(classAttribute);
		return attList;
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
}
