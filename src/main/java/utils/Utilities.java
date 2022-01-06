package utils;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

//import com.sun.xml.internal.bind.v2.runtime.unmarshaller.IntData;

import weka.attributeSelection.InfoGainAttributeEval;
import weka.attributeSelection.Ranker;
import weka.core.Attribute;
import weka.core.AttributeStats;
import weka.core.Instance;
import weka.core.Instances;
import weka.experiment.Stats;
import weka.filters.Filter;
import weka.filters.supervised.attribute.AttributeSelection;
import weka.filters.unsupervised.attribute.Reorder;

public class Utilities {
	/**
	 * This function is used to select a subset of features by using information gain algorithm
	 * Eg. Threshold = 50 means select the first 50 features whose information gain value is highest.
	 * @param dataset
	 * @param threshold: thedhold - 3 = the real number of features selected. Because class, name and bond are not considered as features.
	 * @return
	 * @throws Exception
	 */
	public static Instances selectByInfoGain(Instances dataset, String threshold_string) throws Exception{
		int threshold = Integer.parseInt(threshold_string);
		AttributeSelection filter = new AttributeSelection();
		InfoGainAttributeEval infoGain = new InfoGainAttributeEval();
		Ranker ranker = new Ranker();
	    //System.out.println(newData.numAttributes());
	    /**
	     * Set number of attributes
	     */
	    ranker = new Ranker();
		//ranker.setNumToSelect(100);
		filter.setEvaluator(infoGain);
		//int numFeatures = dataset.numAttributes(); // Get the number of features 
		ranker.setNumToSelect(threshold);
		
		//ranker.getGenerateRanking();
		filter.setSearch(ranker);		
	    filter.setInputFormat(dataset);
	    Instances fianlNewData = Filter.useFilter(dataset, filter);
	    if(dataset.classIndex() == -1){
	    	dataset.setClassIndex(dataset.numAttributes()-1);
	    }
	    return fianlNewData;
	}

	/**
	 * 1st. Existence of the attributes in the testDataSet; 
	 * 2nd. The orders of the attributes in the testDataSet;
	 * @param testInstances
	 * @param trainInstances
	 * @return
	 * @throws Exception
	 */
	public static Instances getSameAttributes_Mp(Instances testInstances, Instances trainInstances) throws Exception{
		Instances allignedInstances = testInstances;
		ArrayList<String> allAtt = new ArrayList();
		ArrayList<Integer> retainAttIdx = new ArrayList<Integer>();
		//This for loop stores all the attributes in the testInstances into ArrayList<Attribute> allAtt
		for(int j = 0; j < testInstances.numAttributes(); j++){
			Attribute oneAtt = testInstances.attribute(j);
			allAtt.add(oneAtt.name());
		}
		/* This for loop will:          Note that attributes in the trainInstances belong to a subset of attributes of testInstances
		 * 1. get the i-th attribute att of the trainInstances
		 * 2. find the index of attribute att in the testInstances
		 * 3. Store i-th as the index of attribute att in the ArrayList<Integer> retainAttIdx. This order will be used in the reorder filter 
		 */
		for(int i = 0; i < trainInstances.numAttributes(); i++){
			Attribute oneTrainAtt = trainInstances.attribute(i);
			int idx = allAtt.indexOf(oneTrainAtt.name());
			if(idx == -1){
				System.out.println(oneTrainAtt.toString() + "idx = " + idx);
			}
			retainAttIdx.add(idx);
		}
		int[] retainList = new int[retainAttIdx.size()];
		for(int k = 0; k < retainAttIdx.size(); k++){
			retainList[k] = retainAttIdx.get(k);
		}
		Reorder reorederFiler = new Reorder(); 
		reorederFiler.setAttributeIndicesArray(retainList);
		reorederFiler.setInputFormat(testInstances);
		allignedInstances = Filter.useFilter(allignedInstances, reorederFiler);

        if(allignedInstances.classIndex()==-1){
        	allignedInstances.setClassIndex(allignedInstances.numAttributes()-1);
        }
		return allignedInstances;
	}
	/**
	 * 1st. Existence of the attributes in the testDataSet; 
	 * 2nd. The orders of the attributes in the testDataSet;
	 * @param testInstances
	 * @param trainInstances
	 * @return
	 * @throws Exception
	 */
	public static Instances getSameAttributes(Instances testInstances, Instances trainInstances) throws Exception{
		Instances allignedInstances = testInstances;
		ArrayList<Attribute> allAtt = new ArrayList<Attribute>();
		ArrayList<Integer> retainAttIdx = new ArrayList<Integer>();
		//This for loop stores all the attributes in the testInstances into ArrayList<Attribute> allAtt
		for(int j = 0; j < testInstances.numAttributes(); j++){
			Attribute oneAtt = testInstances.attribute(j);
			allAtt.add(oneAtt);
		}
		/* This for loop will:          Note that attributes in the trainInstances belong to a subset of attributes of testInstances
		 * 1. get the i-th attribute att of the trainInstances
		 * 2. find the index of attribute att in the testInstances
		 * 3. Store i-th as the index of attribute att in the ArrayList<Integer> retainAttIdx. This order will be used in the reorder filter 
		 */
		for(int i = 0; i < trainInstances.numAttributes(); i++){
			Attribute oneTrainAtt = trainInstances.attribute(i);
			int idx = allAtt.indexOf(oneTrainAtt);
			if(idx == -1){
				System.out.println(oneTrainAtt.toString() + "idx = " + idx);
			}
			retainAttIdx.add(idx);
		}
		int[] retainList = new int[retainAttIdx.size()];
		for(int k = 0; k < retainAttIdx.size(); k++){
			retainList[k] = retainAttIdx.get(k);
		}
		Reorder reorederFiler = new Reorder(); 
		reorederFiler.setAttributeIndicesArray(retainList);
		reorederFiler.setInputFormat(testInstances);
		allignedInstances = Filter.useFilter(allignedInstances, reorederFiler);

        if(allignedInstances.classIndex()==-1){
        	allignedInstances.setClassIndex(allignedInstances.numAttributes()-1);
        }
		return allignedInstances;
	}
	/**This function is used to find indices of name and bond attribute.
	 * Both name and bond attributes need to be removed during the training and testing process.
	 * Note that indices of name and bond are both incremented by 1 because the Remove filter take index of attribute starting at 1.
	 * @param dataset
	 * @return
	 * @throws Exception
	 */
	public static String findIndcesOfNameAndBond(Instances dataset) throws Exception{
		StringBuffer foundIdx = new StringBuffer();
		int nameIdx = -1;
		int bondIdx = -1;
		boolean nameFound = false;
		boolean bondFound = false;
		for(int i = 0; i < dataset.numAttributes(); i++){
			Attribute oneAtt = dataset.attribute(i);
			if(!nameFound && oneAtt.name().equals("name")){
				nameIdx = i + 1;
				nameFound = true;
			}
		}
		for(int i = 0; i < dataset.numAttributes(); i++){
			Attribute oneAtt = dataset.attribute(i);
			if(nameFound && oneAtt.name().equals("bond")){
				bondIdx = i + 1;
				bondFound = true;
			}
		}
		if(!(nameFound&&bondFound)){
			throw new Exception("Name or Bond was removed during feature selection process!");
		}
		if(nameIdx < bondIdx){
			foundIdx.append(String.valueOf(nameIdx) + "," + String.valueOf(bondIdx));
		}
		else{
			foundIdx.append(String.valueOf(bondIdx) + "," + String.valueOf(nameIdx));
		}
		return foundIdx.toString();
		
	}
	/**
	 * Find the index of name attribute. This function is used in stratification of Cross-Validation
	 * @param dataset
	 * @return
	 * @throws Exception
	 */
	public static int findIndexOfName(Instances dataset) throws Exception{
		int nameIdx = -1;
		boolean nameFound = false;
		for(int i = 0; i < dataset.numAttributes(); i++){
			Attribute oneAtt = dataset.attribute(i);
			if(!nameFound && oneAtt.name().equals("name")){
				nameIdx = i;
				nameFound = true;
			}
		}
		if(!nameFound){
			throw new Exception("Name was removed during feature selection process!");
		}
		return nameIdx;
		
	}
	public static int findIndexOfName(Instance dataset) throws Exception{
		int nameIdx = -1;
		boolean nameFound = false;
		for(int i = 0; i < dataset.numAttributes(); i++){
			Attribute oneAtt = dataset.attribute(i);
			if(!nameFound && oneAtt.name().equals("name")){
				nameIdx = i;
				nameFound = true;
			}
		}
		if(!nameFound){
			throw new Exception("Name was removed during feature selection process!");
		}
		return nameIdx;
		
	}
	/**
	 * Find the index of bond attribute. This function is used in stratification of Cross-Validation
	 * @param dataset
	 * @return
	 * @throws Exception
	 */
	public static int findIndexOfBond(Instances dataset) throws Exception{
		int bondIdx = -1;
		boolean bondFound = false;
		for(int i = 0; i < dataset.numAttributes(); i++){
			Attribute oneAtt = dataset.attribute(i);
			if(!bondFound && oneAtt.name().equals("bond")){
				bondIdx = i;
				bondFound = true;
			}
		}
		if(!bondFound){
			throw new Exception("Bond was removed during feature selection process!");
		}
		return bondIdx;
		
	}
	public static int findIndexOfBond(Instance dataset) throws Exception{
		int bondIdx = -1;
		boolean bondFound = false;
		for(int i = 0; i < dataset.numAttributes(); i++){
			Attribute oneAtt = dataset.attribute(i);
			if(!bondFound && oneAtt.name().equals("bond")){
				bondIdx = i;
				bondFound = true;
			}
		}
		if(!bondFound){
			throw new Exception("Bond was removed during feature selection process!");
		}
		return bondIdx;
		
	}
	/**
	 * This function returns the instances whose molecule names are stored in the nameList
	 * @param nameList
	 * @param dataset
	 * @return
	 */
	public Instances getSubInstances(Instances nameList, Instances dataset) throws Exception{
		int nameIdx = Utilities.findIndexOfName(dataset);
		ArrayList<String> name_subList = convertNameList(nameList);
		ArrayList<Attribute> attList = new ArrayList();
		for(int i = 0; i < dataset.numAttributes(); i++){
			Attribute oneAtt = dataset.attribute(i);
			attList.add(oneAtt);
		}
		Instances subInstances = new Instances("Rel", attList, 100000);
		for(int i = 0; i < dataset.numInstances(); i++){
			Instance oneInstance = dataset.get(i);
			String name = oneInstance.stringValue(nameIdx);//The name attribute
			if(name_subList.contains(name)){
				subInstances.add(oneInstance);
			}
		}
		if(subInstances.classIndex()==-1){
			subInstances.setClassIndex(subInstances.numAttributes()-1);
		}
		for(int i = 0; i < subInstances.numAttributes(); i++){
			Stats attStats = subInstances.attributeStats(i).numericStats;
			Double max = attStats.max;
			Double min = attStats.min;
			Double std = attStats.stdDev;
		}
		
		return subInstances;
		
	}
	/**
	 * Convert nameListInstances from "Instances" to "ArrayList<String>"
	 * @param nameListInstances
	 * @return ArrayList<Sting> nameList
	 */
	public ArrayList<String> convertNameList(Instances nameListInstances){
		ArrayList<String> nameList = new ArrayList();
		for(int i = 0; i < nameListInstances.numInstances(); i++){
			Instance nameInst = nameListInstances.get(i);
			String name = nameInst.stringValue(0); // Get the String value of the name attribute of this instance
			nameList.add(name);
		}
		return nameList;
	}
	/**
	 * This function convert a String[] to String, two successive elements are split by ","
	 * @param bondFps
	 * @return
	 */
	public static String convertListToString(String[] bondFps){
		StringBuffer result = new StringBuffer();
		result.append(bondFps[0]);
		for(int i = 1; i < bondFps.length; i++){
			result.append("," + bondFps[i]);
		}
		return result.toString();
	}
	
	
	/**
	 * This function is used to create the title Line of the dataset csv file
	 * @param molecularFeatureTitle
	 * @param bondFeatureitle
	 * @param atomTypeTitle
	 * @param atomicFeatureTitle
	 * @param bondFpTitle
	 * @return
	 */
	public static String creatTypeOneTitleLine(String[] cypList, String[] molecularFeatureTitle, String[] bondFeatureitle, String[] atomTypeTitle,
																String[] atomicFeatureTitle, String[] bondFpTitle, int depth, boolean mergedFp){
		//Molecular + bondfeatures + allfingerprint + atomType + atomicFeature
		// allfingerprint = current + left + right; bondFeatures = currentBond; atomType = left + right; atomicFeature = left + right
		StringBuffer titleLineBuffer = new StringBuffer();
		titleLineBuffer.append("Name" + "," + "Bond");
		/*
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
		/*
		 * Molecular Features for Current Molecule
		 */
		for(int i = 0; i < molecularFeatureTitle.length; i++){
			titleLineBuffer.append("," + molecularFeatureTitle[i]);
		}
		/*
		 * BondFeatures for CurrentBond
		 */
		for(int i = 0; i < bondFeatureitle.length; i++){
			titleLineBuffer.append("," + bondFeatureitle[i]);
		}
		/*
		 * Fingerprints for CurrentBond + LeftBond + RightBond
		 */
		if(!mergedFp){
			for(int i = 0; i < bondFpTitle.length; i++){
				titleLineBuffer.append("," + bondFpTitle[i] + "_current");
			}
			for(int d = 1; d < depth + 1; d++){
				for(int i = 0; i < bondFpTitle.length; i++){
					titleLineBuffer.append("," + "d=" + d +"_" +  bondFpTitle[i] + "_left");
				}
			}
			for(int d = 1; d < depth + 1; d++){
				for(int i = 0; i < bondFpTitle.length; i++){
					titleLineBuffer.append("," + "d=" + d +"_" +  bondFpTitle[i] + "_right");
				}
			}
		}
		else{
			for(int i = 0; i < bondFpTitle.length; i++){
				titleLineBuffer.append("," + bondFpTitle[i] + "_current");
			}
			for(int i = 0; i < bondFpTitle.length; i++){
				titleLineBuffer.append("," + bondFpTitle[i] + "_left");
			}
			for(int i = 0; i < bondFpTitle.length; i++){
				titleLineBuffer.append("," + bondFpTitle[i] + "_right");
			}
		}
		
		/*
		 * AtomType Features for leftAtom + RightAtom
		 */
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + atomTypeTitle[i] + "_left");
		}
		/*
		 * The section below need to be deactivated when using variant of SybylAtomType
		*/
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + atomTypeTitle[i] + "_right");
		}
		/*
		 * Atomic descriptor Features for leftAtom + RightAtom
		 */
		for(int i = 0; i < atomicFeatureTitle.length; i++){
			titleLineBuffer.append("," + atomicFeatureTitle[i] + "_left");
		}
		for(int i = 0; i < atomicFeatureTitle.length; i++){
			titleLineBuffer.append("," + atomicFeatureTitle[i] + "_right");
		}
		return titleLineBuffer.toString();
	}
	
	/**
	 * This function is used to create the title Line of the dataset csv file
	 * @param molecularFeatureTitle
	 * @param bondFeatureitle
	 * @param atomTypeTitle
	 * @param atomicFeatureTitle
	 * @param bondFpTitle
	 * @return
	 */
	public static String creatTypeOneTitleLine_NA(String[] cypList, String[] molecularFeatureTitle, String[] bondFeatureitle, String[] atomTypeTitle,
																String[] atomicFeatureTitle, String[] bondFpTitle, int depth, boolean mergedFp){
		//Molecular + bondfeatures + allfingerprint + atomType + atomicFeature
		// allfingerprint = current + left + right; bondFeatures = currentBond; atomType = left + right; atomicFeature = left + right
		StringBuffer titleLineBuffer = new StringBuffer();
		titleLineBuffer.append("Name" + "," + "Bond");
		/*
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
		/*
		 * Molecular Features for Current Molecule
		 */
		for(int i = 0; i < molecularFeatureTitle.length; i++){
			titleLineBuffer.append("," + molecularFeatureTitle[i]);
		}
		/*
		 * BondFeatures for CurrentBond
		 */
		for(int i = 0; i < bondFeatureitle.length; i++){
			titleLineBuffer.append("," + bondFeatureitle[i]);
		}
		/*
		 * Fingerprints for CurrentBond + LeftBond + RightBond
		 */
		if(!mergedFp){
			for(int i = 0; i < bondFpTitle.length; i++){
				titleLineBuffer.append("," + bondFpTitle[i] + "_current");
			}
			for(int d = 1; d < depth + 1; d++){
				for(int i = 0; i < bondFpTitle.length; i++){
					titleLineBuffer.append("," + "d=" + d +"_" +  bondFpTitle[i] + "_left");
				}
			}
			for(int d = 1; d < depth + 1; d++){
				for(int i = 0; i < bondFpTitle.length; i++){
					titleLineBuffer.append("," + "d=" + d +"_" +  bondFpTitle[i] + "_right");
				}
			}
		}
		else{
			for(int i = 0; i < bondFpTitle.length; i++){
				titleLineBuffer.append("," + bondFpTitle[i] + "_current");
			}
			for(int i = 0; i < bondFpTitle.length; i++){
				titleLineBuffer.append("," + bondFpTitle[i] + "_left");
			}
			for(int i = 0; i < bondFpTitle.length; i++){
				titleLineBuffer.append("," + bondFpTitle[i] + "_right");
			}
		}
		
		/*
		 * AtomType Features for leftAtom + RightAtom
		 */
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + atomTypeTitle[i] + "_left");
		}
		/*
		 * The section below need to be deactivated when using variant of SybylAtomType
		*/
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + atomTypeTitle[i] + "_right");
		}
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + atomTypeTitle[i] + "_left_neighbor");
		}
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + atomTypeTitle[i] + "_right_neighbor");
		}
		/*
		 * Atomic descriptor Features for leftAtom + RightAtom
		 */
		for(int i = 0; i < atomicFeatureTitle.length; i++){
			titleLineBuffer.append("," + atomicFeatureTitle[i] + "_left");
		}
		for(int i = 0; i < atomicFeatureTitle.length; i++){
			titleLineBuffer.append("," + atomicFeatureTitle[i] + "_right");
		}
		return titleLineBuffer.toString();
	}
	/**
	 * This function is used to create the title Line of the dataset csv file
	 * @param molecularFeatureTitle
	 * @param bondFeatureitle
	 * @param atomTypeTitle
	 * @param atomicFeatureTitle
	 * @param bondFpTitle
	 * @return
	 */
	public static String creatTypeOneTitleLine(String[] cypList, String[] molecularFeatureTitle, String[] bondFeatureitle, String[] atomTypeTitle,
																String[] atomicFeatureTitle, String[] bondFpTitle){
		//Molecular + bondfeatures + allfingerprint + atomType + atomicFeature
		// allfingerprint = current + left + right; bondFeatures = currentBond; atomType = left + right; atomicFeature = left + right
		StringBuffer titleLineBuffer = new StringBuffer();
		titleLineBuffer.append("Name" + "," + "Bond");
		/*
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
		/*
		 * Molecular Features for Current Molecule
		 */
		for(int i = 0; i < molecularFeatureTitle.length; i++){
			titleLineBuffer.append("," + molecularFeatureTitle[i]);
		}
		/*
		 * BondFeatures for CurrentBond
		 */
		for(int i = 0; i < bondFeatureitle.length; i++){
			titleLineBuffer.append("," + bondFeatureitle[i]);
		}
		/*
		 * Fingerprints for CurrentBond + LeftBond + RightBond
		 */
		
		for(int i = 0; i < bondFpTitle.length; i++){
			titleLineBuffer.append("," + bondFpTitle[i] + "_current");
		}
		for(int i = 0; i < bondFpTitle.length; i++){
			titleLineBuffer.append("," + bondFpTitle[i] + "_left");
		}
		for(int i = 0; i < bondFpTitle.length; i++){
			titleLineBuffer.append("," + bondFpTitle[i] + "_left");
		}

		
		/*
		 * AtomType Features for leftAtom + RightAtom
		 */
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + atomTypeTitle[i] + "_left");
		}
		/*
		 * The section below need to be deactivated when using variant of SybylAtomType
		*/
		for(int i = 0; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + atomTypeTitle[i] + "_right");
		}
		/*
		 * Atomic descriptor Features for leftAtom + RightAtom
		 */
		for(int i = 0; i < atomicFeatureTitle.length; i++){
			titleLineBuffer.append("," + atomicFeatureTitle[i] + "_left");
		}
		for(int i = 0; i < atomicFeatureTitle.length; i++){
			titleLineBuffer.append("," + atomicFeatureTitle[i] + "_right");
		}
		return titleLineBuffer.toString();
	}
	/**
	 *
	 * @param cypList
	 * @param molecularFeatureTitle
	 * @param atomTypeTitle
	 * @param atomicFeatureTitle
	 * @param bondFpTitle
	 * @return
	 */
	public static String creatTypeTwoTitleLine(String[] cypList, String[] molecularFeatureTitle, String[] atomHydrogenCount, String[] atomTypeTitle,
			String[] atomicFeatureTitle, String[] bondFpTitle){
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
		//molecularDescriptorTitle, atomHydrogenCount, atomTypeTitle, atomicDescriptorTitle, bondFpTitle
		/**
		 * Molecular Features for Current Molecule
		 */
		for(int i = 0; i < molecularFeatureTitle.length; i++){
			titleLineBuffer.append("," + molecularFeatureTitle[i]);
		}
		
		/**
		 *	Atom Hydrogen Count. eg. CH, CH2, etc..
		 */
		for(int i = 0 ; i < atomHydrogenCount.length; i++){
			titleLineBuffer.append("," + atomHydrogenCount[i]);
		}
		/**
		 *	AtomType 
		 */
		for(int i =0 ; i < atomTypeTitle.length; i++){
			titleLineBuffer.append("," + atomTypeTitle[i]);
		}
		
		/**
		 *	Atomic Descriptors
		 */
		for(int i = 0; i < atomicFeatureTitle.length; i++){
			titleLineBuffer.append("," + atomicFeatureTitle[i]);
		}
		
		/**
		 *	Bond fingerprints by radius 
		 */
		for(int i = 0; i < bondFpTitle.length; i++){
			titleLineBuffer.append("," + bondFpTitle[i]);
		}
		return titleLineBuffer.toString();
		
	}
	/**
	 * Read a sdf file and create an IAtomContainerSet that contains all molecules in the sdf file
	 * @param String pathToInputFile---(sdf file)
	 * @return IAtomContainerSet that contains all molecules in the sdf file        
	 * @throws Exception
	 */
	public static IAtomContainerSet readFile(String pathToInputFile)
			throws FileNotFoundException, CDKException {
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(pathToInputFile),
				bldr);
		Properties prop = new Properties();
		prop.setProperty("ForceReadAs3DCoordinates", "true");
		PropertiesListener listener = new PropertiesListener(prop);
		sdfr.addChemObjectIOListener(listener);
		sdfr.customizeJob();
		IAtomContainerSet allMoles = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		while(sdfr.hasNext()){
			IAtomContainer oneMole = sdfr.next();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(oneMole.getBuilder());
			Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
			aromaticity.apply(oneMole);
			adder.addImplicitHydrogens(oneMole);
			AtomContainerManipulator.suppressHydrogens(oneMole);
			allMoles.addAtomContainer(oneMole);
		}
		return allMoles;

	}
	
	/**
	 * Check whether the molecule is a reactant for the cyp or not.
	 * @param oneMole
	 * @param cyp
	 * @return true if the compound has more than one BoM for the given cyp450 enzyme
	 */
	public static boolean isReactant(IAtomContainer oneMole, String cyp){
		boolean isReactant = false;
		String cypBoM = "BOM_" + cyp;
		String bom = (String) oneMole.getProperties().get(cypBoM);
		if(bom!=null&&!bom.equals("")){
			isReactant = true;
		}
		return isReactant;
	}
	/**
	 * This function will check if the molecule exists in the current set
	 * @param newMole
	 * @param currentSet
	 * @return
	 * @throws Exception
	 */
	public static boolean existInSet(IAtomContainer newMole, IAtomContainerSet currentSet) throws Exception{
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		String targetSmile = sg.create(newMole);
		for(int i = 0; i < currentSet.getAtomContainerCount(); i++){
			String oneSmile = sg.create(currentSet.getAtomContainer(i));
			if(oneSmile.equals(targetSmile)) return true;
		}
		return false;
	}
	/**
	 * This function will remove the duplicate structures within the input "molecules" by checking SMILES Strings
	 * @param molecules
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet removeDuplicates(IAtomContainerSet molecules) throws Exception{
		IAtomContainerSet resultSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		ArrayList<String> smilesSet = new ArrayList();
		String currentSmiles;
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			currentSmiles = sg.create(molecules.getAtomContainer(i));
			if(!smilesSet.contains(currentSmiles)){
				smilesSet.add(currentSmiles);
				resultSet.addAtomContainer(molecules.getAtomContainer(i));
			} 
		}
		return resultSet;
	}

	/**
	 * Function that saves an IAtomContainerSet of containers to an SDF file
	 *
	 * @param containers	IAtomContainerSet of containers to save
	 * @param outputFileName	Name of output file to send containers to
	 * @throws CDKException	CDK class exceptions
	 * @throws IOException    Invalid file type
	 */
	public static void saveAtomContainerSetToSDF(IAtomContainerSet containers, String outputFileName) throws CDKException, IOException{
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputFileName));
		sdfWriter.write(containers);
		sdfWriter.close();
	}
	/**
	 * Function that parses SMILES Strings in provided SDF file
	 *
	 * @param sdfFileName	Name of SDF file with SMILES Strings
	 * @return	IAtomContainerSet of parsed SMILES Strings
	 * @throws IOException	Invalid file type
	 */
	public static IAtomContainerSet parseSdf(String sdfFileName) throws IOException {
		IAtomContainerSet containers = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);

		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();

		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(sdfFileName), bldr);

		while (sdfr.hasNext()){
			IAtomContainer mol = sdfr.next();
			containers.addAtomContainer(mol);
		}

		sdfr.close();
		return containers;

	}
	
	
}
