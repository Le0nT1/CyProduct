package instances;

import java.util.ArrayList;

import org.openscience.cdk.interfaces.IAtomContainer;

import weka.core.Instances;

public class GenerateInstance {
	/**
	 * This function generates raw features, not the CYP specific features.
	 * These raw features will be used to generate CYP specific features in another function.
	 * We want to generate CYP specific features from the raw features to save time.
	 * @param molecule
	 * @param type
	 */
	public ArrayList<String> generateRawInstances(IAtomContainer molecule, int type) throws Exception{
		ArrayList<String> tempData;
		//Instances data;
		//long start_time = System.currentTimeMillis();
		//long end_feature_time = System.currentTimeMillis();		
		IAtomContainer oneMole = molecule.clone();
		if(type == 1) {
			CreateTypeOneInstance sts = new CreateTypeOneInstance();	
			tempData = sts.generateTypeOneBondFeatures_OneMole(oneMole,4, 2);
			//end_feature_time = System.currentTimeMillis();	
			//data = sts.generateTypeOneBondInstances_OneMole(oneMole, tempData);		
			//System.out.println("type One features time: " + (end_feature_time - start_time));
		}
		else if(type == 2) {
			CreateTypeTwoThreeInstances twoThree = new CreateTypeTwoThreeInstances();		
			tempData = twoThree.generateTypeTwoAtomBasedFeatures(oneMole,4, false);						
			//end_feature_time = System.currentTimeMillis();				
			//data = twoThree.generateTypeTwoInstances_OneMole(oneMole,tempData,4);							
			//System.out.println("type Two features time: " +  (end_feature_time - start_time));
		}
		else if(type == 3) {
			CreateTypeTwoThreeInstances twoThree = new CreateTypeTwoThreeInstances();
			tempData = twoThree.generateTypeThreeAtomBasedFeatures(oneMole,4, false);
			//end_feature_time = System.currentTimeMillis();			
			//data = twoThree.generateTypeTwoInstances_OneMole(oneMole,tempData,4);			
			//System.out.println("type Three features time: " +  (end_feature_time - start_time));
		}
		else throw new Exception("Unknown Bond of Metabolism type: " + type);
		return tempData;
	}
	/**
	 * This function will generate CYP specific Instances for type One BoMs
	 * @param cyp
	 * @param attributes
	 * @param depth_neighborAtomType
	 * @param depth_neighborAtomDescriptor
	 * @return
	 * @throws Exception
	 */
	public static Instances generateTypeOneInstances_forCYP(String cyp, ArrayList<String> rawFeatureList, ArrayList<String> attributes, int depth_neighborAtomType, int depth_neighborAtomDescriptor) throws Exception{
		CreateTypeOneInstance sts = new CreateTypeOneInstance(cyp, attributes, depth_neighborAtomType, depth_neighborAtomDescriptor);
		Instances data = sts.generateTypeOneBondInstances_OneMole(rawFeatureList);
		return data;
	}
	/**
	 * This function will generate CYP specific Instances for type Two and Three BoMs
	 * @param cyp
	 * @param attributes
	 * @param type
	 * @return
	 * @throws Exception
	 */
	public static Instances generateTypeTwoAndThreeInstances_forCYP(String cyp, ArrayList<String> rawFeatureList, ArrayList<String> attributes, int type) throws Exception{
		CreateTypeTwoThreeInstances twoThree = new CreateTypeTwoThreeInstances(cyp, attributes, 4);
		Instances data;
		if(type == 2) {
			data = twoThree.generateTypeTwoInstances_OneMole(rawFeatureList);
		}
		else if(type == 3) {
			data = twoThree.generateTypeTwoInstances_OneMole(rawFeatureList);
		}
		else throw new Exception("This function only handles type Two and Three Bond of Metabolism. However, the input type is " + type);
		return data;
	}
}
