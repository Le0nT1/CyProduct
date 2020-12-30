 package cyProduct;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import reactantpredictor.BioTransformerAPIs;
import utils.MoleculeExplorer;
import weka.classifiers.Classifier;
import weka.core.Instances;

public class PredictorForAllThreeBoMs {
	public static final String[] cyps = {"1A2","2A6","2B6","2C8","2C9","2C19","2D6","2E1","3A4"};
	public static final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
	Classifier currentClassifier;
	MoleculesToInstances sts;
	ArrayList<String> attributes = new ArrayList<String>();
	Double beta;
	String cyp;
	int depth_neighborAtomType = 4;
	int depth_neighborAtomDescriptor = 2;
	
	public static boolean isReactant(IAtomContainer oneMole, String cyp) throws Exception{
		//ReactantPred reactantPredictor = new ReactantPred();
		BioTransformerAPIs cypReactAPI = new BioTransformerAPIs();
		boolean isReactant = cypReactAPI.predictReactant(oneMole, cyp);
		//boolean isReactant = true;
		return isReactant;
	}
	/**
	 * Set up the model and support file that are used to make prediction
	 * @param cyp
	 * @param inputPath
	 * @param supportFolderPath
	 * @throws Exception
	 * Specifically, to read the file within the jar file, check this link: https://stackoverflow.com/questions/20389255/reading-a-resource-file-from-within-jar
	 */
	public void setup(String cyp, int type) throws Exception{
		this.cyp = cyp;
		String modelPath = "/CYP" + cyp + "/" + cyp + ".model";  
		String supportPath = "/CYP" + cyp + "/" + cyp + "_supportfiles.csv";
		if(type == 1){
			modelPath = "/CYP" + cyp + "/TypeOne/" + cyp + ".model"; 
			supportPath = "/CYP" + cyp + "/TypeOne/" + cyp + "_supportfiles.csv";
		}
		else if(type == 2){
			modelPath = "/CYP" + cyp + "/TypeTwo/" + cyp + ".model"; 
			supportPath = "/CYP" + cyp + "/TypeTwo/" + cyp + "_supportfiles.csv";
		}
		else if(type == 3){
			modelPath = "/CYP" + cyp + "/TypeThree/" + cyp + ".model"; 
			supportPath = "/CYP" + cyp + "/TypeThree/" + cyp + "_supportfiles.csv";
		}
//	    String resourceFile = &quot;resource/test_resource.txt&quot;;	    
//	    InputStream resourceStream = ClassLoader.getSystemClassLoader().getResourceAsStream(resourceFile);
		this.attributes = new ArrayList<String>();
//		ClassLoader classLoader = getClass().getClassLoader();
		//File file = new File(classLoader.getResource(modelPath).getFile());
		//modelPath = file.getPath();
		InputStream model_Stream = getClass().getResourceAsStream(modelPath);
		loadModel(model_Stream);
		//File file = new File(classLoader.getResource(supportPath).getFile());
		//supportPath = file.getPath();
		InputStream supportPath_Stream = getClass().getResourceAsStream(supportPath);
		loadAttributesAndBeta(supportPath_Stream);
		MoleculesToInstances sts = new MoleculesToInstances(cyp, this.attributes, this.depth_neighborAtomType, this.depth_neighborAtomDescriptor);	
		this.sts = sts;
	}
	/**
	 * This is the prediction function that predicts whether each N-N bond is a BoM and save them into the output file.
	 * It returns an ArrayList<String> that contains the moleName, bond, predLabel, predProb String for each N-N bond within each input molecule 
	 * Notice that the name of the molecule is not expected to include comma (",") because we handle String as comma-split String.
	 * @param molecules
	 * @param outputPath
	 * @throws Exception
	 */
	public ArrayList<MoleResultPair> makePrediction(IAtomContainerSet molecules, Integer type) throws Exception{
		//FileWriter fw = new FileWriter(new File(outputPath));
		//fw.write("Name,Bond,Pred,ProbBoM\n");		
		ArrayList<MoleResultPair> predictedResultList = new ArrayList<>();
		for(int moleIdx = 0; moleIdx < molecules.getAtomContainerCount(); moleIdx++){			
			IAtomContainer oneMole = molecules.getAtomContainer(moleIdx);
			if(MoleculeExplorer.isInvalidCandidate(oneMole)){
				//throw new Exception("The molecule " + (moleIdx + 1) + "is not valid");
				System.out.println("The molecule " + (moleIdx + 1) + "is not valid, so it is skipped");
				continue;
			}
			ArrayList<String> resultStringList = new ArrayList<>();
			String name = getMoleName(oneMole).replace(",", "_");
			System.out.println(name);
			ArrayList<String> tempData;
			Instances data;
			if(type == 1){
				tempData = this.sts.generateTypeOneBondFeatures_OneMole(oneMole,depth_neighborAtomType, depth_neighborAtomDescriptor);
				data = this.sts.generateTypeOneBondInstances_OneMole(oneMole, tempData);
			}
			else if(type == 2){
				CreateTypeTwoThreeInstances twoThree = new CreateTypeTwoThreeInstances(this.cyp, this.attributes, 4);
				tempData = twoThree.generateTypeTwoAtomBasedFeatures(oneMole,4, false);
				data = twoThree.generateTypeTwoInstances_OneMole(oneMole,tempData,4);
			}
			else{
				CreateTypeTwoThreeInstances twoThree = new CreateTypeTwoThreeInstances(this.cyp, this.attributes, 4);
				tempData = twoThree.generateTypeThreeAtomBasedFeatures(oneMole,4, false);
				data = twoThree.generateTypeTwoInstances_OneMole(oneMole,tempData,4);
			}
			Double[] probPos = new Double[data.numInstances()];
			for(int i = 0; i<data.numInstances(); i++){
				double[] predDist = currentClassifier.distributionForInstance(data.get(i));
				probPos[i] = predDist[1];
			}
			for(int i = 0 ; i < probPos.length; i++){
				String tempString = tempData.get(i);
				String bond = tempString.split(",")[1];
				String oneResult;
				//Predict as positive
				if(probPos[i] > 1/(this.beta + 1)){
					oneResult = name +"," + bond + "," + "1.0" + "," + probPos[i] + "," + data.get(i).classValue();
				}
				//Predict as Negative
				else {
					oneResult = name +"," + bond + "," + "0.0" + "," + probPos[i] + "," + data.get(i).classValue();
				}
				
				resultStringList.add(oneResult);
			}
			MoleResultPair oneResultPair = new MoleResultPair(oneMole, resultStringList);
			predictedResultList.add(oneResultPair);
			System.out.println("Molecule " + (moleIdx+1) + " has been processed.");
		}
		return predictedResultList;
		
	}
	/**
	 * This function is used to 
	 * @param resutList
	 * @param outputPath
	 * @throws Exception
	 */
	public void generateOutput(ArrayList<MoleResultPair> resultPairList, String outputPath, String cyp, Integer type) throws Exception{
		if(outputPath.contains(".csv")){
			FileWriter fw = new FileWriter(new File(outputPath));
			String titleLine = "Name,Bond, " + cyp + "_Pred,ProbBoM" + ",realLabel" + "\n";
			fw.write(titleLine);
			for(int i = 0; i < resultPairList.size(); i++){
				MoleResultPair oneResultPair = resultPairList.get(i);
				ArrayList<String> resultList = oneResultPair.getResult();
				for(int j = 0; j < resultList.size(); j++){
					fw.write(resultList.get(j) + "\n");
				}
				
			}
			fw.close();
		}
		if(outputPath.contains(".sdf")){
			SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputPath));
			IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			for(int i = 0; i < resultPairList.size(); i++){
				MoleResultPair oneResultPair = resultPairList.get(i);
				IAtomContainer mole = oneResultPair.getMolecule();
				ArrayList<String> resultList = oneResultPair.resultList; 
				String bomString ="";
				boolean first = true;
				for(int j = 0; j < resultList.size(); j++){
					String[] parseResultLine = resultList.get(j).split(",");
					String bom = parseResultLine[2];
					if(bom.equals("1.0")){
						if(first){
							bomString = parseResultLine[1].replace(";",",");
							first = false;
						}
						else bomString += "\n" + parseResultLine[1].replace(";",",");
					}
					
				}
				mole.setProperty("BoM", bomString);
				molecules.addAtomContainer(mole);
			}
			sdfWriter.write(molecules);
			sdfWriter.close();
		}
		
	}
	/**
	 * Load the model obtained from the path to the model file
	 * @param modelPath
	 * @throws Exception
	 */
	public void loadModel(InputStream modelPath) throws Exception{
		Classifier cls = (Classifier) weka.core.SerializationHelper.read(modelPath);//Load the model
		this.currentClassifier = cls;
	}
	/**
	 * Load the selected features that are stored in .csv file in the supportPath
	 * 
	 */
	public void loadAttributesAndBeta(InputStream supportPath) throws Exception{
		BufferedReader br = new BufferedReader(new InputStreamReader(supportPath));
		String oneLine = br.readLine();//Beta
		this.beta = Double.parseDouble(oneLine);
		oneLine = br.readLine();//Should be title "Attributes"
		while((oneLine = br.readLine())!=null){			
			this.attributes.add(oneLine);
		}
		br.close();		
	}
	/**
	 * This function is used to find\generate the name of a molecule
	 * If the molecule doesn't have a name, use the InChiKey instead. 
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public String getMoleName(IAtomContainer oneMole) throws Exception{
		String name = oneMole.getID();
		if(name == null) name = oneMole.getTitle();
		if(name == null) name = oneMole.getProperty("name");
		if(name == null) name = oneMole.getProperty("Name");
		if(name == null) name = oneMole.getProperty("NAME");
		if(name == null) name = oneMole.getProperty("Title");
		if(name == null) name = (String)oneMole.getProperties().get("cdk:Title");
		if(name == null){
	        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
	        InChIGenerator gen = factory.getInChIGenerator(oneMole);
	        String inchi = gen.getInchiKey(); 
	        name = inchi;
		}
		return name;
	}
	
	/**
	 * This function will extract BoMs from the predicted result by converting string to IAtom
	 * @param predResultList
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<ArrayList<IAtom>> getTypeOneBoMList(ArrayList<MoleResultPair> predResultList, IAtomContainer oneMole) throws Exception{
		ArrayList<ArrayList<IAtom>> bomList = new ArrayList<>();
		for(int i = 0; i < predResultList.size(); i++){
			MoleResultPair oneResultPair = predResultList.get(i);
			ArrayList<String> oneResultString = oneResultPair.getResult();
			for(int j = 0; j < oneResultString.size(); j++){
				String[] parseResult = oneResultString.get(j).split(",");
				Double prob_pos = Double.parseDouble(parseResult[3]);
				if(parseResult[2].equals("1.0")){
				//<C.1;O.2>
					String[] bomIdx_string_list = parseResult[1].split(";");
					String bomIdx_left_string = bomIdx_string_list[0];
					String bomIdx_right_string = bomIdx_string_list[1];
					bomIdx_left_string = bomIdx_left_string.replace("<", "").replace(".","")
										.replace("Cl", "").replace("O", "").replace("S", "").replace("N", "").replace("P", "")
										.replace("B", "").replace("F", "").replace("C", "");
					bomIdx_right_string = bomIdx_right_string.replace(">", "").replace(".","")
										.replace("Cl", "").replace("O", "").replace("S", "").replace("N", "").replace("P", "")
										.replace("B", "").replace("F", "").replace("C", "");
					
					Integer bomIdx_left = Integer.parseInt(bomIdx_left_string) - 1;
					Integer bomIdx_right = Integer.parseInt(bomIdx_right_string) - 1;
					IAtom atom_first = oneMole.getAtom(bomIdx_left);
					IAtom atom_second = oneMole.getAtom(bomIdx_right);
					if(bomIdx_left > bomIdx_right){
						atom_first = oneMole.getAtom(bomIdx_right);
						atom_second = oneMole.getAtom(bomIdx_left);
					}
					ArrayList<IAtom> one_typeOne_BoM = new ArrayList<>();
					atom_first.setProperty("Score",prob_pos);
					atom_second.setProperty("Score", prob_pos);
					one_typeOne_BoM.add(atom_first);
					one_typeOne_BoM.add(atom_second);		
					if(!bomList.contains(one_typeOne_BoM)) bomList.add(one_typeOne_BoM);
				}
			}			
		}
		return bomList;
	}
	
	/**
	 * This function will extract BoMs from the predicted result by converting string to IAtom
	 * @param predResultList
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<IAtom> getTypeTwoBoMList(ArrayList<MoleResultPair> predResultList, IAtomContainer oneMole) throws Exception{		
		ArrayList<IAtom> bomList = new ArrayList<>();
		for(int i = 0; i < predResultList.size(); i++){
			MoleResultPair oneResultPair = predResultList.get(i);
			ArrayList<String> oneResultString = oneResultPair.getResult();
			oneResultPair.getMolecule();
			for(int j = 0; j < oneResultString.size(); j++){
				String[] parseResult = oneResultString.get(j).split(",");
				Double prob_pos = Double.parseDouble(parseResult[3]);
				if(parseResult[2].equals("1.0")){
					String bomIdx_string = parseResult[1].replace("<C.", "").replace(";H>", "");
					bomIdx_string = bomIdx_string.replace("<S.", "").replace("<N.", "").replace("<P.", "");
					Integer bomIdx = Integer.parseInt(bomIdx_string) - 1;
					IAtom oneAtom = oneMole.getAtom(bomIdx);
					oneAtom.setProperty("Score", prob_pos);				
					if(!bomList.contains(oneAtom)) bomList.add(oneAtom);
				}
			}
		}
		return bomList;
	}
	
	/**
	 * This function will extract BoMs from the predicted result by converting string to IAtom
	 * @param predResultList
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<IAtom> getTypeThreeBoMList(ArrayList<MoleResultPair> predResultList, IAtomContainer oneMole) throws Exception{
		ArrayList<IAtom> bomList = new ArrayList<>();
		for(int i = 0; i < predResultList.size(); i++){
			MoleResultPair oneResultPair = predResultList.get(i);
			ArrayList<String> oneResultString = oneResultPair.getResult();
			for(int j = 0; j < oneResultString.size(); j++){
				String[] parseResult = oneResultString.get(j).split(",");
				Double prob_pos = Double.parseDouble(parseResult[3]);
				if(parseResult[2].equals("1.0")){
					String bomIdx_string = parseResult[1].replace("<", "").replace(";S>", "").replace(";N>", "").replace(";P>", "");
					Integer bomIdx = Integer.parseInt(bomIdx_string) - 1;
					IAtom oneAtom = oneMole.getAtom(bomIdx);
					oneAtom.setProperty("Score", prob_pos);
					if(!bomList.contains(oneAtom)) bomList.add(oneAtom);
				}
			}
		}
		return bomList;
	}
	
	/**
	 * Write input molecules with predicted results into the target sdf file
	 * @param rawMolecules: The original input compounds
	 * @param outSdfPath: The path to the target sdf file
	 * @param predictedResult: predicted results for all compounds
	 * @throws Exception
	 */

	public static void outputResultIAtomContainerSet(IAtomContainerSet predictedResult, String outSdfPath) throws Exception{
		IAtomContainerSet resultMole = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outSdfPath));	
		if(predictedResult == null || predictedResult.isEmpty()){
			sdfWriter.write(resultMole);
			sdfWriter.close();
			return;
		}
		
		for(int i = 0; i < predictedResult.getAtomContainerCount(); i++){			
			IAtomContainer outMole = predictedResult.getAtomContainer(i);
			AtomContainerManipulator.removeHydrogens(outMole);
			String smile = sg.create(outMole);
			IAtomContainer mole = sp.parseSmiles(smile);
			Map<Object,Object> properties = new HashMap<Object, Object>();
			properties = outMole.getProperties();
			mole.setProperties(properties);
			sdfWriter.write(mole);
			
		}
		
		sdfWriter.close();
	}
}
