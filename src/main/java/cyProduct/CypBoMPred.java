package cyProduct;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.SDFWriter;

import utils.MoleculeExplorer;
import utils.ReadMolecules;
import weka.classifiers.Classifier;
import weka.core.Instances;

public class CypBoMPred {
	Classifier currentClassifier;
	MoleculesToInstances sts;
	ArrayList<String> attributes = new ArrayList<String>();
	Double beta;
	int depth_neighborAtomType = 4;
	int depth_neighborAtomDescriptor = 2;
	public static void main(String[] args) throws Exception{
		//Inputs: Cyp, InputPath, OutputPath
		if(args.length<3){
			System.out.println("Don't have enough arguments");
			return;
		}
		else if(args.length>3){
			System.out.println("# of arguments is more than 4");
			return;			
		}
		String cyp = args[0];
		String inputPath = args[1];
		String outputPath = args[2];
		 
		System.out.println("Running CypBoMPred Software");
		String current_dir = System.getProperty("user.dir");
		String supportFolderPath = String.format("%s/%s/%s/", current_dir,"supportfiles","CYP"+cyp);
		CypBoMPred bomPred = new CypBoMPred();
		bomPred.setup(cyp, supportFolderPath);				
		IAtomContainerSet molecules = ReadMolecules.extractMoleculesFromFile(inputPath);		
		ArrayList<MoleResultPair> predResultList = bomPred.makePrediction(molecules);
		bomPred.generateOutput(predResultList, outputPath, cyp);
		
		
	}
	
	/**
	 * Set up the model and support file that are used to make prediction
	 * @param cyp
	 * @param inputPath
	 * @param supportFolderPath
	 * @throws Exception
	 */
	public void setup(String cyp,String supportFolderPath) throws Exception{
		String modelPath = supportFolderPath +  cyp + ".model"; 
		String supportPath = supportFolderPath + cyp + "_supportfiles.csv";
		loadModel(modelPath);
		loadAttributesAndBeta(supportPath);
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
	public ArrayList<MoleResultPair> makePrediction(IAtomContainerSet molecules) throws Exception{
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
			ArrayList<String> tempData = this.sts.generateTypeOneBondFeatures_OneMole(oneMole,depth_neighborAtomType, depth_neighborAtomDescriptor);
			Instances data = this.sts.generateTypeOneBondInstances_OneMole(oneMole, tempData);
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
					oneResult = name +"," + bond + "," + "1.0" + "," + probPos[i];
				}
				//Predict as Negative
				else {
					oneResult = name +"," + bond + "," + "0.0" + "," + probPos[i];
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
	public void generateOutput(ArrayList<MoleResultPair> resultPairList, String outputPath, String cyp) throws Exception{
		if(outputPath.contains(".csv")){
			FileWriter fw = new FileWriter(new File(outputPath));
			String titleLine = "Name,Bond, " + cyp + "_Pred,ProbBoM\n";
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
	public void loadModel(String modelPath) throws Exception{
		Classifier cls = (Classifier) weka.core.SerializationHelper.read(modelPath);//Load the model
		this.currentClassifier = cls;
	}
	/**
	 * Load the selected features that are stored in .csv file in the supportPath
	 * 
	 */
	public void loadAttributesAndBeta(String supportPath) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(new File(supportPath)));
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

}
