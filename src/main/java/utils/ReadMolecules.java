package utils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Properties;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class ReadMolecules {

	/**
	 * Create IAtomContainerSet containing all molecules in the inputFiles(either sdf or smiles).
	 * @param inputPath
	 * @return IAtomContainerSet containing all molecules with ExplicitHydrogens added.
	 * @throws Exception
	 */	
	public static IAtomContainerSet extractMoleculesFromFile(String inputPath) throws Exception{
		IAtomContainerSet moleculeSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
		SmilesParser sp = new SmilesParser(builder);
		
		String oneLine;
		//If the input is SMILEs
		if(inputPath.contains(".csv")){
			FileReader fr = new FileReader(inputPath);
			BufferedReader br = new BufferedReader(fr);
			while((oneLine = br.readLine())!=null){
				//Create IAtomContainer molecule from smiles
				IAtomContainer mol = sp.parseSmiles(oneLine);
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
				CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
				Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
				aromaticity.apply(mol);
				adder.addImplicitHydrogens(mol);
				IChemObjectBuilder builder_1 = SilentChemObjectBuilder.getInstance();
				ModelBuilder3D mb3d = ModelBuilder3D.getInstance(builder_1);
				IAtomContainer oneMolecule = mb3d.generate3DCoordinates(mol, false);
//				StructureDiagramGenerator sdg = new StructureDiagramGenerator();
//				sdg.setMolecule(mol);
//				sdg.generateCoordinates();
//				IAtomContainer oneMolecule = sdg.getMolecule();
				AtomContainerManipulator.suppressHydrogens(oneMolecule);
				moleculeSet.addAtomContainer(oneMolecule);
			}
			br.close();

		}
		//if the input file is a sdf file
		else if(inputPath.contains(".sdf")){
			 moleculeSet = readFile(inputPath);
		}
		//If the input is a SMILES string
		else if(inputPath.contains("SMILES=")){
			IChemObjectBuilder builder_1 = SilentChemObjectBuilder.getInstance();
			ModelBuilder3D mb3d = ModelBuilder3D.getInstance(builder_1);
			
			String inputSMILE = inputPath.replace("SMILES=", "");
			IAtomContainer mol = sp.parseSmiles(inputSMILE);//The inputPath should be a SMILEs String in this case
			System.out.println("The input SMILES string is: " + inputSMILE);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
			Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
			aromaticity.apply(mol);
			adder.addImplicitHydrogens(mol);
			IAtomContainer oneMolecule = mb3d.generate3DCoordinates(mol, false);
//			StructureDiagramGenerator sdg = new StructureDiagramGenerator();
//			sdg.setMolecule(mol);
//			sdg.generateCoordinates();			
//			IAtomContainer oneMolecule = sdg.getMolecule();
			AtomContainerManipulator.suppressHydrogens(oneMolecule);
			moleculeSet.addAtomContainer(oneMolecule);
		}

		return moleculeSet;
	}
	/**
	 * Read the sdf file and generate a IAtomContainerSet that contains all molecules in it.
	 * @param pathToInputFile
	 * @return
	 * @throws FileNotFoundException
	 * @throws CDKException
	 */
	public static IAtomContainerSet readFile(String pathToInputFile) throws FileNotFoundException, CDKException, Exception {
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
		while (sdfr.hasNext()){
			IAtomContainer oneMole = sdfr.next();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(oneMole.getBuilder());
			Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
			aromaticity.apply(oneMole);
			adder.addImplicitHydrogens(oneMole);
			AtomContainerManipulator.suppressHydrogens(oneMole);
			allMoles.addAtomContainer(oneMole);
		}
		sdfr.close();
		return allMoles;
	}
	
}
