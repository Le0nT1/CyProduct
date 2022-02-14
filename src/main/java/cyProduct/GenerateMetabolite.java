package cyProduct;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.tools.StructureResonanceGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import utils.MergeFragments;

/**
 * The reaction types here are:
 * 1. CarbonOxygenOxidation
 * 2. DeAlkylation
 * 3. EpOxidation
 * 4. Hydroxylation
 * 5. NitrosoReduction
 * 6. RingRearrangement
 * 7. SpecialDealkylaiton
 * 8. SPN-Oxidation
 * 9. Dehydrogenation
 * 10. Special Reaction one: X=C-C -> X-C(=O)-C
 * 11. Special Reaction Two:N=C-N -> N-C(=O)-N
 * 12. Desulfurization
 * 13. Hydrolysis
 * 14. Dehalogon
 * @author Tian
 *
 */
public class GenerateMetabolite {
    public static SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
    public static final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
    static{
	    smrkMan.setFlagApplyStereoTransformation(false);
		smrkMan.setFlagCheckResultStereo(true);
		smrkMan.setFlagFilterEquivalentMappings(true);
		smrkMan.setFlagProcessResultStructures(true);
		smrkMan.setFlagAddImplicitHAtomsOnResultProcess(true);
    }
//	public static void main(String[] args) throws Exception{
//		SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
//		/**
//		 * Test hydroxylation and Oxidation
//		 */
//		String smile = "CCCCCCC";
//		IAtomContainer oneMole = sp.parseSmiles(smile);
//		UniqueIDFunctionSet.assignUniqueID(oneMole);
//		generateHydroxylationMetabolite(oneMole.getAtom(3), oneMole);		
		/**
		 * Test Carbon-Oxygen Oxidation
		 */
//		smile = "CCCCCCO";
//		oneMole = sp.parseSmiles(smile);
//		UniqueIDFunctionSet.assignUniqueID(oneMole);
//		generateCarbonOxygenOxidationMetabolite(oneMole.getAtom(5), oneMole.getAtom(6), oneMole);
		
		
		/**
		 * Test O-,S-,N-Dealylation 
		 */
//		smile = "CCCCOCC";
//		oneMole = sp.parseSmiles(smile);
//		UniqueIDFunctionSet.assignUniqueID(oneMole);
//		generateDeAlkylationMetabolite(oneMole.getAtom(3), oneMole.getAtom(4), oneMole, "O-Dealkylation");
//		
//		smile = "CCCCSCC";
//		oneMole = sp.parseSmiles(smile);
//		UniqueIDFunctionSet.assignUniqueID(oneMole);
//		generateDeAlkylationMetabolite(oneMole.getAtom(3), oneMole.getAtom(4), oneMole, "S-Dealkylation");
//		
//		smile = "CN1C(N)=NC2=C1C=CC1=NC=CC=C21";
//		oneMole = sp.parseSmiles(smile);
//		UniqueIDFunctionSet.assignUniqueID(oneMole);
//		generateDeAlkylationMetabolite(oneMole.getAtom(0), oneMole.getAtom(1), oneMole, "N-Dealkylation");
		
		/**
		 * Test EpOxidation
		 */
//		smile = "CC1=CC2=C3C(C=CC4=C3C(C=C2)=CC=C4)=C1";
//				//"CC1=CC2=C3C(C=CC4=C3C(C=C2)=CC=C4)=C1";//"NC1=CC2=CC=CC=C2C=C1";//"NC1=CC2=CC=CC=C2C=C1";//"NC1=CC=CC=C1";//"c1ccccc1";//"C1=CC=CC=C1";//"CCC\\C=C\\C";
//		//String inputSDFPath = "C:/Users/Tian/Desktop/BioData/SOM-React/BioTransformerDB/WaitToMerge/Merged/TypeTwoFinal/TestSample.sdf";
//		//IAtomContainerSet molecules = Utilities.readFile(inputSDFPath);		
//		oneMole = sp.parseSmiles(smile);
//		//oneMole = molecules.getAtomContainer(0);
//		UniqueIDFunctionSet.assignUniqueID(oneMole);
//		for(int i = 0; i < oneMole.getAtomCount(); i++){
//			if(oneMole.getAtom(i).getSymbol().equals("P")){
//				//System.out.println((Integer) oneMole.getAtom(i).getProperty("UniqueID"));
//			}
//		}
//		generateEpOxidationMetabolite(oneMole.getAtom(3), oneMole.getAtom(4), oneMole, 0);//1,2; 6,9
		
		/**
		 * Test nitroso-reduction reaciton
		 */
//		smile = "C[N+](O)=O";
//		oneMole = sp.parseSmiles(smile);
//		UniqueIDFunctionSet.assignUniqueID(oneMole);
//		ArrayList<ArrayList<IAtom>> bom_List = new ArrayList<>();
//		ArrayList<IAtom> bom_one = new ArrayList<>();
//		ArrayList<IAtom> bom_two = new ArrayList<>();
//		bom_one.add(oneMole.getAtom(1));
//		bom_one.add(oneMole.getAtom(2));
//		bom_two.add(oneMole.getAtom(1));
//		bom_two.add(oneMole.getAtom(3));
//		bom_List.add(bom_one);
//		bom_List.add(bom_two);
//		generateNitrosoReductionMetabolite(bom_List, oneMole);
		/**
		 * Test RingRearrangement
		 */
//		smile = "[CH3:18][O:17][CH2:16][CH2:15][O:14][C:12](=[O:13])[C:11]1=[C:10]([CH3:29])[NH:9][C:8]([CH3:30])=[C:7]([CH:19]1[C:20]1=[CH:21][CH:22]=[CH:23][C:24](=[CH:28]1)[N+:25]([O-:26])=[O:27])[C:2](=[O:1])[O:3][CH:4]([CH3:5])[CH3:6]";		
//		
//		/**
//		 * Test O-C-O dealylation
//		 */
//		
//		
//		/**
//		 * Test SPN-Oxidation
//		 */
//	
//	}
	/**
	 * This function will generate the IAtomContainerSet that contains all the resonance structures for the given molecule
	 * @param molecule
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateResonanceStructure(IAtomContainer molecule) throws Exception{
		 SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric); 
		 StructureResonanceGenerator gRI = new StructureResonanceGenerator();
		 LonePairElectronChecker lpec = new LonePairElectronChecker();
		 AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		 CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		 adder.addImplicitHydrogens(molecule);
		 lpec.saturate(molecule);
		 IAtomContainerSet resonanceStructures = gRI.getStructures(molecule);
		 return resonanceStructures;
	}
	/**
	 * This function will generate the hydroxylated metabolite for the input oneMole on the targetAtom (reactive site)
	 * @param targetAtom
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateHydroxylationMetabolite(IAtom targetAtom, IAtomContainer oneMole) throws Exception{		
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		IAtomContainer copy = oneMole.clone();
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		IAtom targetCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(targetAtom), copy);
		targetCopy.setProperty("ReactiveSite", "Hydroxylation");
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(targetCopy);
		for(int i = 0; i < copy.getConnectedAtomsList(targetCopy).size(); i++){
			if(copy.getConnectedAtomsList(targetCopy).get(i).getSymbol().equalsIgnoreCase("H")){
				fragmentStructure.addAtom(copy.getConnectedAtomsList(targetCopy).get(i));
				fragmentStructure.addBond(copy.getBond(copy.getConnectedAtomsList(targetCopy).get(i), targetCopy));
				
			}
		}
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String hydroxylation_SMIKRS = "[*:1][H]>>[H][#8]-[*:1]";//"[#6;A:1][H]>>[H][#8][#6;A:1]";		
		
		SMIRKSReaction reaction = smrkMan.parse(hydroxylation_SMIKRS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart.getAtomContainerCount() > 1) throw new Exception("There are more than one Hydroxyl metabolite for the same reactive site within the molecule.");		
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		if(resultSet.getAtomContainerCount() > 1) throw new Exception("There are more than one Hydroxyl metabolite for the same reactive site within the molecule after partition.");		
		//System.out.println("SMILES of the changed part: " + sg.create(resultSet.getAtomContainer(0)));		
		IAtomContainer metabolite = mergeFragments(copy, changedPart.getAtomContainer(0));
		Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
		properties.put("ReactionType", "Hydroxylation");
		properties.put("Score", targetAtom.getProperty("Score"));
		properties.put("BoMs", (targetAtom.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(targetAtom)));
		metabolite.addProperties(properties);
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		metabolites.addAtomContainer(metabolite);
		return metabolites;
		
	}

	/**
	 * This function will generate the Oxidated metabolite for the input oneMole on the targetAtom (reactive site)
	 * @param targetAtom
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateSPNOxidationMetabolite(IAtom targetAtom, IAtomContainer oneMole) throws Exception{		
		//SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		IAtomContainer copy = oneMole.clone();
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		IAtom targetCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(targetAtom), copy);
		targetCopy.setProperty("ReactiveSite", "SNP-Oxidation");
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(targetCopy);
		//IAtom newAtom = new Atom(targetAtom.getSymbol());
		//fragmentStructure.addAtom(newAtom);
		for(int i = 0; i < copy.getConnectedAtomsList(targetCopy).size(); i++){
			if(copy.getConnectedAtomsList(targetCopy).get(i).getSymbol().equalsIgnoreCase("H")){
				fragmentStructure.addAtom(copy.getConnectedAtomsList(targetCopy).get(i));
				//IBond newBond = copy.getBond(copy.getConnectedAtomsList(targetCopy).get(i), newAtom);
				fragmentStructure.addBond(copy.getBond(copy.getConnectedAtomsList(targetCopy).get(i), targetCopy));
				
			}
		}
		int a = targetCopy.getValency();
		double b = copy.getBondOrderSum(targetCopy);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String snp_Oxidation_SMIKRS = "[*:1]>>O=[*:1]";		
		if(targetCopy.getSymbol().equals("N")){
			double c = copy.getBondOrderSum(targetCopy);
			if(copy.getBondOrderSum(targetCopy) == 2.0){
				snp_Oxidation_SMIKRS = "[#7:1]>>[#7:1]=O";// "[#7;A:1]>>[#7;A+:1]=O";
			}
			else if(copy.getBondOrderSum(targetCopy)==3.0){
				snp_Oxidation_SMIKRS = "[#7:1]>>[#7+:1]-[#8-]"; //"[#7;A:1]>>[#7;A+:1][#8-]";
			}			
		}
		
		SMIRKSReaction reaction = smrkMan.parse(snp_Oxidation_SMIKRS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		//if(changedPart.getAtomContainerCount() > 1) throw new Exception("There are more than one SNP-Oxidation metabolite for the same reactive site within the molecule.");		
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		//if(resultSet.getAtomContainerCount() > 1) throw new Exception("There are more than one SNP-Oxidation metabolite for the same reactive site within the molecule after partition.");		
		//System.out.println("SMILES of the changed part: " + sg.create(resultSet.getAtomContainer(0)));		
		IAtomContainer metabolite = mergeFragments(copy, changedPart.getAtomContainer(0));
		Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
		properties.put("ReactionType", "SNP-Oxidation");
		properties.put("Score", targetAtom.getProperty("Score"));
		metabolite.addProperties(properties);
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		metabolites.addAtomContainer(metabolite);
		return metabolites;
		
	}
	/**
	 * This function will generate the oxidized metabolite for the input oneMole by Oxidizing the target C-O bond
	 * @param targetAtom
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateCarbonOxygenOxidationMetabolite(IAtom leftAtom, IAtom rightAtom, IAtomContainer oneMole) throws Exception{		
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		IAtomContainer copy = oneMole.clone();		
		
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), copy);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), copy);
		IBond targetBond = copy.getBond(leftCopy, rightCopy);
		if(targetBond.getOrder() != IBond.Order.SINGLE) throw new Exception("This target bond doesn't satisfy the requirement for the C-O oxidation reaction");
		
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		//IBond targetBond = copy.getBond(leftCopy, rightCopy);
		//System.out.println("Indices of the BoM atoms: " + copy.indexOf(leftCopy) + "," + copy.indexOf(rightCopy));
		//System.out.println("UniqueIDs of the BoM atoms: " + UniqueIDFunctionSet.getUniqID(leftAtom) + "," + UniqueIDFunctionSet.getUniqID(rightAtom));
		leftCopy.setProperty("ReactiveSite", "Oxidation");
		rightCopy.setProperty("ReactiveSite", "Oxidation");
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(leftCopy);		
		fragmentStructure.addAtom(rightCopy);
		fragmentStructure.addBond(targetBond);//(copy.getBond(leftCopy, rightCopy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(leftCopy, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(rightCopy, copy));
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String hydroxylation_SMIKRS = "[H][#8:2][#6;A:1]>>[#6;A:1]=[O:2]";		
		
		SMIRKSReaction reaction = smrkMan.parse(hydroxylation_SMIKRS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart.getAtomContainerCount() > 1) throw new Exception("There are more than one Oxidized metabolite for the same reactive site within the molecule.");		
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		if(resultSet.getAtomContainerCount() > 1) throw new Exception("There are more than one Oxidized metabolite for the same reactive site within the molecule after partition.");		
		//System.out.println("SMILES of the changed part: " + sg.create(resultSet.getAtomContainer(0)));		
		IAtomContainer metabolite = mergeFragments_EpOxidation(copy, changedPart.getAtomContainer(0));
		Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
		properties.put("ReactionType", "Oxidation");
		properties.put("Score", leftAtom.getProperty("Score"));
		properties.put("BoMs", (leftAtom.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(leftAtom) + ";" + rightAtom.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(rightAtom)));
		metabolite.addProperties(properties);
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		metabolites.addAtomContainer(metabolite);
		return metabolites;
		
	}
	/**
	 * This function will generate the DeAlkylation metabolite for the Dealkylation reaction
	 * As Dealkylation involves the breakage of a Type One <eta,eta> bond, it input requires both atoms connected by the bond
	 * The dealkylattionType means "N,O,etc" deakylation and will converted to the corresponding SMIRKS String
	 * @param leftAtom
	 * @param rightAtom
	 * @param oneMole
	 * @param dealkylationType
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateDeAlkylationMetabolite(IAtom leftAtom, IAtom rightAtom, IAtomContainer oneMole, String dealkylationType) throws Exception{
		//SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		IAtomContainer copy = oneMole.clone();		
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), copy);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), copy);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		//IBond targetBond = copy.getBond(leftCopy, rightCopy);
		leftCopy.setProperty("ReactiveSite", "Dealkyaltion");
		rightCopy.setProperty("ReactiveSite", "Dealkyaltion");
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		IAtom newLeft = new Atom(leftCopy.getSymbol());
		IAtom newRight = new Atom(rightCopy.getSymbol());
		newLeft.setProperties(leftCopy.getProperties());
		newRight.setProperties(rightCopy.getProperties());
		fragmentStructure.addAtom(newLeft);
		fragmentStructure.addAtom(newRight);
		IAtom[] bondAtoms = {newLeft, newRight}; 
		IBond targetBond = new Bond(bondAtoms, copy.getBond(leftCopy, rightCopy).getOrder());
		fragmentStructure.addBond(targetBond);
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(fragmentStructure.getAtom(0), copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(fragmentStructure.getAtom(1), copy));
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String dealkyltion_SMIRKS;
		if(dealkylationType.equals("O-Dealkylation")) dealkyltion_SMIRKS = "[H][#6:2][#8;A:1]>>[#8;A:1][H].[#6:2]=O";//[H][#6:2](-*)-[#8:1]-*>>[#8:1]-*.*-[#6:2]=O";
		else if(dealkylationType.equals("N-Dealkylation")) dealkyltion_SMIRKS = "[H][#6;A:2][#7;A:1]>>[#7;A:1][H].[#6;A:2]=O";//"[H][#6;A:2][#7;A:1]>>[#7;A:1][H].[#6;A:2]=O";
		else if(dealkylationType.equals("S-Dealkylation")) dealkyltion_SMIRKS = "[H][#6;A:2][#16;A:1]>>[#16;A:1][H].[#6;A:2]=O"; 
		//If none of the dealkylationType matches the input, throws Exception
		else throw new Exception("The input Dealkylation type can not be identified");//"[H][#6:2](-*)-[#8:1]-*>>[#8:1]-*.*-[#6:2]=O";
		
		SMIRKSReaction reaction = smrkMan.parse(dealkyltion_SMIRKS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart.getAtomContainerCount() > 2) throw new Exception("There are more than two dealkylation metabolite for the same reactive site within the molecule.");		
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		if(resultSet.getAtomContainerCount() > 2) throw new Exception("There are more than two dealkylation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		/**
		 * We remove the bond <leftCopy, rightCopy> within the molecule 
		 */
		copy.removeBond(leftCopy, rightCopy);
		IAtomContainerSet subFragments = ConnectivityChecker.partitionIntoMolecules(copy);
		//System.out.println("SMILES of the check part: " + sg.create(copy));	
		//The part contains "O" will be at index = 0. There should be only two subFragment for Dealkylation reaction.
		ArrayList<IAtomContainer> reorderedFragments = reorderFragments(subFragments);
		ArrayList<IAtomContainer> reorderedResultSet = reorderFragments(resultSet);
		for(int i = 0; i < reorderedResultSet.size(); i++){
			//System.out.println("SMILES of the changed part: " + sg.create(reorderedResultSet.get(i)));	
			IAtomContainer metabolite = mergeFragments(reorderedFragments.get(i), reorderedResultSet.get(i));			
			Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
			properties.put("ReactionType", dealkylationType);
			//Use the score of the typeOne BoM when both TypeOne and TypeTwo BoMs are included in one reaction
			if(leftAtom.getSymbol().equalsIgnoreCase("C")) properties.put("Score", rightAtom.getProperty("Score"));
			else properties.put("Score", leftAtom.getProperty("Score"));
			properties.put("BoMs", (leftAtom.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(leftAtom) + ";" + rightAtom.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(rightAtom)));
			metabolite.addProperties(properties);
			metabolites.addAtomContainer(metabolite);
		}				
		return metabolites;
	
	}
	
	/**
	 * An Epoxidation reaction requires an <eta,eta> bond to be modified.
	 * Hence both the atoms connected by the reactive Bond should be used as input
	 * @param leftAtom
	 * @param rightAtom
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateEpOxidationMetabolite(IAtom leftAtom, IAtom rightAtom, IAtomContainer oneMole, int counter) throws Exception{
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		IAtomContainer copy = oneMole.clone();		
		String oringinalSmile = sg.create(oneMole);
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), copy);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), copy);
		IBond targetBond = copy.getBond(leftCopy, rightCopy);
		if(targetBond.getOrder() == IBond.Order.SINGLE){
			IAtomContainerSet resonanceStructureSet = generateResonanceStructure(copy);
			for(int i = 0; i < resonanceStructureSet.getAtomContainerCount(); i++){
				IAtomContainer oneResonaceStruct = resonanceStructureSet.getAtomContainer(i);
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneResonaceStruct);
				leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), oneResonaceStruct);
				rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), oneResonaceStruct);
				targetBond = oneResonaceStruct.getBond(leftCopy, rightCopy);
				if(targetBond.getOrder() == IBond.Order.DOUBLE){
					String smile = sg.create(oneResonaceStruct);
					if(!oringinalSmile.contains("+") && !oringinalSmile.contains("-") && 
						!smile.contains("+") && !smile.contains("-")){
						copy = oneResonaceStruct;
						//System.out.println("selected resonance struct: " + sg.create(oneResonaceStruct));
						break;
					}
					else{
						copy = oneResonaceStruct;
						//System.out.println("selected resonance struct: " + sg.create(oneResonaceStruct));
						break;
					}
				}
			}
		}
		
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		//IBond targetBond = copy.getBond(leftCopy, rightCopy);
		//System.out.println("Indices of the BoM atoms: " + copy.indexOf(leftCopy) + "," + copy.indexOf(rightCopy));
		//System.out.println("UniqueIDs of the BoM atoms: " + UniqueIDFunctionSet.getUniqID(leftAtom) + "," + UniqueIDFunctionSet.getUniqID(rightAtom));
		leftCopy.setProperty("ReactiveSite", "EpOxidation");
		rightCopy.setProperty("ReactiveSite", "EpOxidation");
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(leftCopy);		
		fragmentStructure.addAtom(rightCopy);
		IAtom[] bondAtoms = {leftCopy, rightCopy};
		targetBond = new Bond(bondAtoms, IBond.Order.DOUBLE);
		//IBond tt = copy.getBond(leftCopy, rightCopy);
		//tt = newBond;
		fragmentStructure.addBond(targetBond);//(copy.getBond(leftCopy, rightCopy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(leftCopy, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(rightCopy, copy));
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String epOxidaiton_SMIRS = "[*:1]=[*:2]>>[#8]-1-[*:1]-[*:2]-1";//"[#6;A:2]=[#6;A:1]>>[#6;A:1]1[#6;A:2][#8]1";
	
		SMIRKSReaction reaction = smrkMan.parse(epOxidaiton_SMIRS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart.getAtomContainerCount() > 1) throw new Exception("There are more than one epOxidation metabolite for the same reactive site within the molecule.");
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		if(resultSet.getAtomContainerCount() > 1) throw new Exception("There are more than one epOxidation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		/**
		 * We remove the bond <leftCopy, rightCopy> within the molecule 
		 */
		copy.removeBond(leftCopy, rightCopy);
		//System.out.println("SMILES of the changed part: " + sg.create(resultSet.getAtomContainer(0)));	
		ArrayList<IAtom>  removedAtoms = MergeFragments.getRemovedAtomList(fragmentStructure, resultSet.getAtomContainer(0));
		IAtomContainer metabolite = MergeFragments.mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0), removedAtoms);//mergeFragments_EpOxidation(copy, resultSet.getAtomContainer(0));	
		Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
		properties.put("ReactionType", "EpOxidation");
		properties.put("Score", leftAtom.getProperty("Score"));
		metabolite.addProperties(properties);
		metabolites.addAtomContainer(metabolite);			
		return metabolites;
		
	}
	/**
	 * This function will reduce N+OO- to NH2
	 * @param nitrogen
	 * @param oxygen_one
	 * @param oxygen_two
	 * @param oneMole
	 * @return
	 */
	public static IAtomContainerSet generateNitrosoReductionMetabolite(ArrayList<ArrayList<IAtom>> bom_List, IAtomContainer oneMole) throws Exception{
		if(bom_List == null || bom_List.size() != 2){
			//System.out.println("The input BoMs for the nitrosoReductionReaction doesn't satisfy the corresponding conditions.");
			return null;
		}
		IAtomContainer copy = oneMole.clone();
		ArrayList<IAtom> bom_one = bom_List.get(0);
		ArrayList<IAtom> bom_two = bom_List.get(1);
		
		IAtom atom_nitrogen, atom_two, atom_three;
		
		if(bom_one.get(0).getSymbol().equalsIgnoreCase("N")){
			atom_nitrogen = bom_one.get(0);
			atom_two = bom_one.get(1);
		}
		else{
			atom_nitrogen = bom_one.get(1);
			atom_two = bom_one.get(0);
		}
		if(bom_two.get(0).getSymbol().equalsIgnoreCase("O")) atom_three = bom_two.get(0);
		else atom_three = bom_two.get(1);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		IAtom copy_atomOne = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_nitrogen), copy);
		IAtom copy_atomTwo = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_two), copy);
		IAtom copy_atomThree = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_three), copy);
		copy_atomOne.setProperty("ReactiveSite", "NitrosoReduction");
		copy_atomTwo.setProperty("ReactiveSite", "NitrosoReduction");
		copy_atomThree.setProperty("ReactiveSite", "NitrosoReduction");
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(copy_atomOne);
		fragmentStructure.addAtom(copy_atomTwo);
		fragmentStructure.addAtom(copy_atomThree);
		fragmentStructure.addBond(copy.getBond(copy_atomOne, copy_atomTwo));
		fragmentStructure.addBond(copy.getBond(copy_atomOne, copy_atomThree));
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(copy_atomOne, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(copy_atomTwo, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(copy_atomThree, copy));
	
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String nitroso_reduction_SMARTS = "[#8;A:3][#7;A+:1]=[#8;A:2]>>[#7;A:1].[#8:2].[#8:3]";
		//"[H][#8:3][#7;A+:1]=[O:2]>>[#7;A:1].[#8:2].[#8:3]";
		//"[#8;A:3][#7;A+:1]=[#8;A:2]>>[H:2][#7;A:1][H:3]";//"[#8;A][#7;A+:1]=[#8;A]>>[#7;A:1]";
		//"[#8;A:3][#7;A+:1]=[#8;A:2]>>[H:3][#7;A+:1][H:2]";
		SMIRKSReaction reaction = smrkMan.parse(nitroso_reduction_SMARTS);		
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		IAtomContainer changedPart_Organic = null;// = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		for(int i = 0; i < changedPart.getAtomContainerCount(); i++){
			IAtomContainer oneChangedPart = resultSet.getAtomContainer(i);
			for(int j = 0; j < oneChangedPart.getAtomCount(); j++){
				if(oneChangedPart.getAtom(j).getSymbol().equalsIgnoreCase("N")){
					changedPart_Organic = oneChangedPart;
					break;
				}
			}
		}
		if(changedPart_Organic == null) throw new Exception("No organic metabolite was generated for this Nitroso-Reduction reaction");
		//System.out.println("SMILES of the changed part: " + sg.create(changedPart_Organic));	
		IAtomContainer metabolite = mergeFragmentsForReductionReaction(copy, changedPart_Organic);
		Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
		properties.put("ReactionType", "NitrosoReduction");
		Double score = atom_nitrogen.getProperty("Score");
		if(score == null) score = 0.0;
		else if(score < (Double) atom_two.getProperty("Score")) score = atom_two.getProperty("Score");
		else if(score < (Double) atom_three.getProperty("Score")) score = atom_three.getProperty("Score");
		properties.put("Score", score);
		metabolite.addProperties(properties);
		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		metabolites.addAtomContainer(metabolite);
		return metabolites;
		
	}
	/**
	 * This function will generate metabolite for the R-O-C-O-R within the ring patter of O-Dealkyaltion
	 * @param bom_List
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateOCODealkylationMetabolite(ArrayList<ArrayList<IAtom>> bom_List, IAtomContainer oneMole) throws Exception{
		if(bom_List == null || bom_List.size() != 2){
			//System.out.println("The input BoMs for the specialDealkylaiton reaction doesn't satisfy the corresponding conditions.");
			return null;
		}
		IAtomContainer copy = oneMole.clone();
		ArrayList<IAtom> bom_one = bom_List.get(0);
		ArrayList<IAtom> bom_two = bom_List.get(1);
		
		IAtom atom_carbon, atom_two, atom_three;
		
		if(bom_one.get(0).getSymbol().equalsIgnoreCase("C")){
			atom_carbon = bom_one.get(0);
			atom_two = bom_one.get(1);
		}
		else{
			atom_carbon = bom_one.get(1);
			atom_two = bom_one.get(0);
		}
		if(bom_two.get(0).getSymbol().equalsIgnoreCase("O")) atom_three = bom_two.get(0);
		else atom_three = bom_two.get(1);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		IAtom copy_atomOne = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_carbon), copy);
		IAtom copy_atomTwo = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_two), copy);
		IAtom copy_atomThree = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_three), copy);
		copy_atomOne.setProperty("ReactiveSite", "SpecialDealkyaliton");
		copy_atomTwo.setProperty("ReactiveSite", "SpecialDealkyaliton");
		copy_atomThree.setProperty("ReactiveSite", "SpecialDealkyaliton");
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(copy_atomOne);
		fragmentStructure.addAtom(copy_atomTwo);
		fragmentStructure.addAtom(copy_atomThree);
		fragmentStructure.addBond(copy.getBond(copy_atomOne, copy_atomTwo));
		fragmentStructure.addBond(copy.getBond(copy_atomOne, copy_atomThree));
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(copy_atomOne, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(copy_atomTwo, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(copy_atomThree, copy));
	
	
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		String OCO_match_SMIRKS = "[H][#6;A:1]([#8;A:3])[#8;A:2]";
		boolean matches = MergeFragments.matchSMIRKS(OCO_match_SMIRKS, fragmentStructure);
		if(!matches) return null;
		else{
			String OCO_reaction_SMIRKS = "[H][#6:5]-1-[#8:1]-[#6:2]~[#6:3]-[#8:4]-1>>[H][#8:1]-[#6:2]~[#6:3]-[#8:4][H].[#6:5]=O";
			SMIRKSReaction reaction = smrkMan.parse(OCO_reaction_SMIRKS);
			IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(copy, null, reaction);
			//if(changedPart.getAtomContainerCount() > 2) throw new Exception("There are more than two dealkylation metabolite for the same reactive site within the molecule.");	
			IAtomContainerSet resultSet = null;//ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
			for(int i = 0; i < changedPart.getAtomContainerCount(); i++){
				IAtomContainer onePart = changedPart.getAtomContainer(i);
				if((UniqueIDFunctionSet.containAtom(atom_carbon, onePart)) &&
				   (UniqueIDFunctionSet.containAtom(atom_two, onePart)) &&
				   (UniqueIDFunctionSet.containAtom(atom_three, onePart))){
					resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(i));
				}
			}
			//IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
			IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			Double score = atom_carbon.getProperty("Score");
			if(score < (Double) atom_two.getProperty("Score")) score = atom_two.getProperty("Score");
			if(score < (Double) atom_three.getProperty("Score")) score = atom_three.getProperty("Score");
			for(int i = 0; i < resultSet.getAtomContainerCount(); i++){
				//System.out.println("SMILES of the changed part: " + sg.create(resultSet.getAtomContainer(i)));	
				IAtomContainer metabolite = resultSet.getAtomContainer(i);		
				Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
				properties.put("ReactionType", "SepcialDealkyaltion");
				//Use the score of the typeOne BoM when both TypeOne and TypeTwo BoMs are included in one reaction
				properties.put("Score", score);
				properties.put("BoMs", (atom_carbon.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(atom_carbon) + ";" + 
										atom_two.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(atom_two) + ";" +
										atom_three.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(atom_three)));
				metabolite.addProperties(properties);
				AtomContainerManipulator.suppressHydrogens(metabolite);
				metabolites.addAtomContainer(metabolite);
			}				
			return metabolites;
		}
	}
	/**
	 * This function will generate the metabolites for the ring rearrangement reaction.
	 * The inputs are:
	 * typeOne_BoMList that contains all the BoMs involved in the ring rearrangement reaction.
	 * Note that the all BoMs here are TypeOne BoMs.
	 * The compounds having rearrangement reaction in the training dataset can be found at link: 
	 * https://docs.google.com/presentation/d/1sythOW42ERDXv9BQEl-hdSFmrDK-kS9-6Hb1F2AyOYE/edit#slide=id.g84e799a792_0_46
	 * Special Case: simvastatin. Should create one independent function that generates metabolite for this compound. 
	 * It's not clear that why the hydroxylation and the rearrangement reaction occurs at the same time. 
	 * @param typeOne_BoMList
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateRingRearrangementMetabolite(ArrayList<ArrayList<IAtom>> typeOne_BoMList, IAtomContainer oneMole) throws Exception{
		IAtomContainer copy = oneMole.clone();
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		ArrayList<Integer> checkExist = new ArrayList<>();
		Double score = 0.0;
		for(int i = 0; i < typeOne_BoMList.size(); i++){
			ArrayList<IAtom> one_TypeOne_BoM = typeOne_BoMList.get(i);
			for(int j = 0; j < one_TypeOne_BoM.size(); j++){
				IAtom oneAtom_temp = one_TypeOne_BoM.get(j);
				IAtom newAtom = new Atom(oneAtom_temp.getSymbol());
				newAtom.setProperties(oneAtom_temp.getProperties());
				newAtom.setProperty("ReactiveSite", "Rearrangement");
				Integer uid = newAtom.getProperty("UniqueID");
				if(!checkExist.contains(uid)){
					Double oneScore = newAtom.getProperty("Score");
					if(score < oneScore) score = oneScore;
					fragmentStructure.addAtom(newAtom);
					checkExist.add(uid);
				}				
			}
		}
		for(int i = 0; i < fragmentStructure.getAtomCount(); i++){
			IAtom atom_one = fragmentStructure.getAtom(i);
			Integer uid_one = atom_one.getProperty("UniqueID");
			for(int j = i+1; j < fragmentStructure.getAtomCount(); j++){
				IAtom atom_two = fragmentStructure.getAtom(j);
				Integer uid_two = atom_two.getProperty("UniqueID");
				IBond bond_in_origMole = copy.getBond(UniqueIDFunctionSet.getAtomByUniqueID(uid_one, copy), UniqueIDFunctionSet.getAtomByUniqueID(uid_two, copy));
				//If the bond exists in the original structure (copy here), then we should add them into the fragmentStructure
				if(bond_in_origMole!=null){
					IAtom[] bondAtoms = {atom_one, atom_two};
					IBond newBond = new Bond(bondAtoms, bond_in_origMole.getOrder());
					fragmentStructure.addBond(newBond);
				}
			}					
		}
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		for(int i = 0; i < fragmentStructure.getAtomCount(); i++){
			if(!fragmentStructure.getAtom(i).getSymbol().equals("H")){
				fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(fragmentStructure.getAtom(i), copy));
				IAtom theAtom = fragmentStructure.getAtom(i);
				Double theBondOrderSum = fragmentStructure.getBondOrderSum(theAtom);
				int theValence = theAtom.getValency();
				int missingHydrogen = (int) (theValence - theBondOrderSum);
				theAtom.setImplicitHydrogenCount(missingHydrogen);
			}			
		}
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String ringRearrangement_SMARTS;
		ArrayList<String> ringRearrangement_SMARTS_list = new ArrayList<>();
		String ringRearrangement_paraBenzene_SMARTS = "[#6;A:4]1[#6;A:5]=[#6;A:6][#7;A:1][#6;A:2]=[#6;A:3]1>>[#6;A:4]1=[#6;A:5][#6;A:6]=[#7;A:1][#6;A:2]=[#6;A:3]1";
				//"[H][*:1][#6;A:2]1=[#6;A:3][#6;A:4]=[#6;A:5]([*:8][H])[#6;A:6]=[#6;A:7]1>>[*:1]=[#6;A:2]1[#6;A:7]=[#6;A:6][#6;A:5](=[*:8])[#6;A:4]=[#6;A:3]1";
		String ringRearrangement_NitrogenBenzene_SMARTS = "[H][#6;A:4]1[#6;A:5]=[#6;A:6][#7;A:1]([H])[#6;A:2]=[#6;A:3]1>>[#6;A:4]1=[#6;A:5][#6;A:6]=[#7;A:1][#6;A:2]=[#6;A:3]1";
				//String ringRearrangement_DHP_SMARTS = "[#6;A:1][#7;A:2]1[#6;A:3][#6;A:4]=[#6;A:5][#6;A:6]1>>[#6;A:1][#7;A:2]1[#6;A:3]=[#6;A:4][#6;A:5]=[#6;A:6]1";
		String aromatic_to_nonAromatic_non_attachment_1_SMARTS = "[H][#6;A:4]1=[#6;A:5][#6;A:6]=[#7;A:1][#6;A:2]=[#6;A:3]1>>O=[#6;A:4]1[#6;A:3]=[#6;A:2][#7;A:1][#6;A:6]=[#6;A:5]1";//"[#6;A:4]1=[#6;A:5][#6;A:6]=[#7;A:1][#6;A:2]=[#6;A:3]1>>[#6;A:4]1[#6;A:3]=[#6;A:2][#7;A:1][#6;A:6]=[#6;A:5]1";
		String aromatic_to_nonAromatic_attachement_1_SMARTS = "[H][#7:8]-[#6:1]-1=[#6:2]-[#6:3]=[#6:4](-[#8:7][H])-[#6:5]=[#6:6]-1>>[#7:8]=[#6:1]-1-[#6:6]=[#6:5]-[#6:4](=[O:7])-[#6:3]=[#6:2]-1";
		String nonAromatic_to_aromatic_attachment_1_SMARTS = "[H][#8]-[#6:2]-1-[#6:3]-[#6:4]-[#7:5]-[#6:6]-[#6:1]-1>>[#6:2]-1=[#6:1]-[#6:6]=[#7+:5]-[#6:4]=[#6:3]-1";
		ringRearrangement_SMARTS_list.add(ringRearrangement_paraBenzene_SMARTS);
		ringRearrangement_SMARTS_list.add(ringRearrangement_NitrogenBenzene_SMARTS);
		ringRearrangement_SMARTS_list.add(aromatic_to_nonAromatic_attachement_1_SMARTS);
		ringRearrangement_SMARTS_list.add(aromatic_to_nonAromatic_non_attachment_1_SMARTS);
		ringRearrangement_SMARTS_list.add(nonAromatic_to_aromatic_attachment_1_SMARTS);
		//ringRearrangement_SMARTS_list.add(ringRearrangement_DHP_SMARTS);
		ringRearrangement_SMARTS = ringRearrangement_paraBenzene_SMARTS;
		SMIRKSReaction reaction = smrkMan.parse(ringRearrangement_SMARTS);		
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		int counter = 1;
		/**
		 * Keep running till we found:
		 * 1. metabolites generated OR
		 * 2. all the candidate ringRearrangement SMARTS strings are explored
		 */
		while(changedPart == null && counter < ringRearrangement_SMARTS_list.size()){
			ringRearrangement_SMARTS = ringRearrangement_SMARTS_list.get(counter);
			reaction = smrkMan.parse(ringRearrangement_SMARTS);		
			changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
			counter++;
		}
		if(changedPart.getAtomContainerCount() > 2) throw new Exception("There are more than two dealkylation metabolite for the same reactive site within the molecule.");		
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		//if(resultSet.getAtomContainerCount() > 2) throw new Exception("There are more than two dealkylation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<IAtom>  removedAtoms = MergeFragments.getRemovedAtomList(fragmentStructure, resultSet.getAtomContainer(0));
		IAtomContainer metabolite = MergeFragments.mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0), removedAtoms);
				//mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0));	
		//Set reactionType and Score for the corresponding metabolite
		Map<Object,Object> properties = new HashMap<Object,Object>();
		properties.put("ReactionType", "Rearrangement");
		properties.put("Score", score);
		metabolite.addProperties(properties);
		metabolites.addAtomContainer(metabolite);
		return metabolites;		
	}
	
	/**
	 * An Dehydrogenation reaction requires an <eta,eta> bond to be modified.
	 * Currently, the instances we have seen are C-C => C=C reactions
	 * Hence both the atoms connected by the reactive Bond should be used as input
	 * @param leftAtom
	 * @param rightAtom
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateDehydrogenationMetabolite(IAtom leftAtom, IAtom rightAtom, IAtomContainer oneMole, int counter) throws Exception{		
		IAtomContainer copy = oneMole.clone();		
		String oringinalSmile = sg.create(oneMole);
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), copy);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), copy);
		IBond targetBond = copy.getBond(leftCopy, rightCopy);
		if(targetBond.getOrder() == IBond.Order.SINGLE){
			IAtomContainerSet resonanceStructureSet = generateResonanceStructure(copy);
			for(int i = 0; i < resonanceStructureSet.getAtomContainerCount(); i++){
				IAtomContainer oneResonaceStruct = resonanceStructureSet.getAtomContainer(i);
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneResonaceStruct);
				leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), oneResonaceStruct);
				rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), oneResonaceStruct);
				targetBond = oneResonaceStruct.getBond(leftCopy, rightCopy);
				if(targetBond.getOrder() == IBond.Order.DOUBLE){
					String smile = sg.create(oneResonaceStruct);
					if(!oringinalSmile.contains("+") && !oringinalSmile.contains("-") && 
						!smile.contains("+") && !smile.contains("-")){
						copy = oneResonaceStruct;
						//System.out.println("selected resonance struct: " + sg.create(oneResonaceStruct));
						break;
					}
					else{
						copy = oneResonaceStruct;
						//System.out.println("selected resonance struct: " + sg.create(oneResonaceStruct));
						break;
					}
				}
			}
		}
		
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		//IBond targetBond = copy.getBond(leftCopy, rightCopy);
		//System.out.println("Indices of the BoM atoms: " + copy.indexOf(leftCopy) + "," + copy.indexOf(rightCopy));
		//System.out.println("UniqueIDs of the BoM atoms: " + UniqueIDFunctionSet.getUniqID(leftAtom) + "," + UniqueIDFunctionSet.getUniqID(rightAtom));
		leftCopy.setProperty("ReactiveSite", "Dehydrogenation");
		rightCopy.setProperty("ReactiveSite", "Dehydrogenation");
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(leftCopy);		
		fragmentStructure.addAtom(rightCopy);
		IAtom[] bondAtoms = {leftCopy, rightCopy};
		targetBond = new Bond(bondAtoms, IBond.Order.SINGLE);
		//IBond tt = copy.getBond(leftCopy, rightCopy);
		//tt = newBond;
		fragmentStructure.addBond(targetBond);//(copy.getBond(leftCopy, rightCopy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(leftCopy, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(rightCopy, copy));
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String dehydrogenation_SMIRKS = "[#6:1]-[#6:2]>>[#6:2]=[#6:1]";
	
		SMIRKSReaction reaction = smrkMan.parse(dehydrogenation_SMIRKS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart.getAtomContainerCount() > 1) throw new Exception("There are more than one epOxidation metabolite for the same reactive site within the molecule.");
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		if(resultSet.getAtomContainerCount() > 1) throw new Exception("There are more than one dehydrogenation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		/**
		 * We remove the bond <leftCopy, rightCopy> within the molecule 
		 */
		copy.removeBond(leftCopy, rightCopy);
		//System.out.println("SMILES of the changed part: " + sg.create(resultSet.getAtomContainer(0)));	
		ArrayList<IAtom>  removedAtoms = MergeFragments.getRemovedAtomList(fragmentStructure, resultSet.getAtomContainer(0));
		IAtomContainer metabolite = MergeFragments.mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0), removedAtoms);//mergeFragments_EpOxidation(copy, resultSet.getAtomContainer(0));	
		Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
		properties.put("ReactionType", "Dehydrogenation");
		properties.put("Score", leftAtom.getProperty("Score"));
		metabolite.addProperties(properties);
		metabolites.addAtomContainer(metabolite);			
		return metabolites;
		
	}	
	
	/**
	 * This function will generate the metabolites for the special reaction one: S-C=C -> S-C(=O)-C
	 * @param typeOne_BoMList
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateSpecialReactionOneMetabolite(ArrayList<ArrayList<IAtom>> typeOne_BoMList, IAtomContainer oneMole) throws Exception{
		IAtomContainer copy = oneMole.clone();
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		
		String smirks = "[#6]=,:1[#6]=,:[#6][#16][#6]=,:1";
		ArrayList<ArrayList<IAtom>> matchedList = MergeFragments.findMatchedList(smirks, copy);
		
		ArrayList<Integer> checkExist = new ArrayList<>();
		Double score = 0.0;
		for(int i = 0; i < typeOne_BoMList.size(); i++){
			ArrayList<IAtom> one_TypeOne_BoM = typeOne_BoMList.get(i);
			for(int j = 0; j < one_TypeOne_BoM.size(); j++){
				IAtom oneAtom_temp = one_TypeOne_BoM.get(j);
				IAtom newAtom = new Atom(oneAtom_temp.getSymbol());
				newAtom.setProperties(oneAtom_temp.getProperties());
				newAtom.setProperty("ReactiveSite", "SpecialReactionOne");
				Integer uid = newAtom.getProperty("UniqueID");
				if(!checkExist.contains(uid)){
					Double oneScore = newAtom.getProperty("Score");
					if(score < oneScore) score = oneScore;
					fragmentStructure.addAtom(newAtom);
					checkExist.add(uid);
				}				
			}
		}
		//ArrayList<IAtom> targetMatched = new ArrayList<>();
		for(int i = 0; i < matchedList.size(); i++){
			boolean  istarget = true;
			List<IAtom> oneMatched = matchedList.get(i);
			ArrayList<Integer> matchedUID = new ArrayList<>();
			for(int j = 0; j < oneMatched.size(); j++){
				Integer uid = UniqueIDFunctionSet.getUniqID(oneMatched.get(j));
				if(!matchedUID.contains(uid)){
					matchedUID.add(uid);
				}
			}
			for(int j = 0; j < typeOne_BoMList.size(); j++){
				ArrayList<IAtom> oneBoM = typeOne_BoMList.get(j);
				//If one BoM is not within the matched fragment, then it's not the target candidate, continue to the next one
				if(!matchedUID.contains(UniqueIDFunctionSet.getUniqID(oneBoM.get(0))) || !matchedUID.contains(UniqueIDFunctionSet.getUniqID(oneBoM.get(1)))){
					matchedUID = new ArrayList<>();
					istarget = false;
					break;
				}
			}
			if(istarget){
				for(int j = 0; j < oneMatched.size(); j++){
					IAtom oneAtom_temp = oneMatched.get(j);
					IAtom newAtom = new Atom(oneAtom_temp.getSymbol());
					newAtom.setProperties(oneAtom_temp.getProperties());
					Integer uid = newAtom.getProperty("UniqueID");
					if(!checkExist.contains(uid)){
						fragmentStructure.addAtom(newAtom);
						checkExist.add(uid);
					}
				}
			}			
		}
		for(int i = 0; i < fragmentStructure.getAtomCount(); i++){
			IAtom atom_one = fragmentStructure.getAtom(i);
			Integer uid_one = atom_one.getProperty("UniqueID");
			for(int j = i+1; j < fragmentStructure.getAtomCount(); j++){
				IAtom atom_two = fragmentStructure.getAtom(j);
				Integer uid_two = atom_two.getProperty("UniqueID");
				IBond bond_in_origMole = copy.getBond(UniqueIDFunctionSet.getAtomByUniqueID(uid_one, copy), UniqueIDFunctionSet.getAtomByUniqueID(uid_two, copy));
				//If the bond exists in the original structure (copy here), then we should add them into the fragmentStructure
				if(bond_in_origMole!=null){
					IAtom[] bondAtoms = {atom_one, atom_two};
					IBond newBond = new Bond(bondAtoms, bond_in_origMole.getOrder());
					fragmentStructure.addBond(newBond);
				}
			}					
		}
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		for(int i = 0; i < fragmentStructure.getAtomCount(); i++){
			if(!fragmentStructure.getAtom(i).getSymbol().equals("H")){
				fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(fragmentStructure.getAtom(i), copy));
				IAtom theAtom = fragmentStructure.getAtom(i);
				Double theBondOrderSum = fragmentStructure.getBondOrderSum(theAtom);
				int theValence = theAtom.getValency();
				int missingHydrogen = (int) (theValence - theBondOrderSum);
				theAtom.setImplicitHydrogenCount(missingHydrogen);
			}			
		}
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String specialReactionOne_SMIRKS = "[H][#6:5]1=,:[#6:4][#6:3]=,:[#6:2][#16:1]1>>O=[#6:5]-1-[#16:1]-[#6:2]-[#6:3]=[#6:4]-1";
		SMIRKSReaction reaction = smrkMan.parse(specialReactionOne_SMIRKS);		
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart.getAtomContainerCount() > 2) throw new Exception("There are more than two dealkylation metabolite for the same reactive site within the molecule.");		
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		//if(resultSet.getAtomContainerCount() > 2) throw new Exception("There are more than two dealkylation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<IAtom>  removedAtoms = MergeFragments.getRemovedAtomList(fragmentStructure, resultSet.getAtomContainer(0));
		IAtomContainer metabolite = MergeFragments.mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0), removedAtoms);
				//mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0));	
		//Set reactionType and Score for the corresponding metabolite
		Map<Object,Object> properties = new HashMap<Object,Object>();
		properties.put("ReactionType", "SpecialReactionOne");
		properties.put("Score", score);
		metabolite.addProperties(properties);
		metabolites.addAtomContainer(metabolite);
		return metabolites;
		
	}
	
	/**
	 * This function will generate the metabolites for the special reaction one: N=C-N -> N-C(=O)-N 
	 * @param typeOne_BoMList
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainerSet generateSpecialReactionTwoMetabolite(ArrayList<ArrayList<IAtom>> typeOne_BoMList, IAtomContainer oneMole) throws Exception{
		IAtomContainer copy = oneMole.clone();
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		
		String smirks = "[#6]1~[#6][#6,#7]=,:[#6][#7]1";
		ArrayList<ArrayList<IAtom>> matchedList = MergeFragments.findMatchedList(smirks, copy);
		
		ArrayList<Integer> checkExist = new ArrayList<>();
		Double score = 0.0;
		for(int i = 0; i < typeOne_BoMList.size(); i++){
			ArrayList<IAtom> one_TypeOne_BoM = typeOne_BoMList.get(i);
			for(int j = 0; j < one_TypeOne_BoM.size(); j++){
				IAtom oneAtom_temp = one_TypeOne_BoM.get(j);
				IAtom newAtom = new Atom(oneAtom_temp.getSymbol());
				newAtom.setProperties(oneAtom_temp.getProperties());
				newAtom.setProperty("ReactiveSite", "SpecialReactionTwo");
				Integer uid = newAtom.getProperty("UniqueID");
				if(!checkExist.contains(uid)){
					Double oneScore = newAtom.getProperty("Score");
					if(score < oneScore) score = oneScore;
					fragmentStructure.addAtom(newAtom);
					checkExist.add(uid);
				}				
			}
		}
		//ArrayList<IAtom> targetMatched = new ArrayList<>();
		for(int i = 0; i < matchedList.size(); i++){
			boolean  istarget = true;
			List<IAtom> oneMatched = matchedList.get(i);
			ArrayList<Integer> matchedUID = new ArrayList<>();
			for(int j = 0; j < oneMatched.size(); j++){
				Integer uid = UniqueIDFunctionSet.getUniqID(oneMatched.get(j));
				if(!matchedUID.contains(uid)){
					matchedUID.add(uid);
				}
			}
			for(int j = 0; j < typeOne_BoMList.size(); j++){
				ArrayList<IAtom> oneBoM = typeOne_BoMList.get(j);
				//If one BoM is not within the matched fragment, then it's not the target candidate, continue to the next one
				if(!matchedUID.contains(UniqueIDFunctionSet.getUniqID(oneBoM.get(0))) || !matchedUID.contains(UniqueIDFunctionSet.getUniqID(oneBoM.get(1)))){
					matchedUID = new ArrayList<>();
					istarget = false;
					break;
				}
			}
			if(istarget){
				for(int j = 0; j < oneMatched.size(); j++){
					IAtom oneAtom_temp = oneMatched.get(j);
					IAtom newAtom = new Atom(oneAtom_temp.getSymbol());
					newAtom.setProperties(oneAtom_temp.getProperties());
					Integer uid = newAtom.getProperty("UniqueID");
					if(!checkExist.contains(uid)){
						fragmentStructure.addAtom(newAtom);
						checkExist.add(uid);
					}
				}
			}			
		}
		for(int i = 0; i < fragmentStructure.getAtomCount(); i++){
			IAtom atom_one = fragmentStructure.getAtom(i);
			Integer uid_one = atom_one.getProperty("UniqueID");
			for(int j = i+1; j < fragmentStructure.getAtomCount(); j++){
				IAtom atom_two = fragmentStructure.getAtom(j);
				Integer uid_two = atom_two.getProperty("UniqueID");
				IBond bond_in_origMole = copy.getBond(UniqueIDFunctionSet.getAtomByUniqueID(uid_one, copy), UniqueIDFunctionSet.getAtomByUniqueID(uid_two, copy));
				//If the bond exists in the original structure (copy here), then we should add them into the fragmentStructure
				if(bond_in_origMole!=null){
					IAtom[] bondAtoms = {atom_one, atom_two};
					IBond newBond = new Bond(bondAtoms, bond_in_origMole.getOrder());
					fragmentStructure.addBond(newBond);
				}
			}					
		}
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		for(int i = 0; i < fragmentStructure.getAtomCount(); i++){
			if(!fragmentStructure.getAtom(i).getSymbol().equals("H")){
				fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(fragmentStructure.getAtom(i), copy));
				IAtom theAtom = fragmentStructure.getAtom(i);
				Double theBondOrderSum = fragmentStructure.getBondOrderSum(theAtom);
				int theValence = theAtom.getValency();
				int missingHydrogen = (int) (theValence - theBondOrderSum);
				theAtom.setImplicitHydrogenCount(missingHydrogen);
			}			
		}
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String specialReactionTwo_SMIRKS = "[#6:5]1~[#6:4][#6,#7:3]=,:[#6:2][#7:1]1>>O=[#6:2]-1-[#7:1]-[#6:5]~[#6:4]-[#6,#7:3]-1";
		SMIRKSReaction reaction = smrkMan.parse(specialReactionTwo_SMIRKS);		
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart.getAtomContainerCount() > 2) throw new Exception("There are more than two dealkylation metabolite for the same reactive site within the molecule.");		
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		//if(resultSet.getAtomContainerCount() > 2) throw new Exception("There are more than two dealkylation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<IAtom>  removedAtoms = MergeFragments.getRemovedAtomList(fragmentStructure, resultSet.getAtomContainer(0));
		IAtomContainer metabolite = MergeFragments.mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0), removedAtoms);
				//mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0));	
		//Set reactionType and Score for the corresponding metabolite
		Map<Object,Object> properties = new HashMap<Object,Object>();
		properties.put("ReactionType", "SpecialReactionOne");
		properties.put("Score", score);
		metabolite.addProperties(properties);
		metabolites.addAtomContainer(metabolite);
		return metabolites;
		
	}
	/**
	 * This function will generate the Desulfurization metabolites
	 * @param specialTypeOneReactionsBoMs
	 * @param oneMole
	 * @return
	 */
	public static IAtomContainerSet generateDesulfurizationMetabolite(ArrayList<ArrayList<IAtom>> specialTypeOneReactionsBoMs, IAtomContainer oneMole) throws Exception{
		IAtomContainer copy = oneMole.clone();		
		IAtom leftAtom = specialTypeOneReactionsBoMs.get(0).get(0);
		IAtom rightAtom = specialTypeOneReactionsBoMs.get(0).get(1);
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), copy);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), copy);
		IBond targetBond = copy.getBond(leftCopy, rightCopy);			
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		//IBond targetBond = copy.getBond(leftCopy, rightCopy);
		//System.out.println("Indices of the BoM atoms: " + copy.indexOf(leftCopy) + "," + copy.indexOf(rightCopy));
		//System.out.println("UniqueIDs of the BoM atoms: " + UniqueIDFunctionSet.getUniqID(leftAtom) + "," + UniqueIDFunctionSet.getUniqID(rightAtom));
		leftCopy.setProperty("ReactiveSite", "Desulfurization");
		rightCopy.setProperty("ReactiveSite", "Desulfurization");
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(leftCopy);		
		fragmentStructure.addAtom(rightCopy);
		IAtom[] bondAtoms = {leftCopy, rightCopy};
		targetBond = new Bond(bondAtoms, IBond.Order.DOUBLE);
		//IBond tt = copy.getBond(leftCopy, rightCopy);
		//tt = newBond;
		fragmentStructure.addBond(targetBond);//(copy.getBond(leftCopy, rightCopy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(leftCopy, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(rightCopy, copy));
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String deSulfurization_SMIRKS = "S=[#6,#15:1]>>O=[#6,#15:1]";
	
		SMIRKSReaction reaction = smrkMan.parse(deSulfurization_SMIRKS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart.getAtomContainerCount() > 1) throw new Exception("There are more than one epOxidation metabolite for the same reactive site within the molecule.");
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		if(resultSet.getAtomContainerCount() > 1) throw new Exception("There are more than one dehydrogenation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		/**
		 * We remove the bond <leftCopy, rightCopy> within the molecule 
		 */
		copy.removeBond(leftCopy, rightCopy);
		//System.out.println("SMILES of the changed part: " + sg.create(resultSet.getAtomContainer(0)));	
		ArrayList<IAtom>  removedAtoms = MergeFragments.getRemovedAtomList(fragmentStructure, resultSet.getAtomContainer(0));
		IAtomContainer metabolite = MergeFragments.mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0), removedAtoms);//mergeFragments_EpOxidation(copy, resultSet.getAtomContainer(0));	
		Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
		properties.put("ReactionType", "Desulfurization");
		properties.put("Score", leftAtom.getProperty("Score"));
		metabolite.addProperties(properties);
		metabolites.addAtomContainer(metabolite);			
		return metabolites;
	}
	/**
	 * This function will generate hydrolysis metabolites.
	 * Note that O=C-O-[C,N,S] => O=C-OH , HO-[C,N,S] hydrolysis reactions is covered by dealkylation reaction. So they won't be considered here.
	 * @param specialTypeOneReactionsBoMs
	 * @param oneMole
	 * @return
	 */
	public static IAtomContainerSet generateHydrolysisMetabolite(ArrayList<ArrayList<IAtom>> specialTypeOneReactionsBoMs, IAtomContainer oneMole) throws Exception {
		IAtomContainer copy = oneMole.clone();		
		IAtom leftAtom = specialTypeOneReactionsBoMs.get(0).get(0);
		IAtom rightAtom = specialTypeOneReactionsBoMs.get(0).get(1);
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), copy);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), copy);
		IBond targetBond = copy.getBond(leftCopy, rightCopy);			
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		//IBond targetBond = copy.getBond(leftCopy, rightCopy);
		//System.out.println("Indices of the BoM atoms: " + copy.indexOf(leftCopy) + "," + copy.indexOf(rightCopy));
		//System.out.println("UniqueIDs of the BoM atoms: " + UniqueIDFunctionSet.getUniqID(leftAtom) + "," + UniqueIDFunctionSet.getUniqID(rightAtom));
		leftCopy.setProperty("ReactiveSite", "Hydrolysis");
		rightCopy.setProperty("ReactiveSite", "Hydrolysis");
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(leftCopy);		
		fragmentStructure.addAtom(rightCopy);
		IAtom[] bondAtoms = {leftCopy, rightCopy};
		targetBond = new Bond(bondAtoms, targetBond.getOrder());
		//IBond tt = copy.getBond(leftCopy, rightCopy);
		//tt = newBond;
		fragmentStructure.addBond(targetBond);//(copy.getBond(leftCopy, rightCopy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(leftCopy, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(rightCopy, copy));
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String deSulfurization_SMIRKS = "[#8:2]-[#16,#7,#15:1]>>[#8:2][H].[#8]-[#16,#7,#15:1]";
	
		SMIRKSReaction reaction = smrkMan.parse(deSulfurization_SMIRKS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart.getAtomContainerCount() > 1) throw new Exception("There are more than one epOxidation metabolite for the same reactive site within the molecule.");
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		//if(resultSet.getAtomContainerCount() > 1) throw new Exception("There are more than one dehydrogenation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		/**
		 * We remove the bond <leftCopy, rightCopy> within the molecule 
		 */
		copy.removeBond(leftCopy, rightCopy);
		IAtomContainerSet subFragments = ConnectivityChecker.partitionIntoMolecules(copy);
		//The part contains "O" will be at index = 0. There should be only two subFragment for Dealkylation reaction.
		ArrayList<IAtomContainer> reorderedFragments = reorderFragments(subFragments);
		ArrayList<IAtomContainer> reorderedResultSet = reorderFragments(resultSet);
		for(int i = 0; i < reorderedResultSet.size(); i++){
			//System.out.println("SMILES of the changed part: " + sg.create(reorderedResultSet.get(i)));	
			IAtomContainer metabolite = mergeFragments(reorderedFragments.get(i), reorderedResultSet.get(i));			
			Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
			properties.put("ReactionType", "Hydrolysis");
			//Use the score of the typeOne BoM when both TypeOne and TypeTwo BoMs are included in one reaction
			if(leftCopy.getSymbol().equalsIgnoreCase("C")) properties.put("Score", rightCopy.getProperty("Score"));
			else properties.put("Score", leftAtom.getProperty("Score"));
			metabolite.addProperties(properties);
			metabolites.addAtomContainer(metabolite);
		}
		return metabolites;
	}
	/**
	 * This function will generate the DehalogenMetabolite
	 * @param specialTypeOneReactionsBoMs
	 * @param oneMole
	 * @return
	 */
	public static IAtomContainerSet generateDehalogenMetabolite(ArrayList<ArrayList<IAtom>> specialTypeOneReactionsBoMs, IAtomContainer oneMole) throws Exception{
		IAtomContainer copy = oneMole.clone();		
		IAtom leftAtom = specialTypeOneReactionsBoMs.get(0).get(0);
		IAtom rightAtom = specialTypeOneReactionsBoMs.get(0).get(1);
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), copy);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), copy);
		IBond targetBond = copy.getBond(leftCopy, rightCopy);			
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		//IBond targetBond = copy.getBond(leftCopy, rightCopy);
		//System.out.println("Indices of the BoM atoms: " + copy.indexOf(leftCopy) + "," + copy.indexOf(rightCopy));
		//System.out.println("UniqueIDs of the BoM atoms: " + UniqueIDFunctionSet.getUniqID(leftAtom) + "," + UniqueIDFunctionSet.getUniqID(rightAtom));
		leftCopy.setProperty("ReactiveSite", "Dehalogenation");
		rightCopy.setProperty("ReactiveSite", "Dehalogenation");
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(leftCopy);		
		fragmentStructure.addAtom(rightCopy);
		IAtom[] bondAtoms = {leftCopy, rightCopy};
		targetBond = new Bond(bondAtoms, IBond.Order.SINGLE);
		//IBond tt = copy.getBond(leftCopy, rightCopy);
		//tt = newBond;
		fragmentStructure.addBond(targetBond);//(copy.getBond(leftCopy, rightCopy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(leftCopy, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(rightCopy, copy));
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String deSulfurization_SMIRKS = "[#6:2]-[#9,#17:1]>>[#9,#17:1][H].[#6:2]=O";
	
		SMIRKSReaction reaction = smrkMan.parse(deSulfurization_SMIRKS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		//if(changedPart.getAtomContainerCount() > 1) throw new Exception("There are more than one epOxidation metabolite for the same reactive site within the molecule.");
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		//if(resultSet.getAtomContainerCount() > 1) throw new Exception("There are more than one dehydrogenation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		/**
		 * We remove the bond <leftCopy, rightCopy> within the molecule 
		 */
		copy.removeBond(leftCopy, rightCopy);
		IAtomContainerSet subFragments = ConnectivityChecker.partitionIntoMolecules(copy);
		//The part contains "O" will be at index = 0. There should be only two subFragment for Dealkylation reaction.
		ArrayList<IAtomContainer> reorderedFragments = reorderFragments(subFragments);
		ArrayList<IAtomContainer> reorderedResultSet = reorderFragments(resultSet);
		for(int i = 0; i < reorderedResultSet.size(); i++){
			//System.out.println("SMILES of the changed part: " + sg.create(reorderedResultSet.get(i)));	
			IAtomContainer metabolite = mergeFragments(reorderedFragments.get(i), reorderedResultSet.get(i));			
			Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
			properties.put("ReactionType", "Dehalogenation");
			//Use the score of the typeOne BoM when both TypeOne and TypeTwo BoMs are included in one reaction
			if(leftCopy.getSymbol().equalsIgnoreCase("C")) properties.put("Score", rightCopy.getProperty("Score"));
			else properties.put("Score", leftAtom.getProperty("Score"));
			metabolite.addProperties(properties);
			properties.put("BoMs", (leftAtom.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(leftAtom) + ";" + rightAtom.getSymbol() + "." + UniqueIDFunctionSet.getUniqID(rightAtom)));
			metabolites.addAtomContainer(metabolite);
		}
		return metabolites;
	}
	/**
	 * This function will generate the Reduction metabolites
	 * @param specialTypeOneReactionsBoMs
	 * @param oneMole
	 * @return
	 */
	public static IAtomContainerSet generateReductionMetabolite(ArrayList<ArrayList<IAtom>> specialTypeOneReactionsBoMs, IAtomContainer oneMole) throws Exception{
		IAtomContainer copy = oneMole.clone();		
		IAtom leftAtom = specialTypeOneReactionsBoMs.get(0).get(0);
		IAtom rightAtom = specialTypeOneReactionsBoMs.get(0).get(1);
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(leftAtom), copy);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(rightAtom), copy);
		IBond targetBond = copy.getBond(leftCopy, rightCopy);			
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
		//IBond targetBond = copy.getBond(leftCopy, rightCopy);
		//System.out.println("Indices of the BoM atoms: " + copy.indexOf(leftCopy) + "," + copy.indexOf(rightCopy));
		//System.out.println("UniqueIDs of the BoM atoms: " + UniqueIDFunctionSet.getUniqID(leftAtom) + "," + UniqueIDFunctionSet.getUniqID(rightAtom));
		leftCopy.setProperty("ReactiveSite", "Reduction");
		rightCopy.setProperty("ReactiveSite", "Reduction");
		
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		fragmentStructure.addAtom(leftCopy);		
		fragmentStructure.addAtom(rightCopy);
		IAtom[] bondAtoms = {leftCopy, rightCopy};
		targetBond = new Bond(bondAtoms, IBond.Order.DOUBLE);
		//IBond tt = copy.getBond(leftCopy, rightCopy);
		//tt = newBond;
		fragmentStructure.addBond(targetBond);//(copy.getBond(leftCopy, rightCopy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(leftCopy, copy));
		fragmentStructure.add(generateFragmentWithExplicitHydrogenForReactiveAtom(rightCopy, copy));
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		
		//System.out.println("SMILES of the reacive site(s) : " + sg.create(fragmentStructure));
		
		String reduction_SMIRKS = "O=[#15,#16:1]>>[#15,#16:1]";			
		SMIRKSReaction reaction = smrkMan.parse(reduction_SMIRKS);
		IAtomContainerSet changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		if(changedPart==null){
			reduction_SMIRKS = "[#6:1]=[O:2]>>[#6:1]-[#8:2]";	
			reaction = smrkMan.parse(reduction_SMIRKS);
			changedPart = smrkMan.applyTransformationWithSingleCopyForEachPos(fragmentStructure, null, reaction);
		}
		if(changedPart.getAtomContainerCount() > 1) throw new Exception("There are more than one epOxidation metabolite for the same reactive site within the molecule.");
		IAtomContainerSet resultSet = ConnectivityChecker.partitionIntoMolecules(changedPart.getAtomContainer(0));
		if(resultSet.getAtomContainerCount() > 1) throw new Exception("There are more than one dehydrogenation metabolite for the same reactive site within the molecule after partition.");		
		IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		/**
		 * We remove the bond <leftCopy, rightCopy> within the molecule 
		 */
		copy.removeBond(leftCopy, rightCopy);
		//System.out.println("SMILES of the changed part: " + sg.create(resultSet.getAtomContainer(0)));	
		ArrayList<IAtom>  removedAtoms = MergeFragments.getRemovedAtomList(fragmentStructure, resultSet.getAtomContainer(0));
		IAtomContainer metabolite = MergeFragments.mergeFragmentsForRearrangementReaction(copy, resultSet.getAtomContainer(0), removedAtoms);//mergeFragments_EpOxidation(copy, resultSet.getAtomContainer(0));	
		Map<Object,Object> properties = new HashMap<Object,Object>();//("ReactionType", "Oxidation");
		properties.put("ReactionType", "Reduction");
		properties.put("Score", leftAtom.getProperty("Score"));
		metabolite.addProperties(properties);
		metabolites.addAtomContainer(metabolite);			
		return metabolites;
	}
	/**
	 * This function will merge the structure of the reactant/substrate and the modifed structure based on the reactive site
	 * The reactive sites (bonds) will be removed after the merge 
	 * Then return the merged the Structure as an IAtomContainer 
	 * Note that the reactive sites (atoms here) are labeled with the setting of property <"ReactiveSite","$ReactionType$">
	 * @param reactant
	 * @param modified
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainer mergeFragments(IAtomContainer reactant, IAtomContainer modified) throws Exception{
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		ArrayList<IAtom> atomListConnectedToReactiveSiteInReactant = getAtomListConnectedToReactiveSiteInFragment(reactant);
		ArrayList<IAtom> atomListConnectedToReactiveSiteInModified = getAtomListConnectedToReactiveSiteInFragment(modified);
		IAtom reactiveSite = null;
		IAtom changedSite = null;
		Atom resultAtom = null;
		for(int i = 0; i < reactant.getAtomCount(); i++){
			IAtom atomOne = reactant.getAtom(i);
			if(atomOne.getProperty("ReactiveSite")!=null) reactiveSite = atomOne;
			for(int j = i+1; j < reactant.getAtomCount(); j++){
				IAtom anotherAtom = reactant.getAtom(j);
				if(!fragmentStructure.contains(atomOne)) fragmentStructure.addAtom(atomOne);
				if(!fragmentStructure.contains(anotherAtom)) fragmentStructure.addAtom(anotherAtom);
				if(reactant.getBond(atomOne, anotherAtom)!=null) fragmentStructure.addBond(reactant.getBond(atomOne, anotherAtom));
			}
		}
				
		for(int i = 0; i < modified.getAtomCount(); i++){
			IAtom atomOne = modified.getAtom(i);
			if(atomOne.getProperty("ReactiveSite")!=null) changedSite = atomOne;
			for(int j = i+1; j < modified.getAtomCount(); j++){
				IAtom anotherAtom = modified.getAtom(j);
				if(!fragmentStructure.contains(atomOne)) fragmentStructure.addAtom(atomOne);
				if(!fragmentStructure.contains(anotherAtom)) fragmentStructure.addAtom(anotherAtom);
				if(modified.getBond(atomOne, anotherAtom)!=null) fragmentStructure.addBond(modified.getBond(atomOne, anotherAtom));
			}
		}
		resultAtom = new Atom(changedSite.getSymbol());
		resultAtom.setFormalCharge(changedSite.getFormalCharge());
		
		fragmentStructure.addAtom(resultAtom);
		for(int i = 0; i < atomListConnectedToReactiveSiteInReactant.size(); i++){
			IBond oneBond = reactant.getBond(atomListConnectedToReactiveSiteInReactant.get(i),reactiveSite);
			if(atomListConnectedToReactiveSiteInReactant.get(i).getSymbol().equals("H")){
				fragmentStructure.removeAtom(atomListConnectedToReactiveSiteInReactant.get(i));
				continue;
			}
			if(oneBond == null) continue;
			IAtom[] bondAtoms = {atomListConnectedToReactiveSiteInReactant.get(i), resultAtom};
			IBond newBond = new Bond(bondAtoms, oneBond.getOrder());//changedSite
			fragmentStructure.addBond(newBond);
		}
		
		
		for(int i = 0; i < atomListConnectedToReactiveSiteInModified.size(); i++){
			IBond oneBond = modified.getBond(atomListConnectedToReactiveSiteInModified.get(i),changedSite);
			if(atomListConnectedToReactiveSiteInModified.get(i).getSymbol().equals("H")){
				fragmentStructure.removeAtom(atomListConnectedToReactiveSiteInModified.get(i));
				continue;
			}
			if(oneBond == null) continue;
			IAtom[] bondAtoms = {atomListConnectedToReactiveSiteInModified.get(i), resultAtom};
			IBond newBond = new Bond(bondAtoms, oneBond.getOrder());//changedSite
			fragmentStructure.addBond(newBond);
		}
		fragmentStructure.removeAtom(reactiveSite);
		fragmentStructure.removeAtom(changedSite);
		
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
		adder.addImplicitHydrogens(fragmentStructure);
		AtomContainerManipulator.suppressHydrogens(fragmentStructure);
		//System.out.println(sg.create(fragmentStructure));
		return fragmentStructure;
	}
	/**
	 * This function is used to merge the modified structure into the original structure of the compound to generate the ring-rearrangement metabolite
	 * Note that in this case no non-Hydrogen atom was removed. The bonds are rearranged (electrons).
	 * @param reactant
	 * @param modified
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainer mergeFragmentsForRearrangementReaction(IAtomContainer reactant, IAtomContainer modified) throws Exception{
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		AtomContainerManipulator.suppressHydrogens(reactant);
		AtomContainerManipulator.suppressHydrogens(modified);
		//Check the presence of atoms using its uniqueID. Then create new IAtom and add it into the fragmentStructure
		//Add all atoms and bonds within the reactant into the fragmentStructure. We create atoms and bonds manually instead of using the clone() method.
		for(int i = 0; i < reactant.getAtomCount(); i++){
			IAtom atom_one = new Atom(reactant.getAtom(i).getSymbol());
			if(reactant.getAtom(i).getProperty("ReactiveSite")==null){
				atom_one.setFormalCharge(reactant.getAtom(i).getFormalCharge());
			}
			atom_one.setProperties(reactant.getAtom(i).getProperties());
			fragmentStructure.addAtom(atom_one);
		}
		for(int i = 0; i < fragmentStructure.getAtomCount(); i++){
			IAtom atom_one = fragmentStructure.getAtom(i);
			IAtom atom_one_inReactant = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_one), reactant);
			for(int j = i+1; j < fragmentStructure.getAtomCount(); j++){
				IAtom atom_two = fragmentStructure.getAtom(j);
				IAtom atom_two_inReactant = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_two), reactant);
				//IBond oneBond;
				//If one of the atom is newly added and was not in the original structure, then skip it
				if(atom_one_inReactant == null || atom_two_inReactant == null){
					 continue;//oneBond = fragmentStructure.getBond(atom_one_inReactant, atom_two_inReactant);
				}
				//If both atoms were in the reactant, then add them
				else {
					IBond oneBond = reactant.getBond(atom_one_inReactant, atom_two_inReactant);				
					if(oneBond !=null){
						IAtom[] atoms_for_bond = {atom_one, atom_two};
						IBond newBond = new Bond(atoms_for_bond, oneBond.getOrder());
						fragmentStructure.addBond(newBond);
					}
				}
			}
		}
		for(int i = 0; i < modified.getAtomCount(); i++){
			IAtom atom_one = modified.getAtom(i);
			IAtom atom_one_inFragment = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_one), fragmentStructure);
			atom_one_inFragment.setFormalCharge(atom_one.getFormalCharge());
			for(int j = i+1; j < modified.getAtomCount(); j++){
				IAtom atom_two = modified.getAtom(j);
				IAtom atom_two_inFragment = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_two), fragmentStructure);
				IBond oneBond_inModified = modified.getBond(atom_one, atom_two);
				IBond oneBond_inFragment = fragmentStructure.getBond(atom_one_inFragment, atom_two_inFragment);
				if(oneBond_inModified !=null){
					if(oneBond_inFragment == null){
						IAtom[] atoms_for_bond = {atom_one, atom_two};
						IBond newBond = new Bond(atoms_for_bond, oneBond_inModified.getOrder());
						fragmentStructure.addBond(newBond);
					}
					else{
						oneBond_inFragment.setOrder(oneBond_inModified.getOrder());
					}
				}
			}
		}
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
		adder.addImplicitHydrogens(fragmentStructure);
		AtomContainerManipulator.suppressHydrogens(fragmentStructure);
		//System.out.println(sg.create(fragmentStructure));
		return fragmentStructure;
		
	}
	/**
	 * This function will merge the fragments for the reduction reaction
	 * @param reactant
	 * @param modified
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainer mergeFragmentsForReductionReaction(IAtomContainer reactant, IAtomContainer modified) throws Exception{
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		ArrayList<IAtom> atomListConnectedToReactiveSiteInReactant = getAtomListConnectedToReactiveSiteInFragment(reactant);
		ArrayList<IAtom> atomListConnectedToReactiveSiteInModified = getAtomListConnectedToReactiveSiteInFragment(modified);
		IAtom reactiveSite = null;
		IAtom changedSite = null;
		Atom restulAtom = null;
		for(int i = 0; i < reactant.getAtomCount(); i++){
			IAtom atomOne = reactant.getAtom(i);
			//Here we consider Ntriogen in nitroso reduction reaction. It can be S,P when meet corresponding predicted metabolite
			if(atomOne.getProperty("ReactiveSite")!=null && atomOne.getSymbol().equalsIgnoreCase("N")) reactiveSite = atomOne;
			for(int j = i+1; j < reactant.getAtomCount(); j++){
				IAtom anotherAtom = reactant.getAtom(j);
				if(!fragmentStructure.contains(atomOne)) fragmentStructure.addAtom(atomOne);
				if(!fragmentStructure.contains(anotherAtom)) fragmentStructure.addAtom(anotherAtom);
				if(reactant.getBond(atomOne, anotherAtom)!=null) fragmentStructure.addBond(reactant.getBond(atomOne, anotherAtom));
			}
		}
				
		for(int i = 0; i < modified.getAtomCount(); i++){
			IAtom atomOne = modified.getAtom(i);
			if(atomOne.getProperty("ReactiveSite")!=null) changedSite = atomOne;
			for(int j = i+1; j < modified.getAtomCount(); j++){
				IAtom anotherAtom = modified.getAtom(j);
				if(!fragmentStructure.contains(atomOne)) fragmentStructure.addAtom(atomOne);
				if(!fragmentStructure.contains(anotherAtom)) fragmentStructure.addAtom(anotherAtom);
				if(modified.getBond(atomOne, anotherAtom)!=null) fragmentStructure.addBond(modified.getBond(atomOne, anotherAtom));
			}
		}
		ArrayList<IBond> bondToRemove = new ArrayList<>();
		ArrayList<IAtom> atomsToRemove = new ArrayList<>();
		for(int i = 0; i < fragmentStructure.getAtomCount(); i++){
			if(fragmentStructure.getAtom(i).getProperty("ReactiveSite")!=null){
				IAtom theAtom = fragmentStructure.getAtom(i);
				try{
					IAtom corsAtom = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(theAtom), modified);
					for(int j = 0; j < fragmentStructure.getConnectedBondsList(theAtom).size(); j++){
						bondToRemove.add(fragmentStructure.getConnectedBondsList(theAtom).get(j));						
					}
				}catch (Exception e){		
					atomsToRemove.add(theAtom);					
				}				
			}
		}
		restulAtom = new Atom(changedSite.getSymbol());
		fragmentStructure.addAtom(restulAtom);
		for(int i = 0; i < atomListConnectedToReactiveSiteInReactant.size(); i++){
			IBond oneBond = reactant.getBond(atomListConnectedToReactiveSiteInReactant.get(i),reactiveSite);
			if(atomListConnectedToReactiveSiteInReactant.get(i).getSymbol().equals("H")){
				fragmentStructure.removeAtom(atomListConnectedToReactiveSiteInReactant.get(i));
				continue;
			}
			if(oneBond == null) continue;
			IAtom[] bondAtoms = {atomListConnectedToReactiveSiteInReactant.get(i), restulAtom};
			IBond newBond = new Bond(bondAtoms, oneBond.getOrder());//changedSite
			fragmentStructure.addBond(newBond);
		}
		
		
		for(int i = 0; i < atomListConnectedToReactiveSiteInModified.size(); i++){
			IBond oneBond = modified.getBond(atomListConnectedToReactiveSiteInModified.get(i),changedSite);
			if(atomListConnectedToReactiveSiteInModified.get(i).getSymbol().equals("H")){
				fragmentStructure.removeAtom(atomListConnectedToReactiveSiteInModified.get(i));
				continue;
			}
			if(oneBond == null) continue;
			IAtom[] bondAtoms = {atomListConnectedToReactiveSiteInModified.get(i), restulAtom};
			IBond newBond = new Bond(bondAtoms, oneBond.getOrder());//changedSite
			fragmentStructure.addBond(newBond);
		}
		fragmentStructure.removeAtom(reactiveSite);
		fragmentStructure.removeAtom(changedSite);
		for(int j = 0; j < bondToRemove.size(); j++) fragmentStructure.removeBond(bondToRemove.get(j));
		for(int j = 0; j < atomsToRemove.size(); j++) fragmentStructure.removeAtom(atomsToRemove.get(j));
		
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
		adder.addImplicitHydrogens(fragmentStructure);
		AtomContainerManipulator.suppressHydrogens(fragmentStructure);
		//System.out.println(sg.create(fragmentStructure));
		return fragmentStructure;
	}
	
	/**
	 * This function will merge the structure of the reactant/substrate and the modifed structure based on the reactive site
	 * The reactive sites (bonds) will be removed after the merge 
	 * Then return the merged the Structure as an IAtomContainer 
	 * Note that the reactive sites (atoms here) are labeled with the setting of property <"ReactiveSite","$ReactionType$">
	 * @param reactant
	 * @param modified
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainer mergeFragments_EpOxidation(IAtomContainer reactant, IAtomContainer modified) throws Exception{
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		ArrayList<IAtom> atomListConnectedToReactiveSiteInReactant = getAtomListConnectedToReactiveSiteInFragment(reactant);
		ArrayList<IAtom> atomListConnectedToReactiveSiteInModified = getAtomListConnectedToReactiveSiteInFragment(modified);
		ArrayList<IAtom> reactiveSiteList = new ArrayList<>();
		ArrayList<IAtom> changedSiteList = new ArrayList<>();
		ArrayList<IAtom> resultAtomList = new ArrayList<>();
		//Add all atoms and bonds within the reactant into the fragmentStructure
		for(int i = 0; i < reactant.getAtomCount(); i++){
			IAtom atomOne = reactant.getAtom(i);
			if(atomOne.getProperty("ReactiveSite")!=null){
				if(!reactiveSiteList.contains(atomOne)){
					reactiveSiteList.add(atomOne);
				}
			}
			for(int j = i+1; j < reactant.getAtomCount(); j++){
				IAtom anotherAtom = reactant.getAtom(j);
				if(!fragmentStructure.contains(atomOne)) fragmentStructure.addAtom(atomOne);
				if(!fragmentStructure.contains(anotherAtom)) fragmentStructure.addAtom(anotherAtom);
				if(reactant.getBond(atomOne, anotherAtom)!=null && !fragmentStructure.contains(reactant.getBond(atomOne, anotherAtom))){
					fragmentStructure.addBond(reactant.getBond(atomOne, anotherAtom));
				}
			}
		}
		//Add all atoms and bonds within the modified(metabolite) part into the fragmentStructure		
		for(int i = 0; i < modified.getAtomCount(); i++){
			IAtom atomOne = modified.getAtom(i);
			if(atomOne.getProperty("ReactiveSite")!=null){
				if(!changedSiteList.contains(atomOne)) changedSiteList.add(atomOne);
			}
			for(int j = i+1; j < modified.getAtomCount(); j++){
				IAtom anotherAtom = modified.getAtom(j);
				if(!fragmentStructure.contains(atomOne)) fragmentStructure.addAtom(atomOne);
				if(!fragmentStructure.contains(anotherAtom)) fragmentStructure.addAtom(anotherAtom);
				if(modified.getBond(atomOne, anotherAtom)!=null && !fragmentStructure.contains(reactant.getBond(atomOne, anotherAtom))){
					fragmentStructure.addBond(modified.getBond(atomOne, anotherAtom));
				}
			}
		}
		//Record the bond which connects the two changed atoms within the modified part. The order of the bond should retain.
		IBond temp = modified.getBond(changedSiteList.get(0), changedSiteList.get(1));
		//Create atoms representing the changed sites and add them into the fragmentStructure
		for(int i = 0; i < changedSiteList.size(); i++){
			Atom resultAtom = new Atom(changedSiteList.get(i).getSymbol());
			resultAtom.setProperty("UniqueID", UniqueIDFunctionSet.getUniqID(changedSiteList.get(i)));
			resultAtomList.add(resultAtom);
			if(!fragmentStructure.contains(resultAtom)) fragmentStructure.addAtom(resultAtom);
		}		
		IAtom[] changedPair = {resultAtomList.get(0), resultAtomList.get(1)};
		
		IBond corsBond = new Bond(changedPair, temp.getOrder());//IBond.Order.SINGLE);//changedSite
		if(!fragmentStructure.contains(corsBond)){
			fragmentStructure.addBond(corsBond);
		}
		
		for(int i = 0; i < atomListConnectedToReactiveSiteInReactant.size(); i++){		
			//Remove all explicit hydrogens
			if(atomListConnectedToReactiveSiteInReactant.get(i).getSymbol().equals("H")){
				fragmentStructure.removeAtom(atomListConnectedToReactiveSiteInReactant.get(i));
				continue;
			}
			for(int j = 0; j < reactiveSiteList.size(); j++){
				IBond oneBond = reactant.getBond(atomListConnectedToReactiveSiteInReactant.get(i),reactiveSiteList.get(j));
				if(oneBond == null) continue;
				//Note that both atoms within the resultList are carbon atoms without any other property specified. Hence the order doesn't really matter
				IAtom[] bondAtoms = {atomListConnectedToReactiveSiteInReactant.get(i), resultAtomList.get(j)};
				IBond newBond = new Bond(bondAtoms, oneBond.getOrder());//changedSite
				if(!fragmentStructure.contains(newBond)) fragmentStructure.addBond(newBond);
			}
		}
				
		for(int i = 0; i < atomListConnectedToReactiveSiteInModified.size(); i++){
			//Remove all explicit hydrogens
			if(atomListConnectedToReactiveSiteInModified.get(i).getSymbol().equals("H")){
				fragmentStructure.removeAtom(atomListConnectedToReactiveSiteInModified.get(i));
				continue;
			}
			for(int j = 0; j < changedSiteList.size(); j++){
				IBond oneBond = modified.getBond(atomListConnectedToReactiveSiteInModified.get(i),changedSiteList.get(j));
				if(oneBond == null) continue;
				IAtom[] bondAtoms = {atomListConnectedToReactiveSiteInModified.get(i), resultAtomList.get(j)};
				IBond newBond = new Bond(bondAtoms, oneBond.getOrder());//changedSite
				if(!fragmentStructure.contains(newBond)) fragmentStructure.addBond(newBond);
				
			}
		}
		for(int i = 0; i < reactiveSiteList.size(); i++){
			fragmentStructure.removeAtom(reactiveSiteList.get(i));
			fragmentStructure.removeAtom(changedSiteList.get(i));
		}
		
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(fragmentStructure);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
		adder.addImplicitHydrogens(fragmentStructure);
		AtomContainerManipulator.suppressHydrogens(fragmentStructure);
		//System.out.println(sg.create(fragmentStructure));
		return fragmentStructure;
	}
	/**
	 * This function will find the atoms that are connected to the reactive site in both subStrate and the Metabolite
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<IAtom> getAtomListConnectedToReactiveSiteInFragment(IAtomContainer oneMole) throws Exception{
		ArrayList<IAtom> resultList = new ArrayList<>();
		for(int i = 0; i < oneMole.getAtomCount(); i++){
			IAtom oneAtom = oneMole.getAtom(i);
			List<IAtom> neighbors = oneMole.getConnectedAtomsList(oneAtom);
			for(int j = 0; j < neighbors.size(); j++){
				IAtom tempAtom = neighbors.get(j);
				if(tempAtom.getProperty("ReactiveSite")!=null){
					if(!resultList.contains(oneAtom)) resultList.add(oneAtom);
					break;
				}
						
			}
		}
		return resultList;
	}
	/**
	 * This function will add explicit hydrogens to the "oneAtom" and return them together as IAtomContainer
	 * Note that it is important to add explicit hydrogens when using SMIRKS
	 * @param oneAtom
	 * @param oneMole
	 * @return
	 */
	public static IAtomContainer generateFragmentWithExplicitHydrogenForReactiveAtom(IAtom atom_fragment, IAtomContainer oneMole) throws Exception{
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtomContainer fragmentStructure = builder.newInstance(IAtomContainer.class);
		IAtom oneAtom = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(atom_fragment), oneMole);
		for(int i = 0; i < oneMole.getConnectedAtomsList(oneAtom).size(); i++){
			if(oneMole.getConnectedAtomsList(oneAtom).get(i).getSymbol().equalsIgnoreCase("H")){
				IAtom neighbor_copy = oneMole.getConnectedAtomsList(oneAtom).get(i).clone();
				//fragmentStructure.addAtom(oneMole.getConnectedAtomsList(oneAtom).get(i));
				fragmentStructure.addAtom(neighbor_copy);
				//fragmentStructure.addBond(oneMole.getBond(oneMole.getConnectedAtomsList(oneAtom).get(i), oneAtom));
				IBond oneBond = oneMole.getBond(oneMole.getConnectedAtomsList(oneAtom).get(i), oneAtom);
				IAtom[] bondAtoms = {atom_fragment, neighbor_copy};
				IBond newBond = new Bond(bondAtoms, oneBond.getOrder());
				fragmentStructure.addBond(newBond);
				
			}
		}
		return fragmentStructure;
	}
	/**
	 * This function will be used in those reactions that involves the breakage of <eta,eta> bonds
	 * Here, it is used for Dealkylation reaction specifically.
	 * This means that the IAtomContainerSet molecules should be size = 2
	 * @param molecules
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<IAtomContainer> reorderFragments(IAtomContainerSet molecules) throws Exception{
		ArrayList<IAtomContainer> result = new ArrayList<>();
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			IAtomContainer fragment = molecules.getAtomContainer(i);
			for(int j = 0; j < fragment.getAtomCount(); j++){
				if(fragment.getAtom(j).getProperty("ReactiveSite")!=null){
					if(!fragment.getAtom(j).getSymbol().equalsIgnoreCase("C")){
						result.add(fragment);
						break;
					}
				}
			}
		}
		for(int i = 0; i < molecules.getAtomContainerCount(); i++){
			if(!result.contains(molecules.getAtomContainer(i))) result.add(molecules.getAtomContainer(i));
		}
		return result;
	}
	/**
	 * This function will refine the aromatic part that contains both left and right atoms by altering the SINGLE and DOUBLE bonds one by one
	 * @param leftAtom
	 * @param righAtom
	 * @param oneMole
	 * @throws Exception
	 */
	public static IAtomContainer refineAromaticPart(IAtom leftAtom, IAtom rightAtom, IAtomContainer oneMole) throws Exception{
		AllRingsFinder arf = new AllRingsFinder();
		IRingSet allRings = arf.findAllRings(oneMole);
		//IAtomContainerSet candidateSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<IAtomContainer> candidateSet = new ArrayList<>();
		for(int i = 0; i < allRings.getAtomContainerCount(); i++){
			if(allRings.contains(leftAtom) && allRings.contains(rightAtom)){
				IAtomContainer oneCandidate = allRings.getAtomContainer(i);
				if(!candidateSet.contains(oneCandidate)) candidateSet.add(oneCandidate);
			}
		}
		ArrayList<IAtomContainer> reOrderedCandidateList = new ArrayList<>();
		while(candidateSet!=null && !candidateSet.isEmpty()){
			IAtomContainer bestContainer = getMaxRing(candidateSet);
			reOrderedCandidateList.add(bestContainer);
			candidateSet.remove(bestContainer);
		}		
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		//Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());//new Aromaticity(ElectronDonation.cdkAllowingExocyclic(), Cycles.cdkAromaticSet());
		IAtomContainer largestAromaticRing = null;
		for(int i = 0; i < reOrderedCandidateList.size(); i++){
			IAtomContainer ringToCheck  = reOrderedCandidateList.get(i);
			//if(aromaticity.apply(ringToCheck)){
			if(checkAromatic(ringToCheck, oneMole)){
				if(largestAromaticRing == null || largestAromaticRing.getAtomCount() < ringToCheck.getAtomCount()){			
					largestAromaticRing = ringToCheck;
				}
				break;
			}
		}
		ArrayList<IBond> shardedBondList = new ArrayList<>();
		if(largestAromaticRing==null){
			throw new Exception("Failed to find the largestAromaticRing in the molecule");
		}
		else{
			//System.out.println("Largest Ring to Check: " + sg.create(largestAromaticRing));
			//System.out.println("Indices of atoms within the largest ring: ") ;
			for(int i = 0; i < largestAromaticRing.getBondCount(); i++){
				//System.out.println(largestAromaticRing.getBond(i).getAtom(0).getProperty("UniqueID") + "," + largestAromaticRing.getBond(i).getAtom(1).getProperty("UniqueID")) ;
			}
			//System.out.println() ;
			IBond targetBond = largestAromaticRing.getBond(leftAtom, rightAtom);
			
			if(targetBond.getOrder() == IBond.Order.SINGLE){
				alterSingleDoubleBondsInRing(targetBond,largestAromaticRing);
				//Add the bonds shared by more than one rings into the largestAromaicRing. Note that those bonds were removed when searching cycles
 				for(int t = 0; t < reOrderedCandidateList.size(); t++){
					//Take the rings from the smallest to the largest
					IAtomContainer oneSmallRing = reOrderedCandidateList.get(reOrderedCandidateList.size() - 1 - t);
					AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneSmallRing);
					if(!checkAromatic(oneSmallRing, oneMole)){
						//System.out.println("Invalid Ring to Check: " + sg.create(oneSmallRing));
						continue;
					}
					
					//If the atoms within the small ring do not exist in the largest ring, add them into the largest ring.
					for(int j = 0; j < oneSmallRing.getAtomCount(); j++){
						if(!largestAromaticRing.contains(oneSmallRing.getAtom(j))) largestAromaticRing.addAtom(oneSmallRing.getAtom(j));
					}
					//If the bonds within the smll ring do not exist in the largest ring, add them into the largest ring.
					for(int j = 0; j < oneSmallRing.getBondCount(); j++){
						if(!largestAromaticRing.contains(oneSmallRing.getBond(j))) largestAromaticRing.addBond(oneSmallRing.getBond(j));
					}
				}
 				//System.out.println("ring of interest: " + sg.create(largestAromaticRing));
 				//Add the bonds that are shared by more than one rings into the shardedBondList. 
 				for(int i = 0; i < largestAromaticRing.getBondCount(); i++){
 					boolean isSharedBond = false;
 					IBond oneBond = largestAromaticRing.getBond(i);

 					if(computeBondOrderSum(oneBond.getAtom(0),largestAromaticRing) >= 4 && computeBondOrderSum(oneBond.getAtom(1),largestAromaticRing) >= 4){
 						isSharedBond = true;
 					}
 					if(isSharedBond){
 						shardedBondList.add(oneBond);
 						//System.out.println("Bond: " + oneBond.getAtom(0).getProperty("UniqueID") + " , " + oneBond.getAtom(1).getProperty("UniqueID"));
 					}
 				}
 				////System.out.println("Refined Structure: " + sg.create(oneMole));
				//Add the bonds that are shared by more than one Rings
				for(int i = 0; i < shardedBondList.size(); i++){
					IBond oneBond = shardedBondList.get(i);
					if(oneBond.equals(targetBond)) continue;
					if(computeBondOrderSum(oneBond.getAtom(0),largestAromaticRing) > 4 && computeBondOrderSum(oneBond.getAtom(1),largestAromaticRing) > 4){
						oneBond.setOrder(IBond.Order.SINGLE);
					}
				}
				//Make sure that the target bond is still "DOUBLE"
				
				for(int i = 0; i < shardedBondList.size(); i++){
					IBond oneBond = shardedBondList.get(i);
					if(computeBondOrderSum(oneBond.getAtom(0),largestAromaticRing) <=3 && computeBondOrderSum(oneBond.getAtom(1),largestAromaticRing) <=3){
						oneBond.setOrder(IBond.Order.DOUBLE);
					}
				}
				for(int i = 0; i < shardedBondList.size(); i++){
					IBond oneBond = shardedBondList.get(i);
					if(oneBond.equals(targetBond)) continue;
					if(computeBondOrderSum(oneBond.getAtom(0),largestAromaticRing) > 4 && computeBondOrderSum(oneBond.getAtom(1),largestAromaticRing) > 4){
						oneBond.setOrder(IBond.Order.SINGLE);
					}
				}
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(largestAromaticRing);
				//System.out.println("Ring to Check: " + sg.create(largestAromaticRing));
				////System.out.println("ring of interest: " + sg.create(largestAromaticRing));
								
			}
			
		}
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(largestAromaticRing);
		//System.out.println("Refined ring: " + sg.create(largestAromaticRing));
		//System.out.println("Refined Structure: " + sg.create(oneMole));
		
		SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
		return sp.parseSmiles(sg.create(oneMole));
		
	}
	/**
	 * This function will find the IAtomContainer candidate that has the most atoms 
	 * @param candidateSet
	 * @return
	 * @throws Exception
	 */
	public static IAtomContainer getMaxRing(ArrayList<IAtomContainer> candidateSet) throws Exception{
		IAtomContainer bestContainer = null;
		Integer maxAtoms = 0;
		Integer maxBonds = 0;
		for(int i = 0; i < candidateSet.size(); i++){
//			IAtomContainer completeFormat = candidateSet.get(i).clone();
//			for(int t = 0; t < candidateSet.size(); t++){
//				IAtomContainer oneSmallRing = candidateSet.get(t);
//				for(int j = 0; j < oneSmallRing.getAtomCount(); j++){
//					if(!completeFormat.contains(oneSmallRing.getAtom(j))) completeFormat.addAtom(oneSmallRing.getAtom(j));
//				}
//				//If the bonds within the smll ring do not exist in the largest ring, add them into the largest ring.
//				for(int j = 0; j < oneSmallRing.getBondCount(); j++){
//					if(!completeFormat.contains(oneSmallRing.getBond(j))) completeFormat.addBond(oneSmallRing.getBond(j));
//				}
//			}
			if(maxAtoms < candidateSet.get(i).getAtomCount()){
				maxAtoms = candidateSet.get(i).getAtomCount();
				maxBonds = candidateSet.get(i).getBondCount();
				bestContainer = candidateSet.get(i);
			}
			else if(maxAtoms == candidateSet.get(i).getAtomCount()){
				if(maxBonds > candidateSet.get(i).getBondCount()){
					maxBonds = candidateSet.get(i).getBondCount();
					bestContainer = candidateSet.get(i);
				}
			}
		}
		return bestContainer;
	}
	/**
	 * This function will check if the input structure is aromatic
	 * The structure is aromatic if all the atoms within it are aromatic within the original molecule
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static boolean checkAromatic(IAtomContainer structure, IAtomContainer oneMole) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(oneMole);
		for(int i = 0; i < structure.getAtomCount(); i++){
			IAtom oneAtom = structure.getAtom(i);
			Integer uniqueID = UniqueIDFunctionSet.getUniqID(oneAtom);
			oneAtom = UniqueIDFunctionSet.getAtomByUniqueID(uniqueID, oneMole);
			
			if(!oneAtom.isAromatic()) return false;
		}
		return true;
	}
	
	public static void alterSingleDoubleBondsInRing(IBond targetBond, IAtomContainer largestAromaticRing) throws Exception{
		if(targetBond.getOrder() != IBond.Order.SINGLE) throw new Exception("The target bond is not a SINGLE bond and thus should not be changed to DOUBLE bond");
		targetBond.setOrder(IBond.Order.DOUBLE);
		//ArrayList<IBond> checkedBondList = new ArrayList<>();
		ArrayList<IAtom> checkedAtomList = new ArrayList<>();
		//both terminate_left and terminate_right are used to indicate there is a terminating condition matched
		//They are not necessary and should not affect the performance at all
		boolean terminate_left = false;
		boolean terminate_right = false;
		IAtom atom_left = targetBond.getAtom(0);
		IAtom atom_right = targetBond.getAtom(1);
		//checkedBondList.add(targetBond);
		checkedAtomList.add(atom_left);
		checkedAtomList.add(atom_right);
		
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Canonical);
		
		while(checkedAtomList.size() < largestAromaticRing.getAtomCount() && !terminate_left){		 
			List<IAtom> atoms_toCheck = largestAromaticRing.getConnectedAtomsList(atom_left);
			for(int i = 0; i < atoms_toCheck.size(); i++){
				IAtom checkAtom = atoms_toCheck.get(i);
				if(!checkedAtomList.contains(checkAtom)){
					////System.out.println(atom_left.getBondOrderSum());
					IBond bond_candidate = largestAromaticRing.getBond(atom_left, checkAtom);
					//If atom_left has two DOUBLE bonds connected. Note that the current largestAromaticRing doesn't contain bonds shared by multiple rings.
					if(computeBondOrderSum(atom_left, largestAromaticRing)>=4){						
						if(bond_candidate.getOrder() == IBond.Order.DOUBLE) bond_candidate.setOrder(IBond.Order.SINGLE);
					}
					else if(computeBondOrderSum(atom_left, largestAromaticRing) <=2){
						if(bond_candidate.getOrder() == IBond.Order.SINGLE) bond_candidate.setOrder(IBond.Order.DOUBLE);
					}
					else if(computeBondOrderSum(atom_left, largestAromaticRing) == 3){
						terminate_left = true;
						break;
					}
					checkedAtomList.add(checkAtom);
					atom_left = checkAtom;
					//Check checkAtom, if checkAtom has invalid valence, then alter it.
					//Otherwise set terminate_left to true and break  the while loop										
					break;
				}
				//System.out.println("Alter Ring: " + sg.create(largestAromaticRing));
			}			
		}
		checkedAtomList = new ArrayList<>();
		checkedAtomList.add(targetBond.getAtom(0));
		checkedAtomList.add(targetBond.getAtom(1));
		while(checkedAtomList.size() < largestAromaticRing.getAtomCount() && !terminate_right){		 
			List<IAtom> atoms_toCheck = largestAromaticRing.getConnectedAtomsList(atom_right);
			for(int i = 0; i < atoms_toCheck.size(); i++){
				IAtom checkAtom = atoms_toCheck.get(i);
				if(!checkedAtomList.contains(checkAtom)){
					IBond bond_candidate = largestAromaticRing.getBond(atom_right, checkAtom);
					//If atom_left has two DOUBLE bonds connected. Note that the current largestAromaticRing doesn't contain bonds shared by multiple rings.
					if(computeBondOrderSum(atom_right, largestAromaticRing)>=4){						
						if(bond_candidate.getOrder() == IBond.Order.DOUBLE) bond_candidate.setOrder(IBond.Order.SINGLE);
					}
					else if(computeBondOrderSum(atom_right, largestAromaticRing) <=2){
						if(bond_candidate.getOrder() == IBond.Order.SINGLE) bond_candidate.setOrder(IBond.Order.DOUBLE);
					}
					else if(computeBondOrderSum(atom_right, largestAromaticRing) == 3){
						terminate_right = true;
						break;
					}
					checkedAtomList.add(checkAtom);
					atom_right = checkAtom;
					//Check checkAtom, if checkAtom has invalid valence, then alter it.
					//Otherwise set terminate_left to true and break  the while loop										
					break;
				}
			}			
		}
	}
	/**
	 * This function will find the non-hydrogen atoms connected to the "oneAtom" within the largestAromaticRing and 
	 * compute the sum of orders of those bonds.
	 * The sum is returned as integer
	 * @param oneAtom
	 * @param largestAromaticRing
	 * @return
	 * @throws Exception
	 */
	public static Integer computeBondOrderSum(IAtom oneAtom, IAtomContainer largestAromaticRing) throws Exception{
		Integer sum = 0;
		List<IAtom> connectedAtoms = largestAromaticRing.getConnectedAtomsList(oneAtom);
		for(int i = 0; i < connectedAtoms.size(); i++){
			IAtom nextAtom = connectedAtoms.get(i);
			if(largestAromaticRing.getBond(oneAtom, nextAtom) != null){
				if(largestAromaticRing.getBond(oneAtom,nextAtom).getOrder() == IBond.Order.SINGLE) sum = sum + 1;
				else if(largestAromaticRing.getBond(oneAtom,nextAtom).getOrder() == IBond.Order.DOUBLE) sum = sum + 2;
				else if(largestAromaticRing.getBond(oneAtom,nextAtom).getOrder() == IBond.Order.TRIPLE) sum = sum + 3;
				else if(largestAromaticRing.getBond(oneAtom,nextAtom).getOrder() == IBond.Order.QUADRUPLE) sum = sum + 4;
			}
		}
		return sum;
	}

}		

