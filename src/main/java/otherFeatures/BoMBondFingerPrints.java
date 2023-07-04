package otherFeatures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class BoMBondFingerPrints {
	/**
	 * This function returns an ArrayList<IBond> fingerprints for bonds at given radius
	 * @param radius: The radius to search
	 * @param molecule: The molecule being processed
	 * @param neighAtoms: The neighbor atoms of the given node.
	 * @param node: The center when searching 
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<IBond> FindBondsByRadius(int radius, IAtomContainer molecule, List<IAtom> neighAtoms, IAtom node) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		int length = GeneratePatterns.queriesList.size();
		//String[] fpOfRadius = new String[length];
		ArrayList<IBond> bondsByRadius = new ArrayList<>();
		//Condition section
		for(int atomIdx = 0; atomIdx < neighAtoms.size(); atomIdx++){
			IAtom neighbour = neighAtoms.get(atomIdx);
			IBond oneBond = molecule.getBond(neighbour, node);
			if(radius == 1){//ending condition
				bondsByRadius.add(oneBond);
			}
			else{
				List<IAtom> atom_ConnectList = molecule.getConnectedAtomsList(neighbour);
//				if(atom_ConnectList.isEmpty()||atom_ConnectList==null){
//					return null;
//				}
				atom_ConnectList.remove(node);//The List will be empty, eg. [] if no neighbors are found.
//				if(atom_ConnectList.isEmpty()||atom_ConnectList==null){
//					return null;
//				}
				ArrayList<IBond> tempBonds = FindBondsByRadius(radius-1, molecule, atom_ConnectList, neighbour);//Recursively call this function till depth=1
				bondsByRadius.addAll(tempBonds);
				
			}
		}					
		return bondsByRadius;
	}
	/**
	 * Generate union fingerprints for neighbour bonds. If at least on of the neighbor bonds matches the patter, the value of the tuple in String[] is set to "1".
	 * If the possible radius is smaller than the input depth, the bits for those exceeded parts will all be "0"
	 * @param depth: The radius when searching the neighbor bonds
	 * @param molecule
	 * @param neighAtoms
	 * @param node
	 * @param length
	 * @return A StringList consists of "0" and "1".
	 * @throws Exception
	 */
	public static String[] generateMergedTypeOneNeighbourBondFingerPrints(int depth, IAtomContainer molecule, List<IAtom> neighAtoms, IAtom node, int length) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);		
		String[] neighbourFingerPrints = new String[length];
		//Initialize all finger print to 0
		for(int i = 0; i < neighbourFingerPrints.length; i++){
			neighbourFingerPrints[i] = "0";
		}
		for(int atomIdx = 0; atomIdx < neighAtoms.size(); atomIdx++){
			//List<IBond> neighbourBonds = molecule.getConnectedBondsList(neighAtoms.get(atomIdx));
			IAtom neighbour = neighAtoms.get(atomIdx);
			IBond oneBond = molecule.getBond(neighbour, node);
			GeneratePatterns gn = new GeneratePatterns();
			String[] currentBond = gn.generateClassyfireBondFingeprint(molecule, oneBond);
			for(int i = 0; i < currentBond.length; i++){
				if(currentBond[i].equals("1")){
					neighbourFingerPrints[i] = "1";
				}
			}
			List<IAtom> nextNeighbor = molecule.getConnectedAtomsList(neighbour);
			if((depth>0) && (nextNeighbor!=null)){
				nextNeighbor.remove(node);//Remove the node from the neighbor atoms list of the "neighbor" atom.
				String[] nextNeighborBond = generateMergedTypeOneNeighbourBondFingerPrints(depth-1, molecule, nextNeighbor, neighbour, length);
				for(int i = 0; i < nextNeighborBond.length; i++){
					if(nextNeighborBond[i].equals("1")){
						neighbourFingerPrints[i] = "1";
					}
				}
			}
		}
		return neighbourFingerPrints;
	}
	/**
	 * Generate fingerprints for neighbour bonds. If at least on of the neighbor bonds matches the patter, the value of the tuple in String[] is set to "1".
	 * If the possible radius is smaller than the input depth, the bits for those exceeded parts will all be "0"
	 * @param depth: The radius when searching the neighbor bonds
	 * @param molecule
	 * @param neighAtoms
	 * @param node
	 * @param length
	 * @return A StringList consists of "0" and "1".
	 * @throws Exception
	 */
	public static String[] generateTypeOneNeighbourBondFingerPrints(int depth, IAtomContainer molecule, List<IAtom> neighAtoms, IAtom node, int length) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);		
		String[] neighbourFingerPrints = new String[length*depth];
		//Initialize all finger print to 0
		for(int i = 0; i < neighbourFingerPrints.length; i++){
			neighbourFingerPrints[i] = "0";
		}
		
		for(int radius = 1; radius < (depth+1); radius++){
			//int radius = upToDepth;//Try radius = 1 first.
			ArrayList<IBond> ibondList = new ArrayList<>();		
			ibondList = FindBondsByRadius(radius, molecule, neighAtoms, node);//get all bonds at the given radius 
			String[] fpByRadius = new String[length];
			for(int j = 0; j < fpByRadius.length; j++){
				fpByRadius[j] = "0";
			}
			for(int bondIdx = 0; bondIdx < ibondList.size(); bondIdx++){//Iterate through every bond
				IBond oneBond = ibondList.get(bondIdx);
				GeneratePatterns gn = new GeneratePatterns();
				String[] currentBond = gn.generateClassyfireBondFingeprint(molecule, oneBond);//get the fingerprints of the bond
				for(int i = 0; i < currentBond.length; i++){
					if(currentBond[i].equals("1")){
						neighbourFingerPrints[i + (radius-1)*length] = "1";
						fpByRadius[i] = "1"; 
								
					}
				}
			}
			//System.out.println("Radius = " + radius + ", the size of its finger prints is " + fpByRadius.length);
			//System.out.println(Arrays.toString(fpByRadius));
		}
		/*
		//For every bond at this depth
		for(int atomIdx = 0; atomIdx < neighAtoms.size(); atomIdx++){
			//List<IBond> neighbourBonds = molecule.getConnectedBondsList(neighAtoms.get(atomIdx));
			IAtom neighbour = neighAtoms.get(atomIdx);
			IBond oneBond = molecule.getBond(neighbour, node);
			GeneratePatterns gn = new GeneratePatterns();
			String[] currentBond = gn.generateClassyfireBondFingeprint(molecule, oneBond);
			for(int i = 0; i < currentBond.length; i++){
				if(currentBond[i].equals("1")){
					neighbourFingerPrints[i] = "1";
				}
			}
			List<IAtom> nextNeighbor = molecule.getConnectedAtomsList(neighbour);
			if((depth>1) && (nextNeighbor!=null)){
				nextNeighbor.remove(node);//Remove the node from the neighbor atoms list of the "neighbor" atom.
				String[] nextNeighborBond = generateTypeOneNeighbourBondFingerPrints(depth-1, molecule, nextNeighbor, neighbour, length);
				for(int i = 0; i < nextNeighborBond.length; i++){
					if(nextNeighborBond[i].equals("1")){
						neighbourFingerPrints[i] = "1";
					}
				}
			}
		}
		*/
		return neighbourFingerPrints;
	}
	/**
	 * Generate finger prints for neighbor bonds. If at least on of the neighbor bonds matches the pattern, the value of the tuple in String[] is set to "1".
	 * This function is replaced by generateTypeOneNeighbourBondFingerPrints(int depth, IAtomContainer molecule, List<IAtom> neighAtoms, IAtom node, int length)
	 * @param molecule
	 * @param neighAtoms
	 * @param node
	 * @param length: the size of the fingerprints
	 * @return a String[] contains the merges BondFingerPrints of the neighbour atoms of the "node" atom in the "molecule" with depth = 2
	 * @throws Exception
	 */
	public String[] generateTypeOneNeighbourBondFingerPrints(IAtomContainer molecule, List<IAtom> neighAtoms, IAtom node, int length) throws Exception{
		String[] neighbourFingerPrints = new String[length];
		for(int i = 0; i < neighbourFingerPrints.length; i++){
			neighbourFingerPrints[i] = "0";
		}
		for(int atomIdx = 0; atomIdx < neighAtoms.size(); atomIdx++){
			//List<IBond> neighbourBonds = molecule.getConnectedBondsList(neighAtoms.get(atomIdx));
			IAtom neighbour = neighAtoms.get(atomIdx);
			IBond oneBond = molecule.getBond(neighbour, node);
			if(oneBond==null){
				oneBond = molecule.getBond(node,neighbour);
			}
			GeneratePatterns gn = new GeneratePatterns();
			String[] currentBond = gn.generateClassyfireBondFingeprint(molecule, oneBond);
			for(int i = 0; i < currentBond.length; i++){
				if(currentBond[i].equals("1")){
					neighbourFingerPrints[i] = "1";
				}
			}			
		}
		return neighbourFingerPrints;
	}
	
	/**
	 * Generate fingerprints for neighbour bonds for TypeTwo bonds--<i,H> bond. If at least on of the neighbor bonds matches the patter, the value of the tuple in String[] is set to "1".
	 * @param depth: The radius when searching the neighbor bonds.eg: [detph=1,depth=2,...]
	 * @param molecule
	 * @param neighAtoms
	 * @param node
	 * @param length
	 * @return
	 * @throws Exception
	 */
	public static String[] generateTypeTwoNeighbourBondFingerPrints(IAtomContainer molecule, List<IAtom> neighAtoms, IAtom node, int depth, boolean mergeFP) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		int length = GeneratePatterns.queriesList.size();
		String[] neighbourFingerPrints = new String[length*depth];//*upToDepth	 assign space the size of the neighbourFingerPrints list
		for(int i = 0; i < neighbourFingerPrints.length; i++){//initialize the value of each bit as "0"
			neighbourFingerPrints[i] = "0";
		}
		for(int radius = 1; radius < (depth+1); radius++){
			//int radius = upToDepth;//Try radius = 1 first.
			ArrayList<IBond> ibondList = new ArrayList<>();		
			ibondList = FindBondsByRadius(radius, molecule, neighAtoms, node);//get all bonds at the given radius 
			String[] fpByRadius = new String[length];
			for(int j = 0; j < fpByRadius.length; j++){
				fpByRadius[j] = "0";
			}
			for(int bondIdx = 0; bondIdx < ibondList.size(); bondIdx++){//Iterate through every bond
				IBond oneBond = ibondList.get(bondIdx);
				GeneratePatterns gn = new GeneratePatterns();
				String[] currentBond = gn.generateClassyfireBondFingeprint(molecule, oneBond);//get the fingerprints of the bond
				for(int i = 0; i < currentBond.length; i++){
					if(currentBond[i].equals("1")){
						neighbourFingerPrints[i + (radius-1)*length] = "1";
						fpByRadius[i] = "1"; 
								
					}
				}
			}
			//System.out.println("Radius = " + radius + ", the size of its finger prints is " + fpByRadius.length);
			//System.out.println(Arrays.toString(fpByRadius));
		}
		return neighbourFingerPrints;
	}
}
