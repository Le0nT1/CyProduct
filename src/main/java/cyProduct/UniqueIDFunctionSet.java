package cyProduct;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * This class will assign uniqueIDs for each atom within the molecule.
 * This assignment should only be done once.
 * The purpose of uniqueIDs is to break the limitation of CDK toolkit especially when .clone method is used
 * @author Tian
 *
 */
public class UniqueIDFunctionSet {
	/**
	 * This function will assign an uniqueID for each atom within the molecule
	 * @param oneMole
	 */
	public static void assignUniqueID(IAtomContainer oneMole){
		int uniqueID = 1;
		for(int i = 0; i < oneMole.getAtomCount(); i++){
			oneMole.getAtom(i).setProperty("UniqueID", uniqueID);
			uniqueID++;
		}
	}	
	
	/**
	 * This function will assign the given uid to the given atom within the structure
	 */
	public static void assingUID_toAtom(IAtom oneAtom, Integer uid){
		oneAtom.setProperty("UniqueID", uid);
	}
	
	
	/**
	 * This function will return the uniqueID of the input atom, based on its UniqueID property
	 * @param oneAtom
	 * @return
	 * @throws Exception
	 */
	public static Integer getUniqID(IAtom oneAtom) throws Exception{
		Integer uniqueID = oneAtom.getProperty("UniqueID");
		//if(uniqueID==null) throw new Exception("The Atom doesn't have the UniqueID property");
		return uniqueID;
	}
	/**
	 * This function will get the IAtom based having the input uniqueID within the input molecule
	 * @param uniqueID
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static IAtom getAtomByUniqueID(Integer uniqueID, IAtomContainer oneMole) throws Exception{
		for(int i = 0; i < oneMole.getAtomCount(); i++){
			IAtom oneAtom = oneMole.getAtom(i);
			Integer tempID = oneAtom.getProperty("UniqueID");
			if(tempID == uniqueID) return oneAtom;
		}
		throw new Exception("Cannot find the Atom having the input ID");
	}
	
	/**
	 * This function labels the connected point on the substructure after the "principle" (prior) chain is removed
	 * The bondtype is a property of an neighbor atom, it's very important for those atoms connected to the prior chain
	 * @param visited_PChain
	 * @param oneMole
	 * @throws Exception
	 */
	public static void labelBranchCNAtom(IAtomContainer child, IAtomContainer oneMole) throws Exception{
		ArrayList<IAtom> checkVisited = new ArrayList<>();
		int count = 0;
		for(int i = 0; i < child.getAtomCount(); i++){
			IAtom oneAtom_temp = child.getAtom(i);
			Integer uniqueID = oneAtom_temp.getProperty("UniqueID");
			IAtom oneAtom = getAtomByUniqueID(uniqueID, oneMole);
			for(int j = 0; j < oneMole.getConnectedAtomsList(oneAtom).size(); j++){
				IAtom oneNeighbor = oneMole.getConnectedAtomsList(oneAtom).get(j);
				if(checkVisited.contains(oneNeighbor)){
					continue;
				}
				if(!containAtom(oneNeighbor, child)&&containAtom(oneNeighbor, oneMole)){					
					//Set label to indicate that the atom is connected to the prior chain and record the bondtype that can be used in reconstruction
					String bondType = String.valueOf(oneMole.getBond(oneAtom, oneNeighbor).getOrder());
					oneAtom_temp.setProperty("BondType", bondType);
					oneAtom_temp.setProperty("ConnectPriorChain", "True");
					String cis_trans = oneMole.getBond(oneAtom, oneNeighbor).getProperty("Cis_Trans_Config");
					oneAtom_temp.setProperty("StereoToPrior", cis_trans);
					oneAtom_temp.setProperty("PriorAtom", oneNeighbor);
					if(!oneAtom_temp.getSymbol().equalsIgnoreCase("C")){
						List<IAtom> candidate = child.getConnectedAtomsList(oneAtom_temp);
						if(candidate.size() != 1) throw new Exception("Cannot locate the atom connect to the ehter oxygen on the branch");
						IAtom targetAtom = candidate.get(0);
						bondType = String.valueOf(child.getBond(oneAtom_temp, targetAtom).getOrder());
						targetAtom.setProperty("BondType", bondType);
						targetAtom.setProperty("ConnectPriorChain", "True");
						cis_trans = child.getBond(oneAtom_temp, targetAtom).getProperty("Cis_Trans_Config");
						targetAtom.setProperty("StereoToPrior", cis_trans);
						targetAtom.setProperty("PriorAtom", oneAtom_temp);
					}
										
					count++;
				}
				//checkVisited.add(oneNeighbor);
			}
			//checkVisited.add(oneAtom);
		}
		//System.out.println("# Labeled Atoms: " + count);
	}
	/**
	 * This function will check if the oneAtom is part of the oneMole by checking the UID.
	 * @param oneAtom
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static boolean containAtom(IAtom oneAtom, IAtomContainer oneMole) throws Exception{
		Integer uniqueID = oneAtom.getProperty("UniqueID");
		for(int i = 0; i < oneMole.getAtomCount(); i++){
			Integer tempID = oneMole.getAtom(i).getProperty("UniqueID");
			if(uniqueID == tempID) return true;
		}
		return false;
	}
	
	/**
	 * This function will find the largest Unique ID for the atoms within the molecule.
	 * This can be used when assigning new Unique ID for newly added atoms
	 */
	public static Integer getLargestUID(IAtomContainer oneMole) throws Exception{
		int max_uid = 0;
		for(int i = 0; i < oneMole.getAtomCount(); i++){
			if(UniqueIDFunctionSet.getUniqID(oneMole.getAtom(i))> max_uid){
				max_uid = UniqueIDFunctionSet.getUniqID(oneMole.getAtom(i));
			}
		}
		return max_uid;
	}
	
	
	
	
}
