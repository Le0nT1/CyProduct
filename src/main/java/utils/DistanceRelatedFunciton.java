package utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.vecmath.Point3d;

import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * This class will have some functions that are related to the distance between atoms
 * @author Tian
 *
 */
public class DistanceRelatedFunciton {
	PathTools pathFinder = new PathTools();
	/**
	 * This function will initialize the ClosestAtoms property for each atom.
	 * Input:
	 * IAtomContainer oneMole: the molecule to be initialized
	 * Integer numAtoms: the number of closest atoms by the distance
	 * @param oneMole
	 * @param numAtoms
	 * @throws Exception
	 */
	public void initializeClosestAtomsProperty(IAtomContainer oneMole, Integer numAtoms) throws Exception{		
		for(int i = 0; i < oneMole.getAtomCount(); i++){
			IAtom oneAtom = oneMole.getAtom(i);			
			HashMap<IAtom,Double> atomDisMap = new HashMap<>();
			for(int j = 0; j < oneMole.getAtomCount(); j++){
				if(i == j) continue;
				IAtom nextAtom = oneMole.getAtom(j);
				Double distance = calculateDistance(oneAtom, nextAtom);
				atomDisMap.put(nextAtom, distance);				
			}
			ArrayList<IAtom> closestAtoms = getCloesestAtomListInMap(atomDisMap, numAtoms);			
			oneAtom.setProperty("Closest", closestAtoms);			
		}
	}
	/**
	 * This function will get the numAtoms many closest atoms for the oneAtom within the molecule
	 * @param oneMole
	 * @param oneAtom
	 * @param numAtoms
	 * @return
	 */
	public ArrayList<IAtom> getClosestAtomListInMolecule_3D(IAtomContainer oneMole, IAtom oneAtom, Integer numAtoms){
		HashMap<IAtom,Double> atomDisMap = new HashMap<>();
		for(int i = 0; i < oneMole.getAtomCount(); i++){
			IAtom closeAtom = oneMole.getAtom(i);
			if(!oneAtom.equals(closeAtom)){
				Double distance = calculateDistance(oneAtom, closeAtom);
				atomDisMap.put(closeAtom, distance);	
			}
		}
		ArrayList<IAtom> closestAtoms = getCloesestAtomListInMap(atomDisMap, numAtoms);
		return closestAtoms;
	}
	
	/**
	 * This function will get the numAtoms many closest atoms for the oneAtom within the molecule
	 * @param oneMole
	 * @param oneAtom
	 * @param numAtoms
	 * @return
	 */
	public ArrayList<IAtom> getClosestAtomListInMolecule_NoCoordinates(IAtomContainer oneMole, IAtom oneAtom, Integer numAtoms){
		ArrayList<IAtom> closestAtoms = new ArrayList<>();
		int length = 1;
		//add condition: length < numAtoms, so it won't fall into the infinite loop case 
		while(closestAtoms.size() < numAtoms && length < numAtoms-1){
			if(oneMole.getAtomCount() < numAtoms) break;
			else{
				List<List<IAtom>> allPaths = PathTools.getPathsOfLength(oneMole, oneAtom, length);
				for(int i = 0; i < allPaths.size(); i++){
					IAtom checkAtom = allPaths.get(i).get(allPaths.get(i).size() - 1);
					if(!closestAtoms.contains(checkAtom)) closestAtoms.add(checkAtom);
				}				
				length++;
			}
		}
		return closestAtoms;
	}
	/**
	 * This function will find the closest atoms list based containing "numAtoms" many atoms on the given atomDisMap
	 * @param atomDisMap
	 * @return
	 */
	public ArrayList<IAtom> getCloesestAtomListInMap(HashMap<IAtom, Double> atomDisMap, Integer numAtoms){
		ArrayList<IAtom> resultList = new ArrayList<>();
		Integer counter = 0;
		while(!atomDisMap.keySet().isEmpty() && counter < numAtoms){
			IAtom oneAtom = getClosestAtom(atomDisMap);
			resultList.add(oneAtom);
			atomDisMap.remove(oneAtom);			
			counter++;
		}
		
		return resultList;
		
	}
	/**
	 * This function will find the closest atom in the HashMap<IAtom, Double> atomDisMap
	 * @param atomDisMap
	 * @return
	 */
	public IAtom getClosestAtom(HashMap<IAtom, Double> atomDisMap){
		Double minDis = Double.MAX_VALUE;
		IAtom closestAtom = null;
		for(IAtom oneAtom : atomDisMap.keySet()){
			Double dis = atomDisMap.get(oneAtom);
			if(dis < minDis){
				dis = minDis;
				closestAtom = oneAtom;
			}
		}
		return closestAtom;
	}
	/**
	 * This function will caculate the distance between two atoms by their 3D coordinates
	 * @param atom_one
	 * @param atom_two
	 * @return
	 */
	public Double calculateDistance(IAtom atom_one, IAtom atom_two){
		Point3d atomOneCoor = atom_one.getPoint3d();
		Point3d atomTwoCoor = atom_two.getPoint3d();
		Double dis = Math.sqrt((atomOneCoor.x - atomTwoCoor.x)*(atomOneCoor.x - atomTwoCoor.x) + 
							   (atomOneCoor.y - atomTwoCoor.y)*(atomOneCoor.y - atomTwoCoor.y) +
							   ((atomOneCoor.z - atomTwoCoor.z)*(atomOneCoor.z - atomTwoCoor.z)));
		return dis;
		
	}

}
