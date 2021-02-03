package cyProduct;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smsd.ring.RingFinder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import utils.MergeFragments;

public class ClarifyReactionType {
	/**
	 * 1. Check Ring rearrangement reaction first.
	 * 2. NitrosoReduction
	 * 3. Special Dealkylation reaction
	 * 4. Dealkylation
	 * 5. EpOxidation
	 * 6. Carbon-Oxygen Oxidation
	 * 7.  <i, H> Oxidation
	 * 8. <i, SPN> Oxidation
	 * 9. Dehydrogenation reaciton
	 * 10. Special Reaction one: X=C-C -> X-C(=O)-C
	 * 11. Special Reaction two: N=C-N -> N-C(=O)-N
	 * 12. Desulfurization
	 * 13. Hydrolysis
	 * 14. Dehalogon
	 * 15. Reduction
	 * @param bomList_typeOne
	 * @param bomList_typeTwo
	 * @param bomList_typeThree
	 * @return
	 * @throws Exception
	 */
	
	public static IAtomContainerSet arrangeReactionTypesAndPredictMetabolites(ArrayList<ArrayList<IAtom>> bomList_typeOne, ArrayList<IAtom> bomList_typeTwo, ArrayList<IAtom> bomList_typeThree, IAtomContainer oneMole) throws Exception{
		//<i,j> bond changed type of Reactions
		IAtomContainerSet metabolites_all = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		//1. Ring rearrangement
		ArrayList<ArrayList<IAtom>> specialTypeOneReactionsBoMs = findRingRearrangementReactionType(bomList_typeOne, oneMole);
		if(specialTypeOneReactionsBoMs != null){
			try{
				IAtomContainerSet metabolites = GenerateMetabolite.generateRingRearrangementMetabolite(specialTypeOneReactionsBoMs, oneMole);
				metabolites_all.add(metabolites);
				//See nisoldipine and Acetaminophen
				IAtom deCarbon;
				for(int t = 0; t < specialTypeOneReactionsBoMs.size(); t++){
					ArrayList<IAtom> oneBoM = specialTypeOneReactionsBoMs.get(t);
					for(int k = 0; k < oneBoM.size(); k++){
						deCarbon = oneBoM.get(k);
						for(int j = 0; j < bomList_typeTwo.size(); j++){
							if(UniqueIDFunctionSet.getUniqID(bomList_typeTwo.get(j)) == UniqueIDFunctionSet.getUniqID(deCarbon)){
								boolean toRemove = true;
								for(int x = 0; x < specialTypeOneReactionsBoMs.size(); x++){
									if(specialTypeOneReactionsBoMs.get(x).contains(deCarbon)){
										IAtom temp_one = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(specialTypeOneReactionsBoMs.get(x).get(0)), oneMole);
										IAtom temp_two = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(specialTypeOneReactionsBoMs.get(x).get(1)), oneMole);
										if(oneMole.getBond(temp_one, temp_two).getOrder()!=IBond.Order.SINGLE) toRemove = false;
									}
								}
								if(toRemove){
									bomList_typeTwo.remove(j);
									j = 0;
								}
							}
						}
					}
				}
			}catch (Exception e){
				//Do nothing, but still remove that BoM
			}
			bomList_typeOne.removeAll(specialTypeOneReactionsBoMs);
			System.out.println("---------------------------------");
		}
		//2. NitrosoReduction
		specialTypeOneReactionsBoMs = findNitrosoReductionReactionType(bomList_typeOne);
		//Handle NOO --> NH2 reduction reaction and rearrangement reaction first. These BoMs should be removed after the reactions are applied
		//specialTypeOneReactionsBoMs is used as nitrosoReductionReactionBoMList here
		while(specialTypeOneReactionsBoMs != null && !specialTypeOneReactionsBoMs.isEmpty()){
			specialTypeOneReactionsBoMs = findNitrosoReductionReactionType(bomList_typeOne);
			if(specialTypeOneReactionsBoMs == null) break;
			try{
				IAtomContainerSet metabolites = GenerateMetabolite.generateNitrosoReductionMetabolite(specialTypeOneReactionsBoMs, oneMole);
				metabolites_all.add(metabolites);
			}catch (Exception e){
				//Do nothing, but still remove that BoM
			}
			bomList_typeOne.removeAll(specialTypeOneReactionsBoMs);
			System.out.println("---------------------------------");
		}
		//10. Special reaction one
		specialTypeOneReactionsBoMs = findSpecialReactionOneType(bomList_typeOne, oneMole);
		while(specialTypeOneReactionsBoMs != null && !specialTypeOneReactionsBoMs.isEmpty()){
			try{
				IAtomContainerSet metabolites = GenerateMetabolite.generateSpecialReactionOneMetabolite(specialTypeOneReactionsBoMs, oneMole);
				IAtom deCarbon;
				for(int t = 0; t < specialTypeOneReactionsBoMs.size(); t++){
					ArrayList<IAtom> oneBoM = specialTypeOneReactionsBoMs.get(t);
					for(int k = 0; k < oneBoM.size(); k++){
						if(oneBoM.get(k).getSymbol().equals("C")) deCarbon = oneBoM.get(k);
						else continue;
						for(int j = 0; j < bomList_typeTwo.size(); j++){
							if(UniqueIDFunctionSet.getUniqID(bomList_typeTwo.get(j)) == UniqueIDFunctionSet.getUniqID(deCarbon)){
								bomList_typeTwo.remove(j);
								j = 0;
							}
						}
					}
				}
				metabolites_all.add(metabolites);
			}catch (Exception e){
				//Do nothing, but still remove that BoM
			}			
			bomList_typeOne.removeAll(specialTypeOneReactionsBoMs);
			specialTypeOneReactionsBoMs = findSpecialReactionOneType(bomList_typeOne, oneMole);
		}
		//11.  Special reaction Two
		specialTypeOneReactionsBoMs = findSpecialReactionTwoType(bomList_typeOne, oneMole);
		while(specialTypeOneReactionsBoMs != null && !specialTypeOneReactionsBoMs.isEmpty()){
			for(int i = 0; i < specialTypeOneReactionsBoMs.size(); i++){
				ArrayList<IAtom> oneBoM =  specialTypeOneReactionsBoMs.get(i);
				IAtom atom_left = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(0)), oneMole);
				IAtom atom_right = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(1)), oneMole);
				System.out.println(oneMole.indexOf(atom_left) + ","  + oneMole.indexOf(atom_right));
			}
			try{
				IAtomContainerSet metabolites = GenerateMetabolite.generateSpecialReactionTwoMetabolite(specialTypeOneReactionsBoMs, oneMole);
				IAtom deCarbon;
				for(int t = 0; t < specialTypeOneReactionsBoMs.size(); t++){
					ArrayList<IAtom> oneBoM = specialTypeOneReactionsBoMs.get(t);
					for(int k = 0; k < oneBoM.size(); k++){
						if(oneBoM.get(k).getSymbol().equals("C")) deCarbon = oneBoM.get(k);
						else continue;
						for(int j = 0; j < bomList_typeTwo.size(); j++){
							if(UniqueIDFunctionSet.getUniqID(bomList_typeTwo.get(j)) == UniqueIDFunctionSet.getUniqID(deCarbon)){
								bomList_typeTwo.remove(j);
								j = 0;
							}
						}
					}
				}
				metabolites_all.add(metabolites);
			}catch (Exception e){
				//Do nothing, but still remove that BoM
			}
			bomList_typeOne.removeAll(specialTypeOneReactionsBoMs);
			specialTypeOneReactionsBoMs = findSpecialReactionTwoType(bomList_typeOne, oneMole);
		}
		//specialTypeOneReactionsBoMs is used as 
		//Special Dealkylation reaction
		specialTypeOneReactionsBoMs = findOCODealkylationReactionType(bomList_typeOne, bomList_typeTwo, oneMole);
		while(specialTypeOneReactionsBoMs != null && !specialTypeOneReactionsBoMs.isEmpty()){			
			if(specialTypeOneReactionsBoMs == null) break;
			try{
				IAtomContainerSet metabolites = GenerateMetabolite.generateOCODealkylationMetabolite(specialTypeOneReactionsBoMs, oneMole);
				//Remove the corresponding typeTwo BoM
				IAtom deCarbon;
				for(int t = 0; t < specialTypeOneReactionsBoMs.size(); t++){
					ArrayList<IAtom> oneBoM = specialTypeOneReactionsBoMs.get(t);
					for(int k = 0; k < oneBoM.size(); k++){
						if(oneBoM.get(k).getSymbol().equals("C")) deCarbon = oneBoM.get(k);
						else continue;
						for(int j = 0; j < bomList_typeTwo.size(); j++){
							if(UniqueIDFunctionSet.getUniqID(bomList_typeTwo.get(j)) == UniqueIDFunctionSet.getUniqID(deCarbon)){
								bomList_typeTwo.remove(j);
								j = 0;
							}
						}
					}
				}
				metabolites_all.add(metabolites);
			}catch (Exception e){
				//Do nothing, but still remove that BoM
			}	
			bomList_typeOne.removeAll(specialTypeOneReactionsBoMs);
			System.out.println("---------------------------------");
			specialTypeOneReactionsBoMs = findOCODealkylationReactionType(bomList_typeOne, bomList_typeTwo, oneMole);
			
		}
		//12. Desulfurization
		specialTypeOneReactionsBoMs = findDesulfurizationReactionType(bomList_typeOne, oneMole);
		while(specialTypeOneReactionsBoMs != null && !specialTypeOneReactionsBoMs.isEmpty()){
			try{
				IAtomContainerSet metabolites = GenerateMetabolite.generateDesulfurizationMetabolite(specialTypeOneReactionsBoMs, oneMole);
				metabolites_all.add(metabolites);
			}catch (Exception e){
				//Do nothing, but still remove that BoM
			}	
			bomList_typeOne.removeAll(specialTypeOneReactionsBoMs);
			System.out.println("---------------------------------");
			specialTypeOneReactionsBoMs = findDesulfurizationReactionType(bomList_typeOne, oneMole);
		}
		//13. Hydrolysis
		specialTypeOneReactionsBoMs = findHydrolysisReactionType(bomList_typeOne, oneMole);
		while(specialTypeOneReactionsBoMs != null && !specialTypeOneReactionsBoMs.isEmpty()){
			try{
				IAtomContainerSet metabolites = GenerateMetabolite.generateHydrolysisMetabolite(specialTypeOneReactionsBoMs, oneMole);
				metabolites_all.add(metabolites);						
			}catch (Exception e){
				//Do nothing, but still remove that BoM
			}	
			bomList_typeOne.removeAll(specialTypeOneReactionsBoMs);
			System.out.println("---------------------------------");
			specialTypeOneReactionsBoMs = findHydrolysisReactionType(bomList_typeOne, oneMole);
		}
		//14. Dehalogenation
		specialTypeOneReactionsBoMs = findDehalogenReactionType(bomList_typeOne, oneMole);
		while(specialTypeOneReactionsBoMs != null && !specialTypeOneReactionsBoMs.isEmpty()){
			try{
				IAtomContainerSet metabolites = GenerateMetabolite.generateDehalogenMetabolite(specialTypeOneReactionsBoMs, oneMole);
				IAtom deCarbon;
				for(int t = 0; t < specialTypeOneReactionsBoMs.size(); t++){
					ArrayList<IAtom> oneBoM = specialTypeOneReactionsBoMs.get(t);
					for(int k = 0; k < oneBoM.size(); k++){
						if(oneBoM.get(k).getSymbol().equals("C")) deCarbon = oneBoM.get(k);
						else continue;
						for(int j = 0; j < bomList_typeTwo.size(); j++){
							if(UniqueIDFunctionSet.getUniqID(bomList_typeTwo.get(j)) == UniqueIDFunctionSet.getUniqID(deCarbon)){
								bomList_typeTwo.remove(j);
								j = 0;
							}
						}
					}
				}
			metabolites_all.add(metabolites);
			}catch (Exception e){
				//Do nothing, but still remove that BoM
			}	
			bomList_typeOne.removeAll(specialTypeOneReactionsBoMs);
			System.out.println("---------------------------------");
			specialTypeOneReactionsBoMs = findDehalogenReactionType(bomList_typeOne, oneMole);
		}
		//15. Reduction
		specialTypeOneReactionsBoMs = findReductionReactionType(bomList_typeOne, oneMole);
		while(specialTypeOneReactionsBoMs != null && !specialTypeOneReactionsBoMs.isEmpty()){
			try{
				IAtomContainerSet metabolites = GenerateMetabolite.generateReductionMetabolite(specialTypeOneReactionsBoMs, oneMole);
				metabolites_all.add(metabolites);
				}catch (Exception e){
				//Do nothing, but still remove that BoM
			}	
			bomList_typeOne.removeAll(specialTypeOneReactionsBoMs);
			System.out.println("---------------------------------");
			specialTypeOneReactionsBoMs = findReductionReactionType(bomList_typeOne, oneMole);
		}
		
		//Dealkylation
		for(int i = 0; i < bomList_typeOne.size(); i++){
			ArrayList<IAtom> bom_typeOne = bomList_typeOne.get(i);
			//Carbon-Oxygen Oxidation
			if(bomList_typeOne.isEmpty() || bomList_typeOne == null) break;
			bom_typeOne = bomList_typeOne.get(i);
			String carbonOxygenOxidation = findCarbonOxygenOxidationReaction(bom_typeOne, bomList_typeTwo, oneMole);
			while(carbonOxygenOxidation != null && !bomList_typeOne.isEmpty()){
				try{
					IAtomContainerSet metabolites = GenerateMetabolite.generateCarbonOxygenOxidationMetabolite(bom_typeOne.get(0), bom_typeOne.get(1), oneMole);
					metabolites_all.add(metabolites);
				}catch (Exception e){
					//Do nothing, but still remove that BoM
				}	
				bomList_typeOne.remove(bom_typeOne);
				i = 0;
				System.out.println("---------------------------------");
				if(bomList_typeOne.isEmpty()) break;
				bom_typeOne = bomList_typeOne.get(i);
				carbonOxygenOxidation = findCarbonOxygenOxidationReaction(bom_typeOne, bomList_typeTwo, oneMole);
				//continue;
			}
			//Dealkylation
			if(bomList_typeOne.isEmpty() || bomList_typeOne == null) break;
			bom_typeOne = bomList_typeOne.get(i);
			String dealkylation = findDealkyaltionReactionType(bom_typeOne, bomList_typeTwo, oneMole);
			while(dealkylation != null && !bomList_typeOne.isEmpty()){		
				try{
					IAtomContainerSet metabolites = GenerateMetabolite.generateDeAlkylationMetabolite(bom_typeOne.get(0), bom_typeOne.get(1), oneMole, dealkylation);	
					//Remove the corresponding typeTwo BoM
					IAtom deCarbon;
					if(bom_typeOne.get(0).getSymbol().equals("C")) deCarbon = bom_typeOne.get(0);
					else deCarbon = bom_typeOne.get(1);
					for(int j = 0; j < bomList_typeTwo.size(); j++){
						if(UniqueIDFunctionSet.getUniqID(bomList_typeTwo.get(j)) == UniqueIDFunctionSet.getUniqID(deCarbon)){
							bomList_typeTwo.remove(j);
							break;
						}
					}
					metabolites_all.add(metabolites);
				}catch (Exception e){
					//Do nothing, but still remove that BoM
				}
				bomList_typeOne.remove(bom_typeOne);
				i = 0;
				System.out.println("---------------------------------");
				if(bomList_typeOne.isEmpty()) break;
				bom_typeOne = bomList_typeOne.get(i);
				dealkylation = findDealkyaltionReactionType(bom_typeOne, bomList_typeTwo, oneMole);
				//continue;
			}
			//Dehydrogenation
			boolean shrink = false;
			if(bomList_typeOne.isEmpty() || bomList_typeOne == null) break;
			bom_typeOne = bomList_typeOne.get(i);
			String dehydrogenation = findDehydrogenationReaction(bom_typeOne, bomList_typeTwo, oneMole);
			if(dehydrogenation != null){				
				try{
					IAtomContainerSet metabolites = GenerateMetabolite.generateDehydrogenationMetabolite(bom_typeOne.get(0), bom_typeOne.get(1), oneMole, 0);
					metabolites_all.add(metabolites);
				}catch (Exception e){
					//Do nothing, but still remove that BoM
				}	
				bomList_typeOne.remove(bom_typeOne);
				i = 0;
				shrink = true;
				System.out.println("---------------------------------");
				if(bomList_typeOne.isEmpty() || bomList_typeOne == null) break;
				//continue;
			}
			//EpOxidation
			if(bomList_typeOne.isEmpty() || bomList_typeOne == null) break;
			bom_typeOne = bomList_typeOne.get(i);
//			System.out.println(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(0)) + ";" + UniqueIDFunctionSet.getUniqID(bom_typeOne.get(1)));
			String epOxidation = findEpOxidationReaction(bom_typeOne, bomList_typeTwo, oneMole);
			if(epOxidation != null){			
				try{
					IAtomContainerSet metabolites = GenerateMetabolite.generateEpOxidationMetabolite(bom_typeOne.get(0), bom_typeOne.get(1), oneMole, 0);
					metabolites_all.add(metabolites);
				}catch (Exception e){
					//Do nothing, but still remove that BoM
				}	
				bomList_typeOne.remove(bom_typeOne);
				//i = i-1;
				shrink = true;
				System.out.println("---------------------------------");
				if(bomList_typeOne.isEmpty() || bomList_typeOne == null) break;
				//continue;
			}
			if(shrink) i = i - 1;

		}
		//<i,H> bond changed type of Reactions. Here is Hydroxylation reaction
		for(int i = 0; i < bomList_typeTwo.size(); i++){
			IAtom bom_typeTwo = bomList_typeTwo.get(i);
			String hydroxylation = findHydroxylationReaction(bomList_typeOne, bom_typeTwo);
			if(hydroxylation != null){
				try{
					IAtomContainerSet metabolites = GenerateMetabolite.generateHydroxylationMetabolite(bom_typeTwo, oneMole);
					metabolites_all.add(metabolites);
				}catch (Exception e){
					//Do nothing, but still remove that BoM
				}	
				continue;
			}
		}
		//<i,SPN> bond changed type of reactions. Here is S,N,P Oxidation reactions
		for(int i = 0; i < bomList_typeThree.size(); i++){
			IAtom bom_typeThree = bomList_typeThree.get(i);
			String snp_oxidation = findSNPOxidationReaction(bom_typeThree);
			if(snp_oxidation != null){
				try{
					IAtomContainerSet metabolites = GenerateMetabolite.generateSPNOxidationMetabolite(bom_typeThree, oneMole);
					metabolites_all.add(metabolites);
				}catch (Exception e){
					//Do nothing, but still remove that BoM
				}	
			}
		}
		return metabolites_all;
	}
	/**
	 * If <i,H> is typeTwo boM and i is not part of any <i,j> TypeOne BoM
	 * return null if <i,H> is not Hydroxylation
	 * @param bom_typeOne
	 * @param bomList_typeTwo
	 * @return
	 * @throws Exception
	 */
	public static String findHydroxylationReaction(ArrayList<ArrayList<IAtom>> bomList_typeOne, IAtom bom_typeTwo) throws Exception{
		boolean isValid = true;
		for(int i = 0; i < bomList_typeOne.size(); i++){
			if(bomList_typeOne.get(i).contains(bom_typeTwo)) isValid = false;
		}
		if(isValid) return "Hydroxylation";
		else return null;
	}
	/**
	 * <C.i,O.j> is C-O typeOne BoM and C.i is not a TypeTwo BoM
	 * @param bom_typeOne
	 * @param bomList_typeTwo
	 * @return
	 * @throws Excpetion
	 */
	public static String findCarbonOxygenOxidationReaction(ArrayList<IAtom> bom_typeOne, ArrayList<IAtom> bomList_typeTwo, IAtomContainer oneMole) throws Exception{
		boolean valid = false;
		//If must be a Carbon-Oxygen "SINGLE" bond. The Oxygen must be C-O-H and not connect to other non-Hydrogen atom 
		if(bom_typeOne.get(0).getSymbol().equalsIgnoreCase("O") && bom_typeOne.get(1).getSymbol().equalsIgnoreCase("C")){
			IAtom corsAtom = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(0)), oneMole);
			if(corsAtom.getImplicitHydrogenCount() == 1){
				valid = true;
			}
			List<IAtom> neighborAtomList = oneMole.getConnectedAtomsList(corsAtom);
			int counter = 0;
			for(int i = 0; i < neighborAtomList.size(); i++){
				if(!neighborAtomList.get(i).getSymbol().equalsIgnoreCase("H")) counter++;
			}
			if(counter == 1) valid = true;
		}
		else if(bom_typeOne.get(1).getSymbol().equalsIgnoreCase("O") && bom_typeOne.get(0).getSymbol().equalsIgnoreCase("C")){
			IAtom corsAtom = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(1)), oneMole);
			if(corsAtom.getImplicitHydrogenCount() == 1){
				valid = true;
			}
			List<IAtom> neighborAtomList = oneMole.getConnectedAtomsList(corsAtom);
			int counter = 0;
			for(int i = 0; i < neighborAtomList.size(); i++){
				if(!neighborAtomList.get(i).getSymbol().equalsIgnoreCase("H")) counter++;
			}
			if(counter == 1) valid = true;
		}
		//If the BoM is not a carbon-oxygen bond, then it cannot be a carbonOxygen Oxidation reaction.
		if(!valid) return null;
		IAtom atom_one = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(0)), oneMole);
		IAtom atom_two = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(1)), oneMole);
		if(oneMole.getBond(atom_one, atom_two).getOrder() == IBond.Order.SINGLE) valid = true;
		//The it doesn't satisfy the requirement of Carbon-Oxygen reaction

		for(int i = 0; i < bomList_typeTwo.size(); i++){
			if(bom_typeOne.contains(bomList_typeTwo.get(i))){
				valid = false;
				break;
			}
		}
		if(!valid){
			System.out.println("The input BoM doesn't satify the Carbon-Oxygen reaction");
			return null;
		}
		else return "CarbonOxygenOxidation";
	}
	/**
	 * This function is trivial but add it for consistency.
	 * @param bom_typeThree
	 * @return
	 * @throws Exception
	 */
	public static String findSNPOxidationReaction(IAtom bom_typeThree) throws Exception{
		if(!bom_typeThree.getSymbol().equalsIgnoreCase("S") &&!bom_typeThree.getSymbol().equalsIgnoreCase("N") && !bom_typeThree.getSymbol().equalsIgnoreCase("P")){
			System.out.println("The input typeThree BoM doesn't satisfy the SNP-Oxidation");
			return null;
		}
		else return "SNP-Oxidation";
	}
	/**
	 * If <i,j> is typeOne BoM && bond <i,j> is DOUBLE bond && neither i or j are type Two BoM
	 * return null if the <i,j> is not EpOxidation
	 * @param bom_typeOne
	 * @param bomList_typeTwo
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static String findEpOxidationReaction(ArrayList<IAtom> bom_typeOne, ArrayList<IAtom> bomList_typeTwo, IAtomContainer oneMole) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(oneMole);
		boolean isEpOxidation = true;
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(0)), oneMole);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(1)), oneMole);
		//if()
		if(oneMole.getBond(leftCopy, rightCopy).getOrder() != IBond.Order.DOUBLE && 
				!leftCopy.isAromatic() && !rightCopy.isAromatic()){
			isEpOxidation = false;
		}
		for(int i = 0; i < bomList_typeTwo.size(); i++){
			if(bom_typeOne.contains(bomList_typeTwo.get(i))){
				isEpOxidation = false;
			}
		}
		if(leftCopy.getSymbol().equalsIgnoreCase("O") || rightCopy.getSymbol().equalsIgnoreCase("O")){
			isEpOxidation = false;
		}
		if(isEpOxidation) return "EpOxidation";
		else return null;
	}
	/**
	 * If <i,j> is typeOne BoM and i is typeTwo BoM, then it is a dealkylation reaction 
	 * Return "null" if the reaction is not Dealkylation reaction 
	 * @param bomList_typeOne
	 * @param bomList_typeTwo
	 * @return
	 * @throws Exception
	 */
	public  static String findDealkyaltionReactionType(ArrayList<IAtom> bom_typeOne, ArrayList<IAtom> bomList_typeTwo, IAtomContainer oneMole) throws Exception{
		if(bom_typeOne.size() != 2) throw new Exception("The number of atoms doesn't equal to 2. It must be 2 to match TypeOne BoM");
		IAtom atomToCheck = null;
		if(bom_typeOne.get(0).getSymbol().equalsIgnoreCase("C") && bom_typeOne.get(1).getSymbol().equalsIgnoreCase("C")){
			System.out.println("BoM: " + "<" + 
							    UniqueIDFunctionSet.getUniqID(bom_typeOne.get(0)) + ";" + UniqueIDFunctionSet.getUniqID(bom_typeOne.get(1)) + ">" + 
								"is not Dealkyaltion reaction");
			return null;
		}
		IAtom cors_Atom_One = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(0)), oneMole);
		IAtom cors_Atom_Two = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(1)), oneMole);
		if(oneMole.getBond(cors_Atom_One, cors_Atom_Two).getOrder() != IBond.Order.SINGLE){
			return null;
		}
		else if(!bom_typeOne.get(0).getSymbol().equalsIgnoreCase("C")) atomToCheck = bom_typeOne.get(0);
		else atomToCheck = bom_typeOne.get(1);
		//Because for example, when a C-O bond is predicted as BoM, even though C-H is not predicted as typeTwo BoM, we treat this as Dealkylation reaction
		//This means that the validation check process here is not useful at all.
		//In addition, the hydrolysis reaction is also included by doing this
		boolean valid = true;		
		for(int i = 0; i < bomList_typeTwo.size(); i++){
			Integer uniqueID = UniqueIDFunctionSet.getUniqID(bomList_typeTwo.get(i));
			if(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(0)) == uniqueID ||
			   UniqueIDFunctionSet.getUniqID(bom_typeOne.get(1)) == uniqueID){
				if(!bomList_typeTwo.get(i).getSymbol().equalsIgnoreCase("N")){
					valid = true;
					break;
				}
			}
		}
		if(!valid) return null;
		if(atomToCheck.getSymbol().equalsIgnoreCase("O")){
			return "O-Dealkylation";
		}
		else if(atomToCheck.getSymbol().equalsIgnoreCase("N")){
			return "N-Dealkylation";
		}
		else if(atomToCheck.getSymbol().equalsIgnoreCase("S")){
			return "S-Dealkylation";
		}
		else throw new Exception("Can not find the Dealkylation type for Atom with Symbol: " + atomToCheck.getSymbol());
	}
	/**
	 * This function will explore all typeOne BoMs and identify if there are NOO --> NH2 reduction reactions.
	 * It returns NitrosoReductionReaction if there is at least one such reaction identified
	 * returns null if there is none.
	 * @param typeOne_BoM_list
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<ArrayList<IAtom>> findNitrosoReductionReactionType(ArrayList<ArrayList<IAtom>> typeOne_BoM_list) throws Exception{
		ArrayList<ArrayList<IAtom>> candidteList = new ArrayList<>();
		ArrayList<ArrayList<IAtom>> resultList = new ArrayList<>();
		for(int i = 0; i < typeOne_BoM_list.size(); i++){
			ArrayList<IAtom> oneTypeOneBoM = typeOne_BoM_list.get(i);
			if((oneTypeOneBoM.get(0).getSymbol().equalsIgnoreCase("N") && oneTypeOneBoM.get(1).getSymbol().equalsIgnoreCase("O")) ||
				oneTypeOneBoM.get(1).getSymbol().equalsIgnoreCase("N") && oneTypeOneBoM.get(0).getSymbol().equalsIgnoreCase("O")){
				candidteList.add(oneTypeOneBoM);
			}
		}
		for(int i = 0; i < candidteList.size(); i++){
			ArrayList<IAtom> oneTypeOneBoM = candidteList.get(i);
			IAtom atom_one;// oneTypeOneBoM.get(0);
			
			if(oneTypeOneBoM.get(0).getSymbol().equalsIgnoreCase("N")) atom_one = oneTypeOneBoM.get(0);
			else atom_one = oneTypeOneBoM.get(1);			
			//IAtom atom_two = oneTypeOneBoM.get(1);
			for(int j = i+1; j < candidteList.size(); j++){
				IAtom atom_two;
				ArrayList<IAtom> nextTypeOneBoM = candidteList.get(j);				
				if(nextTypeOneBoM.get(0).getSymbol().equalsIgnoreCase("N")) atom_two = nextTypeOneBoM.get(0);
				else atom_two = nextTypeOneBoM.get(1);
				if(UniqueIDFunctionSet.getUniqID(atom_one) == UniqueIDFunctionSet.getUniqID(atom_two)){
					resultList.add(oneTypeOneBoM);
					resultList.add(nextTypeOneBoM);
					return resultList;
				}
			}
			
		}
		return null;
	}
	/**
	 * This function will find the R-O-C-O-R within the ring patter of O-Dealkyaltion. In this reaction, the ring will open and the carbon will be removed as aldehyde
	 * It will return the ArrayList of the two TypeOne BoMs and the one TypeTwo BoM
	 * @param typeOne_BoM_List
	 * @param typeTwo_BoM_List
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<ArrayList<IAtom>> findOCODealkylationReactionType(ArrayList<ArrayList<IAtom>> typeOne_BoM_List, ArrayList<IAtom> typeTwo_BoM_List, IAtomContainer oneMole) throws Exception{
		for(int i = 0; i < typeTwo_BoM_List.size(); i++){
			ArrayList<ArrayList<IAtom>> oneCandidateResult = new ArrayList<>();
			IAtom oneTypeTwo_BoM = typeTwo_BoM_List.get(i);
			Integer uniqueID = UniqueIDFunctionSet.getUniqID(oneTypeTwo_BoM);
			for(int j = 0; j < typeOne_BoM_List.size(); j++){
				ArrayList<IAtom> one_typeOne_BoM = typeOne_BoM_List.get(j);
				if(UniqueIDFunctionSet.getUniqID(one_typeOne_BoM.get(0)) == uniqueID || 
				   UniqueIDFunctionSet.getUniqID(one_typeOne_BoM.get(1)) == uniqueID		){
					if(one_typeOne_BoM.get(0).getSymbol().equalsIgnoreCase("O") || one_typeOne_BoM.get(1).getSymbol().equalsIgnoreCase("O")){
						oneCandidateResult.add(one_typeOne_BoM);
					}
				}
			}
			if(oneCandidateResult.size() == 2){
				IAtomContainer copy = oneMole.clone();
				IAtom atom_copy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneTypeTwo_BoM), copy);
				copy.removeAtom(atom_copy);
				IAtomContainerSet leftSet = ConnectivityChecker.partitionIntoMolecules(copy);
				//After the shared carbon is removed, if the molecule is still connected as one, then the R-O-C-O-R pattern is within a circle
				if(leftSet.getAtomContainerCount() == 1){
					return oneCandidateResult;
				}
			}
		}
		//If no matched TypeOne BoMs can be found, then there is no such special EpOxidation Pattern.
		return null;

	}
	/**
	 * This function will find the S-C=C within the ring. In this reaction, the ring will open and the carbon will be removed as aldehyde
	 * It will return the ArrayList of the two TypeOne BoMs and the one TypeTwo BoM
	 * @param typeOne_BoM_List
	 * @param typeTwo_BoM_List
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<ArrayList<IAtom>> findSpecialReactionOneType(ArrayList<ArrayList<IAtom>> typeOne_BoM_List, IAtomContainer oneMole) throws Exception{
		ArrayList<ArrayList<IAtom>> result = new ArrayList<>();
		IAtomContainer copy = oneMole.clone();
		String smirks = "[#6]=,:1[#6]=,:[#6][#16][#6]=,:1";//"[#16]-1-[#6]=[#6]-[#6]=[#6]-1";
		ArrayList<ArrayList<IAtom>> matchedList = MergeFragments.findMatchedList(smirks, copy);
		for(int t = 0; t< matchedList.size(); t++){
			ArrayList<IAtom> matchedAtoms = matchedList.get(t);
			ArrayList<Integer> matchedUIDList = new ArrayList<>();
			//ArrayList<IAtom> targetAtom = new ArrayList<>();
			for(int i = 0; i < matchedAtoms.size(); i++){
				IAtom oneAtom = matchedAtoms.get(i);
				Integer oneIdx = UniqueIDFunctionSet.getUniqID(oneAtom);
				if(!matchedUIDList.contains(oneIdx))matchedUIDList.add(oneIdx);			
			}
			//If 3 type one bonds are within the matched fragment, then it's the specifalReactionOne
			int count = 0;
			for(int i = 0; i < typeOne_BoM_List.size(); i++){
				ArrayList<IAtom> oneBoM = typeOne_BoM_List.get(i);
				Integer idx_left = UniqueIDFunctionSet.getUniqID(oneBoM.get(0));
				Integer idx_right = UniqueIDFunctionSet.getUniqID(oneBoM.get(1));
				//If two atoms of the bond is within the the matched fragment
				if(matchedUIDList.contains(idx_left) && matchedUIDList.contains(idx_right)){
					result.add(oneBoM);
					count++;
				}			
			}
			if(count == 3){
				System.out.println("Special reaction one  is classifed: S-C=C => S-C(=O)-S (part of the reaction) and in  a ring");
				return result;
			}
		}
		return null;
		//If no matched TypeOne BoMs can be found, then there is no such special EpOxidation Pattern.

	}
	/**
	 * This function will find the N=C-N -> N-C(=O)-N within the 5-member ring. In this reaction, the ring will open and the carbon will be removed as aldehyde
	 * It will return the ArrayList of the two TypeOne BoMs and the one TypeTwo BoM
	 * @param typeOne_BoM_List
	 * @param typeTwo_BoM_List
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<ArrayList<IAtom>> findSpecialReactionTwoType(ArrayList<ArrayList<IAtom>> typeOne_BoM_List, IAtomContainer oneMole) throws Exception{
		ArrayList<ArrayList<IAtom>> result = new ArrayList<>();
		IAtomContainer copy = oneMole.clone();
		String smirks = "[#6]1~[#6][#6,#7]=,:[#6][#7]1";
		ArrayList<ArrayList<IAtom>> matchedList = MergeFragments.findMatchedList(smirks, copy);
		for(int t = 0; t< matchedList.size(); t++){
			ArrayList<IAtom> matchedAtoms = matchedList.get(t);
			ArrayList<Integer> matchedUIDList = new ArrayList<>();
			//ArrayList<IAtom> targetAtom = new ArrayList<>();
			for(int i = 0; i < matchedAtoms.size(); i++){
				IAtom oneAtom = matchedAtoms.get(i);
				Integer oneIdx = UniqueIDFunctionSet.getUniqID(oneAtom);
				if(!matchedUIDList.contains(oneIdx))matchedUIDList.add(oneIdx);			
			}
			//If 3 type one bonds are within the matched fragment, then it's the specifalReactionOne
			int count = 0;
			for(int i = 0; i < typeOne_BoM_List.size(); i++){
				ArrayList<IAtom> oneBoM = typeOne_BoM_List.get(i);
				Integer idx_left = UniqueIDFunctionSet.getUniqID(oneBoM.get(0));
				Integer idx_right = UniqueIDFunctionSet.getUniqID(oneBoM.get(1));
				//If two atoms of the bond is within the the matched fragment
				if(matchedUIDList.contains(idx_left) && matchedUIDList.contains(idx_right)){
					result.add(oneBoM);
					count++;
				}			
			}
			System.out.println("Special reaction one  is classifed: S-C=C => S-C(=O)-S (part of the reaction) and in  a ring");
			return result;
			
		}
		return null;
		//If no matched TypeOne BoMs can be found, then there is no such special EpOxidation Pattern.

	}
	
	/**
	 * This function will find the Desulfurization reaction by checking typeOne BoMs.
	 * If the BoM is X=S, then it will be oxidiazed into X=O
	 */
	public static ArrayList<ArrayList<IAtom>> findDesulfurizationReactionType(ArrayList<ArrayList<IAtom>> typeOne_BoM_List,IAtomContainer oneMole) throws Exception{
		ArrayList<ArrayList<IAtom>> result = new ArrayList<>();
		for(int i = 0; i < typeOne_BoM_List.size(); i++){
			ArrayList<IAtom> oneBoM = typeOne_BoM_List.get(i);
			IAtom atom_left = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(0)), oneMole);
			IAtom atom_right = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(1)), oneMole);
			if(atom_left.getSymbol().equalsIgnoreCase("S") || atom_right.getSymbol().equalsIgnoreCase("S")){
				IBond oneBond = oneMole.getBond(atom_left, atom_right);
				if(oneBond.getOrder() == IBond.Order.DOUBLE){
					result.add(oneBoM);
					return result;
				}
				
			}
		}
		//If no matched TypeOne BoMs can be found, then there is no such special EpOxidation Pattern.
		return null;
	}	
	/**
	 * This function will find the dehalogen reaction type. For example, Cl-C => HCl, C=O
	 */
	public static ArrayList<ArrayList<IAtom>> findDehalogenReactionType(ArrayList<ArrayList<IAtom>> typeOne_BoM_List,IAtomContainer oneMole) throws Exception{
		ArrayList<ArrayList<IAtom>> result = new ArrayList<>();
		for(int i = 0; i < typeOne_BoM_List.size(); i++){
			ArrayList<IAtom> oneBoM = typeOne_BoM_List.get(i);
			IAtom atom_left = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(0)), oneMole);
			IAtom atom_right = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(1)), oneMole);
			if(atom_left.getSymbol().equalsIgnoreCase("Cl") || atom_left.getSymbol().equalsIgnoreCase("F") ||
			   atom_right.getSymbol().equalsIgnoreCase("Cl") || atom_right.getSymbol().equalsIgnoreCase("F")){
				IBond oneBond = oneMole.getBond(atom_left, atom_right);
				if(oneBond.getOrder() == IBond.Order.SINGLE){
					result.add(oneBoM);
					return result;
				}				
			}
		}
		//If no matched TypeOne BoMs can be found, then there is no such special EpOxidation Pattern.
		return null;
	} 
	/**
	 * This function will find the hydrolysis reaction type.
	 * Note that O=C-[O,S,P]-R => O=C-OH + HO-R type hydrolysis is covered by Dealkylation.
	 * We focus on R1-[O,P,S]-[O,P,S]-R2 => R1-[O,P,S]-[O,P,S]H , HO-[O,P,S]-R2 
	 * @param typeOne_BoM_List
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<ArrayList<IAtom>> findHydrolysisReactionType(ArrayList<ArrayList<IAtom>> typeOne_BoM_List, IAtomContainer oneMole) throws Exception{
		ArrayList<ArrayList<IAtom>> result  = new ArrayList<>();
		for(int i = 0; i < typeOne_BoM_List.size(); i++){
			ArrayList<IAtom> oneBoM = typeOne_BoM_List.get(i);
			IAtom atom_left = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(0)), oneMole);
			IAtom atom_right = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(1)), oneMole);
			if((atom_left.getSymbol().equalsIgnoreCase("S") || atom_left.getSymbol().equalsIgnoreCase("P") || atom_left.getSymbol().equalsIgnoreCase("O")) &&
			   (atom_right.getSymbol().equalsIgnoreCase("S") || atom_right.getSymbol().equalsIgnoreCase("P") ||  atom_right.getSymbol().equalsIgnoreCase("O"))){
				result.add(oneBoM);
				return result;
			}
		}
		return null;
	}
	/**
	 * This function will find the X=O reduction reaction by checking typeOne BoMs.
	 * If the BoM is X=S, then it will be oxidiazed into X=O
	 */
	public static ArrayList<ArrayList<IAtom>> findReductionReactionType(ArrayList<ArrayList<IAtom>> typeOne_BoM_List,IAtomContainer oneMole) throws Exception{
		ArrayList<ArrayList<IAtom>> result = new ArrayList<>();
		for(int i = 0; i < typeOne_BoM_List.size(); i++){
			ArrayList<IAtom> oneBoM = typeOne_BoM_List.get(i);
			IAtom atom_left = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(0)), oneMole);
			IAtom atom_right = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(oneBoM.get(1)), oneMole);
			if(atom_left.getSymbol().equalsIgnoreCase("O") || atom_right.getSymbol().equalsIgnoreCase("O")){
				IBond oneBond = oneMole.getBond(atom_left, atom_right);
				if(oneBond.getOrder() == IBond.Order.DOUBLE){
					result.add(oneBoM);
					return result;
				}
				
			}
		}
		//If no matched TypeOne BoMs can be found, then there is no such special EpOxidation Pattern.
		return null;
	}
	/**
	 * This function will find the ring rearrangement reaction by checking all the TypeOne BoMs.
	 *  
	 * Condition: 
	 * 1. There are at least three TypeOne BoMs that are successive.
	 * 2. Most of them are within the same ring structure.
	 * Algorithm: 
	 * 1. Find all the rings within the molecule, if there is a ring whose bonds are all TypeOne BoMs, then it is a ring arrangement reaction
	 * 2. Find all atoms/Bonds that are BoMs and connected to the rearrangement ring. 
	 * Possible Issue: 
	 * There are fused rings and the bonds of the fused rings are rearranged. 
	 * In this case, both basic/small ring and the fused/large ring will satisfy this condition. This can cause problems. 
	 * However, we didn't see any such substrate-metabolite pairs collected in the training and testing datasets. 
	 * Hence we don't consider/handle that complicated case here for now. 
	 * @param typeOne_BoM_List
	 * @param typeTwo_BoM_List
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<ArrayList<IAtom>> findRingRearrangementReactionType(ArrayList<ArrayList<IAtom>> typeOne_BoM_List,  IAtomContainer oneMole) throws Exception{
		AllRingsFinder arf = new AllRingsFinder();
		IRingSet allRings = arf.findAllRings(oneMole);		
		for(int i = 0; i < allRings.getAtomContainerCount(); i++){
			IAtomContainer oneRing = allRings.getAtomContainer(i);
			ArrayList<ArrayList<IAtom>> ring_BoM_List = new ArrayList<>();
			ArrayList<ArrayList<IAtom>> bom_List = new ArrayList<>();
			boolean isRearrangement = true;
			for(int j = 0; j < oneRing.getBondCount(); j++){
				IBond oneBond = oneRing.getBond(j);
				Integer atom_one_UID = UniqueIDFunctionSet.getUniqID(oneBond.getAtom(0));
				Integer atom_two_UID = UniqueIDFunctionSet.getUniqID(oneBond.getAtom(1));
				boolean findMatchInBoMs = false;
				for(int k = 0; k < typeOne_BoM_List.size(); k++){
					ArrayList<IAtom> candidate = typeOne_BoM_List.get(k);
					Integer atom_one_UID_temp = UniqueIDFunctionSet.getUniqID(candidate.get(0));
					Integer atom_two_UID_temp = UniqueIDFunctionSet.getUniqID(candidate.get(1));
					//Check if every bond within the ring is a TypeOne BoM
					if((atom_one_UID == atom_one_UID_temp && atom_two_UID == atom_two_UID_temp) ||
					   (atom_one_UID == atom_two_UID_temp && atom_two_UID == atom_one_UID_temp)){
						findMatchInBoMs = true;
						ring_BoM_List.add(candidate);
						bom_List.add(candidate);
						break;
					}
				}
				//If the bond within the ring is not a TypeOne BoM, then this ring can't be ringRearrangement reaction
				if(!findMatchInBoMs){
					isRearrangement = false;
					break;
				}
			}
			//if the ring satisfies the condition of the ring rearrangement reaction, find all relevant BoMs including those that are not part of the ring, but connected to the ring.			
			if(isRearrangement){
				for(int j = 0; j < ring_BoM_List.size(); j++){
					ArrayList<IAtom> oneBoM = ring_BoM_List.get(j);
					Integer atom_one_UID = UniqueIDFunctionSet.getUniqID(oneBoM.get(0));
					Integer atom_two_UID = UniqueIDFunctionSet.getUniqID(oneBoM.get(1));
					for(int k = 0; k < typeOne_BoM_List.size(); k++){
						ArrayList<IAtom> candidate = typeOne_BoM_List.get(k);
						Integer atom_one_UID_temp = UniqueIDFunctionSet.getUniqID(candidate.get(0));
						Integer atom_two_UID_temp = UniqueIDFunctionSet.getUniqID(candidate.get(1));
						if((atom_one_UID == atom_one_UID_temp && atom_two_UID != atom_two_UID_temp) || 
						   (atom_one_UID == atom_two_UID_temp && atom_two_UID != atom_one_UID_temp) ||
						   (atom_two_UID == atom_two_UID_temp && atom_one_UID != atom_one_UID_temp) ||
						   (atom_two_UID == atom_one_UID_temp && atom_one_UID != atom_two_UID_temp)){
							if(!bom_List.contains(candidate)) {
								//If it's carbon attached to the nitrogen, then igonor it because it should be dealkylation reaction
								if(candidate.get(0).getSymbol().equals("N") && candidate.get(1).getSymbol().equals("C")){
									IAtom checkAtom = UniqueIDFunctionSet.getAtomByUniqueID(atom_two_UID_temp, oneMole);
									if(oneMole.getBondOrderSum(checkAtom) <= 3) continue;
								}
								else if(candidate.get(0).getSymbol().equals("C") && candidate.get(1).getSymbol().equals("N")){	
									IAtom checkAtom = UniqueIDFunctionSet.getAtomByUniqueID(atom_one_UID_temp, oneMole);
									if(oneMole.getBondOrderSum(checkAtom) <= 3) continue;
								}
								bom_List.add(candidate);
							}
						}
					}
					
				}
				//As the case we are handling now only contains one non-fused ring that has the ring rearrangement reaction, the function can stop onece it finds one.
				return bom_List;
				//break;
			}
		}
		return null;
		
	}
	/**
	 * If <i,j> is typeOne BoM && bond <i,j> is DOUBLE bond && neither i or j are type Two BoM
	 * return null if the <i,j> is not EpOxidation
	 * @param bom_typeOne
	 * @param bomList_typeTwo
	 * @param oneMole
	 * @return
	 * @throws Exception
	 */
	public static String findDehydrogenationReaction(ArrayList<IAtom> bom_typeOne, ArrayList<IAtom> bomList_typeTwo, IAtomContainer oneMole) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(oneMole);
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(oneMole);
		boolean isEpOxidation = false;
		IAtom leftCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(0)), oneMole);
		IAtom rightCopy = UniqueIDFunctionSet.getAtomByUniqueID(UniqueIDFunctionSet.getUniqID(bom_typeOne.get(1)), oneMole);
		//if()
		if(oneMole.getBond(leftCopy, rightCopy).getOrder() == IBond.Order.SINGLE && 
		   !leftCopy.isAromatic() && !rightCopy.isAromatic() &&
		   leftCopy.getSymbol().equals("C") && rightCopy.getSymbol().equals("C")){
			isEpOxidation = true;
		}
		if(isEpOxidation) return "Dehydrogenation";
		else return null;
	}
}
