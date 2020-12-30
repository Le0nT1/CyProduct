package typeOneFeatures;


public class GenerateTitleLine_TypeOne {
	int depth_neighborAtomType;
	int depth_neighborAtomDescriptor;
	final String[] attributes = {"Name", "Bond", "1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2|E1", "3A4","ConnectedAtoms_Left","ConnectedAtoms_Right"};
	final String[] atomicDescriptorTitle = {"AtomDegree", "AtomHybridization", "AtomValence", "EffectiveAtomPolarizability", "PartialSigmaCharge",
			"PartialTChargeMMFF94", "PiElectronegativity", "SigmaElectronegativity", "StabilizationPlusCharge", "Atom_ASA"};
	String[] molecularDescriptorTitle = {"AlogP", "APol",  "MomentOfInertia", "RotatableBondsCount", "TPSA", "Weight", "XLogP"};
	String[] molecularFPTitle = new String[FingerPrintQueries_TypeOne.moleQueryList.size()];
	String[] currentBondFPTitle = new String[FingerPrintQueries_TypeOne.queriesList.size()];
	String[] currentBondAtomTypeTitle = new String[GenerateNeighborFeatures_TypeOne.atomTypeLookupTable.length];
	String[] currentBondAtomDescriptorTitle = new String[10*2];//nine atom descriptors
	String[] neighborAtomTypeTitle;
	String[] neighborAtomDescriptorTitle;
	
	public GenerateTitleLine_TypeOne(int depth_neighborAtomType, int depth_neighborAtomDescriptor){
		this.depth_neighborAtomType = depth_neighborAtomType;
		this.depth_neighborAtomDescriptor = depth_neighborAtomDescriptor;
	}
	

	//String[] atomTypeTitle = BoMAtomFeatures.atomTypeLookupTable;

	//String titleLine = Utilities.creatTypeOneTitleLine(cypList, molecularDescriptorTitle, bondDescriptorTitle, atomTypeTitle, atomicDescriptorTitle, bondFpTitle, depth, mergedFP);//Include left and right features
	
	public void generateMolecularFPTitle(){//2nd
		for(int i = 0; i < molecularFPTitle.length; i++){
			molecularFPTitle[i] = "moleFP_" + i;
		}
	}
	public void generateCurrentBondFPTitle(){//3rd
		for(int i = 0; i < currentBondFPTitle.length; i++){
			currentBondFPTitle[i] = "currentBondFP_" + i;
		}
	}
	
	public void generateCurrentBondAtomTypeTitle(){//4th
		for(int i = 0; i < GenerateNeighborFeatures_TypeOne.atomTypeLookupTable.length; i++){
			currentBondAtomTypeTitle[i] = "CurrentBond_" + GenerateNeighborFeatures_TypeOne.atomTypeLookupTable[i];
		}
	}
	
	public void generateCurrentBondAtomDescriptorTitle(){//5th
		for(int i = 0; i < atomicDescriptorTitle.length; i++){
			currentBondAtomDescriptorTitle[i] = "CurrentBond_" + atomicDescriptorTitle[i] + "_left";
		
							
		}
		for(int i = 0; i < atomicDescriptorTitle.length; i++){
			currentBondAtomDescriptorTitle[i+10] = "CurrentBond_" + atomicDescriptorTitle[i] + "_right";
		}
	}
	
	public void generateNeighborAtomTypeTitle(){
		neighborAtomTypeTitle = new String[GenerateNeighborFeatures_TypeOne.atomTypeLookupTable.length*depth_neighborAtomType];
		for(int d = 1; d < depth_neighborAtomType + 1; d++){
			for(int i = 0; i < GenerateNeighborFeatures_TypeOne.atomTypeLookupTable.length; i++){
				int idx = i + (d-1)*GenerateNeighborFeatures_TypeOne.atomTypeLookupTable.length;
				neighborAtomTypeTitle[idx] = "Neighbor_d=" + d +"_" + GenerateNeighborFeatures_TypeOne.atomTypeLookupTable[i];
			}
		}
	}

	public void generateNeighborAtomDescriptorTitle(){
		neighborAtomDescriptorTitle = new String[GenerateNeighborFeatures_TypeOne.roughAtomTypeTable.length*depth_neighborAtomDescriptor*10];
		for(int d = 1; d < depth_neighborAtomDescriptor + 1; d++){
			for(int i = 0; i < GenerateNeighborFeatures_TypeOne.roughAtomTypeTable.length; i++){
				for(int j = 0; j < atomicDescriptorTitle.length; j++){
					int idx = j + i*atomicDescriptorTitle.length + (d-1)*GenerateNeighborFeatures_TypeOne.roughAtomTypeTable.length*atomicDescriptorTitle.length;
					//System.out.println("d=" + d + "," + "i=" + i + ",j=" + j +",idx="+ idx);
					neighborAtomDescriptorTitle[idx] = "Neighbor_d=" + d + "_" + GenerateNeighborFeatures_TypeOne.roughAtomTypeTable[i] + "_" + atomicDescriptorTitle[j];
					
				}
			}
		}
	}
	
	public String generateTitleLine(){
		generateMolecularFPTitle();
		generateCurrentBondFPTitle();
		generateCurrentBondAtomTypeTitle();
		generateCurrentBondAtomDescriptorTitle();
		generateNeighborAtomTypeTitle();
		generateNeighborAtomDescriptorTitle();
		StringBuilder tempBuffer = new StringBuilder();
		for(int i = 0; i < attributes.length; i++){
			tempBuffer.append(attributes[i] + ",");
		}
		for(int i = 0; i < molecularDescriptorTitle.length; i++){
			tempBuffer.append(molecularDescriptorTitle[i] + ",");
		}
		for(int i = 0; i < molecularFPTitle.length; i++){
			tempBuffer.append(molecularFPTitle[i] + ",");
		}
		for(int i = 0; i < currentBondFPTitle.length; i++){
			tempBuffer.append(currentBondFPTitle[i] + ",");
		}
		for(int i = 0; i < currentBondAtomTypeTitle.length; i++){
			tempBuffer.append(currentBondAtomTypeTitle[i] + ",");
		}
		for(int i = 0; i < currentBondAtomDescriptorTitle.length; i++){
			tempBuffer.append(currentBondAtomDescriptorTitle[i] + ",");
		}

		for(int i = 0; i < neighborAtomTypeTitle.length; i++){
			tempBuffer.append(neighborAtomTypeTitle[i] + ",");
		}
		for(int i = 0; i < neighborAtomDescriptorTitle.length; i++){
			if(i == neighborAtomDescriptorTitle.length-1){
				tempBuffer.append(neighborAtomDescriptorTitle[i]);
				break;
			}
			tempBuffer.append(neighborAtomDescriptorTitle[i] + ",");
		}
		return tempBuffer.toString();
	}





}
