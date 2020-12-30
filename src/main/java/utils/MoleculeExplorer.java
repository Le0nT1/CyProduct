/**
 * Authors: Yannick Djoumbou Feunang
 * Class Description:
 */


package utils;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;

public class MoleculeExplorer {

	
	/**
	 * This function applies some preprocessing operations, such as setting the
	 * flag of atoms from aromatic rings to "ISAROMATIC", and kelulizing
	 * molecules.
	 * 
	 * @param molecule
	 *            : A molecule of interest
	 * @return : A processed molecule (AtomContainer)
	 * @throws CDKException
	 */
	public static IAtomContainer preprocessContainer(IAtomContainer molecule)
			throws CDKException {
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder.getInstance(molecule.getBuilder()).addImplicitHydrogens(molecule);		 		
		for (IBond bond : molecule.bonds()) {
			if (bond.getFlag(CDKConstants.SINGLE_OR_DOUBLE)) {
				bond.setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(0).setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(1).setFlag(CDKConstants.ISAROMATIC, true);

			} 
		}
		Kekulization.kekulize(molecule);				
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(molecule);
		sdg.generateCoordinates();		
		IAtomContainer layedOutMol = sdg.getMolecule();
		return layedOutMol;
	}
	public static boolean isMixture(IAtomContainer molecule) throws CDKException{
		// compound is not a mixture (checkConnectivity returns 2 or more atomContainers)
		boolean mixture = ConnectivityChecker.partitionIntoMolecules(molecule).getAtomContainerCount()>1;
		return mixture;	
	}
	
	public static boolean containsCarbon(IAtomContainer molecule) {
		boolean carbon = false;
		for(IAtom at : molecule.atoms()){
			if(at.getAtomicNumber() == 6){
				carbon = true;
				break;
			}
		}	
		return carbon;
	}	
	public static boolean isEtherLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		String constraints = "[$([#8;X2][#6;A;H2X4]!@-[#6;A;X4](!@-[!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([#8]!@-[#6;A;X4](!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])][#6;A;H2X4]!@-[#6;A;X3](!@=[O;X1])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])"
				+ "]";		
		SmartsPatternCDK pattern = new SmartsPatternCDK(constraints);
		b = pattern.hasSMARTSPattern(molecule)>0;
			
		return b;
		
	}
	public static boolean isGlyceroLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		String constraints = "[$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#6;A;H2X4R0][#8]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])]),$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([#6;A;H2X4R0][OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#8;X2]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])])]"	;
		SmartsPatternCDK pattern = new SmartsPatternCDK(constraints);
		b = pattern.hasSMARTSPattern(molecule)>0;
			
		return b;
		
	}	
	public static boolean isGlycerophosphoLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		String constraints =  				"[#8]P([!#1!#6;OX1-,OX2H1,$([O]-[#6])])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]";
		SmartsPatternCDK pattern = new SmartsPatternCDK(constraints);
		b = pattern.hasSMARTSPattern(molecule)>0;
			
		return b;
		
	}

	public static boolean isSphingoLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		String constraints = "[$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8]),$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X4]=[#6;A;H1X4]-[#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8])]";
		SmartsPatternCDK pattern = new SmartsPatternCDK(constraints);
		b = pattern.hasSMARTSPattern(molecule)>0;
			
		return b;
		
	}

	public static boolean isInvalidCandidate(IAtomContainer molecule) throws CDKException, SMARTSException{
		boolean invalid = false;

		
		// Calculate the weight of specified element type in the supplied
		WeightDescriptor weightD = new WeightDescriptor();
		IDescriptorResult weight = weightD.calculate(molecule).getValue();
		
		
		IAtomContainer pmol = preprocessContainer(molecule);
		
		if(isMixture(pmol)){ //|| containsCarbon(pmol)
			invalid = true;
		} else{
			if(isEtherLipid(pmol) || isGlyceroLipid(pmol) || 
					isGlycerophosphoLipid(pmol)	 || isSphingoLipid(pmol)){
				invalid = true;
			}
		}
		//Double.valueOf(weight.toString()) > 1500.0 || 
		return invalid;
	}	
	
}
