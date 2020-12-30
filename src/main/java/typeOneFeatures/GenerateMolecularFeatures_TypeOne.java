package typeOneFeatures;

import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class GenerateMolecularFeatures_TypeOne {

	/**
	 * Generate molecular features for a given molecule.
	 * The features are molecular descriptors.
	 * It's hard coded as -10.0 so far.
	 * @param molecule
	 * @return
	 * @throws CDKException
	 */
	public static String generateMolecularFeatures_new(IAtomContainer molecule) throws CDKException {
		IAtomContainer mol = molecule;
		String molecular;
		String []res=new String[7];
		ALOGPDescriptor Alogp = new ALOGPDescriptor();
		res[0] =  Alogp.calculate(mol).getValue().toString().split(",")[0];
		APolDescriptor Apol = new APolDescriptor();
		res[1] = Apol.calculate(mol).getValue().toString();
		MomentOfInertiaDescriptor Mo = new MomentOfInertiaDescriptor();
		res[2] = Mo.calculate(mol).getValue().toString().split(",")[6];
		RotatableBondsCountDescriptor Ro = new RotatableBondsCountDescriptor();
		res[3] = Ro.calculate(mol).getValue().toString();
		TPSADescriptor Tp = new TPSADescriptor();
		res[4] = Tp.calculate(mol).getValue().toString();
		WeightDescriptor We = new WeightDescriptor();
		res[5]= We.calculate(mol).getValue().toString();
		XLogPDescriptor Xl = new XLogPDescriptor();
		res[6]= Xl.calculate(mol).getValue().toString();
		StringBuffer tempString=new StringBuffer();
		for(int i=0;i<res.length;i++){
			if(i==res.length-1){
				tempString.append(res[i]);
				break;
			}
			tempString.append(res[i]).append(",");
		}		
		molecular = tempString.toString(); 
		return molecular;
	}
	
	/**
	 * This function generates fingerprint for the molecule to describe the important groups may affect the whole property of the molecule 
	 * @param mole
	 * @param bond
	 * @return
	 * @throws Exception
	 */
	public static String generateMoleculeFP(IAtomContainer molecule) throws Exception {
		IAtomContainer mole = molecule;
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mole);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mole.getBuilder());		
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(mole);
		adder.addImplicitHydrogens(mole);
		String[] moleFP = new String[FingerPrintQueries_TypeOne.moleQueryList.size()];
		for(int i = 0; i < moleFP.length; i++){
			moleFP[i] = "0";
		}
		//Iterate through all queries in order
		Iterator<Entry<String, String>> allQueries = FingerPrintQueries_TypeOne.moleQueryList.entrySet().iterator();
		int counter = 0; //Track the index of a query
		while(allQueries.hasNext()) { // Exploit all querries
			Entry<String,String> onePattern = allQueries.next();
			IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(onePattern.getValue(),builder); //Set the SMART Pattern
			boolean occurrence = smartsPattern.matches(mole);//match the pattern with the molecule
			if(occurrence){
				moleFP[counter] = "1";
			}
			else{
				moleFP[counter] = "0";
			}
			
			//System.out.println(onePattern.getKey());
			counter++;
			
		}
		for(int i = 0; i < moleFP.length; i++){
			if(moleFP[i] == null || moleFP[i].isEmpty()){
				System.out.println(i);
			}
		}
		
		StringBuilder moleFPBuilder = new StringBuilder();
		for(int i = 0; i < moleFP.length; i++){
			if(i == moleFP.length-1){
				moleFPBuilder.append(moleFP[i]);
				break;
			}
			moleFPBuilder.append(moleFP[i] + ",");
		}
		return moleFPBuilder.toString();
	}
}
