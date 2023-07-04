package cyProduct;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smarts.SmartsFragmentExtractor;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class TestSmarts {
	public static void main(String[] args) throws Exception{
		TestSmarts ts = new TestSmarts();
		String smiles = "Nc1ccc(O)cc1";
		//String smiles = "NC1=CC=C(O)C=C1";
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IAtomContainer mole = sp.parseSmiles(smiles);
		ts.checkMatch(mole);
		//ts.checkSMARTSFragmenter(mole);
	}

	public void checkMatch(IAtomContainer molecule) throws Exception{
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);		
		SmilesGenerator sg = new SmilesGenerator(SmiFlavor.Isomeric);
		//String smarts = "[#7]-[#6]-1=[#6]-[#6]=[#6](-[#8])-[#6]=[#6]-1";
		String smarts =  "NC1=CC=C(O)C=C1";
		//String smarts = "[#7]-c1ccc(-[#8])cc1";
		SmartsPattern smp2 = SmartsPattern.create(smarts);
		smp2.setPrepare(false);
		//smp2.setPrepare(true);
		System.out.println(sg.create(molecule));
		molecule = sp.parseSmiles(sg.create(molecule));
		boolean status = smp2.matches(molecule);
		System.out.println(status);
	}
	public void checkSMARTSFragmenter(IAtomContainer molecule) throws Exception{
		SmartsFragmentExtractor sfe = new SmartsFragmentExtractor(molecule);
		sfe.setMode(1);
		String  smarts = sfe.generate(new int[]{0,1,3,4,5,6});
		System.out.println(smarts);
	}
}
