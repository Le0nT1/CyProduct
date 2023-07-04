package otherFeatures;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.qsar.descriptors.bond.AtomicNumberDifferenceDescriptor;
import org.openscience.cdk.qsar.descriptors.bond.BondPartialPiChargeDescriptor;
import org.openscience.cdk.qsar.descriptors.bond.BondPartialSigmaChargeDescriptor;
import org.openscience.cdk.qsar.descriptors.bond.BondPartialTChargeDescriptor;
import org.openscience.cdk.qsar.descriptors.bond.BondSigmaElectronegativityDescriptor;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class BoMBondDescriptors {
	/**
	 * A function generates bond features for a given bond in a molecule
	 * @param molecule
	 * @param oneBond
	 * @return
	 * @throws Exception
	 */
	public static String generateBondFeatures(IAtomContainer molecule, IBond oneBond) throws Exception{
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);
		adder.addImplicitHydrogens(molecule);
		String bondfeatures;
		AtomicNumberDifferenceDescriptor and = new AtomicNumberDifferenceDescriptor();
		BondPartialPiChargeDescriptor bppc = new BondPartialPiChargeDescriptor();
		BondPartialSigmaChargeDescriptor bpsc = new BondPartialSigmaChargeDescriptor();
		BondPartialTChargeDescriptor bptc = new BondPartialTChargeDescriptor();
		BondSigmaElectronegativityDescriptor bse = new BondSigmaElectronegativityDescriptor();
		
		String andValue = and.calculate(oneBond, molecule).getValue().toString();
		String bppcValue = "0";//bppc.calculate(oneBond, molecule).getValue().toString(); //"0";
		String bpscValue = bpsc.calculate(oneBond, molecule).getValue().toString();
		String bptcValue = "0";//bptc.calculate(oneBond, molecule).getValue().toString();//"0";
		String bseValue = bse.calculate(oneBond, molecule).getValue().toString();	
		bondfeatures = andValue + "," + bppcValue + "," + bpscValue + "," + bptcValue + "," + bseValue;		
		return bondfeatures;
	}
	
}
