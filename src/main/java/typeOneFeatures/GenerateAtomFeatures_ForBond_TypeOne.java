package typeOneFeatures;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.NoSuchAtomTypeException;
import org.openscience.cdk.geometry.surface.NumericalSurface;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.atomic.AtomDegreeDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.AtomHybridizationDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.AtomValenceDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.EffectiveAtomPolarizabilityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PartialSigmaChargeDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PartialTChargeMMFF94Descriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PiElectronegativityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.SigmaElectronegativityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.StabilizationPlusChargeDescriptor;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class GenerateAtomFeatures_ForBond_TypeOne {
	/**
	 * This function generates atomic features for a atom in the molecule
	 * @param molecule
	 * @param oneAtom
	 * @return String[] whose each tuple is a real value
	 * @throws CDKException
	 * @throws NoSuchAtomTypeException
	 */
	public static String[] generateAtomDescriptorFeatures(IAtomContainer molecule, IAtom oneAtom) throws CDKException, NoSuchAtomTypeException {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		//HydrogenAdder adder_extra = new HydrogenAdder();
		
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(molecule);	
		adder.addImplicitHydrogens(molecule);
		
		
		AtomDegreeDescriptor degree = new AtomDegreeDescriptor();
		AtomHybridizationDescriptor hy = new AtomHybridizationDescriptor();
		AtomValenceDescriptor va = new AtomValenceDescriptor();
		EffectiveAtomPolarizabilityDescriptor ep = new EffectiveAtomPolarizabilityDescriptor();
		
		PartialSigmaChargeDescriptor psc = new PartialSigmaChargeDescriptor();
		//This PartialTChargeMMFF94Descriptor() may generate NaN for some Atoms in some Molecules 
		PartialTChargeMMFF94Descriptor ptc = new PartialTChargeMMFF94Descriptor();
		PiElectronegativityDescriptor pen = new PiElectronegativityDescriptor();
		SigmaElectronegativityDescriptor sen = new SigmaElectronegativityDescriptor();
		StabilizationPlusChargeDescriptor spc = new StabilizationPlusChargeDescriptor();
		
		
		String[] atomicFeatures =new String[10];		
		DescriptorValue d = degree.calculate(oneAtom, molecule);
		atomicFeatures[0] = d.getValue().toString();
		DescriptorValue h = hy.calculate(oneAtom, molecule);
		atomicFeatures[1] = h.getValue().toString();
		DescriptorValue v = va.calculate(oneAtom, molecule);
		atomicFeatures[2] = v.getValue().toString();
		DescriptorValue e = ep.calculate(oneAtom, molecule);
		atomicFeatures[3] = e.getValue().toString();
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		DescriptorValue ps = psc.calculate(oneAtom, molecule);
		atomicFeatures[4] = ps.getValue().toString();
		DescriptorValue pt = ptc.calculate(oneAtom, molecule);
		atomicFeatures[5] = pt.getValue().toString();
		DescriptorValue pe = pen.calculate(oneAtom, molecule);
		atomicFeatures[6]= pe.getValue().toString();
		DescriptorValue se = sen.calculate(oneAtom, molecule);
		atomicFeatures[7] = se.getValue().toString();
		//System.out.println(res[7]);
		DescriptorValue sp = spc.calculate(oneAtom, molecule);
		atomicFeatures[8] = sp.getValue().toString();
		Integer atomIdx = molecule.indexOf(oneAtom);
		AtomContainerManipulator.suppressHydrogens(molecule);
		Double asa, full_asa;
		try{
			NumericalSurface ns = new NumericalSurface(molecule);
			ns.calculateSurface();
			asa = ns.getSurfaceArea(atomIdx);
			full_asa = ns.getTotalSurfaceArea();
		}catch(Exception no3DException){
			asa = 0.0;
			full_asa = 0.0;
			System.out.println("Warning: the input molecule does not have 3D coordinates. The area ratios of atoms are set as -1.0");
		}
		if(full_asa < 0.0001){
			atomicFeatures[9] = Double.toString(-1.0);
		}
		else{
			atomicFeatures[9] = Double.toString(asa/full_asa);
		}
		if(oneAtom.getSymbol().contains("H")){
			System.out.println("Problematic: " + molecule.getProperties().get("cdk:Title") + "," + oneAtom.getAtomTypeName());
		}
		for(int i = 0; i < atomicFeatures.length; i++){
			if(atomicFeatures[i].equals("NaN")){
				throw new CDKException("PartialTChargeMMFF94Descriptor is not valid for molecule" + molecule.getProperties().get("cdk:Title") + "," + oneAtom.getAtomTypeName());
			}
		}
		return atomicFeatures;
	}
	


}
