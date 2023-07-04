package cyProduct.predictionhelpers;

import java.util.ArrayList;

import org.openscience.cdk.interfaces.IAtomContainer;

public class MoleResultPair {
	public final IAtomContainer molecule;
	public final ArrayList<String> resultList;
	public MoleResultPair(IAtomContainer molecule, ArrayList<String> resultList){
		this.molecule = molecule;
		this.resultList = resultList;
	}
	public IAtomContainer getMolecule() {
		return molecule;
	}
	public ArrayList<String> getResult() {
		return resultList;
	}
}
