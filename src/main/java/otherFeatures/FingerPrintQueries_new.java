package otherFeatures;

import java.util.LinkedHashMap;
import java.util.Map.Entry;

import otherFeatures.GeneratePatterns;

/**
 * This class will collect all fingerprints that will be used 
 * @author Tian
 *
 */
public class FingerPrintQueries_new {
	public static LinkedHashMap<String, String> newQueriesList;
	public static LinkedHashMap<String, String> queriesList; 
	public static LinkedHashMap<String, String> moleQueryList;
	static {
		LinkedHashMap<String, String> queries_new = new LinkedHashMap<String, String>();
		LinkedHashMap<String, String> queries = new LinkedHashMap<>();
		LinkedHashMap<String, String> molequeries = new LinkedHashMap<String, String>();
		queries_new.put("N-N=O", "[#6;A][#7]-[#7]=O");
		queries_new.put("C-NH2", "[#7;AH2]*");
		queries_new.put("C-N-POOH", "[#8]P(*)(=O)[#7]-*");
		//queries.put("benzene-OCH3", "[#6;AH3][#8]-[#6]1=[#6;AH1][#6](-*)=[#6](-*)-[#6](-*)=[#6;AH1]1 |c:9,t:2,5|");
		queries_new.put("benzene-OCH3", "[#6;AH3][#8]-[#6]1=[#6;AH1][#6]=[#6]-[#6]=[#6;AH1]1 |c:4,6,t:2|");
		queries_new.put("C-N-C-CO-C", "*-[#7]-[#6]-[#6]=O");
		queries_new.put("C-N-CO-C", "[#6;A][#7]-[#6]=O");
		queries_new.put("C-CO", "*-[#6]=O");
		queries_new.put("C-O-C", "*-[#8]-*");
		//queries.put("C-O-C-O-C", "[#6]-[#8]-[#6]-[#8]-[#6]");
		queries_new.put("O-C-O", "[#8;A][#6][#8;A]");
		queries_new.put("3-fused-Ring", "[#6;AH3][#6]-1=[#6]-[#6]=[#6]-2-[#6]-[#8]-[#6]-3=[#6]-[#6]=[#6]-[#6]=[#6]-3-[#6]-2=[#6]-1 |c:9,11,15,t:1,3,7|");
		queries_new.put("benzene-N", "*-[#6]-1=[#6]-[#7]=[#6](-*)-[#6]=[#6]-1 |c:6,t:1,3|");
		queries_new.put("(P=O)(O)-O", "*-[#8]P(*)(=O)[#8]-*");
		queries_new.put("CH3-NCO(CO)", "[#6;AH3][#7](-[#6]=O)-[#6]=O");
		queries_new.put("A-O-C-C-OH","[#8]-[#6]-[#6]-[#8]-*");
		queries_new.put("5-member-ring-N-side-chain","[#6,#1]-[#6]-[#6]-[#6]-1=[#7]-[#6]-[#6]-[#6,#7]-1 |t:3|");
		queries_new.put("A-N-C-C-C=C", "*-[#7]-[#6]-[#6]-[#6]=[#6]");
		queries_new.put("A-C(OH)-C-N", "[#8]-[#6](-*)-[#6](-*)-[#7]-*");
		queries_new.put("A-COOH", "[#8]-[#6](-*)=O");
		
		
		molequeries.put("ben-O-CH3", "c-[OX2]-[CH3]");
		molequeries.put("N-C", "[N]-[C!H0]"); //N-C in a chain
		molequeries.put("S=C-N-", "[C!H0]-[N!R]-[C!R]=[SX1]");//test mole_9
		molequeries.put("NNX-P=O", "P(-N-[C!H0])(-[N,O])(-[N,O])(=O)");//test mole_13
		molequeries.put("R-N-N-CH3", "N-N-[C!H0]");
		molequeries.put("(CH3-N)-S(=O)2-R", "S(-N-[C!H0])(=O)(=O)(-*)");//test mole_15
		molequeries.put("R-N-C-CF3", "C(-[F,Cl])(-[F,Cl])(-[F,Cl])(-C-N)"); // test mole 16
		molequeries.put("C-P=O(OCH3)2", "P(-O-[C!H0])(-[N,O])(-*)(=O)");//test mole 17
		molequeries.put("(CH3-O)2-P=S(S-R)", "P(-O-[C!H0])(-O-[C!H0])(=[SX1])(-S-*)");//test mole 19
		molequeries.put("Fused benzene rings_6to3", "a12aaaaa1aaa2");// test mole_6,7
		molequeries.put("Fused benzene rings_6to4", "a12aaaaa1aaaa2");
		molequeries.put("Fused benzene rings_6to5", "a12aaaaa1aaaaa2");// test mole_6,7
		molequeries.put("rearrangement", "[NR]1-[CR]=[CR]-[CR]-[CR]=[CR]1");
		//molequeries.put("sulfuric acid",
		//		"[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[([OX2H]),$([OX1-])]");	// 5
//		queries.put("aryl aldehyde", "[#6;X3H1](-[#6;a])=[O;v2X1]");	//7
		molequeries.put("indole",
				"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");	//8
		
		
		//Check ring. whether is in a ring, aromatic ring
		//queries.put("carboxyl_carbon","[$([CX3](=[OX1])(-[OX2])),$([CX3+]-[OX1-])]");
		/**
		 * 30 bond fingerprint queries
		 */
		queries.put("ben-O-CH3", "c-[OX2]-[CH3]"); //test Mole_3; 
		queries.put("connected to aromatic_l1", "*-[a;R]"); //test Mole_3 ; check the bond connected to benzene
		queries.put("connected to aromatic_l2", "*-*-[a;R]"); //test Mole_3; check the bond is one bond away from the benzene
		
		
		/**
		 * Useful
		 */
		//May not need theses fused/unfused ring patterns.
		queries.put("Fused benzene rings_6to3", "a12aaaaa1aaa2");// test mole_6,7
		queries.put("Fused benzene rings_6to4", "a12aaaaa1aaaa2");
		queries.put("Fused benzene rings_6to5", "a12aaaaa1aaaaa2");// test mole_6,7
		queries.put("rearrangement", "[NR]1-[CR]=[CR]-[CR]-[CR]=[CR]1");
		/**
		 * C-N
		 */
		queries.put("N-C", "[N!R]-[C!H0!R]"); //N-C in a chain
		queries.put("N_r-C", "[NR]-[C!H0!R]"); //N-C that N is in a ring but C is not
		queries.put("N-C_r", "[N!R]-[C!H0R]"); //N-C that C is in a ring but N is not
		queries.put("N_r-C_r", "N-[C!H0R]");//N-C that both N and C are in a ring
		queries.put("ben-C-N-R", "c-[C!r!H0]-[N!-r]");// test molecule_5
		queries.put("C-N-C=O", "[C!H0]-[N!R]-[C!R]=[OX1]"); //test mole_7,8
		queries.put("C-N-C=O_ring", "[C!H0]-[NR]-[CR]=[OX1]"); //test mole_7,8
		queries.put("S=C-N-", "[C!H0]-[N!R]-[C!R]=[SX1]");//test mole_9
		queries.put("C-N-C=C-","[C!H0]-[N!R]-C=C");//test mole_10
		queries.put("N-C-C-OH", "N-C-C-[OH1]");//test mole_12
		queries.put("-N-C-C=O", "[C!H0]-N-C-C=[OX1]");//test mole_11
		queries.put("NNX-P=O", "P(-N-[C!H0])(-[N,O])(-[N,O])(=O)");//test mole_13
		queries.put("R-N-N-CH3", "N-N-[C!H0]");
		queries.put("(CH3-N)-S(=O)2-R", "S(-N-[C!H0])(=O)(=O)(-*)");//test mole_15
		queries.put("R-N-C-CF3", "C(-[F,Cl])(-[F,Cl])(-[F,Cl])(-C-N)"); // test mole 16


//		/**
//		 * C-O
//		 */
		queries.put("C-P=O(OCH3)2", "P(-O-[C!H0])(-[N,O])(-*)(=O)");//test mole 17
		queries.put("O-C-O in a ring", "[O;R]-[C]-[O;R]");
		queries.put("C-O-Alkyl", "[C!R!H0]-[OX2]-[C!R]");//test mole 18
		queries.put("(CH3-O)2-P=S(S-R)", "P(-O-[C!H0])(-O-[C!H0])(=[SX1])(-S-*)");//test mole 19	
		queries.put("ester", "[CX3](=[OX1])(-O-*)");//ester will not hydrolysis because it requires hydrolase; test mole 20
		queries.put("carboxyl_", "[#8;A;X2H1,X1-][#6]([#6,#1;A])=O");  // 0
		queries.put("hydroxyl", "[#6][OX2H1]");						  // 1
		queries.put("indole",
				"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");	//8
		
		
		newQueriesList = queries_new;
		newQueriesList.putAll(GeneratePatterns.queriesList);
		queriesList = queries;
		moleQueryList = molequeries;
		moleQueryList.putAll(queries_new);
	}
}
