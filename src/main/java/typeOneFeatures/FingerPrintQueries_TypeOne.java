package typeOneFeatures;

import java.util.LinkedHashMap;

public class FingerPrintQueries_TypeOne {
	public static LinkedHashMap<String, String> queriesList; 
	public static LinkedHashMap<String, String> moleQueryList;
	static {
		LinkedHashMap<String, String> queries = new LinkedHashMap<String, String>();
		LinkedHashMap<String, String> molequeries = new LinkedHashMap<String, String>();
		/**
		 * 14 molecular fingerprint queries
		 */
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
		
		moleQueryList = molequeries;
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
		//queries.put("c-N(C)-c in a ring", "N(-[CR])(-[CR])");//test mole_14 trivial
		//queries.put("C-NH-C in a ring", value); //suspicious

//		/**
//		 * C-O
//		 */
		queries.put("C-P=O(OCH3)2", "P(-O-[C!H0])(-[N,O])(-*)(=O)");//test mole 17
		queries.put("O-C-O in a ring", "[O;R]-[C]-[O;R]");
		queries.put("C-O-Alkyl", "[C!R!H0]-[OX2]-[C!R]");//test mole 18
		queries.put("(CH3-O)2-P=S(S-R)", "P(-O-[C!H0])(-O-[C!H0])(=[SX1])(-S-*)");//test mole 19
		//queries.put("-C-OH-C-O-benzene", value);//Trivial		
		queries.put("ester", "[CX3](=[OX1])(-O-*)");//ester will not hydrolysis because it requires hydrolase; test mole 20

//		
		queries.put("carboxyl_", "[#8;A;X2H1,X1-][#6]([#6,#1;A])=O");  // 0
		queries.put("hydroxyl", "[#6][OX2H1]");						  // 1
//		queries.put("aromatic_4", "a1aaa1");						  // 2 Added by ST
//		queries.put("aromatic_5", "a1aaaa1");						  // 3 Added by ST
//		queries.put("aromatic_6", "a1aaaaa1");						  // 4 Added by ST. The original [*,a] is removed.
//		queries.put("sulfuric acid",
//				"[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[([OX2H]),$([OX1-])]");	// 5
//		queries.put("aryl aldehyde", "[#6;X3H1](-[#6;a])=[O;v2X1]");	//7
		queries.put("indole",
				"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");	//8
		queriesList = queries;
	}
}
