package multibinning.business;

import java.util.*;

public class Constants 
{
	public static String FILE_INPUT; 					
	public static String FILE_CP_OUTPUT; 				
	public static String FILE_RUNTIME_OUTPUT;
	public static String FILE_DATA_OUTPUT;
	public static String FILE_ORG_INPUT;
	public static String FILE_ORG_NORMAL_INPUT;
	public static String FILE_IDIST = null;
	
	public static int NUM_ROWS;
	public static int NUM_MEASURE_COLS;
	public static int NUM_CAT_CONTEXT_COLS = 0;	// number of categorical context dimensions, currently unused
	public static int NUM_ORG_CAT_CONTEXT_COLS = 0;
	public static String FIELD_DELIMITER = ";";
	public static double MAX_VAL;
	public static double DELTA;
	public static double GAMMA = 0.01;
	public static int QUANTILE = 3; // Added by kailash: t_q parameter from the paper (roughly 3) 
	public static int INIT_BIN_COUNT = -1;
	public static int MAX_BINS = 50;
	public static boolean USE_THRESHOLD = false;
	public static boolean USE_CE = false;
	public static boolean USE_CLUSTERING = false;
	public static int NUM_CENTROIDS = 5;
	public static double ALPHA;
	public static int NUM_SUBSAMPLING;
	public static Random RAN_INDEX_BLOCK;
	
	public static final int MB_MDL_DP = 0;	// MDL dynamic programming
	public static final int MB_MDL_GD = 1;	// MDL greedy
	public static final int MB_MDL_EGD = 2;	// normal greedy
	public static final int EW = 3;			// equal-width
	public static final int EF = 4;			// equal-frequency
	public static final int CAIM = 5;		// equal-frequency
	public static final int MB_NM_GD = 6;	// normal greedy
	public static final int J_MDL = 9;		// normal greedy
	public static final int J_MDL_HM = 10;	// normal greedy
	public static final int J_GD_HM = 11;	// normal greedy
	public static final int M_GD_CAIM = 7;	// normal greedy
	public static final int M_GD = 8;		// normal greedy
	public static final int MMIC = 12;		// normal greedy
	public static final int FayyadMDL = 14;	// normal greedy
	public static final int KoMDL = 15;		// normal greedy
	public static final int DP_MEAN = 16;	// MDL dynamic programming
	public static final int DP_MIN = 17;	// MDL dynamic programming
	public static final int DP_MAX = 18;	// MDL dynamic programming
	public static final int MVD = 19;	// MDL dynamic programming
	public static final int PCABINNING = 20;
	public static int METHOD = 3;			// method used
	
	public static int LP_NORM = 2;							
	public static int LOG_BASE = 2;
	//public static int MAX_BINS;
	
	public static final double MAX_ERROR = 0.00001;
	
	public static ArrayList<String> CLASS_LABELS;
	
	public static double[] MIN_COLS;
	public static double[] MAX_COLS;
	
	public static final double C0 = 2.865064;
	
	public static double[] DATA_MEANS;
	public static double[] DATA_DEVS;
	public static double[] CRES;
	public static int[][] DISC_DATA;
	public static int[] ODIMS;
}
