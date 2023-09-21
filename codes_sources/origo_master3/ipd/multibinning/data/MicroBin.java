package multibinning.data;

import java.util.*;

public class MicroBin 
{
	public double lowerBound;
	public double upperBound;
	public ArrayList<Double> distinctValues;
	public ArrayList<Integer> distinctValueSupports;
	public ArrayList<Double> allValues;
	public ArrayList<Integer> pointIDs;
	public DataMatrix dataMatrix;
	public IndexMatrix indexMatrix;
	public ArrayList<Integer> dims;
	public double[] means;
	public double[] devs;
	//public int[] mdh;
	//public int numMDHs;
	
	public MicroBin(double lowerBound, double upperBound)
	{
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		pointIDs = new ArrayList<Integer>();
		distinctValues = new ArrayList<Double>();
		distinctValueSupports = new ArrayList<Integer>();
		allValues = new ArrayList<Double>();
		dims = new ArrayList<Integer>();
	}
}