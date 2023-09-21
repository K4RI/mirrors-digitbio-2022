package multibinning.data;

import java.util.ArrayList;

public class MacroBin 
{
	public double lowerBound;
	public double upperBound;
	public int numPoints;
	public ArrayList<Integer> microBinIDs;
	public ArrayList<Integer> pointIDs;
	public int[] mdh;
	public int numMDHs;
	public DataMatrix dataMatrix;
	public double[] means;
	public double[] devs;
	public String name;
	
	public MacroBin(double lowerBound, double upperBound)
	{
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		microBinIDs = new ArrayList<Integer>();
		pointIDs = new ArrayList<Integer>();
	}
}