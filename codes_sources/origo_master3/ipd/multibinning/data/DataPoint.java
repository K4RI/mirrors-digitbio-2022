package multibinning.data;

import java.util.*;

public class DataPoint 
{
	public int universalID;									
	public int localID;										
	public double[] measures;
	public String[] catContexts;	// contains categorical context dimensions
	
	public ArrayList<Integer> neighbors;					
	public ArrayList<Double> distToNeighbors;				
	public double lrd;										
	public String classID;										
	
	public DataPoint(int numMeasures, int numCatContexts)
	{
		if (numMeasures > 0)
			measures = new double[numMeasures];
		
		if (numCatContexts > 0)
			catContexts = new String[numCatContexts];
	}
	
	public DataPoint(DataPoint p, int numMeasures, int numCatContexts)
	{
		this.universalID = p.universalID;
		this.localID = p.localID;
		
		if (numMeasures > 0)
			measures = new double[numMeasures];
		
		if (numCatContexts > 0)
			catContexts = new String[numCatContexts];
		
		for (int i = 0; i < numMeasures; i++)
			measures[i] = p.measures[i];
		
		for (int i = 0; i < numCatContexts; i++)
			catContexts[i] = p.catContexts[i];
	}
	
	public static double distanceLNorm(int expo, ArrayList<Integer> dims, DataPoint p, DataPoint q)
	{
		int numMeasuresCols = dims.size();
		double dist = 0;
		int curCol;
		for (int i = 0; i < numMeasuresCols; i++)
		{
			curCol = dims.get(i).intValue();
			dist += Math.pow(Math.abs(p.measures[curCol] - q.measures[i]), expo);
		}
		
		return Math.pow(dist, 1.0 / expo);
	}
	
	public static double distanceLNorm(int expo, DataPoint p, DataPoint q)
	{
		int numMeasuresCols = p.measures.length;
		double dist = 0;
		for (int i = 0; i < numMeasuresCols; i++)
			dist += Math.pow(Math.abs(p.measures[i] - q.measures[i]), expo);
		
		return Math.pow(dist, 1.0 / expo);
	}
	
	public static double distanceLNorm(int expo, DataPoint p, DataPoint q, int start, int end)
	{
		double dist = 0;
		for (int i = start; i <= end; i++)
			dist += Math.pow(Math.abs(p.measures[i] - q.measures[i]), expo);
		
		return Math.pow(dist, 1.0 / expo);
	}
}
