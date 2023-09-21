package multibinning.business;

import java.util.*;

import multibinning.data.*;

public class MVD 
{
	@SuppressWarnings("unchecked")
	public static ArrayList<MacroBin>[] discretizeData(DataMatrix dataMatrix) throws Exception
	{
		int numCols = dataMatrix.cols;
		ArrayList<MacroBin>[] ret = new ArrayList[numCols];
		
		// for each dimension, find the optimal binning strategy
		for (int dim = 0; dim < numCols; dim++)
			ret[dim] = findBinning(dim, dataMatrix);
		
		return ret;
	}
	
	// find binning strategy for each dimension
	public static ArrayList<MacroBin> findBinning(int curDim, DataMatrix dataMatrix) throws Exception
	{
		 ArrayList<MacroBin> ret = null;
		 
		 // get the total numbers of rows and columns
		 int numRows = dataMatrix.rows;
		 int numDims = dataMatrix.cols;
		 
		 // sort the data along the curDim dimension in ascending order
		 DataPoint curPoint = null;
		 SortedObjectComparator c = new SortedObjectComparator();
		 SortedObject[] tmp = new SortedObject[numRows];
		 for (int i = 0; i < numRows; i++)
		 {
			 curPoint = dataMatrix.data.get(i);
			 tmp[i] = new SortedObject(i, curPoint.measures[curDim]);
		 }
		 Arrays.sort(tmp, c);
		 
		 // divide data into equal-frequency micro bins
		 int initBinCount = Constants.INIT_BIN_COUNT;
		 ArrayList<MicroBin> initialMicroBins = findEqualFrequencyBinning(curDim, dataMatrix, initBinCount, tmp);
		 
		 // get the number of micro bins
		 int numMicroBins = initialMicroBins.size();
		 System.out.println("dim " + curDim + " --- num of initial bins = " + numMicroBins);
		 
		 // get the inter-distance of micro bins
		 int[] interBinsDiffBinary = new int[numMicroBins - 1];
		 for (int i = 0; i < interBinsDiffBinary.length; i++)
			 interBinsDiffBinary[i] = calculateBinsDiff(initialMicroBins.get(i), initialMicroBins.get(i + 1));
		 
		 // for each micro bin, create a macro bin containing it
		 ArrayList<MacroBin> candidateMacroBins = new ArrayList<MacroBin>();
		 MacroBin tmpMacroBin = null;
		 MicroBin tmpMicroBin = null;
		 for (int i = 0; i < numMicroBins; i++)
		 {
			 tmpMicroBin = initialMicroBins.get(i);
			 tmpMacroBin = new MacroBin(tmpMicroBin.lowerBound, tmpMicroBin.upperBound);
			 tmpMacroBin.numPoints = tmpMicroBin.allValues.size();
			 tmpMacroBin.microBinIDs.add(new Integer(i));
			 tmpMacroBin.numMDHs = 0;
			 tmpMacroBin.mdh = new int[0];
			 candidateMacroBins.add(tmpMacroBin);
		 }
		 
		 if (numMicroBins == 0)
			 throw new Exception("Dimension with no candidate bin!");
			 
		 // mine the optimal macro bins by means of dynamic programming
		 ret = greedyNormalBinning(initialMicroBins, numRows, candidateMacroBins.size(), initialMicroBins, interBinsDiffBinary, numDims - 1);
		 return ret;
	}
	
	public static ArrayList<MacroBin> greedyNormalBinning(ArrayList<MicroBin> candidateMicroBins, int N, int T, ArrayList<MicroBin> initialBins, int[] interBinsDiffBinary, int numDims) throws Exception
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		
		int numCanMicroBins = candidateMicroBins.size();
		MicroBin microBin1 = null;
		MicroBin microBin2 = null;
		MicroBin newMicroBin = null;
		int minSupport;
		int tmpSupport;
		int minIndex;
		
		// initialize fixed boundary status
		boolean[] fixedBoundary = new boolean[numCanMicroBins - 1];
		for (int i = 0; i < numCanMicroBins - 1; i++)
		{
			if (interBinsDiffBinary[i] == 1)
				fixedBoundary[i] = true;
			else
				fixedBoundary[i] = false;
		}
			
		boolean[] newFixedBoundary = null;
		int diff;
		while (numCanMicroBins > 1) // loop till the number of bins drops to 1
		{
			minSupport = Integer.MAX_VALUE;
			minIndex = -1;
			
			// loop through all the bins and pick two that have no fixed boundary and have minimum combined support
			for (int binIndex = 0; binIndex < numCanMicroBins - 1; binIndex++)
			{
				if (fixedBoundary[binIndex] == false)
				{
					microBin1 = candidateMicroBins.get(binIndex);
					microBin2 = candidateMicroBins.get(binIndex + 1);
					tmpSupport = microBin1.allValues.size() + microBin2.allValues.size();
					if (tmpSupport < minSupport)
					{
						minSupport = tmpSupport;
						minIndex = binIndex;
					} // end if
				} // end if
			} // end for
			
			// merge the two bins
			if (minIndex != -1)
			{
				// merge two bins
				newMicroBin = mergeMicroBins(candidateMicroBins.get(minIndex), candidateMicroBins.get(minIndex + 1));
				
				// update boundary
				newFixedBoundary = new boolean[numCanMicroBins - 2];
				for (int i = 0; i < minIndex; i++)
					newFixedBoundary[i] = fixedBoundary[i];
				for (int i = minIndex + 1; i < numCanMicroBins - 1; i++)
					newFixedBoundary[i - 1] = fixedBoundary[i];
				fixedBoundary = newFixedBoundary;
				
				// remove two old bins and add the newly created bin
				candidateMicroBins.remove(minIndex);
				candidateMicroBins.remove(minIndex);
				candidateMicroBins.add(minIndex, newMicroBin);
				
				// reduce the number of bins
				numCanMicroBins--;
				
				// adjust boundary
				if (minIndex < numCanMicroBins - 1)
				{
					diff = calculateBinsDiff(candidateMicroBins.get(minIndex), candidateMicroBins.get(minIndex + 1));
					if (diff == 1)
						fixedBoundary[minIndex] = true;
					else
						fixedBoundary[minIndex] = false;
				}
				if (minIndex > 0)
				{
					diff = calculateBinsDiff(candidateMicroBins.get(minIndex - 1), candidateMicroBins.get(minIndex));
					if (diff == 1)
						fixedBoundary[minIndex - 1] = true;
					else
						fixedBoundary[minIndex - 1] = false;
				}
			}
			else // terminate since no two bins without a fixed boundary
				break;
		}
		
		// for each micro bin, create a macro bin containing it
		numCanMicroBins = candidateMicroBins.size();
		MacroBin tmpMacroBin = null;
		MicroBin tmpMicroBin = null;
		for (int i = 0; i < numCanMicroBins; i++)
		{
			tmpMicroBin = candidateMicroBins.get(i);
			tmpMacroBin = new MacroBin(tmpMicroBin.lowerBound, tmpMicroBin.upperBound);
			tmpMacroBin.numPoints = tmpMicroBin.allValues.size();
			tmpMacroBin.microBinIDs.add(new Integer(i));
			tmpMacroBin.numMDHs = 0;
			tmpMacroBin.mdh = new int[0];
			ret.add(tmpMacroBin);
		}
		
		return ret;
	}
	
	public static MicroBin mergeMicroBins(MicroBin a, MicroBin b)
	{	
		// init the new macro bin
		MicroBin ret = new MicroBin(a.lowerBound, b.upperBound);
		
		// combine the distinct values
		int numDistinctValuesA = a.distinctValues.size();
		int numDistinctValuesB = b.distinctValues.size();
		for (int i = 0; i < numDistinctValuesA; i++)
		{
			ret.distinctValues.add(new Double(a.distinctValues.get(i).doubleValue()));
			ret.distinctValueSupports.add(new Integer(a.distinctValueSupports.get(i).intValue()));
		}
		for (int i = 0; i < numDistinctValuesB; i++)
		{
			ret.distinctValues.add(new Double(b.distinctValues.get(i).doubleValue()));
			ret.distinctValueSupports.add(new Integer(b.distinctValueSupports.get(i).intValue()));
		}
		
		// combine all values of two bins
		int numValuesA = a.allValues.size();
		int numValuesB = b.allValues.size();
		int numRowsA = numValuesA;
		int numRowsB = numValuesB;
		ret.dataMatrix = new DataMatrix(numRowsA + numRowsB, a.dataMatrix.cols, 0);
		int numCols = ret.dataMatrix.cols;
		DataPoint tmpPoint = null;
		for (int i = 0; i < numValuesA; i++)
		{
			ret.allValues.add(new Double(a.allValues.get(i).intValue()));
			ret.pointIDs.add(new Integer(a.pointIDs.get(i).intValue()));
			
			// add points to the data matrix
			tmpPoint = a.dataMatrix.data.get(i);
			ret.dataMatrix.data.add(tmpPoint);
		}
		for (int i = 0; i < numValuesB; i++)
		{
			ret.allValues.add(new Double(b.allValues.get(i).intValue()));
			ret.pointIDs.add(new Integer(b.pointIDs.get(i).intValue()));
			
			// add points to the data matrix
			tmpPoint = b.dataMatrix.data.get(i);
			ret.dataMatrix.data.add(tmpPoint);
		}
		
		for (int i = 0; i < numCols; i++)
			ret.dims.add(new Integer(i));
		
		return ret;
	}
	
	// divide a dimension into equal-frequency bins
	public static ArrayList<MicroBin> findEqualFrequencyBinning(int curDim, DataMatrix dataMatrix, int numBins, SortedObject[] tmp)
	{
		ArrayList<MicroBin> ret = new  ArrayList<MicroBin>();
		
		int numRows = dataMatrix.rows;
		int numCols = dataMatrix.cols;
		
		// number of distinct values per bin
		int binCount = (int)Math.floor(numRows * 1.0 / numBins);
		
		int lowerIndex;
		int upperIndex;
		MicroBin tmpBin = null;
		DataPoint curPoint = null;
		DataPoint tmpPoint = null;
		int curTmpCol;
		int curNumRows;
		int curNumCols;
		double curVal;
		int curCount;
		int curNumPoints;
		for (int i = 0; i < numBins; i++)
		{
			// lower and upper indices of distinct values
			lowerIndex = i * binCount;
			if (i < numBins - 1)
				upperIndex = (i + 1) * binCount - 1;
			else
				upperIndex = numRows - 1;
			curNumPoints = upperIndex - lowerIndex + 1;
			
			// get the bin's lower and upper bounds
			if (i == 0)
				tmpBin = new MicroBin(tmp[lowerIndex].value - 1, tmp[upperIndex].value);
			else
				tmpBin = new MicroBin(tmp[lowerIndex].value, tmp[upperIndex].value);
			
			// calculate the bin's data matrix, index matrix, means, and devs
			tmpBin.dataMatrix = new DataMatrix();
			tmpBin.dataMatrix.data = new ArrayList<DataPoint>();
			for (int j = 0; j < numCols - 1; j++)	// set up the bin's multivariate dimensions
				tmpBin.dims.add(new Integer(j));
			tmpBin.means = new double[numCols - 1];
			tmpBin.devs = new double[numCols - 1];
			
			// get current value and add the corresponding points to the bin
			curVal = tmp[lowerIndex].value;
			curPoint = dataMatrix.data.get(tmp[lowerIndex].index);
			tmpPoint = new DataPoint(numCols - 1, 0);
			curTmpCol = -1;
			for (int l = 0; l < numCols; l++)
			{
				if (l != curDim)
				{
					curTmpCol++;
					tmpPoint.measures[curTmpCol] = curPoint.measures[l];
					tmpBin.means[curTmpCol] = tmpPoint.measures[curTmpCol];
				}
			}
			tmpBin.dataMatrix.data.add(tmpPoint);
			tmpBin.pointIDs.add(new Integer(tmp[lowerIndex].index));
			curCount = 1;
			
			// loop through the remaining value
			for (int j = lowerIndex + 1; j <= upperIndex; j++)
			{	
				if (tmp[j].value == curVal)
				{
					curCount++;
					if (j == curNumPoints - 1)
						tmpBin.distinctValueSupports.add(new Integer(curCount));
				}
				else
				{
					tmpBin.distinctValueSupports.add(new Integer(curCount));
					curVal = tmp[j].value;
					tmpBin.distinctValues.add(new Double(curVal));
					curCount = 1;
					if (j == curNumPoints - 1)
						 tmpBin.distinctValueSupports.add(new Integer(curCount));
				}
				
				// add the respective data point to the bin's data matrix
				curPoint = dataMatrix.data.get(tmp[j].index);
				tmpPoint = new DataPoint(numCols - 1, 0);
				curTmpCol = -1;
				for (int l = 0; l < numCols; l++)
				{
					if (l != curDim)
					{
						curTmpCol++;
						tmpPoint.measures[curTmpCol] = curPoint.measures[l];
						tmpBin.means[curTmpCol] += tmpPoint.measures[curTmpCol];
					}
				}
				tmpBin.dataMatrix.data.add(tmpPoint);
				tmpBin.pointIDs.add(new Integer(tmp[j].index));
			}
			
			// calculate the bin's index matrix
			tmpBin.dataMatrix.rows = tmpBin.pointIDs.size();
			tmpBin.dataMatrix.cols = numCols - 1;
			tmpBin.indexMatrix = new IndexMatrix();
			tmpBin.indexMatrix.rows = tmpBin.dataMatrix.rows;
			tmpBin.indexMatrix.cols = tmpBin.dataMatrix.cols;
			tmpBin.indexMatrix.data = new int[tmpBin.dataMatrix.rows][tmpBin.dataMatrix.cols];
			DataProcess.calculateIndices(tmpBin.dataMatrix, tmpBin.indexMatrix);
			
			// calculate the bin's means and deviations
			curNumRows = tmpBin.dataMatrix.rows;
			curNumCols = tmpBin.dataMatrix.cols;
			for (int j = 0; j < curNumCols; j++)
			{
				tmpBin.means[j] = tmpBin.means[j] / curNumRows;
				if (curNumRows == 1)
					tmpBin.devs[j] = 0;
				else
				{
					for (int k = 0; k < curNumRows; k++)
					{
						curPoint = tmpBin.dataMatrix.data.get(k);
						tmpBin.devs[j] += (curPoint.measures[j] - tmpBin.means[j]) * (curPoint.measures[j] - tmpBin.means[j]) / (curNumRows - 1);
					}
					tmpBin.devs[j] = Math.sqrt(tmpBin.devs[j]);
				}
			}
			
			// add the bin to the list of bins
			ret.add(tmpBin);
		}
		
		return ret;
	}
	
	// calculate the MCKL distance of two micro bins
	public static int calculateBinsDiff(MicroBin a, MicroBin b) throws Exception
	{
		int numPointsA = a.dataMatrix.rows;
		int numPointsB = b.dataMatrix.rows;
		int cols = a.dataMatrix.cols;
		int expo = Constants.LP_NORM;
		DataPoint curPointA = null;
		DataPoint curPointB = null;
		double dist;
		
		int supportA;
		int supportB;
		double delta = Constants.GAMMA * Math.min(numPointsA, numPointsB);
		for (int i = 0; i < numPointsA; i++)
		{
			supportA = 0;
			curPointA = a.dataMatrix.data.get(i);
			
			for (int start = 0; start < cols; start++)
				for (int end = start; end < cols; end++)
				{
					for (int j = 0; j < numPointsA; j++)
					{
						curPointB = a.dataMatrix.data.get(j);
						dist = DataPoint.distanceLNorm(expo, curPointA, curPointB, start, end);
						if (dist == 0)
							supportA++;
					}
					
					supportB = 0;
					for (int j = 0; j < numPointsB; j++)
					{
						curPointB = b.dataMatrix.data.get(j);
						dist = DataPoint.distanceLNorm(expo, curPointA, curPointB, start, end);
						if (dist == 0)
							supportB++;
					}
					
					if (Math.abs(supportA - supportB) > delta)
						return 1;
				}
		}
		
		for (int i = 0; i < numPointsB; i++)
		{
			supportB = 0;
			curPointB = b.dataMatrix.data.get(i);
			
			for (int start = 0; start < cols; start++)
				for (int end = start; end < cols; end++)
				{
					for (int j = 0; j < numPointsB; j++)
					{
						curPointA = b.dataMatrix.data.get(j);
						dist = DataPoint.distanceLNorm(expo, curPointB, curPointA, start, end);
						if (dist == 0)
							supportB++;
					}
					
					supportA = 0;
					for (int j = 0; j < numPointsA; j++)
					{
						curPointA = a.dataMatrix.data.get(j);
						dist = DataPoint.distanceLNorm(expo, curPointB, curPointA, start, end);
						if (dist == 0)
							supportA++;
					}
					
					if (Math.abs(supportB - supportA) > delta)
						return 1;
				}
		}
		
		return 0;
	}
}
