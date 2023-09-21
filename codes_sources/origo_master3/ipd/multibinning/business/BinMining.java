package multibinning.business;

//import java.lang.reflect.Array;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;
import multibinning.data.*;

class DataWithIndex
{
	public DataMatrix data;
	public IndexMatrix index;
	public ArrayList<double[][]> sortedData;
}

public class BinMining 
{
	@SuppressWarnings("unchecked")
	public static ArrayList<MacroBin>[] discretizeData(DataMatrix dataMatrix) throws Exception
	{
		int numCols = dataMatrix.cols;
		ArrayList<MacroBin>[] ret = new ArrayList[numCols];
		
		// for each dimension, find the optimal binning strategy
		int method = Constants.METHOD;
		if (method == Constants.MB_MDL_DP || method == Constants.MB_MDL_GD  || method == Constants.MB_MDL_EGD || method == Constants.MB_NM_GD || Constants.METHOD == Constants.DP_MEAN || Constants.METHOD == Constants.DP_MIN || Constants.METHOD == Constants.DP_MAX)
		{
			Constants.INIT_BIN_COUNT = (int)Math.ceil(Math.sqrt(dataMatrix.rows));
			for (int dim = 0; dim < numCols; dim++)
				ret[dim] = findBinning(dim, dataMatrix, method);
		}
		else if (method == Constants.EF)
		{
			for (int dim = 0; dim < numCols; dim++)
				ret[dim] = findEqualFrequencyBinning(dim, dataMatrix);
		}
		else if (method == Constants.EW)
		{
			for (int dim = 0; dim < numCols; dim++)
				ret[dim] = findEqualWidthBinning(dim, dataMatrix);
		}
		
		return ret;
	}
	
	public static ArrayList<MacroBin> findEqualFrequencyBinning(int curDim, DataMatrix dataMatrix) throws Exception
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		 
		// get the total numbers of rows and columns
		int numRows = dataMatrix.rows;
		 
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
		 
		// get distinct values, the supports of distinct values, and IDs of member points of distinct values
		/*int curCount;
		ArrayList<Double> distinctValues = new ArrayList<Double>();
		ArrayList<Integer> distinctValueSupports = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> memberPointIDs = new ArrayList<ArrayList<Integer>>();
		double curVal = tmp[0].value;
		distinctValues.add(new Double(curVal));
		curCount = 1;
		ArrayList<Integer> tmpPointIDs = new ArrayList<Integer>();
		tmpPointIDs.add(new Integer(tmp[0].index));
		for (int i = 1; i < numRows; i++)
		{
			if (tmp[i].value == curVal)
			{
				curCount++;
				tmpPointIDs.add(new Integer(tmp[i].index));
				if (i == numRows - 1)
				{
					distinctValueSupports.add(new Integer(curCount));
					memberPointIDs.add(tmpPointIDs);
				}
			}
			else
			{
				distinctValueSupports.add(new Integer(curCount));
				memberPointIDs.add(tmpPointIDs);
				curVal = tmp[i].value;
				distinctValues.add(new Double(curVal));
				curCount = 1;
				tmpPointIDs = new ArrayList<Integer>();
				tmpPointIDs.add(new Integer(tmp[i].index));
				if (i == numRows - 1)
				{
					distinctValueSupports.add(new Integer(curCount));
					memberPointIDs.add(tmpPointIDs);
				}
			}
		}*/
		 
		// divide data into equal-frequency micro bins
		int initBinCount = Constants.INIT_BIN_COUNT;
		ArrayList<MicroBin> initialMicroBins = findEqualFrequencyBinning(curDim, dataMatrix, initBinCount, tmp);
		 
		// for each micro bin, create a macro bin containing it
		int numMicroBins = initialMicroBins.size();
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
			ret.add(tmpMacroBin);
		}
		
		return ret;
	}
	
	public static ArrayList<MacroBin> findEqualWidthBinning(int curDim, DataMatrix dataMatrix) throws Exception
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		 
		// get the total numbers of rows and columns
		int numRows = dataMatrix.rows;
		 
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
		 
		// get distinct values, the supports of distinct values, and IDs of member points of distinct values
		int curCount;
		ArrayList<Double> distinctValues = new ArrayList<Double>();
		ArrayList<Integer> distinctValueSupports = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> memberPointIDs = new ArrayList<ArrayList<Integer>>();
		double curVal = tmp[0].value;
		distinctValues.add(new Double(curVal));
		curCount = 1;
		ArrayList<Integer> tmpPointIDs = new ArrayList<Integer>();
		tmpPointIDs.add(new Integer(tmp[0].index));
		for (int i = 1; i < numRows; i++)
		{
			if (tmp[i].value == curVal)
			{
				curCount++;
				tmpPointIDs.add(new Integer(tmp[i].index));
				if (i == numRows - 1)
				{
					distinctValueSupports.add(new Integer(curCount));
					memberPointIDs.add(tmpPointIDs);
				}
			}
			else
			{
				distinctValueSupports.add(new Integer(curCount));
				memberPointIDs.add(tmpPointIDs);
				curVal = tmp[i].value;
				distinctValues.add(new Double(curVal));
				curCount = 1;
				tmpPointIDs = new ArrayList<Integer>();
				tmpPointIDs.add(new Integer(tmp[i].index));
				if (i == numRows - 1)
				{
					distinctValueSupports.add(new Integer(curCount));
					memberPointIDs.add(tmpPointIDs);
				}
			}
		}
		 
		// divide data into equal-frequency micro bins
		ArrayList<MicroBin> initialMicroBins = findEqualWidthBinning(curDim, dataMatrix, Constants.INIT_BIN_COUNT, distinctValues, distinctValueSupports, memberPointIDs);
		 
		// for each micro bin, create a macro bin containing it
		int numMicroBins = initialMicroBins.size();
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
			ret.add(tmpMacroBin);
		}
		
		return ret;
	}
	
	// find binning strategy for each dimension
	public static ArrayList<MacroBin> findBinning(int curDim, DataMatrix dataMatrix, int method) throws Exception
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
		 //System.out.println(curDim + " --- num of distinct values = " + distinctValues.size());
		 int initBinCount = Constants.INIT_BIN_COUNT;
		 //ArrayList<MicroBin> initialMicroBins = findEqualFrequencyBinning(curDim, dataMatrix, initBinCount, tmp);
		 ArrayList<MicroBin> initialMicroBins = findEqualFrequencyBinningNew(curDim, Constants.ODIMS, dataMatrix, initBinCount, tmp);
		 //ArrayList<MicroBin> initialMicroBins = findEqualWidthBinning(curDim, dataMatrix, initBinCount, distinctValues, distinctValueSupports, memberPointIDs);
		 
		 // get the number of micro bins
		 int numMicroBins = initialMicroBins.size();
		 System.out.println("dim " + curDim + " --- num of initial bins = " + numMicroBins);
		 
		 // get the inter-distance of micro bins
		 double[] interBinsDiff = new double[numMicroBins - 1];
		 double[] tmpBinsDiff = new double[numMicroBins - 1];
		 int[] interBinsDiffBinary = new int[numMicroBins - 1];
		 //double meanDiff = 0;
		 //double devDiff = 0;
		 for (int i = 0; i < interBinsDiff.length; i++)
		 {
			 // compute bin difference
			 if (Constants.METHOD != Constants.DP_MEAN && Constants.METHOD != Constants.DP_MIN && Constants.METHOD != Constants.DP_MAX)
				 interBinsDiff[i] = calculateBinsDiffNewNew(initialMicroBins.get(i), initialMicroBins.get(i + 1), Constants.MAX_VAL);
			 else if (Constants.METHOD == Constants.DP_MEAN)
				 interBinsDiff[i] = calculateBinsDiffByMean(initialMicroBins.get(i), initialMicroBins.get(i + 1));
			 else if (Constants.METHOD == Constants.DP_MIN)
				 interBinsDiff[i] = calculateBinsDiffByMinDist(initialMicroBins.get(i), initialMicroBins.get(i + 1));
			 else
				 interBinsDiff[i] = calculateBinsDiffByMaxDist(initialMicroBins.get(i), initialMicroBins.get(i + 1));
			 
			 // copy bin difference
			 tmpBinsDiff[i] = interBinsDiff[i];
			 //meanDiff += interBinsDiff[i];
			 
			 System.out.println("end MCKL inter-bin " + i);
		 }
		 Arrays.sort(tmpBinsDiff);
//		 BufferedWriter writerIDIST = new BufferedWriter(new FileWriter(Constants.FILE_IDIST, true));
//		 for (int i = 0; i < interBinsDiff.length; i++)
//			 writerIDIST.write(Double.toString(tmpBinsDiff[i]) + ",");
//		 writerIDIST.newLine();
//		 writerIDIST.flush();
//		 writerIDIST.close();
		 double thyThreshold = tmpBinsDiff[(int)Math.ceil(tmpBinsDiff.length / Constants.QUANTILE) - 1];
		 
		 // binarize the bin difference
		 for (int i = 0; i < interBinsDiff.length; i++)
		 {
			 if (interBinsDiff[i] > thyThreshold)
				 interBinsDiffBinary[i] = 1;
			 else
				 interBinsDiffBinary[i] = 0;
		 }
		 
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
		 if (method == Constants.MB_MDL_DP || Constants.METHOD == Constants.DP_MEAN || Constants.METHOD == Constants.DP_MIN || Constants.METHOD == Constants.DP_MAX)
			 ret = dynamicProgrammingMDLBinning(candidateMacroBins, numRows, candidateMacroBins.size(), initialMicroBins, interBinsDiffBinary);
		 else if (method == Constants.MB_MDL_GD || method == Constants.MB_MDL_EGD)
		 {
			 if (numMicroBins > 2)
				 ret = greedyMDLBinning(candidateMacroBins, numRows, candidateMacroBins.size(), initialMicroBins, interBinsDiffBinary);
			 else
				 ret = candidateMacroBins;
		 }
		 else
			 ret = greedyNormalBinning(initialMicroBins, numRows, candidateMacroBins.size(), initialMicroBins, interBinsDiff, interBinsDiffBinary, numDims - 1);
		 return ret;
	}
	
	public static ArrayList<MacroBin> greedyMDLBinning(ArrayList<MacroBin> candidateMacroBins, int N, int T, ArrayList<MicroBin> initialBins, int[] interBinsDiffBinary) throws Exception
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		
		int numCanMacroBins = T;
		//System.out.println(numCanMacroBins);
		MacroBin macroBin1 = null;
		MacroBin macroBin2 = null;
		MacroBin newMacroBin = null;
		double minOverallCost = Double.MAX_VALUE;
		double minCost;
		double tmpCost;
		int optimalBinIndex;
		MacroBin optimalNewBin = null;
		while (numCanMacroBins > 1)
		{
			//System.out.println(numCanMacroBins);
			optimalBinIndex = -1;
			optimalNewBin = null;
			minCost = Double.MAX_VALUE;
			//System.out.println(numCanMacroBins);
			for (int binIndex = 0; binIndex < numCanMacroBins - 1; binIndex++)
			{
				// pick the two consecutive bins
				macroBin1 = candidateMacroBins.get(binIndex);
				macroBin2 = candidateMacroBins.get(binIndex + 1);
				
				// merge these two bins
				newMacroBin = mergeMacroBins(macroBin1, macroBin2, initialBins, interBinsDiffBinary);
				
				// remove the two old bins and insert the new bins
				candidateMacroBins.remove(binIndex);
				candidateMacroBins.remove(binIndex);
				candidateMacroBins.add(binIndex, newMacroBin);
				
				// compute the coding cost of the new binning
				tmpCost = computeEncodingCost(candidateMacroBins, N, T, initialBins);
				//System.out.println(tmpCost);
				
				// if the coding cost is less than the minimum coding cost seen so far, then record the respective binIndex
				if (tmpCost < minCost)
				{
					minCost = tmpCost;
					optimalBinIndex = binIndex;
					optimalNewBin = newMacroBin;
				}
				
				// add back the two old bins
				candidateMacroBins.remove(binIndex);
				candidateMacroBins.add(binIndex, macroBin1);
				candidateMacroBins.add(binIndex + 1, macroBin2);
			}
			
			// update the binning with the best bin (hence greedy)
			candidateMacroBins.remove(optimalBinIndex);
			candidateMacroBins.remove(optimalBinIndex);
			candidateMacroBins.add(optimalBinIndex, optimalNewBin);
			numCanMacroBins--;
			if (minCost < minOverallCost)
			{
				minOverallCost = minCost;
				ret.clear();
				for (int binIndex = 0; binIndex < numCanMacroBins; binIndex++)
					ret.add(cloneBin(candidateMacroBins.get(binIndex)));
			}
			else if (Constants.METHOD == Constants.MB_MDL_EGD)
			{
				break;
			} // end if
			//System.out.println(minOverallCost);
		} // end while
		
		System.out.println("Cost = " + minOverallCost);
		System.out.println("-----------");
		
		return ret;
	}
	
	public static ArrayList<MacroBin> greedyNormalBinning(ArrayList<MicroBin> candidateMicroBins, int N, int T, ArrayList<MicroBin> initialBins, double[] interBinsDiff, int[] interBinsDiffBinary, int numDims) throws Exception
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		
		// threshold for micro bins' similarity
		int numCanMicroBins = candidateMicroBins.size();
		MicroBin newMicroBin = null;
		double minDiff;
		double tmpDiff;
		int minIndex;
		
		// initialize fixed boundary status and inter-bin distances
		boolean[] fixedBoundary = new boolean[numCanMicroBins - 1];
		double[] binsDiff = new double[numCanMicroBins - 1];
		for (int i = 0; i < numCanMicroBins - 1; i++)
		{
			if (interBinsDiffBinary[i] == 1)
				fixedBoundary[i] = true;
			else
				fixedBoundary[i] = false;
			
			binsDiff[i] = interBinsDiff[i];
		}
		
		boolean[] newFixedBoundary = null;
		double[] newBinsDiff = null;
		double threshold = quantileThreshold(binsDiff);;
		while (numCanMicroBins > 1) // loop till the number of bins drops to 1
		{
			minDiff = Double.MAX_VALUE;
			minIndex = -1;
			
			// loop through all the bins and pick two that have no fixed boundary and are most similar
			for (int binIndex = 0; binIndex < numCanMicroBins - 1; binIndex++)
			{
				if (fixedBoundary[binIndex] == false)
				{
					tmpDiff = binsDiff[binIndex];
					if (tmpDiff < minDiff)
					{
						minDiff = tmpDiff;
						minIndex = binIndex;
					}
				}
			}
			
			// if the difference is less than threshold, then merge the two bins
			if (minIndex != -1 && minDiff <= threshold)
			{
				// merge two bins
				newMicroBin = mergeMicroBins(candidateMicroBins.get(minIndex), candidateMicroBins.get(minIndex + 1));
				
				// update boundary
				newFixedBoundary = new boolean[numCanMicroBins - 2];
				newBinsDiff = new double[numCanMicroBins - 2];
				for (int i = 0; i < minIndex; i++)
				{
					newFixedBoundary[i] = fixedBoundary[i];
					newBinsDiff[i] = binsDiff[i];
				}
				for (int i = minIndex + 1; i < numCanMicroBins - 1; i++)
				{
					newFixedBoundary[i - 1] = fixedBoundary[i];
					newBinsDiff[i - 1] = binsDiff[i];
				}
				fixedBoundary = newFixedBoundary;
				binsDiff = newBinsDiff;
				
				// remove two old bins and add the newly created bin
				candidateMicroBins.remove(minIndex);
				candidateMicroBins.remove(minIndex);
				candidateMicroBins.add(minIndex, newMicroBin);
				
				// reduce the number of bins
				numCanMicroBins--;
				
				// adjust bin distances
				if (minIndex < numCanMicroBins - 1)
					binsDiff[minIndex] = calculateBinsDiffNew(candidateMicroBins.get(minIndex), candidateMicroBins.get(minIndex + 1));
				if (minIndex > 0)
					binsDiff[minIndex - 1] = calculateBinsDiffNew(candidateMicroBins.get(minIndex - 1), candidateMicroBins.get(minIndex));
				
				// adjust fixed boundary
				threshold = quantileThreshold(binsDiff);
				if (minIndex < numCanMicroBins - 1)
				{
					if (binsDiff[minIndex] > threshold)
						fixedBoundary[minIndex] = true;
					else
						fixedBoundary[minIndex] = false;
				}
				
				if (minIndex > 0)
				{
					if (binsDiff[minIndex - 1] > threshold)
						fixedBoundary[minIndex - 1] = true;
					else
						fixedBoundary[minIndex - 1] = false;
				}
			}
			else // terminate since no two bins without a fixed boundary is similar enough for merging
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
	
	@SuppressWarnings("unchecked")
	public static ArrayList<MacroBin> dynamicProgrammingMDLBinning(ArrayList<MacroBin> candidateMacroBins, double N, double T, ArrayList<MicroBin> initialBins, int[] interBinsDiffBinary) throws Exception
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		
		double logBase = Math.log(Constants.LOG_BASE);
		
		// create the table for dynamic programming to store already solved sub-problems
		int numTotalMacroBins = candidateMacroBins.size();
		ArrayList<MacroBin>[][] dynamicBinnings = new ArrayList[numTotalMacroBins + 1][numTotalMacroBins + 1];
		double[][] dynamicBinningCosts = new double[numTotalMacroBins + 1][numTotalMacroBins + 1];
		int[][] dynamicNumPoints = new int[numTotalMacroBins + 1][numTotalMacroBins + 1];
		
		// compute binning of size 1
		MacroBin tmpMergeMacroBin = null;
		MacroBin tmpMergeMacroBinRight = null;
		double tmpCost;
		ArrayList<MacroBin> tmpMacroBins = new ArrayList<MacroBin>();
		tmpMergeMacroBin = cloneBin(candidateMacroBins.get(0));
		for (int maxIndex = 1; maxIndex < numTotalMacroBins; maxIndex++)
			tmpMergeMacroBin = mergeMacroBins(tmpMergeMacroBin, candidateMacroBins.get(maxIndex), initialBins, interBinsDiffBinary);
		tmpMacroBins.add(tmpMergeMacroBin);
		dynamicBinningCosts[numTotalMacroBins][1] = computeEncodingCost(tmpMacroBins, tmpMergeMacroBin.numPoints, tmpMergeMacroBin.microBinIDs.size(), initialBins);
		dynamicBinnings[numTotalMacroBins][1] = new ArrayList<MacroBin>();
		dynamicBinnings[numTotalMacroBins][1].add(tmpMergeMacroBin);
		
		// find optimal binning of size 2
		for (int numFirstMacroBins = 2; numFirstMacroBins <= numTotalMacroBins; numFirstMacroBins++)
		{
			dynamicBinningCosts[numFirstMacroBins][2] = Double.MAX_VALUE;
			
			// loop through each possible bin position
			// maxIndex contains the minimum bin index of the right macro bin
			for (int maxIndex = 1; maxIndex < numFirstMacroBins; maxIndex++)
			{
				// clear the temporary binning
				tmpMacroBins.clear();
				
				// merge all macro bins from position 0 to (maxIndex - 1)
				tmpMergeMacroBin = cloneBin(candidateMacroBins.get(0));
				for (int i = 1; i < maxIndex; i++)
					tmpMergeMacroBin = mergeMacroBins(tmpMergeMacroBin, candidateMacroBins.get(i), initialBins, interBinsDiffBinary);
				
				// merge all macro bins from position maxIndex to (numFirstMacroBins - 1)
				tmpMergeMacroBinRight = cloneBin(candidateMacroBins.get(maxIndex));
				for (int i = maxIndex + 1; i < numFirstMacroBins; i++)
					tmpMergeMacroBinRight = mergeMacroBins(tmpMergeMacroBinRight, candidateMacroBins.get(i), initialBins, interBinsDiffBinary);
				
				// add two new macro bins to the temporary binning
				tmpMacroBins.add(tmpMergeMacroBin);
				tmpMacroBins.add(tmpMergeMacroBinRight);
				
				// compute the coding cost of the temporary binning
				tmpCost = computeEncodingCost(tmpMacroBins, tmpMergeMacroBin.numPoints + tmpMergeMacroBinRight.numPoints, tmpMergeMacroBin.microBinIDs.size() + tmpMergeMacroBinRight.microBinIDs.size(), initialBins);
				
				// if the cost is less than the current cost, then update the cost, the binning, and the total number of points of the binning
				if (tmpCost < dynamicBinningCosts[numFirstMacroBins][2])
				{
					dynamicBinningCosts[numFirstMacroBins][2] = tmpCost;
					//dynamicNumPoints[numFirstMacroBins][2] = tmpMergeMacroBin.numPoints + tmpMergeMacroBinRight.numPoints;
					dynamicNumPoints[numFirstMacroBins][2] = tmpMergeMacroBin.microBinIDs.size() + tmpMergeMacroBinRight.microBinIDs.size();
					if (dynamicBinnings[numFirstMacroBins][2] == null)
						dynamicBinnings[numFirstMacroBins][2] = new ArrayList<MacroBin>();
					else
						dynamicBinnings[numFirstMacroBins][2].clear();
					dynamicBinnings[numFirstMacroBins][2].add(tmpMergeMacroBin);
					dynamicBinnings[numFirstMacroBins][2].add(tmpMergeMacroBinRight);
				}
			}
		}
		
		// build the rest of the table
		// for each (numMacroBins, numFirstMacroBins) pairs, find the binning that yields the smallest cost
		//int tmpNumPoints;
		int tmpNumMicroBins;
		int optimalPosition;
		int tmpTotalNumMicroBins;
		MacroBin optimalNewBin = null;
		for (int numMacroBins = 3; numMacroBins <= numTotalMacroBins; numMacroBins++)
		{
			for (int numFirstMacroBins = numMacroBins; numFirstMacroBins <= numTotalMacroBins; numFirstMacroBins++)
			{
				dynamicBinningCosts[numFirstMacroBins][numMacroBins] = Double.MAX_VALUE;
				optimalPosition = -1;
				optimalNewBin = null;
				
				// loop through each possible bin position
				for (int maxIndex = numMacroBins - 1; maxIndex < numFirstMacroBins; maxIndex++)
				{
					// retrieve the stored results
					tmpCost = dynamicBinningCosts[maxIndex][numMacroBins - 1];
					
					// merge all the macro bins from position maxIndex to (numFirstMacroBins - 1)
					tmpMergeMacroBinRight = cloneBin(candidateMacroBins.get(maxIndex));
					for (int i = maxIndex + 1; i < numFirstMacroBins; i++)
						tmpMergeMacroBinRight = mergeMacroBins(tmpMergeMacroBinRight, candidateMacroBins.get(i), initialBins, interBinsDiffBinary);
					
					// compute the new binning's cost based on stored information
					//tmpNumPoints = dynamicNumPoints[maxIndex][numMacroBins - 1] + tmpMergeMacroBinRight.numPoints;
					tmpNumMicroBins = tmpMergeMacroBinRight.microBinIDs.size();
					tmpTotalNumMicroBins = dynamicNumPoints[maxIndex][numMacroBins - 1] + tmpNumMicroBins;
					
					// increase LcpM
					tmpCost += LInteger(numMacroBins - 1) + logOfCombination(numFirstMacroBins - 1, numMacroBins - 1);
					tmpCost -= LInteger(numMacroBins - 2) +  logOfCombination(maxIndex - 1, numMacroBins - 2);
					
					// increase coding cost of the number of micro bins
					//tmpCost += LInteger(numFirstMacroBins) - LInteger(maxIndex);
					
					// increase LbidM
					//if (tmpMergeMacroBinRight.numPoints != 0)
					//	tmpCost -=  Math.log(tmpMergeMacroBinRight.numPoints * 1.0 / tmpNumPoints) / logBase;
					if (tmpNumMicroBins != 0)
						tmpCost -=  Math.log(tmpNumMicroBins * 1.0 / tmpTotalNumMicroBins) / logBase;
					
					// increase LbidM
					//if (tmpNumPoints - tmpMergeMacroBinRight.numPoints != 0)
					//	tmpCost -=  (numMacroBins - 1) * Math.log((tmpNumPoints - tmpMergeMacroBinRight.numPoints) * 1.0 / tmpNumPoints) / logBase;
					if (dynamicNumPoints[maxIndex][numMacroBins - 1] != 0)
						tmpCost -=  (numMacroBins - 1) * Math.log(dynamicNumPoints[maxIndex][numMacroBins - 1] * 1.0 / tmpTotalNumMicroBins) / logBase;
					
					// increase LicM
					if (Constants.USE_CE == true)
						tmpCost += computeBinInformationContentWithCE(tmpMergeMacroBinRight, initialBins);
					else
						tmpCost += computeBinInformationContentWithoutCE(tmpMergeMacroBinRight, initialBins);
					
					// increase L(DS | M)
					//if (tmpMergeMacroBinRight.numPoints != 0)
						//tmpCost -= tmpMergeMacroBinRight.numPoints * Math.log(tmpMergeMacroBinRight.numPoints * 1.0 / tmpNumPoints) / logBase;
					if (tmpNumMicroBins != 0)
						tmpCost -= tmpNumMicroBins * Math.log(tmpNumMicroBins * 1.0 / tmpTotalNumMicroBins) / logBase;
					
					// increase L(DS | M)
					//if (tmpNumPoints - tmpMergeMacroBinRight.numPoints != 0)
					//	tmpCost -= (tmpNumPoints - tmpMergeMacroBinRight.numPoints) * Math.log((tmpNumPoints - tmpMergeMacroBinRight.numPoints) * 1.0 / tmpNumPoints) / logBase;
					if (dynamicNumPoints[maxIndex][numMacroBins - 1] != 0)
						tmpCost -= dynamicNumPoints[maxIndex][numMacroBins - 1] * Math.log(dynamicNumPoints[maxIndex][numMacroBins - 1] * 1.0 / tmpTotalNumMicroBins) / logBase;
					
					// increase L(DS | M)
					tmpCost += LInteger(tmpNumMicroBins);
					
					// increase LmdhM
					if (tmpMergeMacroBinRight.numMDHs != 0)
						tmpCost += LInteger(tmpMergeMacroBinRight.numMDHs);
					if (tmpNumMicroBins - 1 != 0)
						tmpCost += tmpMergeMacroBinRight.numMDHs * Math.log(tmpNumMicroBins - 1) / logBase;
					/*if (tmpNumMicroBins - 1 != 0)
					{
						if (tmpMergeMacroBinRight.numMDHs != 0)
							tmpCost -= tmpMergeMacroBinRight.numMDHs * Math.log(tmpMergeMacroBinRight.numMDHs * 1.0 / (tmpNumMicroBins - 1)) / logBase;
						if (tmpNumMicroBins - 1 - tmpMergeMacroBinRight.numMDHs != 0)
							tmpCost -= (tmpNumMicroBins - 1 - tmpMergeMacroBinRight.numMDHs) * Math.log((tmpNumMicroBins - 1 - tmpMergeMacroBinRight.numMDHs) * 1.0 / (tmpNumMicroBins - 1)) / logBase;
					}*/
					
					//System.out.println(tmpCost);
					// if the new cost is less than the current cost, then update the binning
					if (tmpCost < dynamicBinningCosts[numFirstMacroBins][numMacroBins])
					{
						dynamicBinningCosts[numFirstMacroBins][numMacroBins] = tmpCost;
						dynamicNumPoints[numFirstMacroBins][numMacroBins] = tmpTotalNumMicroBins;
						optimalPosition = maxIndex;
						optimalNewBin = tmpMergeMacroBinRight;
						//System.out.println("here");
					} // end if
				} // end for
				
				// update the table's entry with the binning yielding the smallest cost
				dynamicBinnings[numFirstMacroBins][numMacroBins] = new ArrayList<MacroBin>();
				for (int i = 0; i < numMacroBins - 1; i++)
				{
					//System.out.println(optimalPosition);
					tmpMergeMacroBin = cloneBin(dynamicBinnings[optimalPosition][numMacroBins - 1].get(i));
					dynamicBinnings[numFirstMacroBins][numMacroBins].add(tmpMergeMacroBin);
					//dynamicNumPoints[numFirstMacroBins][numMacroBins] += tmpMergeMacroBin.numPoints;
				}
				dynamicBinnings[numFirstMacroBins][numMacroBins].add(optimalNewBin);
				//dynamicNumPoints[numFirstMacroBins][numMacroBins] += tmpMergeMacroBinRight.numPoints;
			} // end for
			
			System.out.println("numMacroBins " + numMacroBins);
		} // end for
		
		double minCost = Double.MAX_VALUE;
		int minNumMacroBins = -1;
		for (int numMacroBins = 1; numMacroBins <= numTotalMacroBins; numMacroBins++)
		{
			if (dynamicBinningCosts[numTotalMacroBins][numMacroBins] < minCost)
			{
				minCost = dynamicBinningCosts[numTotalMacroBins][numMacroBins];
				minNumMacroBins = numMacroBins;
			}
		}
		
		MacroBin curMacroBin = null;
		for (int i = 0; i < minNumMacroBins; i++)
		{
			curMacroBin = dynamicBinnings[numTotalMacroBins][minNumMacroBins].get(i);
			curMacroBin.pointIDs.clear();
			curMacroBin.microBinIDs.clear();
			curMacroBin.means = null;
			curMacroBin.devs = null;
			curMacroBin.mdh = null;
			//if (curMacroBin.dataMatrix.data != null)
			//	curMacroBin.dataMatrix.data.clear();
			//ret.add(dynamicBinnings[numTotalMacroBins][minNumMacroBins].get(i));
			ret.add(curMacroBin);
		}
		
		System.out.println("Cost = " + minCost);
		System.out.println("-----------");
		return ret;
	}
	
	// clone a macro bin
	public static MacroBin cloneBin(MacroBin a)
	{
		// init the new macro bin
		MacroBin ret = new MacroBin(a.lowerBound, a.upperBound);
		
		// get the total number of data points
		ret.numPoints = a.numPoints;
		
		// get the IDs of micro bins
		int numMicroBinsA = a.microBinIDs.size();
		for (int i = 0; i < numMicroBinsA; i++)
			ret.microBinIDs.add(new Integer(a.microBinIDs.get(i).intValue()));
		
		// adjust MDH accordingly
		ret.numMDHs = a.numMDHs;
		ret.mdh = new int[a.mdh.length];
		for (int i = 0; i < a.mdh.length; i++)
			ret.mdh[i] = a.mdh[i];
		
		return ret;
	}
	
	// merge two macro bins
	public static MacroBin mergeMacroBins(MacroBin a, MacroBin b, ArrayList<MicroBin> initialBins, int[] interBinsDiffBinary) throws Exception
	{	
		// init the new macro bin
		if (a.lowerBound > b.upperBound)
			throw new Exception("Invalid bin merge");
		
		MacroBin ret = new MacroBin(a.lowerBound, b.upperBound);
		
		// get the total number of data points
		ret.numPoints = a.numPoints + b.numPoints;
		
		// get the IDs of micro bins
		int numMicroBinsA = a.microBinIDs.size();
		int numMicroBinsB = b.microBinIDs.size();
		for (int i = 0; i < numMicroBinsA; i++)
			ret.microBinIDs.add(new Integer(a.microBinIDs.get(i).intValue()));
		for (int i = 0; i < numMicroBinsB; i++)
			ret.microBinIDs.add(new Integer(b.microBinIDs.get(i).intValue()));
		
		// adjust MDH accordingly
		int maxMicroBinIDA = a.microBinIDs.get(numMicroBinsA - 1).intValue();
		//int minMicroBinIDB = b.microBinIDs.get(0).intValue();
		ret.numMDHs = a.numMDHs + b.numMDHs;
		if (a.mdh.length == 0 && b.mdh.length == 0)
		{
			ret.mdh = new int[1];
			ret.mdh[0] = interBinsDiffBinary[maxMicroBinIDA];
			if (ret.mdh[0] == 1)
				ret.numMDHs++;
			//else
			//	ret.numMDHs = 0;
		}
		else if (a.mdh.length == 0)
		{
			ret.mdh = new int[b.mdh.length + 1];
			for (int i = 0; i < b.mdh.length; i++)
				ret.mdh[i + 1] = b.mdh[i];
			
			ret.mdh[0] = interBinsDiffBinary[maxMicroBinIDA];
			if (ret.mdh[0] == 1)
				ret.numMDHs++;
			//else
			//	ret.numMDHs = b.numMDHs;
		}
		else if (b.mdh.length == 0)
		{
			ret.mdh = new int[a.mdh.length + 1];
			for (int i = 0; i < a.mdh.length; i++)
				ret.mdh[i] = a.mdh[i];
			
			ret.mdh[a.mdh.length] = interBinsDiffBinary[maxMicroBinIDA];
			if (ret.mdh[a.mdh.length] == 1)
				ret.numMDHs++;
			//else
			//	ret.numMDHs = a.numMDHs;
		}
		else
		{
			ret.mdh = new int[a.mdh.length + b.mdh.length + 1];
			for (int i = 0; i < a.mdh.length; i++)
				ret.mdh[i] = a.mdh[i];
			for (int i = 0; i < b.mdh.length; i++)
				ret.mdh[i + 1 + a.mdh.length] = b.mdh[i];
			
			ret.mdh[a.mdh.length] = interBinsDiffBinary[maxMicroBinIDA];
			if (ret.mdh[a.mdh.length] == 1)
				ret.numMDHs++;
			//else
			//	ret.numMDHs = a.numMDHs + b.numMDHs;
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
		ret.means = new double[numCols];
		ret.devs = new double[numCols];
		DataPoint tmpPoint = null;
		for (int i = 0; i < numValuesA; i++)
		{
			ret.allValues.add(new Double(a.allValues.get(i).intValue()));
			ret.pointIDs.add(new Integer(a.pointIDs.get(i).intValue()));
			
			// add points to the data matrix
			tmpPoint = a.dataMatrix.data.get(i);
			for (int j = 0; j < numCols; j++)
			{
				if (i == 0)
					ret.means[j] = tmpPoint.measures[j];
				else
					ret.means[j] += tmpPoint.measures[j];
			}
			ret.dataMatrix.data.add(a.dataMatrix.data.get(i));
		}
		for (int i = 0; i < numValuesB; i++)
		{
			ret.allValues.add(new Double(b.allValues.get(i).intValue()));
			ret.pointIDs.add(new Integer(b.pointIDs.get(i).intValue()));
			
			// add points to the data matrix
			tmpPoint = b.dataMatrix.data.get(i);
			for (int j = 0; j < numCols; j++)
				ret.means[j] += tmpPoint.measures[j];
			ret.dataMatrix.data.add(b.dataMatrix.data.get(i));
		}
		
		// combine index matrices
		ret.indexMatrix = new IndexMatrix();
		ret.indexMatrix.rows = ret.dataMatrix.rows;
		ret.indexMatrix.cols = ret.dataMatrix.cols;
		ret.indexMatrix.data = new int[ret.dataMatrix.rows][ret.dataMatrix.cols];
		DataProcess.calculateIndices(ret.dataMatrix, ret.indexMatrix);
		
		// compute means and devs
		int numRows = ret.dataMatrix.rows;
		for (int i = 0; i < numCols; i++)
		{
			ret.dims.add(new Integer(i));
			ret.means[i] = ret.means[i] / numRows;
			if (numRows == 1)
				ret.devs[i] = 0;
			else
			{
				for (int j = 0; j < numRows; j++)
				{
					tmpPoint = ret.dataMatrix.data.get(j);
					ret.devs[i] += (tmpPoint.measures[i] - ret.means[i]) * (tmpPoint.measures[i] - ret.means[i]) / (numRows - 1);
				}
				ret.devs[i] = Math.sqrt(ret.devs[i]);
			}
		}
		
		return ret;
	}
	
	// T = total number of micro bins contained in all macro bins
	// N = total number of data points contained in DS
	public static double computeEncodingCost(ArrayList<MacroBin> bins, int N, int T, ArrayList<MicroBin> initialBins)
	{
		// the total coding cost
		double ret = 0;
		
		double logBase = Math.log(Constants.LOG_BASE);
		
		// get the total number of macro bins
		int numBins = bins.size();
		
		// L(DS | M)
		double LDSM = 0;
		
		// L_{cp}(M)
		double LcpM = LInteger(numBins - 1) + logOfCombination(T - 1, numBins - 1);
		
		// L_{bid}(M)
		double LbidM = 0;
		
		// L_{ic}(M)
		double LicM = 0;
		
		// L_{mdh}(M)
		double LmdhM = 0;
		
		
		// iterate through each macro bins
		MacroBin tmpBin = null;
		//int tmpBinSupport;
		int tmpNumMicroBins;
		double logTmpBinSupport;
		for (int i = 0; i < numBins; i++)
		{
			// get the macro bin
			tmpBin = bins.get(i);
			
			// get the total number of micro bins contained in the macro bin
			tmpNumMicroBins = tmpBin.microBinIDs.size();
			if (tmpNumMicroBins != 0)
				logTmpBinSupport = Math.log(tmpNumMicroBins * 1.0 / T) / logBase;
			else
				logTmpBinSupport = 0;
			
			// increase L(DS | M)
			LDSM -= tmpNumMicroBins * logTmpBinSupport;
			LDSM += LInteger(tmpNumMicroBins);
			
			// increase L_{bid}(M)
			LbidM -= logTmpBinSupport;
			
			// increase L_{ic}(M)
			if (Constants.USE_CE == true)
				LicM += computeBinInformationContentWithCE(tmpBin, initialBins);
			else
				LicM += computeBinInformationContentWithoutCE(tmpBin, initialBins);
			
			// increase L_{mdh}(M)
			if (tmpBin.numMDHs != 0)
				LmdhM += LInteger(tmpBin.numMDHs);
			if (tmpNumMicroBins != 1)
				LmdhM += tmpBin.numMDHs * Math.log(tmpNumMicroBins - 1) / logBase;
			/*if (tmpNumMicroBins != 1)
			{
				if (tmpBin.numMDHs != 0)
					LmdhM -= tmpBin.numMDHs * Math.log(tmpBin.numMDHs * 1.0 / (tmpNumMicroBins - 1)) / logBase;
				if (tmpNumMicroBins - 1 - tmpBin.numMDHs != 0)
					LmdhM -= (tmpNumMicroBins - 1 - tmpBin.numMDHs) * Math.log((tmpNumMicroBins - 1 - tmpBin.numMDHs) * 1.0 / (tmpNumMicroBins - 1)) / logBase;
			}*/
		}
		
		// get the total coding cost
		ret = LDSM + LcpM + LbidM + LicM + LmdhM;
		
		return ret;
	}
	
	// compute bin's IC using normal entropy
	public static double computeBinInformationContentWithoutCE(MacroBin macroBin, ArrayList<MicroBin> initialBins)
	{
		// the total IC of the macro bin
		double ret = 0;
		
		double logBase = Math.log(Constants.LOG_BASE);
		
		// get the number of micro bins
		int numMicroBins = macroBin.microBinIDs.size();
		
		// get the total number of data points covered by the macro bins
		//int numPoints = macroBin.numPoints;
		
		// iterate through all micro bins contained in the macro bins
		/*double tmpValueSupport;
		MicroBin tmpMicroBin = null;
		for (int i = 0; i < numMicroBins; i++)
		{
			// get the actual micro bins
			tmpMicroBin = initialBins.get(macroBin.microBinIDs.get(i).intValue());
			
			// get the total number of data points (including duplicates) contained in the micro bin
			tmpValueSupport = tmpMicroBin.allValues.size();
			
			// add the IC of the micro bin into the total IC of the macro bin
			if (tmpValueSupport != 0)
				ret -= tmpValueSupport * Math.log(tmpValueSupport * 1.0 / numPoints) / logBase;
		}*/
		ret = numMicroBins * Math.log(numMicroBins * 1.0) / logBase;
		
		return ret;
	}
	
	// compute bin's IC using CE
	public static double computeBinInformationContentWithCE(MacroBin macroBin, ArrayList<MicroBin> initialBins)
	{
		// the total IC of the macro bin
		double ret = 0;
		
		double logBase = Math.log(Constants.LOG_BASE);
		
		// get the number of micro bins
		int numMicroBins = macroBin.microBinIDs.size();
		
		// get the total number of data points covered by the macro bins
		int numPoints = macroBin.numPoints;
		
		// iterate through all micro bins contained in the macro bins
		MicroBin tmpMicroBin1;
		MicroBin tmpMicroBin2;
		int curPointCount = initialBins.get(macroBin.microBinIDs.get(0).intValue()).allValues.size();
		double maxVal = Constants.MAX_VAL;
		for (int i = 0; i < numMicroBins - 1; i++)
		{
			// get two consecutive micro bins
			tmpMicroBin1 = initialBins.get(macroBin.microBinIDs.get(i).intValue());
			tmpMicroBin2 = initialBins.get(macroBin.microBinIDs.get(i + 1).intValue());
			
			// add the IC of the micro bin into the total IC of the macro bin
			if (curPointCount != 0)
				ret -= (tmpMicroBin2.lowerBound - tmpMicroBin1.lowerBound) * (curPointCount / (1.0 * numPoints)) * Math.log(curPointCount / (1.0 * numPoints)) / (2 * maxVal * logBase);
			
			// increase the current point count
			curPointCount += tmpMicroBin2.allValues.size();
		}
		
		ret = ret * numPoints;
		return ret;
	}
	
	// divide a dimension into equal-frequency bins
	public static ArrayList<MicroBin> findEqualFrequencyBinningNew(int curDim, int[] ODIMS, DataMatrix dataMatrix, int numBins, SortedObject[] tmp)
	{
		ArrayList<MicroBin> ret = new  ArrayList<MicroBin>();
		
		int numRows = dataMatrix.rows;
		int numCols = dataMatrix.cols;
		
		// number of distinct values per bin
		int binCount = (int)Math.floor(numRows * 1.0 / numBins);
		//System.out.println(T + " --- " + numBins);
		
		int lowerIndex;
		int upperIndex;
		MicroBin tmpBin = null;
		//ArrayList<Integer> tmpPointIDs = null;
		//int tmpPointCount;
		DataPoint curPoint = null;
		DataPoint tmpPoint = null;
		int curTmpCol;
		//int tmpPointID;
		int curNumRows;
		int curNumCols;
		//int curSupport;
		double curVal;
		int curCount;
		int curNumPoints;
		//System.out.println("num of distinct values = " + distinctValues.size());
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
			//System.out.println(tmpBin.lowerBound + " --- " + tmpBin.upperBound);
			
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
				if (ODIMS[l] != curDim)
				{
					curTmpCol++;
					tmpPoint.measures[curTmpCol] = curPoint.measures[ODIMS[l]];
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
					if (ODIMS[l] != curDim)
					{
						curTmpCol++;
						tmpPoint.measures[curTmpCol] = curPoint.measures[ODIMS[l]];
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
		
		//System.out.println("---end---");
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
		//System.out.println(T + " --- " + numBins);
		
		int lowerIndex;
		int upperIndex;
		MicroBin tmpBin = null;
		//ArrayList<Integer> tmpPointIDs = null;
		//int tmpPointCount;
		DataPoint curPoint = null;
		DataPoint tmpPoint = null;
		int curTmpCol;
		//int tmpPointID;
		int curNumRows;
		int curNumCols;
		//int curSupport;
		double curVal;
		int curCount;
		int curNumPoints;
		//System.out.println("num of distinct values = " + distinctValues.size());
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
			//System.out.println(tmpBin.lowerBound + " --- " + tmpBin.upperBound);
			
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
		
		//System.out.println("---end---");
		return ret;
	}
	
	// divide a dimension into equal-frequency bins
	public static ArrayList<MicroBin> findEqualWidthBinning(int curDim, DataMatrix dataMatrix, int numDesiredBins, ArrayList<Double> distinctValues, ArrayList<Integer> distinctValueSupports, ArrayList<ArrayList<Integer>> memberpointIDs)
	{
		ArrayList<MicroBin> ret = new  ArrayList<MicroBin>();
		
		int T = distinctValues.size();
		double minValue = distinctValues.get(0).doubleValue();
		double maxValue = distinctValues.get(T - 1).doubleValue();
		int numCols = dataMatrix.cols;
		
		// number of distinct values per bin
		double binWidth = (maxValue - minValue) / numDesiredBins;
		
		double lowerBound;
		double upperBound;
		MicroBin tmpBin = null;
		ArrayList<Integer> tmpPointIDs = null;
		int tmpPointCount;
		DataPoint curPoint = null;
		DataPoint tmpPoint = null;
		int curTmpCol;
		int tmpPointID;
		int curNumRows;
		int curNumCols;
		int curSupport;
		double curVal;
		for (int i = 0; i < numDesiredBins; i++)
		{
			// lower and upper indices of distinct values
			lowerBound = minValue + i * binWidth;
			if (i == numDesiredBins - 1)
				upperBound = maxValue + 1;
			else
				upperBound = minValue + (i + 1) * binWidth;
			
			// get the bin's lower and upper bounds
			tmpBin = new MicroBin(lowerBound, upperBound);
			//System.out.println(tmpBin.lowerBound + " --- " + tmpBin.upperBound);
			
			// calculate the bin's data matrix, index matrix, means, and devs
			tmpBin.dataMatrix = new DataMatrix();
			tmpBin.dataMatrix.data = new ArrayList<DataPoint>();
			for (int j = 0; j < numCols - 1; j++)	// set up the bin's multivariate dimensions
				tmpBin.dims.add(new Integer(j));
			tmpBin.means = new double[numCols - 1];
			tmpBin.devs = new double[numCols - 1];
			for (int j = 0; j < T; j++)
			{
				// get current value
				curVal = distinctValues.get(j).doubleValue();
				
				if (curVal >= lowerBound && curVal < upperBound)
				{
					// add the current value into the list of bin's distinct values
					tmpBin.distinctValues.add(new Double(curVal));
					
					// add the support of the new distinct value into the bin's list of supports
					tmpBin.distinctValueSupports.add(new Integer(distinctValueSupports.get(j).intValue()));
					
					// add values (including duplicates) into the bin's list of all values
					curSupport = distinctValueSupports.get(j).intValue();
					for (int k = 0; k < curSupport; k++)
						tmpBin.allValues.add(new Double(curVal));
					
					// add all respective data points to the bin's data matrix
					tmpPointIDs = memberpointIDs.get(j);
					tmpPointCount = tmpPointIDs.size();
					for (int k = 0; k < tmpPointCount; k++)
					{
						tmpPointID = tmpPointIDs.get(k).intValue();
						curPoint = dataMatrix.data.get(tmpPointID);
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
						tmpBin.pointIDs.add(new Integer(tmpPointID));
					}
				}
			}
			
			if (tmpBin.distinctValues.size() > 0)
			{
				if (tmpBin.distinctValues.size() > 0)
				{
					tmpBin.lowerBound = tmpBin.distinctValues.get(0).doubleValue();
					tmpBin.upperBound = tmpBin.distinctValues.get(tmpBin.distinctValues.size() - 1).doubleValue();
				}
				else
				{
					tmpBin.lowerBound = lowerBound;
					tmpBin.upperBound = upperBound;
				}
				//System.out.println(tmpBin.upperBound);
				
				// calculate the bin's index matrix
				tmpBin.dataMatrix.rows = tmpBin.pointIDs.size();
				tmpBin.dataMatrix.cols = numCols - 1;
				tmpBin.indexMatrix = new IndexMatrix();
				tmpBin.indexMatrix.rows = tmpBin.dataMatrix.rows;
				tmpBin.indexMatrix.cols = tmpBin.dataMatrix.cols;
				tmpBin.indexMatrix.data = new int[tmpBin.dataMatrix.rows][tmpBin.dataMatrix.cols];
				if (tmpBin.distinctValues.size() > 0)
					DataProcess.calculateIndices(tmpBin.dataMatrix, tmpBin.indexMatrix);
				
				// calculate the bin's means and deviations
				if (tmpBin.distinctValues.size() > 0)
				{
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
					} // end for
				} // end if
				
				// add the bin to the list of bins
				ret.add(tmpBin);
			}
		}
		
		//System.out.println("---end---");
		return ret;
	}
	
	// calculate the MCKL distance of two micro bins
	public static double calculateBinsDiffNewNew(MicroBin a, MicroBin b, double v) throws Exception
	{
		double diff = MCKLComputer.computeDistanceNewNew(a.dataMatrix, b.dataMatrix, v);
		return diff;
	}
	
	// calculate the MCKL distance of two micro bins
	public static double calculateBinsDiffNew(MicroBin a, MicroBin b) throws Exception
	{
		double t0 = 2 * Constants.MAX_VAL * a.dims.size() * (2 + 1 / Math.log(2));
		double diffAB = MCKLComputer.computeDistanceNew(a.dataMatrix, a.indexMatrix, b.dataMatrix, b.indexMatrix, a.dims, a.means, a.devs);
		if (diffAB < 0 || diffAB >= t0)
		{
			System.out.println(diffAB + " --- " + t0);
			throw new Exception("Invalid MCKL calculation");
		}
		
		double diffBA = MCKLComputer.computeDistanceNew(b.dataMatrix, b.indexMatrix, a.dataMatrix, a.indexMatrix, b.dims, b.means, b.devs);
		if (diffBA < 0 || diffBA >= t0)
		{
			System.out.println(diffBA + " --- " + t0);
			throw new Exception("Invalid MCKL calculation");
		}
		
		return diffAB + diffBA;
	}
	
	// calculate the MCKL distance of two micro bins
	public static double calculateBinsDiffByMean(MicroBin a, MicroBin b) throws Exception
	{
		int numPointsA = a.dataMatrix.rows;
		int numPointsB = b.dataMatrix.rows;
		int numCols = a.dataMatrix.cols;
		DataPoint curPoint = null;
		
		// compute meanA
		double[] meanA = new double[numCols];
		for (int i = 0; i < numPointsA; i++)
		{
			curPoint = a.dataMatrix.data.get(i);
			for (int j = 0; j < numCols; j++)
			{
				if (i == 0)
					meanA[j] = curPoint.measures[j] / numPointsA;
				else
					meanA[j] += curPoint.measures[j] / numPointsA;
			}
		}
		
		// compute meanB
		double[] meanB = new double[numCols];
		for (int i = 0; i < numPointsB; i++)
		{
			curPoint = b.dataMatrix.data.get(i);
			for (int j = 0; j < numCols; j++)
			{
				if (i == 0)
					meanB[j] = curPoint.measures[j] / numPointsB;
				else
					meanB[j] += curPoint.measures[j] / numPointsB;
			}
		}
		
		// compute the distance between two means
		int expo = Constants.LP_NORM;
		double dist = 0;
		for (int i = 0; i < numCols; i++)
			dist += Math.pow(Math.abs(meanA[i] - meanB[i]), expo);
		
		return Math.pow(dist, 1.0 / expo);
	}
	
	// calculate the MCKL distance of two micro bins
	public static double calculateBinsDiffByMinDist(MicroBin a, MicroBin b) throws Exception
	{
		int numPointsA = a.dataMatrix.rows;
		int numPointsB = b.dataMatrix.rows;
		int expo = Constants.LP_NORM;
		DataPoint curPointA = null;
		DataPoint curPointB = null;
		double minDist = Double.MAX_VALUE;
		double dist;
		
		for (int i = 0; i < numPointsA; i++)
		{
			curPointA = a.dataMatrix.data.get(i);
			for (int j = 0; j < numPointsB; j++)
			{
				curPointB = b.dataMatrix.data.get(j);
				dist = DataPoint.distanceLNorm(expo, curPointA, curPointB);
				if (dist < minDist)
					minDist = dist;
			}
		}
		
		return minDist;
	}
	
	// calculate the MCKL distance of two micro bins
	public static double calculateBinsDiffByMaxDist(MicroBin a, MicroBin b) throws Exception
	{
		int numPointsA = a.dataMatrix.rows;
		int numPointsB = b.dataMatrix.rows;
		int expo = Constants.LP_NORM;
		DataPoint curPointA = null;
		DataPoint curPointB = null;
		double maxDist = -Double.MAX_VALUE;
		double dist;
		
		for (int i = 0; i < numPointsA; i++)
		{
			curPointA = a.dataMatrix.data.get(i);
			for (int j = 0; j < numPointsB; j++)
			{
				curPointB = b.dataMatrix.data.get(j);
				dist = DataPoint.distanceLNorm(expo, curPointA, curPointB);
				if (dist > maxDist)
					maxDist = dist;
			}
		}
		
		return maxDist;
	}
	
	public static double LInteger(int N)
	{
		double logBase = Math.log(Constants.LOG_BASE);
		double ret = Math.log(Constants.C0) / logBase;
		
		double tmpVal = Math.log(N) / logBase;
		while (tmpVal > 0)
		{
			ret += tmpVal;
			tmpVal = Math.log(tmpVal) / logBase;
		}
		
		return ret;
	}
	
	public static double logOfCombination(int N, int k)
	{
		double ret = 0;
		
		int lowerBound = N - k + 1;
		for (int i = lowerBound; i <= N; i++)
			ret += Math.log(i);
		
		for (int i = 1; i <= k; i++)
			ret -= Math.log(i);
		
		ret = ret / Math.log(Constants.LOG_BASE);
		
		return ret;
	}
	
	public static int[] binarizeDoubleArray(double[] vals)
	{
		int[] ret = new int[vals.length];
		
		double[] tmp = new double[vals.length];
		for (int i = 0; i < vals.length; i++)
			tmp[i] = vals[i];
		
		Arrays.sort(tmp);
		double thyThreshold = tmp[(int)Math.ceil(tmp.length / Constants.QUANTILE) - 1];
		for (int i = 0; i < vals.length; i++)
		{
			if (vals[i] > thyThreshold)
				ret[i] = 1;
			else
				ret[i] = 0;
		}
		
		return ret;
	}
	
	public static double quantileThreshold(double[] vals)
	{
		double[] tmp = new double[vals.length];
		for (int i = 0; i < vals.length; i++)
			tmp[i] = vals[i];
		
		Arrays.sort(tmp);
		double thyThreshold = tmp[(int)Math.ceil(tmp.length / Constants.QUANTILE) - 1];
		return thyThreshold;
	}
}