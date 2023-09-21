package multibinning.business;

import java.util.*;

import multibinning.data.*;

public class MMIC 
{
	@SuppressWarnings("unchecked")
	public static ArrayList<MacroBin>[] discretizeData(DataMatrix dataMatrix, double[] CRES, double[] means, double[] devs) throws Exception
	{
		// initialize discretized data
		Constants.DISC_DATA = new int[dataMatrix.rows][dataMatrix.cols];
		
		// get all dimension
		int numCols = dataMatrix.cols;
		ArrayList<Integer> dims = new ArrayList<Integer>();
		for (int i = 0; i < numCols; i++)
			dims.add(new Integer(i));
		
		// initialize the list of cells and respective member points
		ArrayList<ArrayList<Integer>> cells = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> cellPointIDs = new ArrayList<ArrayList<Integer>>();
		
		// choose two most correlated dimensions
		int numDims = dims.size();
		ArrayList<Integer> tmpDims = new ArrayList<Integer>();
		ArrayList<Integer> tmpFrontDims = new ArrayList<Integer>();
		double maxPairContrast = -Double.MAX_VALUE;
		ArrayList<Integer> maxPair = new ArrayList<Integer>();
		double tmpContrast;
		double[][] tmpContrasts;
		int dim1;
		int dim2;
		for (int i = 0; i < numDims; i++)
			for (int j = i + 1; j < numDims; j++)
			{
				dim1 = dims.get(i);
				dim2 = dims.get(j);
				
				tmpDims.clear();
				tmpFrontDims.clear();
				tmpDims.add(new Integer(dim1));
				tmpFrontDims.add(new Integer(dim2));
				tmpContrasts = computeContrasts(dataMatrix, tmpDims, tmpFrontDims, CRES, Constants.NUM_CENTROIDS, means, devs);
				tmpContrast = tmpContrasts[1][0];
				
				tmpDims.clear();
				tmpFrontDims.clear();
				tmpDims.add(new Integer(dim2));
				tmpFrontDims.add(new Integer(dim1));
				tmpContrasts = computeContrasts(dataMatrix, tmpDims, tmpFrontDims, CRES, Constants.NUM_CENTROIDS, means, devs);
				tmpContrast += tmpContrasts[1][0];
				
				if (tmpContrast > maxPairContrast)
				{
					maxPairContrast = tmpContrast;
					maxPair.clear();
					maxPair.add(new Integer(dim1));
					maxPair.add(new Integer(dim2));
				}
			}
		
		// discretize two best dimensions
		double[] maxCorrelation = new double[1];
		ArrayList<MacroBin>[] existingMacroBins = new ArrayList[numCols];
		ArrayList<Integer> processedDims = new ArrayList<Integer>();
		int secondDim = maxPair.get(1);
		int firstDim = maxPair.get(0);
		System.out.println(firstDim + " --- " + secondDim);
		ArrayList<MacroBin>[] firstSecondBins = MICBinningForTwoDims(firstDim, secondDim, dataMatrix, existingMacroBins, Constants.DISC_DATA, maxCorrelation);
		
		existingMacroBins[firstDim] = firstSecondBins[0];
		updateDiscreteData(firstDim, existingMacroBins[firstDim], Constants.DISC_DATA);
		updateCells(existingMacroBins[firstDim], cells, cellPointIDs);
		processedDims.add(new Integer(firstDim));
		
		existingMacroBins[secondDim] = firstSecondBins[1];
		updateDiscreteData(secondDim, existingMacroBins[secondDim], Constants.DISC_DATA);
		updateCells(existingMacroBins[secondDim], cells, cellPointIDs);
		processedDims.add(new Integer(secondDim));
		
		System.out.println("dim " + firstDim + ": number of bins = " + existingMacroBins[firstDim].size());
		System.out.println("dim " + secondDim + ": number of bins = " + existingMacroBins[secondDim].size());
		System.out.println("---------------");
		
		// get the remaining dimensions
		ArrayList<Integer> remainingDims = new ArrayList<Integer>();
		int curDim;
		for (int i = 0; i < numDims; i++)
		{
			curDim = dims.get(i);
			if (curDim != firstDim && curDim != secondDim)
				remainingDims.add(new Integer(i));
		}
		
		// discretize each of the remaining dimensions
		int curLength = 2;
		int tmpDim;
		ArrayList<Double> tmpResult;
		ArrayList<MicroBin> initialMicroBins;
		ArrayList<MacroBin> candidateMacroBins;
		ArrayList<Double> distinctValues = new ArrayList<Double>();
		ArrayList<Integer> distinctValueSupports = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> memberPointIDs = new ArrayList<ArrayList<Integer>>();
		int B = (int)Math.pow(dataMatrix.rows, 0.6);
		int MAX_BINS;
		int maxNumBinsSoFar = (int)Math.max(existingMacroBins[firstDim].size(), existingMacroBins[secondDim].size());
		SortedObject[] sos = null;
		while (curLength < numDims)
		{
			distinctValues.clear();
			distinctValueSupports.clear();
			memberPointIDs.clear();
			
			tmpResult = computeContrastSpecial(dataMatrix, maxPair, remainingDims, CRES, cellPointIDs);
			tmpDim = tmpResult.get(0).intValue();
			
			sos = getDistinctValues(tmpDim, dataMatrix, distinctValues, distinctValueSupports, memberPointIDs);
			initialMicroBins = findEqualFrequencyBinning(tmpDim, dataMatrix, (int)Math.min(Constants.INIT_BIN_COUNT, Constants.MAX_BINS), sos);
			candidateMacroBins = convertMicroToMacroBins(initialMicroBins);
			
			// discretize tmpDim
			MAX_BINS = B / maxNumBinsSoFar;
			existingMacroBins[tmpDim] = MMICBinning(tmpDim, (int)Math.min(Math.min(Constants.INIT_BIN_COUNT, MAX_BINS), Constants.MAX_BINS), dataMatrix, processedDims, existingMacroBins, candidateMacroBins, initialMicroBins, cells, cellPointIDs, maxCorrelation, Constants.DISC_DATA);
			updateDiscreteData(tmpDim, existingMacroBins[tmpDim], Constants.DISC_DATA);
			updateCells(existingMacroBins[tmpDim], cells, cellPointIDs);
			processedDims.add(new Integer(tmpDim));
			maxNumBinsSoFar = (int)Math.max(maxNumBinsSoFar, existingMacroBins[tmpDim].size());
			
			curLength++;
		}
		
		return existingMacroBins;
	}
	
	public static void updateCells(ArrayList<MacroBin> macroBins, ArrayList<ArrayList<Integer>> cells, ArrayList<ArrayList<Integer>> cellPointIDs)
	{
		int numBins = macroBins.size();
		int curNumCells = cells.size();
		int curCellSize;
		int binNumPoints;
		int count;
		int pointID;
		int cellPointID;
		ArrayList<Integer> newCell = null;
		ArrayList<Integer> newCellPointIDs = null;
		ArrayList<Integer> curCell = null;
		ArrayList<Integer> curCellPointIDs = null;
		ArrayList<ArrayList<Integer>> newListCells = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> newListCellPointIDs = new ArrayList<ArrayList<Integer>>();
		MacroBin bin = null;
		for (int i = 0; i < numBins; i++)
		{
			bin = macroBins.get(i);
			binNumPoints = bin.numPoints;
			if (curNumCells > 0)
			{
				for (int j = 0; j < curNumCells; j++)
				{
					// create the new cell
					newCell = new ArrayList<Integer>();
					curCell = cells.get(j);
					curCellSize = curCell.size();
					for (int z = 0; z < curCellSize; z++)
						newCell.add(new Integer(curCell.get(z)));
					newCell.add(new Integer(i));
					
					// update point IDs of the new cell
					curCellPointIDs = cellPointIDs.get(j);
					count = curCellPointIDs.size();
					newCellPointIDs = new ArrayList<Integer>();
					for (int h = 0; h < binNumPoints; h++)
					{
						pointID = bin.pointIDs.get(h);
						for (int k = 0; k < count; k++)
						{
							cellPointID = curCellPointIDs.get(k);
							if (pointID == cellPointID)
							{
								newCellPointIDs.add(new Integer(pointID));
								break;
							}
						} // end for
					} // end for
					
					// add the new cell only if it is non-empty
					if (newCellPointIDs.size() > 0)
					{
						newListCellPointIDs.add(newCellPointIDs);
						newListCells.add(newCell);
					}
				} // end for
			}
			else
			{
				newCell = new ArrayList<Integer>();
				newCell.add(new Integer(i));
				cells.add(newCell);
				
				binNumPoints = bin.numPoints;
				newCellPointIDs = new ArrayList<Integer>();
				for (int h = 0; h < binNumPoints; h++)
				{
					pointID = bin.pointIDs.get(h);
					newCellPointIDs.add(new Integer(pointID));
				}
				cellPointIDs.add(newCellPointIDs);
			}
		} // end for
		
		// copy the information to cells
		if (curNumCells > 0)
		{
			cells.clear();
			cellPointIDs.clear();
			curNumCells = newListCells.size();
			for (int i = 0; i < curNumCells; i++)
			{
				cells.add(newListCells.get(i));
				cellPointIDs.add(newListCellPointIDs.get(i));
			}
		}
	}
	
	public static void updateDiscreteData(int curDim, ArrayList<MacroBin> macroBins, int[][] discreteData)
	{
		int numMacroBins = macroBins.size();
		MacroBin macroBin = null;
		int numPoints;
		int pointID;
		for (int i = 0; i < numMacroBins; i++)
		{
			macroBin = macroBins.get(i);
			numPoints = macroBin.pointIDs.size();
			for (int j = 0; j < numPoints; j++)
			{
				pointID = macroBin.pointIDs.get(j);
				discreteData[pointID][curDim] = i;
			}
		}
	}
	
	@SuppressWarnings("unchecked")
	public static ArrayList<MacroBin>[] MICBinningForTwoDims(int firstDim, int secondDim, DataMatrix dataMatrix, ArrayList<MacroBin>[] existingMacroBins, int[][] discreteData, double[] maxCorrelation) throws Exception
	{
		ArrayList<Integer> processedDims = new ArrayList<Integer>();
		ArrayList<MacroBin> candidateMacroBinsFirst;
		ArrayList<MacroBin> candidateMacroBinsSecond;
		ArrayList<MicroBin> initialMicroBinsFirst;
		ArrayList<MicroBin> initialMicroBinsSecond;
		ArrayList<MicroBin> initialMicroBins;
		ArrayList<ArrayList<Integer>> cells = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> cellPointIDs = new ArrayList<ArrayList<Integer>>();
		ArrayList<Double> distinctValuesFirstDim = new ArrayList<Double>();
		ArrayList<Integer> distinctValueSupportsFirstDim = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> memberPointIDsFirstDim = new ArrayList<ArrayList<Integer>>();
		ArrayList<Double> distinctValuesSecondDim = new ArrayList<Double>();
		ArrayList<Integer> distinctValueSupportsSecondDim = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> memberPointIDsSecondDim = new ArrayList<ArrayList<Integer>>();
		
		// get distinct values
		SortedObject[] sos1 = getDistinctValues(firstDim, dataMatrix, distinctValuesFirstDim, distinctValueSupportsFirstDim, memberPointIDsFirstDim);
		SortedObject[] sos2 = getDistinctValues(secondDim, dataMatrix, distinctValuesSecondDim, distinctValueSupportsSecondDim, memberPointIDsSecondDim);
		
		// get initial bins
		initialMicroBinsFirst = findEqualFrequencyBinning(firstDim, dataMatrix, (int)Math.min(Constants.INIT_BIN_COUNT, Constants.MAX_BINS), sos1);
		candidateMacroBinsFirst = convertMicroToMacroBins(initialMicroBinsFirst);
		initialMicroBinsSecond = findEqualFrequencyBinning(secondDim, dataMatrix, (int)Math.min(Constants.INIT_BIN_COUNT, Constants.MAX_BINS), sos2);
		candidateMacroBinsSecond = convertMicroToMacroBins(initialMicroBinsSecond);
		
		int B = (int)Math.pow(dataMatrix.rows, 0.6);
		int BHalf = (int)Math.min(B / 2, Constants.MAX_BINS);
		int x, y;
		double MAX_CORRELATION = -Double.MAX_VALUE;
		ArrayList<MacroBin> OPTIMAL_BINS_FIRST = new ArrayList<MacroBin>();
		ArrayList<MacroBin> OPTIMAL_BINS_SECOND = new ArrayList<MacroBin>();
		
		// process first dimension
		processedDims.clear();
		processedDims.add(new Integer(firstDim));
		ArrayList<MacroBin> tmpOptimal;
		double[] tmpOptimalCost = new double[1];
		int numBins;
		for (x = 2; x <= BHalf; x++)
		{
			cells.clear();
			cellPointIDs.clear();
			initialMicroBins = findEqualFrequencyBinning(firstDim, dataMatrix, (int)Math.min(Constants.MAX_BINS, x), sos1);
			existingMacroBins[firstDim] = convertMicroToMacroBins(initialMicroBins);
			updateCells(existingMacroBins[firstDim], cells, cellPointIDs);
			updateDiscreteData(firstDim, existingMacroBins[firstDim], discreteData);
			
			y = B / x;
			tmpOptimal = MMICBinning(secondDim, (int)Math.min(Constants.MAX_BINS, y), dataMatrix, processedDims, existingMacroBins, candidateMacroBinsSecond, initialMicroBinsSecond, cells, cellPointIDs, tmpOptimalCost, discreteData);
			if (tmpOptimalCost[0] > MAX_CORRELATION)
			{
				MAX_CORRELATION = tmpOptimalCost[0];
				
				OPTIMAL_BINS_FIRST.clear();
				numBins = existingMacroBins[firstDim].size();
				for (int i = 0; i < numBins; i++)
					OPTIMAL_BINS_FIRST.add(existingMacroBins[firstDim].get(i));
				
				OPTIMAL_BINS_SECOND.clear();
				numBins = tmpOptimal.size();
				for (int i = 0; i < numBins; i++)
					OPTIMAL_BINS_SECOND.add(tmpOptimal.get(i));
			}
		}
		
		// process second dimension
		processedDims.clear();
		processedDims.add(new Integer(secondDim));
		for (y = 2; y <= BHalf; y++)
		{
			cells.clear();
			cellPointIDs.clear();
			initialMicroBins = findEqualFrequencyBinning(secondDim, dataMatrix, (int)Math.min(Constants.MAX_BINS, y), sos2);
			existingMacroBins[secondDim] = convertMicroToMacroBins(initialMicroBins);
			updateCells(existingMacroBins[secondDim], cells, cellPointIDs);
			updateDiscreteData(secondDim, existingMacroBins[secondDim], discreteData);
			
			x = B / y;
			tmpOptimal = MMICBinning(firstDim, (int)Math.min(Constants.MAX_BINS, x), dataMatrix, processedDims, existingMacroBins, candidateMacroBinsFirst, initialMicroBinsFirst, cells, cellPointIDs, tmpOptimalCost, discreteData);
			if (tmpOptimalCost[0] > MAX_CORRELATION)
			{
				MAX_CORRELATION = tmpOptimalCost[0];
				
				OPTIMAL_BINS_FIRST.clear();
				numBins = tmpOptimal.size();
				for (int i = 0; i < numBins; i++)
					OPTIMAL_BINS_FIRST.add(tmpOptimal.get(i));
				
				OPTIMAL_BINS_SECOND.clear();
				numBins = existingMacroBins[secondDim].size();
				for (int i = 0; i < numBins; i++)
					OPTIMAL_BINS_SECOND.add(existingMacroBins[secondDim].get(i));
			}
		}
		
		maxCorrelation[0] = MAX_CORRELATION;
		ArrayList<MacroBin>[] OPTIMAL_BINS = new ArrayList[2];
		OPTIMAL_BINS[0] = OPTIMAL_BINS_FIRST;
		OPTIMAL_BINS[1] = OPTIMAL_BINS_SECOND;
		return OPTIMAL_BINS;
	}
	
	public static ArrayList<MacroBin> convertMicroToMacroBins(ArrayList<MicroBin> initialMicroBins)
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		
		// for each micro bin, create a macro bin containing it
		int numMicroBins = initialMicroBins.size();
		int pointID;
		MacroBin tmpMacroBin = null;
		MicroBin tmpMicroBin = null;
		for (int i = 0; i < numMicroBins; i++)
		{
			tmpMicroBin = initialMicroBins.get(i);
			tmpMacroBin = new MacroBin(tmpMicroBin.lowerBound, tmpMicroBin.upperBound);
			//tmpMacroBin.numPoints = tmpMicroBin.allValues.size();
			tmpMacroBin.numPoints = tmpMicroBin.pointIDs.size();
			tmpMacroBin.microBinIDs.add(new Integer(i));
			for (int j = 0; j < tmpMacroBin.numPoints; j++)
			{
				pointID = tmpMicroBin.pointIDs.get(j);
				tmpMacroBin.pointIDs.add(new Integer(pointID));
			}
			ret.add(tmpMacroBin);
		}
		
		return ret;
	}
	
	public static SortedObject[] getDistinctValues(int curDim, DataMatrix dataMatrix, ArrayList<Double> distinctValues, ArrayList<Integer> distinctValueSupports, ArrayList<ArrayList<Integer>> memberPointIDs)
	{
		// sort the data along the curDim dimension in ascending order
		int numRows = dataMatrix.rows;
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
		} // end for
		
		return tmp;
	}
	
	@SuppressWarnings("unchecked")
	public static ArrayList<MacroBin> MMICBinning(int curDim, int MAX_BINS, DataMatrix dataMatrix, ArrayList<Integer> processedDims, ArrayList<MacroBin>[] existingMacroBins, ArrayList<MacroBin> candidateMacroBins, ArrayList<MicroBin> initialBins, ArrayList<ArrayList<Integer>> cells, ArrayList<ArrayList<Integer>> cellPointIDs, double[] cost, int[][] discreteData) throws Exception
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		
		double logBase = Math.log(Constants.LOG_BASE);
		int numProcessedDims = processedDims.size();
		int[] procDims = new int[numProcessedDims];
		double sumLog = 0;
		double maxLog = -Double.MAX_VALUE;
		double val;
		for (int i = 0; i < numProcessedDims; i++)
		{
			procDims[i] = processedDims.get(i);
			val = Math.log(existingMacroBins[procDims[i]].size()) / logBase;
			sumLog += val;
			if (val > maxLog)
				maxLog = val;
		}
		
		// create the table for dynamic programming to store already solved sub-problems
		int numTotalMacroBins = candidateMacroBins.size();
		ArrayList<MacroBin>[][] dynamicBinnings = new ArrayList[numTotalMacroBins + 1][MAX_BINS + 1];
		double[][] dynamicBinningCosts = new double[numTotalMacroBins + 1][MAX_BINS + 1];
		int[][] dynamicNumPoints = new int[numTotalMacroBins + 1][MAX_BINS + 1];
		
		// find optimal binning of size 2
		MacroBin tmpMergeMacroBin = null;
		MacroBin tmpMergeMacroBinRight = null;
		//MacroBin macroBin = null;
		int numBins;
		int numPoints = 0;
		double tmpCost;
		ArrayList<MacroBin> tmpMacroBins = new ArrayList<MacroBin>();
		double[] entropies = new double[numTotalMacroBins + 1];
		ArrayList<Integer> pointIDs = new ArrayList<Integer>();
		//int countJoin;
		//int count;
		int pointID;
		//int cellPointID;
		int[] counts = null;
		for (int numFirstMacroBins = 2; numFirstMacroBins <= numTotalMacroBins; numFirstMacroBins++)
		{
			dynamicBinningCosts[numFirstMacroBins][2] = -Double.MAX_VALUE;
			pointIDs.clear();
			
			// loop through each possible bin position
			// maxIndex contains the minimum bin index of the right macro bin
			for (int maxIndex = 1; maxIndex < numFirstMacroBins; maxIndex++)
			{
				// clear the temporary binning
				tmpMacroBins.clear();
				
				// merge all macro bins from position 0 to (maxIndex - 1)
				tmpMergeMacroBin = cloneBin(candidateMacroBins.get(0));
				for (int i = 1; i < maxIndex; i++)
					tmpMergeMacroBin = mergeMacroBins(tmpMergeMacroBin, candidateMacroBins.get(i), initialBins);
				
				// merge all macro bins from position maxIndex to (numFirstMacroBins - 1)
				tmpMergeMacroBinRight = cloneBin(candidateMacroBins.get(maxIndex));
				for (int i = maxIndex + 1; i < numFirstMacroBins; i++)
					tmpMergeMacroBinRight = mergeMacroBins(tmpMergeMacroBinRight, candidateMacroBins.get(i), initialBins);
				
				// populate pointIDs
				if (maxIndex == 1)
				{
					for (int i = 0; i < tmpMergeMacroBin.numPoints; i++)
						pointIDs.add(new Integer(tmpMergeMacroBin.pointIDs.get(i)));
					
					for (int i = 0; i < tmpMergeMacroBinRight.numPoints; i++)
						pointIDs.add(new Integer(tmpMergeMacroBinRight.pointIDs.get(i)));
					
					numPoints = tmpMergeMacroBin.numPoints + tmpMergeMacroBinRight.numPoints;
				}
				
				// add two new macro bins to the temporary binning
				tmpMacroBins.add(tmpMergeMacroBin);
				tmpMacroBins.add(tmpMergeMacroBinRight);
				
				// compute the coding cost of the temporary binning
				tmpCost = computeEncodingCost(dataMatrix, procDims, discreteData, curDim, tmpMacroBins, tmpMergeMacroBin.numPoints + tmpMergeMacroBinRight.numPoints, dataMatrix.rows, initialBins, cells, cellPointIDs);
				
				// if the cost is less than the current cost, then update the cost, the binning, and the total number of points of the binning
				if (tmpCost > dynamicBinningCosts[numFirstMacroBins][2])
				{
					dynamicBinningCosts[numFirstMacroBins][2] = tmpCost;
					dynamicNumPoints[numFirstMacroBins][2] = tmpMergeMacroBin.numPoints + tmpMergeMacroBinRight.numPoints;
					if (dynamicBinnings[numFirstMacroBins][2] == null)
						dynamicBinnings[numFirstMacroBins][2] = new ArrayList<MacroBin>();
					else
						dynamicBinnings[numFirstMacroBins][2].clear();
					dynamicBinnings[numFirstMacroBins][2].add(tmpMergeMacroBin);
					dynamicBinnings[numFirstMacroBins][2].add(tmpMergeMacroBinRight);
				}
			} // end for
			
			// compute summation entropies of all dimensions
			entropies[numFirstMacroBins] = 0;
			for (int i = 0; i < numProcessedDims; i++)
			{
				numBins = existingMacroBins[procDims[i]].size();
				counts = new int[numBins];
				for (int j = 0; j < numPoints; j++)
				{
					pointID = pointIDs.get(j);
					counts[discreteData[pointID][procDims[i]]]++;
				} // end for
				
				for (int j = 0; j < numBins; j++)
				{
					if (counts[j] > 0)
						entropies[numFirstMacroBins] -= (counts[j] * 1.0 / numPoints) * Math.log(counts[j] * 1.0 / numPoints) / logBase;
				} // end for
			} // end for
			
			//System.out.println("numFirstMacroBins = " + numFirstMacroBins);
		} // end for
		
		for (int numFirstMacroBins = 2; numFirstMacroBins <= numTotalMacroBins; numFirstMacroBins++)
		{
			if (entropies[numFirstMacroBins] < 0)
				throw new Exception("negative single entropy = " + entropies[numFirstMacroBins]);
			
			dynamicBinningCosts[numFirstMacroBins][2] += entropies[numFirstMacroBins];
			if (Math.abs(dynamicBinningCosts[numFirstMacroBins][2]) < Constants.MAX_ERROR)
				dynamicBinningCosts[numFirstMacroBins][2] = 0;
			if (dynamicBinningCosts[numFirstMacroBins][2] < 0)
				throw new Exception("negative entropy = " + dynamicBinningCosts[numFirstMacroBins][2]);
		}
		
		// build the rest of the table
		// for each (numMacroBins, numFirstMacroBins) pairs, find the binning that yields the smallest cost
		int tmpTotalNumPoints;
		int optimalPosition;
		MacroBin optimalNewBin = null;
		int numCells = cellPointIDs.size();
		counts = new int[numCells];
		//ArrayList<Integer> curCellPointIDs = null;
		ArrayList<Integer> curCell = null;
		ArrayList<Integer> cellIDs = new ArrayList<Integer>();
		int curBinID;
		int cellBinID;
		int totalPoints;
		boolean match;
		for (int numMacroBins = 3; numMacroBins <= MAX_BINS; numMacroBins++)
		{
			for (int numFirstMacroBins = numMacroBins; numFirstMacroBins <= numTotalMacroBins; numFirstMacroBins++)
			{
				dynamicBinningCosts[numFirstMacroBins][numMacroBins] = -Double.MAX_VALUE;
				optimalPosition = -1;
				optimalNewBin = null;
				
				// loop through each possible bin position
				for (int maxIndex = numMacroBins - 1; maxIndex < numFirstMacroBins; maxIndex++)
				{
					// merge all the macro bins from position maxIndex to (numFirstMacroBins - 1)
					tmpMergeMacroBinRight = cloneBin(candidateMacroBins.get(maxIndex));
					for (int i = maxIndex + 1; i < numFirstMacroBins; i++)
						tmpMergeMacroBinRight = mergeMacroBins(tmpMergeMacroBinRight, candidateMacroBins.get(i), initialBins);
					
					// compute the total number of data points
					tmpTotalNumPoints = dynamicNumPoints[maxIndex][numMacroBins - 1] + tmpMergeMacroBinRight.numPoints;
					
					// retrieve the stored results
					tmpCost = dynamicNumPoints[maxIndex][numMacroBins - 1] * (dynamicBinningCosts[maxIndex][numMacroBins - 1] - entropies[maxIndex]) / tmpTotalNumPoints;
					
					// clear point counts
					for (int j = 0; j < numCells; j++)
						counts[j] = 0;
					
					// compute the entropy of tmpMergeMacroBinRight
					totalPoints = 0;
					for (int h = 0; h < tmpMergeMacroBinRight.numPoints; h++)
					{
						pointID = tmpMergeMacroBinRight.pointIDs.get(h);
						cellIDs.clear();
						for (int k = 0; k < numProcessedDims; k++)
							cellIDs.add(new Integer(discreteData[pointID][procDims[k]]));
						
						for (int j = 0; j < numCells; j++)
						{
							if (h == 0)
								counts[j] = 0;
							
							curCell = cells.get(j);
							match = true;
							for (int k = 0; k < numProcessedDims; k++)
							{
								curBinID = cellIDs.get(k);
								cellBinID = curCell.get(k);
								if (curBinID != cellBinID)
								{
									match = false;
									break;
								}
							} // end for
							
							if (match == true)
							{
								counts[j]++;
								break;
							}
						} // end for
					} // end for
					
					// increase cost
					for (int j = 0; j < numCells; j++)
					{
						totalPoints += counts[j];
						if (counts[j] > 0)
							tmpCost += (tmpMergeMacroBinRight.numPoints * 1.0 / tmpTotalNumPoints) * (counts[j] * 1.0 / tmpMergeMacroBinRight.numPoints) * Math.log(counts[j] * 1.0 / tmpMergeMacroBinRight.numPoints) / logBase;
					}
					
					// this serves as control check
					if (totalPoints != tmpMergeMacroBinRight.numPoints)
					{
						System.out.println(totalPoints + " --- " + tmpMergeMacroBinRight.numPoints);
						throw new Exception("total points not match");
					}
						
					// compute the entropy of tmpMergeMacroBinRight
					/*for (int j = 0; j < numCells; j++)
					{
						curCellPointIDs = cellPointIDs.get(j);
						count = curCellPointIDs.size();
						countJoin = 0;
						for (int h = 0; h < tmpMergeMacroBinRight.numPoints; h++)
						{
							pointID = tmpMergeMacroBinRight.pointIDs.get(h);
							for (int k = 0; k < count; k++)
							{
								cellPointID = curCellPointIDs.get(k);
								if (pointID == cellPointID)
								{
									countJoin++;
									break;
								}
							} // end for
						} // end for
						
						if (countJoin > 0)
							tmpCost += (tmpMergeMacroBinRight.numPoints * 1.0 / tmpTotalNumPoints) * (countJoin * 1.0 / tmpMergeMacroBinRight.numPoints) * Math.log(countJoin * 1.0 / tmpMergeMacroBinRight.numPoints) / logBase;
					}*/
					
					// if the new cost is larger than the current cost, then update the binning
					tmpCost += entropies[numFirstMacroBins];
					if (tmpCost < 0)
						throw new Exception("negative entropy = " + tmpCost);
					
					if (tmpCost > dynamicBinningCosts[numFirstMacroBins][numMacroBins])
					{
						dynamicBinningCosts[numFirstMacroBins][numMacroBins] = tmpCost;
						dynamicNumPoints[numFirstMacroBins][numMacroBins] = tmpTotalNumPoints;
						optimalPosition = maxIndex;
						optimalNewBin = tmpMergeMacroBinRight;
					} // end if
				} // end for
				
				// update the table's entry with the binning yielding the smallesy cost
				dynamicBinnings[numFirstMacroBins][numMacroBins] = new ArrayList<MacroBin>();
				for (int i = 0; i < numMacroBins - 1; i++)
				{
					tmpMergeMacroBin = cloneBin(dynamicBinnings[optimalPosition][numMacroBins - 1].get(i));
					dynamicBinnings[numFirstMacroBins][numMacroBins].add(tmpMergeMacroBin);
				}
				dynamicBinnings[numFirstMacroBins][numMacroBins].add(optimalNewBin);
				
				//System.out.println("numMacroBins = " + numMacroBins + " --- numFirstMacroBins = " + numFirstMacroBins);
			} // end for
		} // end for
		
		double maxCost = -Double.MAX_VALUE;
		int maxNumMacroBins = -1;
		double normalizedCost;
		for (int numMacroBins = 2; numMacroBins <= MAX_BINS; numMacroBins++)
		{
			val = Math.log(numMacroBins) / logBase;
			normalizedCost = dynamicBinningCosts[numTotalMacroBins][numMacroBins] / (sumLog + val - Math.max(maxLog, val));
			if (normalizedCost < 0 || normalizedCost > 1)
				throw new Exception("out of range cost = " + normalizedCost);
			
			if (normalizedCost > maxCost)
			{
				maxCost = normalizedCost;
				maxNumMacroBins = numMacroBins;
			}
		}
		for (int i = 0; i < maxNumMacroBins; i++)
			ret.add(dynamicBinnings[numTotalMacroBins][maxNumMacroBins].get(i));
		cost[0] = maxCost;
		
		System.out.println("dim " + curDim + ": number of bins = " + maxNumMacroBins);
		System.out.println("Cost = " + maxCost);
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
		for (int i = 0; i < a.numPoints; i++)
			ret.pointIDs.add(new Integer(a.pointIDs.get(i)));
		
		// get the IDs of micro bins
		int numMicroBinsA = a.microBinIDs.size();
		for (int i = 0; i < numMicroBinsA; i++)
			ret.microBinIDs.add(new Integer(a.microBinIDs.get(i)));
		
		return ret;
	}
	
	// merge two macro bins
	public static MacroBin mergeMacroBins(MacroBin a, MacroBin b, ArrayList<MicroBin> initialBins) throws Exception
	{	
		// init the new macro bin
		if (a.lowerBound > b.upperBound)
			throw new Exception("Invalid bin merge");
		
		MacroBin ret = new MacroBin(a.lowerBound, b.upperBound);
		
		// get the total number of data points
		ret.numPoints = a.numPoints + b.numPoints;
		ret.pointIDs = a.pointIDs;
		//for (int i = 0; i < a.numPoints; i++)
		//	ret.pointIDs.add(new Integer(a.pointIDs.get(i).intValue()));
		for (int i = 0; i < b.numPoints; i++)
			ret.pointIDs.add(new Integer(b.pointIDs.get(i)));
		
		// get the IDs of micro bins
		//int numMicroBinsA = a.microBinIDs.size();
		int numMicroBinsB = b.microBinIDs.size();
		ret.microBinIDs = a.microBinIDs;
		//for (int i = 0; i < numMicroBinsA; i++)
		//	ret.microBinIDs.add(new Integer(a.microBinIDs.get(i).intValue()));
		for (int i = 0; i < numMicroBinsB; i++)
			ret.microBinIDs.add(new Integer(b.microBinIDs.get(i)));
		
		return ret;
	}
	
	public static double computeEncodingCost(DataMatrix dataMatrix, int[] procDims, int[][] discreteData, int curDim, ArrayList<MacroBin> bins, int numPointsSoFar, int N, ArrayList<MicroBin> initialBins, ArrayList<ArrayList<Integer>> cells, ArrayList<ArrayList<Integer>> cellPointIDs) throws Exception
	{
		// the total coding cost
		double ret = 0;
		double logBase = Math.log(Constants.LOG_BASE);
		
		// get the total number of macro bins
		int numBins = bins.size();
		
		// iterate through each macro bins
		MacroBin tmpBin = null;
		int numCells = cellPointIDs.size();
		//ArrayList<Integer> curCellPointIDs = null;
		int binNumPoints;
		int pointID;
		//nt cellPointID;
		//int countJoin;
		//int count;
		int totalPoints;
		int[] counts = new int[numCells];
		//ArrayList<Integer> curCellPointIDs = null;
		ArrayList<Integer> curCell = null;
		ArrayList<Integer> cellIDs = new ArrayList<Integer>();
		int curBinID;
		int cellBinID;
		boolean match;
		int numProcessedDims = procDims.length;
		double correlation = 0;
		for (int i = 0; i < numBins; i++)
		{
			// get the macro bin
			tmpBin = bins.get(i);
			binNumPoints = tmpBin.numPoints;
			correlation -= (binNumPoints * 1.0 / numPointsSoFar) * Math.log((binNumPoints * 1.0 / numPointsSoFar)) / logBase;
			
			// clear point counts
			for (int j = 0; j < numCells; j++)
				counts[j] = 0;
			
			// compute the entropy of tmpMergeMacroBinRight
			totalPoints = 0;
			for (int h = 0; h < binNumPoints; h++)
			{
				pointID = tmpBin.pointIDs.get(h);
				cellIDs.clear();
				for (int k = 0; k < numProcessedDims; k++)
					cellIDs.add(new Integer(discreteData[pointID][procDims[k]]));
				
				for (int j = 0; j < numCells; j++)
				{
					curCell = cells.get(j);
					match = true;
					for (int k = 0; k < numProcessedDims; k++)
					{
						curBinID = cellIDs.get(k);
						cellBinID = curCell.get(k);
						if (curBinID != cellBinID)
						{
							match = false;
							break;
						}
					} // end for
					
					if (match == true)
					{
						counts[j]++;
						break;
					}
				} // end for
			} // end for
			
			// increase correlation
			for (int j = 0; j < numCells; j++)
			{
				totalPoints += counts[j];
				if (counts[j] > 0)
					correlation += (counts[j] * 1.0 / numPointsSoFar) * Math.log(counts[j] * 1.0 / numPointsSoFar) / logBase;
			}
			
			// this serves as control check
			if (totalPoints != binNumPoints)
			{
				System.out.println(totalPoints + " --- " + binNumPoints);
				throw new Exception("total points not match");
			}
			
			// increase correlation
			/*numCells = cellPointIDs.size();
			totalPoints = 0;
			for (int j = 0; j < numCells; j++)
			{
				curCellPointIDs = cellPointIDs.get(j);
				count = curCellPointIDs.size();
				totalPoints += count;
				countJoin = 0;
				for (int h = 0; h < binNumPoints; h++)
				{
					pointID = tmpBin.pointIDs.get(h);
					for (int k = 0; k < count; k++)
					{
						cellPointID = curCellPointIDs.get(k);
						if (pointID == cellPointID)
						{
							countJoin++;
							break;
						}
					} // end for
				} // end for
				
				if (countJoin > 0)
					correlation += (countJoin * 1.0 / numPointsSoFar) * Math.log(countJoin * 1.0 / numPointsSoFar) / logBase;
			} // end for*/
			
			/*if (totalPoints != N)
			{
				System.out.println(totalPoints + " --- " + N);
				throw new Exception("total points not match");
			}*/
		} // end for
		
		// get the total coding cost
		//System.out.println(LB + " --- " + LBinSize + " --- " + LBinEntropy + " --- " + LBinCorrelation);
		ret = correlation;
		return ret;
	}
	
	// divide a dimension into equal-frequency bins
	public static ArrayList<MicroBin> findEqualFrequencyBinning(int curDim, DataMatrix dataMatrix, int numDesiredBins, SortedObject[] tmp)
	{
		ArrayList<MicroBin> ret = new  ArrayList<MicroBin>();
		int numRows = dataMatrix.rows;
		int numCols = dataMatrix.cols;
		//int T = distinctValues.size();
		
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
		
		// number of distinct values per bin
		int totalNumBins = numDesiredBins;
		if (numDesiredBins > numRows)
			totalNumBins = numRows;
		int binCount = (int)Math.floor(numRows * 1.0 / totalNumBins);
		
		for (int i = 0; i < totalNumBins; i++)
		{
			// lower and upper indices of distinct values
			lowerIndex = i * binCount;
			if (i < totalNumBins - 1)
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
		} // end for
		
		return ret;
	}
	
	public static ArrayList<Double> computeContrastSpecial(DataMatrix dataMatrix, ArrayList<Integer> dims, ArrayList<Integer> frontDims, double[] CRES, ArrayList<ArrayList<Integer>> cellPointIDs) throws Exception
	{
		int rows = dataMatrix.rows;
		int curNumPoints;
		double[] vals = null;
		int curDim;
		int numFrontDims = frontDims.size();
		int numCells = cellPointIDs.size();
		ArrayList<Integer> curCell = null;
		double ce;
		double maxCE = -Double.MAX_VALUE;
		int maxPos = -1;
		for (int i = 0; i < numFrontDims; i++)
		{
			curDim = frontDims.get(i);
			ce = 0;
			for (int j = 0; j < numCells; j++)
			{
				curCell = cellPointIDs.get(j);
				curNumPoints = curCell.size();
				vals = new double[curNumPoints];
				for (int k = 0; k < curNumPoints; k++)
					vals[k] = dataMatrix.data.get(curCell.get(k)).measures[curDim];
			
				ce += curNumPoints * computeCRE(vals, false) / rows;
			} // end for
			
			ce = Math.abs(CRES[curDim] - ce);
			if (ce > maxCE)
			{
				maxCE = ce;
				maxPos = i;
			}
		}
		
		int selectedDim = frontDims.get(maxPos);
		dims.add(new Integer(selectedDim));
		frontDims.remove(maxPos);
		ArrayList<Double> ret = new ArrayList<Double>();
		ret.add(new Double(selectedDim));
		ret.add(new Double(Math.max(maxCE, 0)));
		
		return ret;
	}
	
	public static ArrayList<Double> computeContrastSpecial(DataMatrix dataMatrix, ArrayList<Integer> dims, ArrayList<Integer> frontDims, double[] CRES, int NumSeeds, double[] means, double[] devs) throws Exception
	{
		int num_dims = dims.size();
		int rows = dataMatrix.rows;
		
		ArrayList<Integer> dimsForClustering = new ArrayList<Integer>();
		ArrayList<Integer>[] clusters = null;
		ArrayList<Integer> curCluster = null;
		int curNumPoints;
		double[] vals = null;
		for (int i = 0; i < num_dims; i++)
			dimsForClustering.add(new Integer(dims.get(i)));
		
		clusters = incrementalClusteringNewSeed(dataMatrix, dimsForClustering, NumSeeds, means, devs);
		
		int curDim;
		int numFrontDims = frontDims.size();
		double ce;
		double maxCE = -Double.MAX_VALUE;
		int maxPos = -1;
		for (int i = 0; i < numFrontDims; i++)
		{
			curDim = frontDims.get(i);
			ce = 0;
			for (int j = 0; j < NumSeeds; j++)
			{
				curCluster = clusters[j];
				curNumPoints = curCluster.size();
				vals = new double[curNumPoints];
				for (int k = 0; k < curNumPoints; k++)
					vals[k] = dataMatrix.data.get(curCluster.get(k)).measures[curDim];
			
				ce += curNumPoints * computeCRE(vals, false) / rows;
			} // end for
			
			ce = Math.abs(CRES[curDim] - ce);
			if (ce > maxCE)
			{
				maxCE = ce;
				maxPos = i;
			}
		}
		
		int selectedDim = frontDims.get(maxPos);
		dims.add(new Integer(selectedDim));
		frontDims.remove(maxPos);
		ArrayList<Double> ret = new ArrayList<Double>();
		ret.add(new Double(selectedDim));
		ret.add(new Double(Math.max(maxCE, 0)));
		
		return ret;
	}
	
	public static double[][] computeContrasts(DataMatrix dataMatrix, ArrayList<Integer> dims, ArrayList<Integer> frontDims, double[] CRES, int NumSeeds, double[] means, double[] devs) throws Exception
	{
		int num_dims = dims.size();
		int rows = dataMatrix.rows;
		
		ArrayList<Integer> dimsForClustering = new ArrayList<Integer>();
		ArrayList<Integer>[] clusters = null;
		ArrayList<Integer> curCluster = null;
		int curNumPoints;
		double[] vals = null;
		for (int i = 0; i < num_dims; i++)
			dimsForClustering.add(new Integer(dims.get(i)));
		
		clusters = incrementalClusteringNewSeed(dataMatrix, dimsForClustering, NumSeeds, means, devs);
		
		int curDim;
		int numFrontDims = frontDims.size();
		double ce;
		double[][] ret = new double[2][numFrontDims];
		for (int i = 0; i < numFrontDims; i++)
		{
			curDim = frontDims.get(i);
			ce = 0;
			for (int j = 0; j < NumSeeds; j++)
			{
				curCluster = clusters[j];
				curNumPoints = curCluster.size();
				vals = new double[curNumPoints];
				for (int k = 0; k < curNumPoints; k++)
					vals[k] = dataMatrix.data.get(curCluster.get(k)).measures[curDim];
			
				ce += curNumPoints * computeCRE(vals, false) / rows;
			} // end for
			
			ce = Math.abs(CRES[curDim] - ce);
			ret[0][i] = curDim;
			ret[1][i] = ce;
		}
		
		return ret;
	}
	
	public static double computeCRE(double[] vals, boolean hasSorted)
	{
		int num_items = vals.length;
		if (hasSorted == false)
			Arrays.sort(vals);
		
		double cre = 0;
		double logBase = Math.log(Constants.LOG_BASE);
		for (int i = 0; i < num_items - 1; i++)
			if (vals[i + 1] != vals[i])
				cre += (vals[i + 1] - vals[i]) * ((i + 1) / (1.0 * num_items)) * Math.log((i + 1) / (1.0 * num_items)) / logBase;
		
		return -cre;
	}
	
	@SuppressWarnings("unchecked")
	public static ArrayList<Integer>[] incrementalClusteringNewSeed(DataMatrix dataMatrix, ArrayList<Integer> clusteringDims, int NumSeeds, double[] means, double[] devs)
	{
		ArrayList<Integer>[] tmpRet = new ArrayList[NumSeeds];
		for (int i = 0; i < NumSeeds; i++)
			tmpRet[i] = new ArrayList<Integer>();
		
		DataPoint[] tmpSeeds = new DataPoint[NumSeeds];
		int rows = dataMatrix.rows;
		int numDimsClustering = clusteringDims.size();
		DataPoint curPoint = null;
		int curDim;
		for (int i = 0; i < NumSeeds; i++)
		{
			tmpSeeds[i] = new DataPoint(numDimsClustering, 0);
			for (int j = 0; j < numDimsClustering; j++)
			{
				curDim = clusteringDims.get(j);
				if (j % 2 == 0)
					tmpSeeds[i].measures[j] = means[curDim] + i * devs[curDim] / (numDimsClustering * (NumSeeds - 1));
				else
					tmpSeeds[i].measures[j] = means[curDim] - i * devs[curDim] / (numDimsClustering * (NumSeeds - 1));
			}
		}
		
		DataPoint curSeed = null;
		double min;
		int minSeedId;
		double dist;
		int numMSteps = (int)Math.sqrt(rows);
		for (int iteration = 0; iteration < 3; iteration++)
		{
			for (int i = 0; i < NumSeeds; i++)
				tmpRet[i].clear();
			
			for (int i = 0; i < rows; i++)
			{
				curPoint = dataMatrix.data.get(i);
				min = Double.MAX_VALUE;
				minSeedId = -1;
				for (int j = 0; j < NumSeeds; j++)
				{
					curSeed = tmpSeeds[j];
					dist = DataPoint.distanceLNorm(Constants.LP_NORM, clusteringDims, curPoint, curSeed);
					if (dist < min)
					{
						min = dist;
						minSeedId = j;
					} // end if
				} // end for
				
				tmpRet[minSeedId].add(new Integer(i));
				
				if (iteration == 0 && i % numMSteps == 0)
					performMStep(dataMatrix, clusteringDims, tmpSeeds, tmpRet);
			} // end for
			
			performMStep(dataMatrix, clusteringDims, tmpSeeds, tmpRet);
		} // end for
		
		return tmpRet;
	}
	
	public static void performMStep(DataMatrix dataMatrix, ArrayList<Integer> clusteringDims, DataPoint[] seeds, ArrayList<Integer>[] members)
	{
		int numSeeds = seeds.length;
		int numPoints;
		DataPoint curPoint  = null;
		int numDimsClustering = clusteringDims.size();
		int curDim;
		for (int i = 0; i < numSeeds; i++)
		{
			numPoints = members[i].size();
			for (int j = 0; j < numPoints; j++)
			{
				curPoint = dataMatrix.data.get(members[i].get(j));
				for (int k = 0; k < numDimsClustering; k++)
				{
					curDim = clusteringDims.get(k);
					if (j == 0)
						seeds[i].measures[k] = curPoint.measures[curDim] / numPoints;
					else
						seeds[i].measures[k] += curPoint.measures[curDim] / numPoints;
				}
			} // end for
		} // end for
	}
}
