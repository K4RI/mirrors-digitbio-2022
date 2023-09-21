package multibinning.business;

import java.util.*;

import multibinning.data.*;

public class JBinning 
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
		
		// choose two most correlated dimensions
		int numDims = dims.size();
		ArrayList<Integer> tmpDims = new ArrayList<Integer>();
		ArrayList<Integer> tmpFrontDims = new ArrayList<Integer>();
		double maxPairContrast = -Double.MAX_VALUE;
		ArrayList<Integer> maxPair = new ArrayList<Integer>();
		ArrayList<Double> tmpResult;
		double tmpContrast;
		for (int i = 0; i < numDims; i++)
			for (int j = 0; j < numDims; j++)
				if (j != i)
				{
					tmpDims.clear();
					tmpFrontDims.clear();
					tmpDims.add(new Integer(dims.get(i)));
					tmpFrontDims.add(new Integer(dims.get(j)));
					tmpResult = computeContrastSpecial(dataMatrix, tmpDims, tmpFrontDims, CRES, Constants.NUM_CENTROIDS, means, devs);
					tmpContrast = tmpResult.get(1);
					if (tmpContrast > maxPairContrast)
					{
						maxPairContrast = tmpContrast;
						maxPair.clear();
						maxPair.add(new Integer(tmpDims.get(0)));
						maxPair.add(new Integer(tmpDims.get(1)));
					}
				}
		
		// discretize firstDim
		ArrayList<Integer> processedDims = new ArrayList<Integer>();
		ArrayList<MacroBin>[] existingMacroBins = new ArrayList[numCols];
		int firstDim = maxPair.get(0);
		existingMacroBins[firstDim] = findBinning(firstDim, dataMatrix, Constants.METHOD, processedDims, existingMacroBins, Constants.DISC_DATA);
		updateDiscreteData(firstDim, existingMacroBins[firstDim], Constants.DISC_DATA);
		processedDims.add(new Integer(firstDim));
		
		// discretize secondDim
		int secondDim = maxPair.get(1);
		existingMacroBins[secondDim] = findBinning(secondDim, dataMatrix, Constants.METHOD, processedDims, existingMacroBins, Constants.DISC_DATA);
		updateDiscreteData(secondDim, existingMacroBins[secondDim], Constants.DISC_DATA);
		processedDims.add(new Integer(secondDim));
		
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
		while (curLength < numDims)
		{
			tmpResult = computeContrastSpecial(dataMatrix, maxPair, remainingDims, CRES, Constants.NUM_CENTROIDS, means, devs);
			tmpDim = tmpResult.get(0).intValue();
			
			// discretize tmpDim
			existingMacroBins[tmpDim] = findBinning(tmpDim, dataMatrix, Constants.METHOD, processedDims, existingMacroBins, Constants.DISC_DATA);
			updateDiscreteData(tmpDim, existingMacroBins[tmpDim], Constants.DISC_DATA);
			processedDims.add(new Integer(tmpDim));
			
			curLength++;
		}
		
		return existingMacroBins;
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
	
	// find binning strategy for each dimension
	public static ArrayList<MacroBin> findBinning(int curDim, DataMatrix dataMatrix, int method, ArrayList<Integer> processedDims, ArrayList<MacroBin>[] existingMacroBins, int[][] discreteData) throws Exception
	{
		System.out.println("Find binning for dimension " + curDim);
		 ArrayList<MacroBin> ret = null;
		 
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
		 } // end for
		 
		 // divide data into equal-frequency micro bins
		 int initBinCount = Constants.INIT_BIN_COUNT;
		 ArrayList<MicroBin> initialMicroBins = findEqualWidthBinning(curDim, dataMatrix, initBinCount, distinctValues, distinctValueSupports, memberPointIDs);
		 
		 // get the number of micro bins
		 int numMicroBins = initialMicroBins.size();
		 System.out.println("num of initial bins = " + numMicroBins);
		 
		 // for each micro bin, create a macro bin containing it
		 ArrayList<MacroBin> candidateMacroBins = new ArrayList<MacroBin>();
		 MacroBin tmpMacroBin = null;
		 MicroBin tmpMicroBin = null;
		 for (int i = 0; i < numMicroBins; i++)
		 {
			 tmpMicroBin = initialMicroBins.get(i);
			 tmpMacroBin = new MacroBin(tmpMicroBin.lowerBound, tmpMicroBin.upperBound);
			 tmpMacroBin.numPoints = tmpMicroBin.allValues.size();
			 for (int j = 0; j < tmpMacroBin.numPoints; j++)
				 tmpMacroBin.pointIDs.add(new Integer(tmpMicroBin.pointIDs.get(j)));
			 tmpMacroBin.microBinIDs.add(new Integer(i));
			 tmpMacroBin.numMDHs = 0;
			 tmpMacroBin.mdh = new int[0];
			 candidateMacroBins.add(tmpMacroBin);
		 }
		 
		 if (numMicroBins == 0)
			 throw new Exception("Dimension with no candidate bin!");
			 
		 // mine the optimal macro bins by means of dynamic programming
		 if (Constants.METHOD != Constants.J_GD_HM)
			 ret = dynamicProgrammingMDLBinning(dataMatrix, curDim, processedDims, existingMacroBins, candidateMacroBins, numRows, candidateMacroBins.size(), initialMicroBins, discreteData);
		 else
			 ret = greedyMDLBinning(dataMatrix, curDim, processedDims, existingMacroBins, numRows, numMicroBins, initialMicroBins, discreteData);
		 return ret;
	}
	
	public static ArrayList<MacroBin> greedyMDLBinning(DataMatrix dataMatrix, int curDim, ArrayList<Integer> processedDims, ArrayList<MacroBin>[] existingMacroBins, int N, int T, ArrayList<MicroBin> initialBins, int[][] discreteData) throws Exception
	{
		int numProcessedDims = processedDims.size();
		int[] procDims = new int[numProcessedDims];
		for (int i = 0; i < numProcessedDims; i++)
			procDims[i] = processedDims.get(i);
		
		ArrayList<Integer> remainingCPs = new ArrayList<Integer>();
		for (int i = 1; i < T - 1; i++)
			remainingCPs.add(new Integer(i));
		
		int numRemainingValues = T - 2;
		double GlobalMinCost = Double.MAX_VALUE;
		double tmpCost;
		double tmpMinCost;
		int tmpMinIndex;
		int tmpMinPos;
		int tmpNumBins;
		int candidateCP;
		int pos;
		int tmpLowerBound;
		MacroBin tmpBin = null;
		ArrayList<MacroBin> GlobalBins = new ArrayList<MacroBin>();
		ArrayList<MacroBin> tmpBins = new ArrayList<MacroBin>();
		ArrayList<Integer> cutPoints = new ArrayList<Integer>();
		cutPoints.add(new Integer(T - 1));
		int numCPs = 1;
		int curCP;
		MicroBin tmpMicroBin = null;
		int tmpNumPoints;
		while (numRemainingValues > 0)
		{
			tmpMinCost = Double.MAX_VALUE;
			tmpMinIndex = -1;
			tmpMinPos = -1;
			for (int i = 0; i < numRemainingValues; i++)
			{
				candidateCP = remainingCPs.get(i);
				pos = numCPs;
				for (int j = numCPs - 1; j >= 0; j--)
				{
					//System.out.println(j + " -- " + cutPoints.size() + " --- " + cutPoints.get(j).doubleValue());
					if (candidateCP >= cutPoints.get(j))
						break;
					else
						pos--;
				}
				
				cutPoints.add(pos, new Integer(candidateCP));
				numCPs++;
				tmpBins.clear();
				tmpLowerBound = 0;
				for (int j = 0; j < numCPs; j++)
				{
					curCP = cutPoints.get(j);
					tmpBin = new MacroBin(initialBins.get(tmpLowerBound).lowerBound, initialBins.get(curCP).upperBound);
					
					for (int k = tmpLowerBound; k <= curCP; k++)
					{
						tmpBin.microBinIDs.add(new Integer(k));
						tmpMicroBin = initialBins.get(k);
						tmpNumPoints = tmpMicroBin.pointIDs.size();
						for (int h = 0; h < tmpNumPoints; h++)
							tmpBin.pointIDs.add(new Integer(tmpMicroBin.pointIDs.get(h)));
					}
					
					tmpBin.numPoints = tmpBin.pointIDs.size();
					tmpBins.add(tmpBin);
					tmpLowerBound = curCP + 1;
				}
				
				tmpCost = computeEncodingCost(dataMatrix, curDim, procDims, tmpBins, existingMacroBins, N, T, initialBins, discreteData);
				if (tmpCost < tmpMinCost)
				{
					tmpMinCost = tmpCost;
					tmpMinIndex = i;
					tmpMinPos = pos;
				}
				
				cutPoints.remove(pos);
				numCPs--;
			}
			cutPoints.add(tmpMinPos, remainingCPs.get(tmpMinIndex).intValue());
			numCPs++;
			remainingCPs.remove(tmpMinIndex);
			//System.out.println("------");
			
			tmpNumBins = tmpBins.size();
			if (tmpMinCost < GlobalMinCost)
			{
				GlobalMinCost = tmpMinCost;
				GlobalBins.clear();
				for (int i = 0; i < tmpNumBins; i++)
					GlobalBins.add(tmpBins.get(i));
			}
			else
				break;
			
			numRemainingValues--;
		}
		
		return GlobalBins;
	}
	
	@SuppressWarnings("unchecked")
	public static ArrayList<MacroBin> dynamicProgrammingMDLBinning(DataMatrix dataMatrix, int curDim, ArrayList<Integer> processedDims, ArrayList<MacroBin>[] existingMacroBins, ArrayList<MacroBin> candidateMacroBins, double N, double T, ArrayList<MicroBin> initialBins, int[][] discreteData) throws Exception
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		
		double logBase = Math.log(Constants.LOG_BASE);
		int numProcessedDims = processedDims.size();
		int[] procDims = new int[numProcessedDims];
		for (int i = 0; i < numProcessedDims; i++)
			procDims[i] = processedDims.get(i);
		
		// create the table for dynamic programming to store already solved sub-problems
		int numTotalMacroBins = candidateMacroBins.size();
		ArrayList<MacroBin>[][] dynamicBinnings = new ArrayList[numTotalMacroBins + 1][numTotalMacroBins + 1];
		double[][] dynamicBinningCosts = new double[numTotalMacroBins + 1][numTotalMacroBins + 1];
		int[][] dynamicNumPoints = new int[numTotalMacroBins + 1][numTotalMacroBins + 1];
		
		// find optimal binning of size 2
		MacroBin tmpMergeMacroBin = null;
		MacroBin tmpMergeMacroBinRight = null;
		double tmpCost;
		ArrayList<MacroBin> tmpMacroBins = new ArrayList<MacroBin>();
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
					tmpMergeMacroBin = mergeMacroBins(tmpMergeMacroBin, candidateMacroBins.get(i), initialBins);
				
				// merge all macro bins from position maxIndex to (numFirstMacroBins - 1)
				tmpMergeMacroBinRight = cloneBin(candidateMacroBins.get(maxIndex));
				for (int i = maxIndex + 1; i < numFirstMacroBins; i++)
					tmpMergeMacroBinRight = mergeMacroBins(tmpMergeMacroBinRight, candidateMacroBins.get(i), initialBins);
				
				// add two new macro bins to the temporary binning
				tmpMacroBins.add(tmpMergeMacroBin);
				tmpMacroBins.add(tmpMergeMacroBinRight);
				
				// compute the coding cost of the temporary binning
				tmpCost = computeEncodingCost(dataMatrix, curDim, procDims, tmpMacroBins, existingMacroBins, tmpMergeMacroBin.numPoints + tmpMergeMacroBinRight.numPoints, tmpMergeMacroBin.microBinIDs.size() + tmpMergeMacroBinRight.microBinIDs.size(), initialBins, discreteData);
				
				// if the cost is less than the current cost, then update the cost, the binning, and the total number of points of the binning
				if (tmpCost < dynamicBinningCosts[numFirstMacroBins][2])
				{
					dynamicBinningCosts[numFirstMacroBins][2] = tmpCost;
					//dynamicNumPoints[numFirstMacroBins][2] = tmpMergeMacroBin.numPoints + tmpMergeMacroBinRight.numPoints;
					dynamicNumPoints[numFirstMacroBins][2] = tmpMergeMacroBin.numPoints + tmpMergeMacroBinRight.numPoints;
					if (dynamicBinnings[numFirstMacroBins][2] == null)
						dynamicBinnings[numFirstMacroBins][2] = new ArrayList<MacroBin>();
					else
						dynamicBinnings[numFirstMacroBins][2].clear();
					dynamicBinnings[numFirstMacroBins][2].add(tmpMergeMacroBin);
					dynamicBinnings[numFirstMacroBins][2].add(tmpMergeMacroBinRight);
				}
			}
			
			System.out.println(numFirstMacroBins);
		}
		
		// build the rest of the table
		// for each (numMacroBins, numFirstMacroBins) pairs, find the binning that yields the smallest cost
		int tmpTotalNumPoints;
		int tmpNumPoints;
		//int tmpNumMacroBins;
		int optimalPosition;
		MacroBin optimalNewBin = null;
		int tmpPointID;
		int tmpPointIDOther;
		int tmpCol;
		//MacroBin tmpMacroBin = null;
		int countJoin;
		int count;
		//int countWithoutNew;
		DataPoint tmpPoint = null;
		boolean eligible;
		//ArrayList<Integer> allPointIDs = new ArrayList<Integer>();
		int totalNumPoints;
		//for (int numMacroBins = 3; numMacroBins <= numTotalMacroBins; numMacroBins++)
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
						tmpMergeMacroBinRight = mergeMacroBins(tmpMergeMacroBinRight, candidateMacroBins.get(i), initialBins);
					
					// get all point IDs
					/*allPointIDs.clear();
					tmpNumMacroBins = dynamicBinnings[maxIndex][numMacroBins - 1].size();
					for (int i = 0; i < tmpNumMacroBins; i++)
					{
						tmpMacroBin = dynamicBinnings[maxIndex][numMacroBins - 1].get(i);
						tmpNumPoints = tmpMacroBin.pointIDs.size();
						for (int j = 0; j < tmpNumPoints; j++)
							allPointIDs.add(new Integer(tmpMacroBin.pointIDs.get(j)));
					}
					for (int i = 0; i < tmpMergeMacroBinRight.numPoints; i++)
						allPointIDs.add(new Integer(tmpMergeMacroBinRight.pointIDs.get(i)));*/
					
					// compute the new binning's cost based on stored information
					tmpTotalNumPoints = dynamicNumPoints[maxIndex][numMacroBins - 1] + tmpMergeMacroBinRight.numPoints;
					
					// adjust LB
					tmpCost += LInteger(numMacroBins - 1) + logOfCombination(numFirstMacroBins - 1, numMacroBins - 1);
					tmpCost -= LInteger(numMacroBins - 2) + logOfCombination(maxIndex - 1, numMacroBins - 2);
					
					// adjust LBinSize
					tmpCost += LInteger(tmpMergeMacroBinRight.numPoints);
					
					// adjust LBinEntropy
					if (Constants.USE_CE == true)
						tmpCost += computeBinInformationContentWithCE(tmpMergeMacroBinRight, initialBins);
					else
						tmpCost += computeBinInformationContentWithoutCE(tmpMergeMacroBinRight, initialBins);
					
					// adjust LBinCorrelation
					if (numProcessedDims > 0)
					{
						//totalNumPoints = allPointIDs.size();
						totalNumPoints = dataMatrix.rows;
						//tmpNumMacroBins = dynamicBinnings[maxIndex][numMacroBins - 1].size();
						
						// for each macro bin in the discretization of the first maxIndex points, recompute its LBinCorrelation
						/*for (int i = 0; i < tmpNumMacroBins; i++)
						{
							tmpMacroBin = dynamicBinnings[maxIndex][numMacroBins - 1].get(i);
							tmpNumPoints = tmpMacroBin.pointIDs.size();
							
							// loop through each point of the bin
							for (int j = 0; j < tmpNumPoints; j++)
							{
								tmpPointID = tmpMacroBin.pointIDs.get(j);
								count = 0;
								countWithoutNew = 0;
								
								// loop through each point of the points considered so far
								for (int row = 0; row < totalNumPoints; row++)
								{
									eligible = true;
									tmpPointIDOther = allPointIDs.get(row);
									
									// process depending on if hidden markov is used
									if (Constants.METHOD == Constants.J_MDL)
									{
										for (int col = 0; col < numProcessedDims; col++)
										{
											tmpCol = procDims[col];
											if (discreteData[tmpPointIDOther][tmpCol] != discreteData[tmpPointID][tmpCol])
											{
												eligible = false;
												break;
											}
										}
									}
									else
									{
										tmpCol = procDims[numProcessedDims - 1];
										if (discreteData[tmpPointIDOther][tmpCol] != discreteData[tmpPointID][tmpCol])
											eligible = false;
									}
									
									// if the point satisfies all conditions then increase the counts
									if (eligible)
									{
										count++;
										if (row < totalNumPoints - tmpMergeMacroBinRight.numPoints)
											countWithoutNew++;
									}
								} // end for
								
								tmpCost -= Math.log(countWithoutNew * 1.0 / count) / logBase;
							} // end for
						} // end for
						*/
						
						// loop through each point of the new macro bin
						tmpNumPoints = tmpMergeMacroBinRight.pointIDs.size();
						for (int k = 0; k < tmpNumPoints; k++)
						{
							tmpPointID = tmpMergeMacroBinRight.pointIDs.get(k);
							countJoin = 0;
							count = 0;
							
							// loop through each point of the points considered so far
							if (Constants.METHOD == Constants.J_MDL)
							{
								for (int row = 0; row < totalNumPoints; row++)
								{
									eligible = true;
									//tmpPointIDOther = allPointIDs.get(row);
									tmpPointIDOther = row;
									
									// process each already discretized dimensions
									for (int col = 0; col < numProcessedDims; col++)
									{
										tmpCol = procDims[col];
										if (discreteData[tmpPointIDOther][tmpCol] != discreteData[tmpPointID][tmpCol])
										{
											eligible = false;
											break;
										}
									}
									
									// if the point satisfies all conditions then increase the counts
									if (eligible)
									{
										count++;
										tmpPoint = dataMatrix.data.get(tmpPointIDOther);
										if (tmpPoint.measures[curDim] >= tmpMergeMacroBinRight.lowerBound && tmpPoint.measures[curDim] <= tmpMergeMacroBinRight.upperBound)
											countJoin++;
									}
								} // end for
							}
							else
							{
								tmpCol = procDims[numProcessedDims - 1];
								for (int row = 0; row < tmpNumPoints; row++)
								{
									tmpPointIDOther = tmpMergeMacroBinRight.pointIDs.get(row);
									if (discreteData[tmpPointIDOther][tmpCol] == discreteData[tmpPointID][tmpCol])
										countJoin++;
								} // end for
								
								count = existingMacroBins[tmpCol].get(discreteData[tmpPointID][tmpCol]).numPoints;
							}
							
							tmpCost -= Math.log(countJoin * 1.0 / count) / logBase;
						} // end for
					}
					else
					{
						tmpCost -= tmpMergeMacroBinRight.numPoints * Math.log(tmpMergeMacroBinRight.numPoints * 1.0 / tmpTotalNumPoints) / logBase;
						tmpCost -= dynamicNumPoints[maxIndex][numMacroBins - 1] * Math.log(dynamicNumPoints[maxIndex][numMacroBins - 1] * 1.0 / tmpTotalNumPoints) / logBase;
					}
					
					//System.out.println(tmpCost);
					// if the new cost is less than the current cost, then update the binning
					if (tmpCost < dynamicBinningCosts[numFirstMacroBins][numMacroBins])
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
			} // end for
			
			System.out.println(numMacroBins);
		} // end for
		
		double minCost = Double.MAX_VALUE;
		int minNumMacroBins = -1;
		//for (int numMacroBins = 2; numMacroBins <= numTotalMacroBins; numMacroBins++)
		for (int numMacroBins = 2; numMacroBins <= numTotalMacroBins; numMacroBins++)
		{
			if (dynamicBinningCosts[numTotalMacroBins][numMacroBins] < minCost)
			{
				minCost = dynamicBinningCosts[numTotalMacroBins][numMacroBins];
				minNumMacroBins = numMacroBins;
			}
		}
		
		for (int i = 0; i < minNumMacroBins; i++)
			ret.add(dynamicBinnings[numTotalMacroBins][minNumMacroBins].get(i));
		
		System.out.println("Number of bins = " + minNumMacroBins);
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
		if (a.lowerBound >= b.upperBound)
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
	
	public static MicroBin mergeMicroBins(MicroBin a, MicroBin b)
	{	
		// init the new macro bin
		MicroBin ret = new MicroBin(a.lowerBound, b.upperBound);
		
		// combine the distinct values
		int numDistinctValuesA = a.distinctValues.size();
		int numDistinctValuesB = b.distinctValues.size();
		for (int i = 0; i < numDistinctValuesA; i++)
		{
			ret.distinctValues.add(new Double(a.distinctValues.get(i)));
			ret.distinctValueSupports.add(new Integer(a.distinctValueSupports.get(i)));
		}
		for (int i = 0; i < numDistinctValuesB; i++)
		{
			ret.distinctValues.add(new Double(b.distinctValues.get(i)));
			ret.distinctValueSupports.add(new Integer(b.distinctValueSupports.get(i)));
		}
		
		// combine all values of two bins
		int numValuesA = a.allValues.size();
		int numValuesB = b.allValues.size();
		for (int i = 0; i < numValuesA; i++)
		{
			ret.allValues.add(new Double(a.allValues.get(i)));
			ret.pointIDs.add(new Integer(a.pointIDs.get(i)));
		}
		for (int i = 0; i < numValuesB; i++)
		{
			ret.allValues.add(new Double(b.allValues.get(i)));
			ret.pointIDs.add(new Integer(b.pointIDs.get(i)));
		}
		
		return ret;
	}
	
	// T = total number of micro bins contained in all macro bins
	// N = total number of data points contained in DS
	public static double computeEncodingCost(DataMatrix dataMatrix, int curDim, int[] processedDims, ArrayList<MacroBin> bins, ArrayList<MacroBin>[] existingMacroBins, int N, int T, ArrayList<MicroBin> initialBins, int[][] discreteData) throws Exception
	{
		// the total coding cost
		double ret = 0;
		
		double logBase = Math.log(Constants.LOG_BASE);
		int numProcessedDims = processedDims.length;
		
		// get the total number of macro bins
		int numBins = bins.size();
		
		// LB
		double LB = LInteger(numBins - 1);
		if (T != 1)
			LB += logOfCombination(T - 1, numBins - 1);
		
		// LBinSize
		double LBinSize = 0;
		
		// LBinEntropy
		double LBinEntropy = 0;
		
		// LBinCorrelation
		double LBinCorrelation = 0;
		
		// iterate through each macro bins
		MacroBin tmpBin = null;
		int tmpNumPoints;
		int tmpPointID;
		int tmpPointIDOther;
		int tmpCol;
		DataPoint tmpPoint = null;
		int countJoin;
		int count;
		boolean eligible;
		
		// get all point IDs
		/*ArrayList<Integer> allPointIDs = new ArrayList<Integer>();
		for (int i = 0; i < numBins; i++)
		{
			// get the macro bin
			tmpBin = bins.get(i);
			tmpNumPoints = tmpBin.pointIDs.size();
			for (int k = 0; k < tmpNumPoints; k++)
				allPointIDs.add(new Integer(tmpBin.pointIDs.get(k)));
		}*/
		
		//int totalNumPoints = allPointIDs.size();
		int totalNumPoints = dataMatrix.rows;
		for (int i = 0; i < numBins; i++)
		{
			// get the macro bin
			tmpBin = bins.get(i);
			
			// increase LBinSize
			LBinSize += LInteger(tmpBin.numPoints);
			
			// increase LBinEntropy
			if (Constants.USE_CE == true)
				LBinEntropy += computeBinInformationContentWithCE(tmpBin, initialBins);
			else
				LBinEntropy += computeBinInformationContentWithoutCE(tmpBin, initialBins);
			
			// increase LBinCorrelation
			if (numProcessedDims > 0)
			{
				tmpNumPoints = tmpBin.pointIDs.size();
				for (int k = 0; k < tmpNumPoints; k++)
				{
					tmpPointID = tmpBin.pointIDs.get(k);
					countJoin = 0;
					count = 0;
					
					// loop through each point of the points considered so far
					if (Constants.METHOD == Constants.J_MDL)
					{
						for (int row = 0; row < totalNumPoints; row++)
						{
							eligible = true;
							//tmpPointIDOther = allPointIDs.get(row);
							tmpPointIDOther = row;
							
							// process each already discretized dimensions
							for (int col = 0; col < numProcessedDims; col++)
							{
								tmpCol = processedDims[col];
								if (discreteData[tmpPointIDOther][tmpCol] != discreteData[tmpPointID][tmpCol])
								{
									eligible = false;
									break;
								}
							}
							
							// if the point satisfies all conditions then increase the counts
							if (eligible)
							{
								count++;
								tmpPoint = dataMatrix.data.get(tmpPointIDOther);
								if (tmpPoint.measures[curDim] >= tmpBin.lowerBound && tmpPoint.measures[curDim] <= tmpBin.upperBound)
									countJoin++;
							}
						} // end for
					}
					else
					{
						tmpCol = processedDims[numProcessedDims - 1];
						for (int row = 0; row < tmpNumPoints; row++)
						{
							tmpPointIDOther = tmpBin.pointIDs.get(row);
							if (discreteData[tmpPointIDOther][tmpCol] == discreteData[tmpPointID][tmpCol])
								countJoin++;
						} // end for
						
						count = existingMacroBins[tmpCol].get(discreteData[tmpPointID][tmpCol]).numPoints;
					}
					
					LBinCorrelation -= Math.log(countJoin * 1.0 / count) / logBase;
				} // end for
			}
			else
				LBinCorrelation -= tmpBin.numPoints * Math.log(tmpBin.numPoints * 1.0 / N) / logBase;
		} // end for
		
		// get the total coding cost
		ret = LB + LBinSize + LBinEntropy + LBinCorrelation;
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
		int numPoints = macroBin.numPoints;
		MicroBin tmpMicroBin = null;
		int tmpNumDistinctValues;
		int tmpSupport;
		for (int i = 0; i < numMicroBins; i++)
		{
			tmpMicroBin = initialBins.get(macroBin.microBinIDs.get(i));
			tmpNumDistinctValues = tmpMicroBin.distinctValues.size();
			for (int j = 0; j < tmpNumDistinctValues; j++)
			{
				tmpSupport = tmpMicroBin.distinctValueSupports.get(j);
				ret -= tmpSupport * Math.log(tmpSupport * 1.0 / numPoints) / logBase;
			}
		}
		
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
		MicroBin tmpMicroBin = null;
		int tmpNumDistinctValues;
		int tmpSupport;
		int curPointCount = 0;
		for (int i = 0; i < numMicroBins; i++)
		{
			tmpMicroBin = initialBins.get(macroBin.microBinIDs.get(i));
			tmpNumDistinctValues = tmpMicroBin.distinctValues.size();
			for (int j = 0; j < tmpNumDistinctValues; j++)
			{
				tmpSupport = tmpMicroBin.distinctValueSupports.get(j);
				curPointCount += tmpSupport;
				ret -= curPointCount * Math.log(curPointCount * 1.0 / numPoints) / logBase;
			}
		}
		
		return ret;
	}
	
	// divide a dimension into equal-frequency bins
	public static ArrayList<MicroBin> findEqualWidthBinning(int curDim, DataMatrix dataMatrix, int numDesiredBins, ArrayList<Double> distinctValues, ArrayList<Integer> distinctValueSupports, ArrayList<ArrayList<Integer>> memberpointIDs)
	{
		ArrayList<MicroBin> ret = new  ArrayList<MicroBin>();
		
		int T = distinctValues.size();
		double minValue = distinctValues.get(0);
		double maxValue = distinctValues.get(T - 1);
		
		// number of distinct values per bin
		double binWidth = (maxValue - minValue) / numDesiredBins;
		
		double lowerBound;
		double upperBound;
		MicroBin tmpBin = null;
		ArrayList<Integer> tmpPointIDs = null;
		int tmpPointCount;
		int tmpPointID;
		int curSupport;
		double curVal;
		if (numDesiredBins > 0)
		{
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
				
				// calculate the bin
				for (int j = 0; j < T; j++)
				{
					// get current value
					curVal = distinctValues.get(j);
					
					if (curVal >= lowerBound && curVal < upperBound)
					{
						// add the current value into the list of bin's distinct values
						tmpBin.distinctValues.add(new Double(curVal));
						
						// add the support of the new distinct value into the bin's list of supports
						tmpBin.distinctValueSupports.add(new Integer(distinctValueSupports.get(j)));
						
						// add values (including duplicates) into the bin's list of all values
						curSupport = distinctValueSupports.get(j);
						for (int k = 0; k < curSupport; k++)
							tmpBin.allValues.add(new Double(curVal));
						
						// add all respective data points to the bin's data matrix
						tmpPointIDs = memberpointIDs.get(j);
						tmpPointCount = tmpPointIDs.size();
						for (int k = 0; k < tmpPointCount; k++)
						{
							tmpPointID = tmpPointIDs.get(k);
							tmpBin.pointIDs.add(new Integer(tmpPointID));
						}
					}
				}
				
				if (tmpBin.distinctValues.size() > 0)
				{
					if (tmpBin.distinctValues.size() > 0)
					{
						tmpBin.lowerBound = tmpBin.distinctValues.get(0);
						tmpBin.upperBound = tmpBin.distinctValues.get(tmpBin.distinctValues.size() - 1);
					}
					else
					{
						tmpBin.lowerBound = lowerBound;
						tmpBin.upperBound = upperBound;
					}
					//System.out.println(tmpBin.upperBound);
					
					// add the bin to the list of bins
					ret.add(tmpBin);
				}
			}
		}
		else
		{
			// calculate the bin
			for (int i = 0; i < T; i++)
			{
				// get current value
				curVal = distinctValues.get(i);
				tmpBin = new MicroBin(curVal, curVal);
				
				// add the current value into the list of bin's distinct values
				tmpBin.distinctValues.add(new Double(curVal));
				
				// add the support of the new distinct value into the bin's list of supports
				tmpBin.distinctValueSupports.add(new Integer(distinctValueSupports.get(i)));
				
				// add values (including duplicates) into the bin's list of all values
				curSupport = distinctValueSupports.get(i);
				for (int k = 0; k < curSupport; k++)
					tmpBin.allValues.add(new Double(curVal));
				
				// add all respective data points to the bin's data matrix
				tmpPointIDs = memberpointIDs.get(i);
				tmpPointCount = tmpPointIDs.size();
				for (int k = 0; k < tmpPointCount; k++)
				{
					tmpPointID = tmpPointIDs.get(k);
					tmpBin.pointIDs.add(new Integer(tmpPointID));
				}
				
				// add the bin to the list of bins
				ret.add(tmpBin);
			}
		}
		
		//System.out.println("---end---");
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
}
