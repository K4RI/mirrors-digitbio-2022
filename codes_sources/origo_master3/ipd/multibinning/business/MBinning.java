package multibinning.business;

import java.util.*;

import multibinning.data.*;

public class MBinning 
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
		existingMacroBins[firstDim] = findBinning(firstDim, dataMatrix, Constants.METHOD, processedDims, existingMacroBins, Constants.DISC_DATA, cells, cellPointIDs);
		updateDiscreteData(firstDim, existingMacroBins[firstDim], Constants.DISC_DATA);
		updateCells(existingMacroBins[firstDim], cells, cellPointIDs);
		processedDims.add(new Integer(firstDim));
		
		// discretize secondDim
		int secondDim = maxPair.get(1);
		existingMacroBins[secondDim] = findBinning(secondDim, dataMatrix, Constants.METHOD, processedDims, existingMacroBins, Constants.DISC_DATA, cells, cellPointIDs);
		updateDiscreteData(secondDim, existingMacroBins[secondDim], Constants.DISC_DATA);
		updateCells(existingMacroBins[secondDim], cells, cellPointIDs);
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
			existingMacroBins[tmpDim] = findBinning(tmpDim, dataMatrix, Constants.METHOD, processedDims, existingMacroBins, Constants.DISC_DATA, cells, cellPointIDs);
			updateDiscreteData(tmpDim, existingMacroBins[tmpDim], Constants.DISC_DATA);
			updateCells(existingMacroBins[tmpDim], cells, cellPointIDs);
			processedDims.add(new Integer(tmpDim));
			
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
	
	// find binning strategy for each dimension
	public static ArrayList<MacroBin> findBinning(int curDim, DataMatrix dataMatrix, int method, ArrayList<Integer> processedDims, ArrayList<MacroBin>[] existingMacroBins, int[][] discreteData, ArrayList<ArrayList<Integer>> cells, ArrayList<ArrayList<Integer>> cellPointIDs) throws Exception
	{
		System.out.println("Find binning for dimension " + curDim);
		ArrayList<MacroBin> ret = null;
		 
		int numProcessedDims = processedDims.size();
		int[] procDims = new int[numProcessedDims];
		for (int i = 0; i < numProcessedDims; i++)
			procDims[i] = processedDims.get(i);
		 
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
		ArrayList<MicroBin> initialMicroBins = null;
		if (numProcessedDims == 0)
			initialMicroBins = findEqualWidthBinning(curDim, dataMatrix, 5, distinctValues, distinctValueSupports, memberPointIDs);
		else
			initialMicroBins = findEqualWidthBinning(curDim, dataMatrix, initBinCount, distinctValues, distinctValueSupports, memberPointIDs);
		 
		// get the number of micro bins
		int numMicroBins = initialMicroBins.size();
		System.out.println("num of initial bins = " + numMicroBins);
		 
		// for each micro bin, create a macro bin containing it
		ArrayList<MacroBin> candidateMacroBins = new ArrayList<MacroBin>();
		MacroBin tmpMacroBin = null;
		MicroBin tmpMicroBin = null;
		int pointID;
		DataPoint tmpPoint = null;
		DataPoint newPoint = null;
		for (int i = 0; i < numMicroBins; i++)
		{
			tmpMicroBin = initialMicroBins.get(i);
			tmpMacroBin = new MacroBin(tmpMicroBin.lowerBound, tmpMicroBin.upperBound);
			tmpMacroBin.numPoints = tmpMicroBin.allValues.size();
			tmpMacroBin.dataMatrix = new DataMatrix();
			tmpMacroBin.dataMatrix.data = new ArrayList<DataPoint>();
			tmpMacroBin.means = new double[numProcessedDims];
			tmpMacroBin.devs = new double[numProcessedDims];
			for (int j = 0; j < tmpMacroBin.numPoints; j++)
			{
				pointID = tmpMicroBin.pointIDs.get(j);
				tmpPoint = dataMatrix.data.get(pointID);
				newPoint = new DataPoint(numProcessedDims, 0);
				for (int k = 0; k < numProcessedDims; k++)
				{
					newPoint.measures[k] = tmpPoint.measures[procDims[k]];
					tmpMacroBin.means[k] += newPoint.measures[k];
				}
				tmpMacroBin.dataMatrix.data.add(newPoint);
				tmpMacroBin.pointIDs.add(new Integer(pointID));
			}
			for (int j = 0; j < numProcessedDims; j++)
			{
				tmpMacroBin.means[j] = tmpMacroBin.means[j] / tmpMacroBin.numPoints;
				if (tmpMacroBin.numPoints == 1)
					tmpMacroBin.devs[j] = 0;
				else
				{
					for (int k = 0; k < tmpMacroBin.numPoints; k++)
					{
						tmpPoint = tmpMacroBin.dataMatrix.data.get(k);
						tmpMacroBin.devs[j] += (tmpPoint.measures[j] - tmpMacroBin.means[j]) * (tmpPoint.measures[j] - tmpMacroBin.means[j]) / (tmpMacroBin.numPoints - 1);
					}
					tmpMacroBin.devs[j] = Math.sqrt(tmpMacroBin.devs[j]);
				}
			} // end for
			tmpMacroBin.microBinIDs.add(new Integer(i));
			candidateMacroBins.add(tmpMacroBin);
		}
		 
		if (numMicroBins == 0)
			throw new Exception("Dimension with no candidate bin!");
			 
		// mine the optimal macro bins by means of dynamic programming
		if (numProcessedDims > 0)
		{
			if (Constants.METHOD == Constants.M_GD_CAIM)
				ret = greedyMDLCAIMBinning(dataMatrix, curDim, existingMacroBins, numRows, numMicroBins, initialMicroBins, discreteData, cells, cellPointIDs);
			//else
			//	ret = greedyMDLBinning(dataMatrix, curDim, procDims, existingMacroBins, candidateMacroBins, numRows, numMicroBins, initialMicroBins, discreteData, cells, cellPointIDs);
		}
		else
			ret = candidateMacroBins;
		
		System.out.println("number of output bins = " + ret.size());
		System.out.println("------------");
		return ret;
	}
	
	/*public static ArrayList<MacroBin> greedyMDLBinning(DataMatrix dataMatrix, int curDim, int[] procDims, ArrayList<MacroBin>[] existingMacroBins, ArrayList<MacroBin> candidateMacroBins, int N, int T, ArrayList<MicroBin> initialBins, int[][] discreteData, ArrayList<ArrayList<Integer>> cells, ArrayList<ArrayList<Integer>> cellPointIDs) throws Exception
	{
		ArrayList<Integer> dims = new ArrayList<Integer>();
		for (int i = 0; i < procDims.length; i++)
			dims.add(new Integer(i));
		double t0 = 2 * Constants.MAX_VAL * dims.size() * (2 + 1 / Math.log(2));
		
		int numOfBins = candidateMacroBins.size();
		ArrayList<MacroBin> macroBins = new ArrayList<MacroBin>();
		for (int i = 0; i < numOfBins; i++)
			macroBins.add(candidateMacroBins.get(i));
		
		MacroBin bin1 = null;
		MacroBin bin2 = null;
		MacroBin newBin = null;
		//boolean found;
		double diffAB;
		double diffBA;
		double tmpCost;
		int pos;
		double GlobalMinCost = Double.MAX_VALUE;
		ArrayList<MacroBin> GlobalMacroBins = new ArrayList<MacroBin>();
		SortedObjectComparator c = new SortedObjectComparator();
		SortedObject[] pairwiseMCKLs = null;
		while (numOfBins > 1)
		{
			pairwiseMCKLs = new SortedObject[numOfBins - 1];
			for (int i = 0; i < numOfBins - 1; i++)
			{
				bin1 = macroBins.get(i);
				bin2 = macroBins.get(i + 1);
				diffAB = MCKLComputer.computeDistance(bin1.dataMatrix, null, bin2.dataMatrix, null, dims, Constants.ALPHA, Constants.RAN_INDEX_BLOCK, bin1.means, bin1.devs);
				if (diffAB < 0 || diffAB >= t0)
				{
					System.out.println(diffAB + " --- " + t0);
					throw new Exception("Invalid MCKL calculation");
				}
				
				diffBA = MCKLComputer.computeDistance(bin2.dataMatrix, null, bin1.dataMatrix, null, dims, Constants.ALPHA, Constants.RAN_INDEX_BLOCK, bin2.means, bin2.devs);
				if (diffBA < 0 || diffBA >= t0)
				{
					System.out.println(diffBA + " --- " + t0);
					throw new Exception("Invalid MCKL calculation");
				}
				
				pairwiseMCKLs[i] = new SortedObject(i, diffAB + diffBA);
			}
			
			Arrays.sort(pairwiseMCKLs, c);
			//found = false;
			for (int i = 0; i < numOfBins - 1; i++)
			{
				pos = pairwiseMCKLs[i].index;
				bin1 = macroBins.get(pos);
				bin2 = macroBins.get(pos + 1);
				newBin = mergeMacroBins(dataMatrix, bin1, bin2, initialBins, procDims);
				macroBins.remove(pos);
				macroBins.remove(pos);
				macroBins.add(pos, newBin);
				numOfBins--;
				tmpCost = computeEncodingCost(dataMatrix, curDim, macroBins, N, T, initialBins, discreteData, cellPointIDs);
				if (tmpCost < GlobalMinCost)
				{
					GlobalMinCost = tmpCost;
					GlobalMacroBins.clear();
					for (int j = 0; j < numOfBins; j++)
						GlobalMacroBins.add(cloneBin(macroBins.get(j)));
					
					//found = true;
				}
				macroBins.remove(pos);
				macroBins.add(pos, bin1);
				macroBins.add(pos + 1, bin2);
				numOfBins++;
			}
			
			//if (found == false)
			//	break;
		} // end while
		
		return GlobalMacroBins;
	}*/
	
	public static ArrayList<MacroBin> greedyMDLCAIMBinning(DataMatrix dataMatrix, int curDim, ArrayList<MacroBin>[] existingMacroBins, int N, int T, ArrayList<MicroBin> initialBins, int[][] discreteData, ArrayList<ArrayList<Integer>> cells, ArrayList<ArrayList<Integer>> cellPointIDs) throws Exception
	{	
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
				
				tmpCost = computeEncodingCost(dataMatrix, curDim, tmpBins, N, T, initialBins, discreteData, cellPointIDs);
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
	
	// clone a macro bin
	public static MacroBin cloneBin(MacroBin a)
	{
		// init the new macro bin
		MacroBin ret = new MacroBin(a.lowerBound, a.upperBound);
		
		// get the total number of data points
		ret.numPoints = a.numPoints;
		ret.dataMatrix = new DataMatrix();
		ret.dataMatrix.data = new ArrayList<DataPoint>();
		int numCols = a.dataMatrix.cols;
		ret.means = new double[numCols];
		ret.devs = new double[numCols];
		int pointID;
		DataPoint tmpPoint = null;
		DataPoint newPoint = null;
		for (int i = 0; i < a.numPoints; i++)
		{
			pointID = a.pointIDs.get(i);
			tmpPoint = a.dataMatrix.data.get(i);
			newPoint = new DataPoint(numCols, 0);
			for (int j = 0; j < numCols; j++)
				newPoint.measures[j] = tmpPoint.measures[j];
			ret.dataMatrix.data.add(newPoint);
			ret.pointIDs.add(new Integer(pointID));
		}
		for (int i = 0; i < numCols; i++)
		{
			ret.means[i] = a.means[i];
			ret.devs[i] = a.devs[i];
		}
		
		// get the IDs of micro bins
		int numMicroBinsA = a.microBinIDs.size();
		for (int i = 0; i < numMicroBinsA; i++)
			ret.microBinIDs.add(new Integer(a.microBinIDs.get(i)));
		
		return ret;
	}
	
	// merge two macro bins
	public static MacroBin mergeMacroBins(DataMatrix dataMatrix, MacroBin a, MacroBin b, ArrayList<MicroBin> initialBins, int[] procDims) throws Exception
	{	
		// init the new macro bin
		if (a.lowerBound >= b.upperBound)
			throw new Exception("Invalid bin merge");
		
		MacroBin ret = new MacroBin(a.lowerBound, b.upperBound);
		
		// get the total number of data points
		int numCols = procDims.length;
		ret.numPoints = a.numPoints + b.numPoints;
		ret.dataMatrix = new DataMatrix();
		ret.dataMatrix.data = new ArrayList<DataPoint>();
		ret.means = new double[numCols];
		ret.devs = new double[numCols];
		int pointID;
		DataPoint tmpPoint = null;
		DataPoint newPoint = null;
		for (int i = 0; i < a.numPoints; i++)
		{
			pointID = a.pointIDs.get(i);
			tmpPoint = dataMatrix.data.get(pointID);
			newPoint = new DataPoint(numCols, 0);
			for (int j = 0; j < numCols; j++)
			{
				newPoint.measures[j] = tmpPoint.measures[procDims[j]];
				ret.means[j] += newPoint.measures[j];
			}
			ret.dataMatrix.data.add(newPoint);
			ret.pointIDs.add(new Integer(pointID));
		}
		for (int i = 0; i < b.numPoints; i++)
		{
			pointID = b.pointIDs.get(i);
			tmpPoint = dataMatrix.data.get(pointID);
			newPoint = new DataPoint(numCols, 0);
			for (int j = 0; j < numCols; j++)
			{
				newPoint.measures[j] = tmpPoint.measures[procDims[j]];
				ret.means[j] += newPoint.measures[j];
			}
			ret.dataMatrix.data.add(newPoint);
			ret.pointIDs.add(new Integer(pointID));
		}
		for (int i = 0; i < numCols; i++)
		{
			ret.means[i] = ret.means[i] / ret.numPoints;
			if (ret.numPoints == 1)
				ret.devs[i] = 0;
			else
			{
				for (int j = 0; j < ret.numPoints; j++)
				{
					tmpPoint = ret.dataMatrix.data.get(j);
					ret.devs[i] += (tmpPoint.measures[i] - ret.means[i]) * (tmpPoint.measures[i] - ret.means[i]) / (ret.numPoints - 1);
				}
				ret.devs[i] = Math.sqrt(ret.devs[i]);
			}
		} // end for
		
		// get the IDs of micro bins
		int numMicroBinsA = a.microBinIDs.size();
		int numMicroBinsB = b.microBinIDs.size();
		for (int i = 0; i < numMicroBinsA; i++)
			ret.microBinIDs.add(new Integer(a.microBinIDs.get(i).intValue()));
		for (int i = 0; i < numMicroBinsB; i++)
			ret.microBinIDs.add(new Integer(b.microBinIDs.get(i)));
		
		return ret;
	}
	
	// T = total number of micro bins contained in all macro bins
	// N = total number of data points contained in DS
	public static double computeEncodingCost(DataMatrix dataMatrix, int curDim, ArrayList<MacroBin> bins, int N, int T, ArrayList<MicroBin> initialBins, int[][] discreteData, ArrayList<ArrayList<Integer>> cellPointIDs) throws Exception
	{
		// the total coding cost
		double ret = 0;
		double logBase = Math.log(Constants.LOG_BASE);
		
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
		int numCells;
		ArrayList<Integer> curCellPointIDs = null;
		int binNumPoints;
		int pointID;
		int cellPointID;
		int countJoin;
		int count;
		int totalPoints;
		for (int i = 0; i < numBins; i++)
		{
			// get the macro bin
			tmpBin = bins.get(i);
			binNumPoints = tmpBin.numPoints;
			
			// increase LBinSize
			LBinSize += LInteger(tmpBin.numPoints);
			
			// increase LBinEntropy
			if (Constants.USE_CE == true)
				LBinEntropy += computeBinInformationContentWithCE(tmpBin, initialBins);
			else
				LBinEntropy += computeBinInformationContentWithoutCE(tmpBin, initialBins);
			
			// increase LBinCorrelation
			numCells = cellPointIDs.size();
			if (numCells > 0)
			{
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
						LBinCorrelation -= countJoin * Math.log(countJoin * 1.0 / count) / logBase;
				} // end for
				
				if (totalPoints != N)
				{
					System.out.println(totalPoints + " --- " + N);
					throw new Exception("total points not match");
				}
			}
			else
				LBinCorrelation -= binNumPoints * Math.log(binNumPoints * 1.0 / N) / logBase;
		} // end for
		
		// get the total coding cost
		//System.out.println(LB + " --- " + LBinSize + " --- " + LBinEntropy + " --- " + LBinCorrelation);
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
