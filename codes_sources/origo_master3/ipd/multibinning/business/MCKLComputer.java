package multibinning.business;

import java.util.*;
import multibinning.data.*;

public class MCKLComputer 
{
	public static double computeDistanceNewNew(DataMatrix dataMatrixA, DataMatrix dataMatrixB, double v) throws Exception
	{
		int numDims = dataMatrixA.cols;
		int k = dataMatrixA.rows;
		int l = dataMatrixB.rows;
		
		// pick the dimension with largest difference between two sets of points
		double diff = 0;
		
		double sum1 = 0;
		double prod;
		for (int j1 = 0; j1 < k; j1++)
			for (int j2 = 0; j2 < k; j2++)
			{
				prod = 1;
				for (int i = 0; i < numDims; i++)
					prod = prod * (v - Math.max(dataMatrixA.data.get(j1).measures[i], dataMatrixA.data.get(j2).measures[i]));
				sum1 += prod;
			}
		
		double sum2 = 0;
		for (int j1 = 0; j1 < k; j1++)
			for (int j2 = 0; j2 < l; j2++)
			{
				prod = 1;
				for (int i = 0; i < numDims; i++)
					prod = prod * (v - Math.max(dataMatrixA.data.get(j1).measures[i], dataMatrixB.data.get(j2).measures[i]));
				sum2 += prod;
			}
		
		double sum3 = 0;
		for (int j1 = 0; j1 < l; j1++)
			for (int j2 = 0; j2 < l; j2++)
			{
				prod = 1;
				for (int i = 0; i < numDims; i++)
					prod = prod * (v - Math.max(dataMatrixB.data.get(j1).measures[i], dataMatrixB.data.get(j2).measures[i]));
				sum3 += prod;
			}
		
		diff = sum1 / (k * 1.0 * k) - 2 * sum2 / (k * 1.0 * l) + sum3 / (l * 1.0 * l); 
		
		return diff;
	}
	
	public static double computeDistanceNew(DataMatrix dataMatrixA, IndexMatrix indexMatrixA, DataMatrix dataMatrixB, IndexMatrix indexMatrixB, ArrayList<Integer> dims, double[] means, double[] devs) throws Exception
	{
		int numDims = dims.size();
		int numRowsA = dataMatrixA.rows;
		int numRowsB = dataMatrixB.rows;
		
		// pick the dimension with largest difference between two sets of points
		double maxDimDiff = -Double.MAX_VALUE;
		int maxDim = -1;
		int tmpDim;
		double tmpDiff;
		double[] a = new double[numRowsA];
		double[] b = new double[numRowsB];
		tmpDim = 0;
		for (int j = 0; j < numRowsA; j++)
			a[j] = dataMatrixA.data.get(indexMatrixA.data[j][tmpDim]).measures[tmpDim];
		
		for (int j = 0; j < numRowsB; j++)
			b[j] = dataMatrixB.data.get(indexMatrixB.data[j][tmpDim]).measures[tmpDim];
		
		tmpDiff = MCKL(a, true, b, true);
		maxDim = tmpDim;
		maxDimDiff = tmpDiff;
		
		// get the remaining dimensions
		//System.out.println(maxDim);
		ArrayList<Integer> processedDims = new ArrayList<Integer>();
		processedDims.add(new Integer(maxDim));
		ArrayList<Integer> remainingDims = new ArrayList<Integer>();
		int curDim;
		for (int i = 0; i < numDims; i++)
		{
			curDim = dims.get(i).intValue();
			if (curDim != maxDim)
				remainingDims.add(new Integer(i));
		}
		
		// compute the difference in the remaining dimensions
		int curLength = 1;
		double diff = maxDimDiff;
		while (curLength < numDims)
		{
			tmpDiff = computeDiffSpecialNew(dataMatrixA, dataMatrixB, processedDims, remainingDims, Constants.NUM_CENTROIDS, means, devs);
			diff += tmpDiff;
			curLength++;
		}
		
		return diff;
	}
	
	@SuppressWarnings({ "unchecked" })
	public static double computeDiffSpecialNew(DataMatrix dataMatrixA, DataMatrix dataMatrixB, ArrayList<Integer> processedDims, ArrayList<Integer> remainingDims, int NumSeeds, double[] means, double[] devs) throws Exception
	{
		int numProcessedDims = processedDims.size();
		int numRowsA = dataMatrixA.rows;
		int numRowsB = dataMatrixB.rows;
		
		if (numRowsA > 0)
		{
			ArrayList<Integer> dimsForClustering = new ArrayList<Integer>();
			ArrayList<Integer>[] clustersA = null;
			ArrayList<Integer>[] clustersB = null;
			ArrayList<Integer> curClusterA = null;
			ArrayList<Integer> curClusterB = null;
			int curClusterSizeA;
			int curClusterSizeB;
			for (int i = 0; i < numProcessedDims; i++)
				dimsForClustering.add(new Integer(processedDims.get(i).intValue()));
			
			// perform data clustering
			ArrayList<DataPoint> tmpSeeds = new ArrayList<DataPoint>();
			ArrayList<Double> tmpClusterRadius = new ArrayList<Double>();
			//System.out.println("start clustering");
			ArrayList<ArrayList<Integer>> tmpClusters = incrementalClusteringNewSeed(dataMatrixA, dimsForClustering, NumSeeds, means, devs, tmpClusterRadius, tmpSeeds); 
			//System.out.println("end clustering");
			int numClusters = tmpClusters.size();
			DataPoint[] seeds = new DataPoint[numClusters];
			double[] clusterRadius = new double[numClusters];
			clustersA = new ArrayList[numClusters];
			clustersB = new ArrayList[numClusters];
			for (int i = 0; i < numClusters; i++)
			{
				clustersA[i] = tmpClusters.get(i);
				seeds[i] = tmpSeeds.get(i);
				clusterRadius[i] = tmpClusterRadius.get(i).doubleValue();
				clustersB[i] = new ArrayList<Integer>();
			}
			
			// assign the data points of b to nearest clusters without breaking clusters' radius
			DataPoint curPoint;
			DataPoint curSeed;
			double min;
			int minSeedId;
			double tmpDist;
			for (int i = 0; i < numRowsB; i++)
			{
				min = Double.MAX_VALUE;
				minSeedId = -1;
				curPoint = dataMatrixB.data.get(i);
				for (int j = 0; j < numClusters; j++)
				{
					curSeed = seeds[j];
					tmpDist = DataPoint.distanceLNorm(Constants.LP_NORM, dimsForClustering, curPoint, curSeed);
					if (tmpDist < min)
					{
						min = tmpDist;
						minSeedId = j;
					} // end if
				}
				
				if (min <= clusterRadius[minSeedId])
					clustersB[minSeedId].add(new Integer(i));
			}
			
			// for each remaining dimension, compute its gain in overall difference and pick the one with largest gain
			double[] a = null;
			double[] b = null;
			int curDim;
			double diff;
			double maxDiff = -Double.MAX_VALUE;
			int maxPos = -1;
			//System.out.println(numRemainingDims);
			curDim = remainingDims.get(0).intValue();
			diff = 0;
			for (int j = 0; j < numClusters; j++)
			{
				// pick two clusters
				curClusterA = clustersA[j];
				curClusterSizeA = curClusterA.size();
				curClusterB = clustersB[j];
				curClusterSizeB = curClusterB.size();
				//System.out.println(curClusterSizeA + " --- " + curClusterSizeB);
				
				// get the data points of two clusters
				a = new double[curClusterSizeA];
				b = new double[curClusterSizeB];
				for (int k = 0; k < curClusterSizeA; k++)
					a[k] = dataMatrixA.data.get(curClusterA.get(k).intValue()).measures[curDim];
				for (int k = 0; k < curClusterSizeB; k++)
					b[k] = dataMatrixB.data.get(curClusterB.get(k).intValue()).measures[curDim];
			
				// compute their difference
				diff += curClusterSizeA * MCKL(a, false, b, false) / numRowsA;
			} // end for
			maxDiff = diff;
			maxPos = 0;
			
			// remove the picked dimension
			processedDims.add(new Integer(remainingDims.get(maxPos).intValue()));
			remainingDims.remove(maxPos);
			return maxDiff;
		}
		else if (numRowsB > 0)
		{
			// for each remaining dimension, compute its gain in overall difference and pick the one with largest gain
			double[] a = null;
			double[] b = null;
			int curDim;
			double diff;
			double maxDiff = -Double.MAX_VALUE;
			int maxPos = -1;
			//System.out.println(numRemainingDims);
			curDim = remainingDims.get(0).intValue();
			a = new double[0];
			b = new double[numRowsB];
			for (int j = 0; j < numRowsB; j++)
				b[j] = dataMatrixB.data.get(j).measures[curDim];
		
			// compute their difference
			diff = MCKL(a, false, b, false);
			maxDiff = diff;
			maxPos = 0;
			
			// remove the picked dimension
			processedDims.add(new Integer(remainingDims.get(maxPos).intValue()));
			remainingDims.remove(maxPos);
			return maxDiff;
		}
		else
			return 0;
	}
	
	public static double computeDiffSpecial(DataMatrix dataMatrixA, IndexMatrix indexMatrixA, DataMatrix dataMatrixB, IndexMatrix indexMatrixB, ArrayList<Integer> processedDims, ArrayList<Integer> remainingDims, double alpha, Random ranIndexBlock) throws Exception
	{
		int numProcessedDims = processedDims.size();
		int numRemainingDims = remainingDims.size();
		int numRowsA = dataMatrixA.rows;
		int numRowsB = dataMatrixB.rows;
		
		int[] indexBlockA = null;
		int currentCols;
		int startPosition;
		ArrayList<Integer> keepPointsA = new ArrayList<Integer>();
		ArrayList<Integer> keepPointsB = new ArrayList<Integer>();
		int blockSize = (int)Math.ceil(dataMatrixA.rows * Math.pow(alpha, 1.0 / numProcessedDims));
		int maxBlockPosition = dataMatrixA.rows - blockSize;
		double[] minIndexBlockVal = new double[numProcessedDims];
		double[] maxIndexBlockVal = new double[numProcessedDims];
		double[] a = null;
		double[] b = null;
		int countA;
		int countB;
		DataPoint tmpPoint = null;
		double tmpDiff;
		
		// initialize all difference per dimension to be 0
		double[] diffs = new double[numRemainingDims];
		for (int i = 0; i < numRemainingDims; i++)
			diffs[i] = 0;
		
		// perform sampling NUM_SUBSAMPLING times
		for (int h = 0; h < Constants.NUM_SUBSAMPLING; h++)
		{
			// clear two sets of points
			keepPointsA.clear();
			keepPointsB.clear();
			
			// get index blocks for all dimensions
			for (int i = 0; i < numProcessedDims; i++)
			{
				currentCols = processedDims.get(i).intValue();
				
				startPosition = DataProcess.generateRandomNumber(ranIndexBlock, 0, maxBlockPosition, false);
				indexBlockA = DataProcess.getIndexBlock(indexMatrixA, currentCols, startPosition, blockSize);
				
				minIndexBlockVal[i] = dataMatrixA.data.get(indexBlockA[0]).measures[currentCols];
				maxIndexBlockVal[i] = dataMatrixA.data.get(indexBlockA[indexBlockA.length - 1]).measures[currentCols];
			}
			
			// extract points of a
			for (int i = 0; i < numRowsA; i++)
			{
				tmpPoint = dataMatrixA.data.get(i);
				for (int j = 0; j < numProcessedDims; j++)
				{
					currentCols = processedDims.get(j).intValue();
					if (tmpPoint.measures[currentCols] >= minIndexBlockVal[j] && tmpPoint.measures[currentCols] <= maxIndexBlockVal[j])
						keepPointsA.add(new Integer(i));
				}
			}
			
			// extract points of b
			for (int i = 0; i < numRowsB; i++)
			{
				tmpPoint = dataMatrixB.data.get(i);
				for (int j = 0; j < numProcessedDims; j++)
				{
					currentCols = processedDims.get(j).intValue();
					if (tmpPoint.measures[currentCols] >= minIndexBlockVal[j] && tmpPoint.measures[currentCols] <= maxIndexBlockVal[j])
						keepPointsB.add(new Integer(i));
				}
			}
			
			// compute the difference between two sets of points
			countA = keepPointsA.size();
			countB = keepPointsB.size();
			a = new double[countA];
			b = new double[countB];
			for (int i = 0; i < numRemainingDims; i++)
			{
				currentCols = remainingDims.get(i).intValue();
				for (int j = 0; j < countA; j++)
					a[j] = dataMatrixA.data.get(keepPointsA.get(j).intValue()).measures[currentCols];
				for (int j = 0; j < countB; j++)
					b[j] = dataMatrixB.data.get(keepPointsB.get(j).intValue()).measures[currentCols];
				
				tmpDiff = MCKL(a, false, b, false);
				diffs[i] += 1.0 * tmpDiff / Constants.NUM_SUBSAMPLING;
			}
		}
		
		// pick the dimension that yields largest gain in overall difference
		double maxDiff = -Double.MAX_VALUE;
		int maxPos = -1;
		for (int i = 0; i < numRemainingDims; i++)
		{
			if (diffs[i] > maxDiff)
			{
				maxDiff = diffs[i];
				maxPos = i;
			}
		}
		
		// remove the picked dimension
		processedDims.add(new Integer(remainingDims.get(maxPos).intValue()));
		remainingDims.remove(maxPos);
		
		return maxDiff;
	}
	
	@SuppressWarnings({ "unchecked" })
	public static double computeDiffSpecial(DataMatrix dataMatrixA, DataMatrix dataMatrixB, ArrayList<Integer> processedDims, ArrayList<Integer> remainingDims, int NumSeeds, double[] means, double[] devs) throws Exception
	{
		int numProcessedDims = processedDims.size();
		int numRemainingDims = remainingDims.size();
		int numRowsA = dataMatrixA.rows;
		int numRowsB = dataMatrixB.rows;
		
		if (numRowsA > 0)
		{
			ArrayList<Integer> dimsForClustering = new ArrayList<Integer>();
			ArrayList<Integer>[] clustersA = null;
			ArrayList<Integer>[] clustersB = null;
			ArrayList<Integer> curClusterA = null;
			ArrayList<Integer> curClusterB = null;
			int curClusterSizeA;
			int curClusterSizeB;
			for (int i = 0; i < numProcessedDims; i++)
				dimsForClustering.add(new Integer(processedDims.get(i).intValue()));
			
			// perform data clustering
			ArrayList<DataPoint> tmpSeeds = new ArrayList<DataPoint>();
			ArrayList<Double> tmpClusterRadius = new ArrayList<Double>();
			//System.out.println("start clustering");
			ArrayList<ArrayList<Integer>> tmpClusters = incrementalClusteringNewSeed(dataMatrixA, dimsForClustering, NumSeeds, means, devs, tmpClusterRadius, tmpSeeds); 
			//System.out.println("end clustering");
			int numClusters = tmpClusters.size();
			DataPoint[] seeds = new DataPoint[numClusters];
			double[] clusterRadius = new double[numClusters];
			clustersA = new ArrayList[numClusters];
			clustersB = new ArrayList[numClusters];
			for (int i = 0; i < numClusters; i++)
			{
				clustersA[i] = tmpClusters.get(i);
				seeds[i] = tmpSeeds.get(i);
				clusterRadius[i] = tmpClusterRadius.get(i).doubleValue();
				clustersB[i] = new ArrayList<Integer>();
			}
			
			// assign the data points of b to nearest clusters without breaking clusters' radius
			DataPoint curPoint;
			DataPoint curSeed;
			double min;
			int minSeedId;
			double tmpDist;
			for (int i = 0; i < numRowsB; i++)
			{
				min = Double.MAX_VALUE;
				minSeedId = -1;
				curPoint = dataMatrixB.data.get(i);
				for (int j = 0; j < numClusters; j++)
				{
					curSeed = seeds[j];
					tmpDist = DataPoint.distanceLNorm(Constants.LP_NORM, dimsForClustering, curPoint, curSeed);
					if (tmpDist < min)
					{
						min = tmpDist;
						minSeedId = j;
					} // end if
				}
				
				if (min <= clusterRadius[minSeedId])
					clustersB[minSeedId].add(new Integer(i));
			}
			
			// for each remaining dimension, compute its gain in overall difference and pick the one with largest gain
			double[] a = null;
			double[] b = null;
			int curDim;
			double diff;
			double maxDiff = -Double.MAX_VALUE;
			int maxPos = -1;
			//System.out.println(numRemainingDims);
			for (int i = 0; i < numRemainingDims; i++)
			{
				curDim = remainingDims.get(i).intValue();
				diff = 0;
				for (int j = 0; j < numClusters; j++)
				{
					// pick two clusters
					curClusterA = clustersA[j];
					curClusterSizeA = curClusterA.size();
					curClusterB = clustersB[j];
					curClusterSizeB = curClusterB.size();
					//System.out.println(curClusterSizeA + " --- " + curClusterSizeB);
					
					// get the data points of two clusters
					a = new double[curClusterSizeA];
					b = new double[curClusterSizeB];
					for (int k = 0; k < curClusterSizeA; k++)
						a[k] = dataMatrixA.data.get(curClusterA.get(k).intValue()).measures[curDim];
					for (int k = 0; k < curClusterSizeB; k++)
						b[k] = dataMatrixB.data.get(curClusterB.get(k).intValue()).measures[curDim];
				
					// compute their difference
					diff += curClusterSizeA * MCKL(a, false, b, false) / numRowsA;
				} // end for
				
				// update the maximum gain
				if (diff > maxDiff)
				{
					maxDiff = diff;
					maxPos = i;
				}
				//System.out.println("dim " + i);
			}
			
			// remove the picked dimension
			processedDims.add(new Integer(remainingDims.get(maxPos).intValue()));
			remainingDims.remove(maxPos);
			return maxDiff;
		}
		else if (numRowsB > 0)
		{
			// for each remaining dimension, compute its gain in overall difference and pick the one with largest gain
			double[] a = null;
			double[] b = null;
			int curDim;
			double diff;
			double maxDiff = -Double.MAX_VALUE;
			int maxPos = -1;
			//System.out.println(numRemainingDims);
			for (int i = 0; i < numRemainingDims; i++)
			{
				curDim = remainingDims.get(i).intValue();
				a = new double[0];
				b = new double[numRowsB];
				for (int j = 0; j < numRowsB; j++)
					b[j] = dataMatrixB.data.get(j).measures[curDim];
			
				// compute their difference
				diff = MCKL(a, false, b, false);
				
				// update the maximum gain
				if (diff > maxDiff)
				{
					maxDiff = diff;
					maxPos = i;
				}
			}
			
			// remove the picked dimension
			processedDims.add(new Integer(remainingDims.get(maxPos).intValue()));
			remainingDims.remove(maxPos);
			return maxDiff;
		}
		else
			return 0;
	}
	
	@SuppressWarnings("unchecked")
	public static ArrayList<ArrayList<Integer>> incrementalClusteringNewSeed(DataMatrix dataMatrix, ArrayList<Integer> clusteringDims, int NumSeeds, double[] means, double[] devs, ArrayList<Double> clusterRadius, ArrayList<DataPoint> seeds)
	{
		/*int totalRows = Constants.NUM_ROWS;
		int totalColumns = Constants.NUM_MEASURE_COLS;
		boolean isLargeData = false;
		if (totalRows > 1200000)
			isLargeData = true;
		else if (totalColumns > 500)
			isLargeData = true;
		else if (totalRows > 5000)
		{
			if (totalColumns > 150)
				isLargeData = true;
		}*/
		//System.out.println("isLargeData " + isLargeData);
		
		ArrayList<ArrayList<Integer>> ret = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer>[] tmpRet = new ArrayList[NumSeeds];
		DataPoint[] tmpSeeds = new DataPoint[NumSeeds];
		double[] tmpClusterRadius =  new double[NumSeeds];
		//for (int i = 0; i < NumSeeds; i++)
		//	tmpClusterRadius[i] = 0;
		for (int i = 0; i < NumSeeds; i++)
			tmpRet[i] = new ArrayList<Integer>();
		
		int rows = dataMatrix.rows;
		int numDimsClustering = clusteringDims.size();
		DataPoint curPoint = null;
		int curDim;
		for (int i = 0; i < NumSeeds; i++)
		{
			tmpSeeds[i] = new DataPoint(numDimsClustering, 0);
			for (int j = 0; j < numDimsClustering; j++)
			{
				curDim = clusteringDims.get(j).intValue();
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
		int clusterSize;
		double tmpDist;
		double[][] sums = new double[NumSeeds][numDimsClustering];
		for (int iteration = 0; iteration < 3; iteration++)
		{
			for (int i = 0; i < NumSeeds; i++)
			{
				for (int c = 0; c < numDimsClustering; c++)
					sums[i][c] = 0;
				tmpRet[i].clear();
			}
			
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
				for (int c = 0; c < numDimsClustering; c++)
					sums[minSeedId][c] += curPoint.measures[clusteringDims.get(c).intValue()];
				
				if (iteration == 0 && i % numMSteps == 0)
				{
					//System.out.println("perform frequent M steps");
					performMStep(dataMatrix, clusteringDims, tmpSeeds, tmpRet, sums);
				}
			} // end for
			
			performMStep(dataMatrix, clusteringDims, tmpSeeds, tmpRet, sums);
			
			// if the last iteration, then compute the clusters' radius
			if (iteration == 2)
			{
				for (int i = 0; i < NumSeeds; i++)
				{
					tmpClusterRadius[i] = -Double.MAX_VALUE;
					curSeed = tmpSeeds[i];
					clusterSize = tmpRet[i].size();
					for (int j = 0; j < clusterSize; j++)
					{
						curPoint = dataMatrix.data.get(tmpRet[i].get(j).intValue());
						tmpDist = DataPoint.distanceLNorm(Constants.LP_NORM, clusteringDims, curPoint, curSeed);
						if (tmpDist > tmpClusterRadius[i])
							tmpClusterRadius[i] = tmpDist;
					}
				}
			} // end if
			//System.out.println("iteration " + iteration);
		} // end for
		
		// return only non-empty clusters
		for (int i = 0; i < NumSeeds; i++)
		{
			if (tmpRet[i].size() > 0)
			{
				ret.add(tmpRet[i]);
				seeds.add(tmpSeeds[i]);
				clusterRadius.add(new Double(tmpClusterRadius[i]));
			}
		}
		
		return ret;
	}
	
	public static void performMStep(DataMatrix dataMatrix, ArrayList<Integer> clusteringDims, DataPoint[] seeds, ArrayList<Integer>[] members, double[][] sums)
	{
		int numSeeds = seeds.length;
		int numPoints;
		//DataPoint curPoint  = null;
		int numDimsClustering = clusteringDims.size();
		//int curDim;
		for (int i = 0; i < numSeeds; i++)
		{
			numPoints = members[i].size();
			if (numPoints > 0)
			{
				for (int k = 0; k < numDimsClustering; k++)
					seeds[i].measures[k] = sums[i][k] / numPoints;
			}
			/*for (int j = 0; j < numPoints; j++)
			{
				curPoint = dataMatrix.data.get(members[i].get(j).intValue());
				for (int k = 0; k < numDimsClustering; k++)
				{
					curDim = clusteringDims.get(k).intValue();
					if (j == 0)
						seeds[i].measures[k] = curPoint.measures[curDim] / numPoints;
					else
						seeds[i].measures[k] += curPoint.measures[curDim] / numPoints;
				}
			} // end for
			*/
		} // end for
	}
	
	public static double MCKL(double[] a, boolean hasSortedA, double[] b, boolean hasSortedB) throws Exception
	{
		//System.out.println("----");
		double ret = 0;
		double cp1 = 0;
		if (a.length > 0)
			cp1 = computeOneSequenceWithLog(a, hasSortedA);
		ret += cp1;
		//System.out.println(ret);
		
		// if b is empty, there's no need to compute this value
		double cp2 = 0;
		double cp3 = 0;
		//System.out.println("b.length = " + b.length);
		if (b.length > 0)
		{
			if (a.length > 0)
				cp2 = computeTwistedSequencesWithLog(a, hasSortedA, b, hasSortedB);
			else
				cp2 = computeLogOneSequence(b, hasSortedB);
			ret -= cp2;
			//System.out.println(ret);
			
			cp3 = computeOneSequenceWithoutLog(b, hasSortedB) / Math.log(Constants.LOG_BASE);;
			ret += cp3;
			//System.out.println(ret);
		}
		
		double cp4 = 0;
		if (a.length > 0)
			cp4 = computeOneSequenceWithoutLog(a, hasSortedA) / Math.log(Constants.LOG_BASE);
		ret -= cp4;
		//System.out.println(ret);
		
		if (Math.abs(ret) < Constants.MAX_ERROR)
			ret = 0;
		double MAX_VALUE = 2 * Constants.MAX_VAL * (2 + 1 / Math.log(2));
		if (ret < 0 || ret >= MAX_VALUE)
		{
			System.out.println(cp1 + " --- " + cp2 + " --- " + cp3 + " --- " + cp4);
			System.out.println(MAX_VALUE + " --- " + ret);
			throw new Exception("Invalid calculation of MCKL");
		}
			
		return ret;
	}
	
	public static double computeOneSequenceWithoutLog(double[] a, boolean hasSorted)
	{
		double ret = 0;
		
		int numA = a.length;
		if (hasSorted == false)
			Arrays.sort(a);
		
		for (int i = 0; i < numA - 1; i++)
			if (a[i + 1] != a[i])
				ret += (a[i + 1] - a[i]) * (i + 1) / (1.0 * numA);
		
		ret += (Constants.MAX_VAL - a[numA - 1]);
		
		return ret;
	}
	
	public static double computeOneSequenceWithLog(double[] a, boolean hasSorted)
	{
		double ret = 0;
		
		double logBase = Math.log(Constants.LOG_BASE);
		int numA = a.length;
		if (hasSorted == false)
			Arrays.sort(a);
		
		for (int i = 0; i < numA - 1; i++)
			if (a[i + 1] != a[i])
				ret += (a[i + 1] - a[i]) * ((i + 1) / (1.0 * numA) + 1) * Math.log((i + 1) / (1.0 * numA) + 1) / logBase;
		
		ret += (Constants.MAX_VAL - a[numA - 1]) * 2;
		
		return ret;
	}
	
	public static double computeLogOneSequence(double[] a, boolean hasSorted)
	{
		double ret = 0;
		
		double logBase = Math.log(Constants.LOG_BASE);
		int numA = a.length;
		if (hasSorted == false)
			Arrays.sort(a);
		
		for (int i = 0; i < numA - 1; i++)
			if (a[i + 1] != a[i])
				ret += (a[i + 1] - a[i]) * Math.log((i + 1) / (1.0 * numA) + 1) / logBase;
		
		ret += (Constants.MAX_VAL - a[numA - 1]);
		
		return ret;
	}
	
	public static double computeTwistedSequencesWithLog(double[] a, boolean hasSortedA, double[] b, boolean hasSortedB) throws Exception
	{
		double ret = 0;
		
		double logBase = Math.log(Constants.LOG_BASE);
		
		// sort the a array
		int numA = a.length;
		if (hasSortedA == false)
			Arrays.sort(a);
		
		// sort the b array
		int numB = b.length;
		if (hasSortedB == false)
			Arrays.sort(b);
		
		// get the total number of items in both arrays
		int numTotal = numA + numB;
		
		// form the sorting of all values
		SortedObjectComparator c = new SortedObjectComparator();
		SortedObject[] tmp = new SortedObject[numTotal];
		
		// copy data to tmp
		for (int i = 0; i < numA; i++)
			tmp[i] = new SortedObject(i, a[i]);
		for (int i = 0; i < numB; i++)
			tmp[i + numA] = new SortedObject(i + numA, b[i]);
		
		// sort tmp
		Arrays.sort(tmp, c);
		
		// compute twisted result
		int curIndexA = 0;
		int curIndexB = 0;
		double part1;
		double part2;
		for (int i = 0; i < numTotal - 1; i++)
			if (tmp[i + 1].value != tmp[i].value)
			{
				// identify the types of two bounds
				if (tmp[i + 1].index < numA && tmp[i].index < numA)
				{
					curIndexA = tmp[i + 1].index;
				}
				else if (tmp[i + 1].index >= numA && tmp[i].index >= numA)
				{
					curIndexB = tmp[i + 1].index - numA;
				}
				else if (tmp[i + 1].index >= numA && tmp[i].index < numA)
				{
					curIndexA = tmp[i].index + 1;
					curIndexB = tmp[i + 1].index - numA;
				}
				else if (tmp[i + 1].index < numA && tmp[i].index >= numA)
				{
					curIndexA = tmp[i + 1].index;
					curIndexB = tmp[i].index + 1 - numA;
				}
				
				// get a robust computation
				if (numA != 0)
					part1 = curIndexA / (1.0 * numA) + 1;
				else
				{
					if (curIndexA != 0)
						throw new Exception("Division by zero");
					else
						part1 = 1;
				}
				
				// get a robust computation
				if (numB != 0)
					part2 = curIndexB / (1.0 * numB) + 1;
				else
				{
					if (curIndexB != 0)
						throw new Exception("Division by zero");
					else
						part2 = 1;
				}
				
				// increase the twisted value
				ret += (tmp[i + 1].value - tmp[i].value) * part1 * Math.log(part2) / logBase;
			}
		
		//if (numTotal > 0)
		ret += (Constants.MAX_VAL - tmp[numTotal - 1].value) * 2;
		
		return ret;
	}
}
