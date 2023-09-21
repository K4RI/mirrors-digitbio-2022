package multibinning.business;

import java.util.*;

import multibinning.data.*;

public class CAIM 
{
	@SuppressWarnings("unchecked")
	public static ArrayList<MacroBin>[] discretizeData(DataMatrix dataMatrix) throws Exception
	{
		int numCols = dataMatrix.cols;
		ArrayList<MacroBin>[] ret = new ArrayList[numCols];
		
		for (int dim = 0; dim < numCols; dim++)
			ret[dim] = caimBinning(dim, dataMatrix, Constants.CLASS_LABELS);
		
		return ret;
	}
	
	public static ArrayList<MacroBin> caimBinning(int curDim, DataMatrix dataMatrix, ArrayList<String> classLabels) throws Exception
	{
		ArrayList<MacroBin> ret = new ArrayList<MacroBin>();
		
		int numRows = dataMatrix.rows;
		//int numCols = dataMatrix.cols;
		int numClasses = classLabels.size();
		
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
		ArrayList<Double> distinctValues = new ArrayList<Double>();
		double curVal = tmp[0].value;
		distinctValues.add(new Double(curVal));
		for (int i = 1; i < numRows; i++)
		{
			if (tmp[i].value != curVal)
			{
				curVal = tmp[i].value;
				distinctValues.add(new Double(curVal));
			}
		}
		
		int numRemainingValues = distinctValues.size();
		double minValue = distinctValues.get(0);
		double maxValue = distinctValues.get(numRemainingValues - 1);
		
		ArrayList<Double> remainingCPs = new ArrayList<Double>();
		for (int i = 1; i < numRemainingValues - 1; i++)
			remainingCPs.add(new Double(distinctValues.get(i)));
		numRemainingValues -= 2;
			
		double GlobalCAIM = 0;
		double tmpCAIM;
		
		ArrayList<Double> cutPoints = new ArrayList<Double>();
		cutPoints.add(new Double(maxValue));
		int numCPs = 1;
		ArrayList<MicroBin> GlobalBins = new ArrayList<MicroBin>();
		GlobalBins.add(new MicroBin(minValue - 1, maxValue));
		
		double candidateCP;
		int pos;
		ArrayList<MicroBin> tmpBins = new ArrayList<MicroBin>();
		MicroBin tmpBin = null;
		double tmpLowerBound;
		double tmpMaxCAIM;
		int tmpMaxIndex;
		int tmpMaxPos;
		int tmpNumBins;
		while (numRemainingValues > 0)
		{
			tmpMaxCAIM = -Double.MAX_VALUE;
			tmpMaxIndex = -1;
			tmpMaxPos = -1;
			for (int i = 0; i < numRemainingValues; i++)
			{
				candidateCP = remainingCPs.get(i).doubleValue();
				pos = numCPs;
				for (int j = numCPs - 1; j >= 0; j--)
				{
					//System.out.println(j + " -- " + cutPoints.size() + " --- " + cutPoints.get(j).doubleValue());
					if (candidateCP >= cutPoints.get(j).doubleValue())
						break;
					else
						pos--;
				}
				
				cutPoints.add(pos, new Double(candidateCP));
				numCPs++;
				tmpBins.clear();
				tmpLowerBound = minValue - 1;
				for (int j = 0; j < numCPs; j++)
				{
					tmpBin = new MicroBin(tmpLowerBound, cutPoints.get(j).doubleValue());
					for (int r = 0; r < numRows; r++)
					{
						curPoint = dataMatrix.data.get(r);
						if (curPoint.measures[curDim] > tmpBin.lowerBound && curPoint.measures[curDim] <= tmpBin.upperBound)
							tmpBin.pointIDs.add(new Integer(r));
					}
					tmpBins.add(tmpBin);
					tmpLowerBound = cutPoints.get(j).doubleValue();
				}
				
				tmpCAIM = computeCAIM(tmpBins, dataMatrix, classLabels);
				if (tmpCAIM > tmpMaxCAIM)
				{
					tmpMaxCAIM = tmpCAIM;
					tmpMaxIndex = i;
					tmpMaxPos = pos;
				}
				
				cutPoints.remove(pos);
				numCPs--;
			}
			cutPoints.add(tmpMaxPos, remainingCPs.get(tmpMaxIndex).doubleValue());
			numCPs++;
			remainingCPs.remove(tmpMaxIndex);
			//System.out.println("------");
			
			tmpNumBins = tmpBins.size();
			if (tmpMaxCAIM > GlobalCAIM || tmpNumBins < numClasses)
			{
				GlobalCAIM = tmpMaxCAIM;
				GlobalBins.clear();
				for (int i = 0; i < tmpNumBins; i++)
					GlobalBins.add(tmpBins.get(i));
			}
			else
				break;
			
			numRemainingValues--;
		}
		
		// for each micro bin, create a macro bin containing it
		int numMicroBins = GlobalBins.size();
		MacroBin tmpMacroBin = null;
		MicroBin tmpMicroBin = null;
		for (int i = 0; i < numMicroBins; i++)
		{
			tmpMicroBin = GlobalBins.get(i);
			tmpMacroBin = new MacroBin(tmpMicroBin.lowerBound, tmpMicroBin.upperBound);
			ret.add(tmpMacroBin);
		}
		
		return ret;
	}
	
	public static double computeCAIM(ArrayList<MicroBin> bins, DataMatrix dataMatrix, ArrayList<String> classLabels) throws Exception
	{
		int numBins = bins.size();
		MicroBin tmpBin = null;
		int tmpNumPoints;
		DataPoint tmpPoint = null;
		int numClasses = classLabels.size();
		int[] classSupports = new int[numClasses];
		boolean exists;
		int maxClassSupport;
		double caim = 0;
		for (int i = 0; i < numBins; i++)
		{
			for (int j = 0; j < numClasses; j++)
				classSupports[j] = 0;
				
			tmpBin = bins.get(i);
			tmpNumPoints = tmpBin.pointIDs.size();
			for (int j = 0; j < tmpNumPoints; j++)
			{
				tmpPoint = dataMatrix.data.get(tmpBin.pointIDs.get(j).intValue());
				exists = false;
				for (int k = 0; k < numClasses; k++)
				{
					if (tmpPoint.classID.equals(classLabels.get(k)))
					{
						classSupports[k]++;
						exists = true;
						break;
					}
				}
				
				if (exists == false)
					throw new Exception("CAIM: Class not found!");
			}
			
			maxClassSupport = -1;
			for (int j = 0; j < numClasses; j++)
				if (classSupports[j] > maxClassSupport)
					maxClassSupport = classSupports[j];
			
			caim += maxClassSupport * maxClassSupport / (1.0 * tmpNumPoints);
		}
		
		caim = caim / numBins;
		
		return caim;
	}
}
