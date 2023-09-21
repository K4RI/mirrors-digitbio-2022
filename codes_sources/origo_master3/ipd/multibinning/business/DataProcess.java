package multibinning.business;

import java.util.*;
import multibinning.data.*;

public class DataProcess 
{
	public static void calculateIndices(DataMatrix dataMatrix, IndexMatrix indexMatrix)
	{
		int numRows = dataMatrix.rows;
		int numCols = dataMatrix.cols;
		SortedObjectComparator c = new SortedObjectComparator();
		
		SortedObject[] tmp = new SortedObject[numRows];
		for (int i = 0; i < numRows; i++)
			tmp[i] = new SortedObject(0, 0);
		
		double[] vals = new double[numRows];
		for (int j = 0; j < numCols; j++)
		{
			for (int i = 0; i < numRows; i++)
			{
				tmp[i].index = i;
				tmp[i].value = dataMatrix.data.get(i).measures[j];
			}
			
			Arrays.sort(tmp, c);
			for (int i = 0; i < numRows; i++)
			{
				vals[i] = tmp[i].value;
				indexMatrix.data[i][j] = tmp[i].index;
				tmp[i].index = 0;
				tmp[i].value = 0;
			}
		}
	}
	
	public static int[] getIndexBlock(IndexMatrix indexMatrix, int colID, int startPosition, int size)
	{
		int[] ret = new int[size];
		int endPosition = startPosition + size - 1;
		for (int i = startPosition; i <= endPosition; i++)
			ret[i - startPosition] = indexMatrix.data[i][colID];
			
		return ret;
	}
	
	public static int[] getIndexBlock(int[] indexBlock, int startPosition, int size)
	{
		int[] ret = new int[size];
		int endPosition = startPosition + size - 1;
		for (int i = startPosition; i <= endPosition; i++)
			ret[i - startPosition] = indexBlock[i];
			
		return ret;
	}
	
	public static int generateRandomNumber(Random generator, int start, int end, boolean useFraction)
	{
		int retVal;
		if (useFraction == false)
			retVal = start + generator.nextInt(end - start + 1);
		else
		{
			long range = (long)end - (long)start + 1;
		    long fraction = (long)(range * generator.nextDouble());
		    retVal =  (int)(fraction + start);    
		}
		
		return retVal;
	}
}
