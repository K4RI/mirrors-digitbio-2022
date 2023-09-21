package multibinning.business;

import java.io.*;
import java.util.*;
import multibinning.data.*;

public class DataReader 
{
	public static void readData(String fileName, DataMatrix dataMatrix, ArrayList<String>[] singleContexts) throws Exception
	{
		try
		{
			Scanner scanner = new Scanner(new FileInputStream(fileName));
			
			String s = null;
			int rows = dataMatrix.rows;
			int cols = dataMatrix.cols;
			int catCols = dataMatrix.catCols;
			String[] vals = new String[catCols + cols + 1];
			DataPoint newPoint = null;
			int count = 0;
			
			double[] maxDims = new double[cols];
			double[] minDims = new double[cols];
			Constants.DATA_MEANS = new double[cols];
			Constants.DATA_DEVS = new double[cols];
			Constants.MAX_COLS = new double[cols];
			Constants.MIN_COLS = new double[cols];
			Constants.CRES  = new double[cols];
			for (int i = 0; i < cols; i++)
			{
				maxDims[i] = -Double.MAX_VALUE;
				minDims[i] = Double.MAX_VALUE;
				Constants.DATA_MEANS[i] = 0;
				Constants.DATA_DEVS[i] = 0;
			}
			
			Constants.CLASS_LABELS = new ArrayList<String>();
			int numClasses;
			boolean labelExists;
			String tmpLabel = null;
			boolean exists;
			int curNumSingleContexts;
			//scanner = scanner.useDelimiter(Constants.FIELD_DELIMITER + "|\\n");
			
			//while (scanner.hasNextLine())
			for (int r = 0; r < rows; r++)
			{
				s = scanner.nextLine();
				if (s.contains("?"))
					continue;
				
	        	vals = s.split(Constants.FIELD_DELIMITER);
	        	//for (int i = 0; i < vals.length; i++)
				//	vals[i] = scanner.next();
	        	
	        	if (vals.length != catCols + cols + 1)
	        	{
	        		System.out.println(s);
	        		throw new Exception("Invalid input!");
	        	}
	        	s = null;
	        	
		        newPoint = new DataPoint(cols, catCols);
		        newPoint.universalID = count;
		        
		        // update the list of class labels
		        newPoint.classID = new String(vals[catCols + cols]);
		        numClasses = Constants.CLASS_LABELS.size();
		        labelExists = false;
		        for (int j = 0; j < numClasses; j++)
		        {
		        	tmpLabel = Constants.CLASS_LABELS.get(j);
		        	if (tmpLabel.equals(newPoint.classID))
		        	{
		        		labelExists = true;
		        		break;
		        	}
		        }
		        if (!labelExists)
		        	Constants.CLASS_LABELS.add(new String(newPoint.classID));
		        
		        for (int i = 0; i < cols; i++)
		        {
		        	//vals[i] = vals[i].replace(',', '.');
		        	newPoint.measures[i] = Double.parseDouble(vals[i]);
		        	
		        	if (newPoint.measures[i] > maxDims[i])
		        		maxDims[i] = newPoint.measures[i];
		        	
		        	if (newPoint.measures[i] < minDims[i])
		        		minDims[i] = newPoint.measures[i];
		        }
		        
		        for (int i = 0; i < catCols; i++)
		        {
		        	//newPoint.catContexts[i] = new String(vals[i + cols]);
		        	newPoint.catContexts[i] = vals[i + cols];
		        	
		        	curNumSingleContexts = singleContexts[i].size();
		        	exists = false;
		        	for (int j = 0; j < curNumSingleContexts; j++)
		        		if (vals[i + cols].equals(singleContexts[i].get(j)))
		        		{
		        			exists = true;
		        			break;
		        		}
		        	if (exists == false)
		        		singleContexts[i].add(new String(vals[i + cols]));
		        } // end for
		        
		        dataMatrix.data.add(newPoint);
		        vals = null;
		        
		        //System.out.println("row " + r);
		        count++;
		        //if (count == rows)
		        //	break;
		    } // end while
			rows = count;
			Constants.NUM_ROWS = count;
			
			for (int i = 0; i < cols; i++)
			{
				Constants.MAX_COLS[i] = maxDims[i];
				Constants.MIN_COLS[i] = minDims[i];
			}
			
			for (int i = 0; i < rows; i++)
			{
				newPoint = dataMatrix.data.get(i);
				for (int j = 0; j < cols; j++)
				{
					if (maxDims[j] == minDims[j])
						newPoint.measures[j] = -Constants.MAX_VAL;
					else
						newPoint.measures[j] = 2 * Constants.MAX_VAL * (newPoint.measures[j] - minDims[j]) / (maxDims[j] - minDims[j]) - Constants.MAX_VAL;
					
					if (newPoint.measures[j] > Constants.MAX_VAL)
						newPoint.measures[j] = Constants.MAX_VAL;
					else if (newPoint.measures[j] < -Constants.MAX_VAL)
						newPoint.measures[j] = -Constants.MAX_VAL;
					
					Constants.DATA_MEANS[j] += newPoint.measures[j] / rows;
				}
			}
			
			double[] dimData = new double[rows];
			SortedObject[] so = new SortedObject[cols];
			for (int i = 0; i < cols; i++)
			{
				//dimData = new double[rows];
				for (int j = 0; j < rows; j++)
				{
					dimData[j] = dataMatrix.data.get(j).measures[i];
					Constants.DATA_DEVS[i] += (dimData[j] - Constants.DATA_MEANS[i]) * (dimData[j] - Constants.DATA_MEANS[i]) / (rows - 1);
				}
				Constants.DATA_DEVS[i] = Math.sqrt(Constants.DATA_DEVS[i]);
				
				Constants.CRES[i] = JBinning.computeCRE(dimData, false);
				so[i] = new SortedObject(i, Constants.CRES[i]);
			}
			
			// sort data according to CRES
			SortedObjectComparator c = new SortedObjectComparator();
			Arrays.sort(so, c);
			Constants.ODIMS = new int[cols];
			for (int i = 0; i < cols; i++)
				Constants.ODIMS[i] = so[i].index;
			
			scanner.close();
			
			dataMatrix.rows = dataMatrix.data.size();
		}
		catch (Exception ex)
		{
			throw ex;
		}
	}	
}