import java.util.*;
import java.io.*;
import multibinning.business.*;
import multibinning.data.*;

public class MainClass 
{
	public static long totalTimeSortingSubspaces = 0;
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args)
	{
		try
		{
			//Get the jvm heap size.
	        long heapSize = Runtime.getRuntime().totalMemory();
	        System.out.println("Heap Size = " + heapSize);
	        
			if (args.length == 6)
			{
				SupUniMDL.Discretize(args);
			}
			else
			{
				readInputString(args);
				DataMatrix dataMatrix = new DataMatrix(Constants.NUM_ROWS, Constants.NUM_MEASURE_COLS, Constants.NUM_CAT_CONTEXT_COLS);
				ArrayList<String>[] singleContexts = new ArrayList[Constants.NUM_CAT_CONTEXT_COLS];
				for (int i = 0; i < Constants.NUM_CAT_CONTEXT_COLS; i++)
					singleContexts[i] = new ArrayList<String>();
				
				System.out.println("start reading input data...");
				DataReader.readData(Constants.FILE_INPUT, dataMatrix, singleContexts);
				System.out.println("end reading input data...");
				
				System.out.println("start discretization...");
				BufferedWriter writerCP = new BufferedWriter(new FileWriter(new File(Constants.FILE_CP_OUTPUT)));
				BufferedWriter writerRuntime = new BufferedWriter(new FileWriter(new File(Constants.FILE_RUNTIME_OUTPUT)));
				
				long start = System.currentTimeMillis();
				ArrayList<MacroBin>[] finalBins = null;
				if (Constants.METHOD != Constants.CAIM && Constants.METHOD != Constants.J_MDL && Constants.METHOD != Constants.J_MDL_HM && Constants.METHOD != Constants.J_GD_HM && Constants.METHOD != Constants.M_GD_CAIM && Constants.METHOD != Constants.M_GD && Constants.METHOD != Constants.MMIC && Constants.METHOD != Constants.MVD && Constants.METHOD != Constants.PCABINNING)
					finalBins = BinMining.discretizeData(dataMatrix);
				else if (Constants.METHOD == Constants.J_MDL || Constants.METHOD == Constants.J_MDL_HM || Constants.METHOD == Constants.J_GD_HM)
					finalBins = JBinning.discretizeData(dataMatrix, Constants.CRES, Constants.DATA_MEANS, Constants.DATA_DEVS);
				else if (Constants.METHOD == Constants.M_GD_CAIM || Constants.METHOD == Constants.M_GD)
					finalBins = MBinning.discretizeData(dataMatrix, Constants.CRES, Constants.DATA_MEANS, Constants.DATA_DEVS);
				else if (Constants.METHOD == Constants.MMIC)
					finalBins = MMIC.discretizeData(dataMatrix, Constants.CRES, Constants.DATA_MEANS, Constants.DATA_DEVS);
				else if (Constants.METHOD == Constants.MVD)
					finalBins = MVD.discretizeData(dataMatrix);
				else if (Constants.METHOD == Constants.PCABINNING)
				{
					finalBins = PCABinning.discretizeData(dataMatrix, Constants.FILE_ORG_INPUT);
					
					dataMatrix = new DataMatrix(Constants.NUM_ROWS, Constants.NUM_MEASURE_COLS, Constants.NUM_CAT_CONTEXT_COLS);
					
					singleContexts = new ArrayList[Constants.NUM_CAT_CONTEXT_COLS];
					for (int i = 0; i < Constants.NUM_CAT_CONTEXT_COLS; i++)
						singleContexts[i] = new ArrayList<String>();
					
					DataReader.readData(Constants.FILE_ORG_NORMAL_INPUT, dataMatrix, singleContexts);
				}
				else
					finalBins = CAIM.discretizeData(dataMatrix);
				long end = System.currentTimeMillis();
				writerRuntime.write(Long.toString(end - start));
				writerRuntime.flush();
				writerRuntime.close();
				
				int numCols = Constants.NUM_MEASURE_COLS;
				int numCatCols = Constants.NUM_CAT_CONTEXT_COLS;
				int tmpNumMacroBins;
				MacroBin tmpMacroBin = null;
				MacroBin tmpPrevMacroBin = null;
				double maxCol;
				double minCol;
				ArrayList<MacroBin>[] outputBins = new ArrayList[numCols];
				for (int dim = 0; dim < numCols; dim++)
				{
					outputBins[dim] = new ArrayList<MacroBin>();
					maxCol = Constants.MAX_COLS[dim];
					minCol = Constants.MIN_COLS[dim];
					tmpNumMacroBins = finalBins[dim].size();
					tmpPrevMacroBin = finalBins[dim].get(0);
					outputBins[dim].add(tmpPrevMacroBin);
					for (int binIndex = 1; binIndex < tmpNumMacroBins; binIndex++)
					{
						tmpMacroBin = finalBins[dim].get(binIndex);
						if (tmpMacroBin.upperBound > tmpPrevMacroBin.upperBound)
							outputBins[dim].add(tmpMacroBin);
						tmpPrevMacroBin = tmpMacroBin;
					}
					tmpNumMacroBins = outputBins[dim].size();
					
					writerCP.write("dimension " + Integer.toString(dim) + " (" + tmpNumMacroBins + " bins)");
					writerCP.newLine();
					for (int binIndex = 0; binIndex < tmpNumMacroBins; binIndex++)
					{
						tmpMacroBin = outputBins[dim].get(binIndex);
						//writerCP.write(Double.toString(tmpMacroBin.lowerBound) + " --- " + Double.toString(tmpMacroBin.upperBound));
						writerCP.write(Double.toString(getOriginalValue(tmpMacroBin.upperBound, minCol, maxCol)));
						writerCP.newLine();
					}
					writerCP.write("-------------------------------------");
					writerCP.newLine();
				}
				System.out.println("end discretization...");
				
				writerCP.flush();
				writerCP.close();
				
				// write discretized data
				BufferedWriter writerData = new BufferedWriter(new FileWriter(new File(Constants.FILE_DATA_OUTPUT)));
				writerData.write("@relation DB");
				writerData.newLine();
				writerData.newLine();
				
				// write discretized continuous data
				int numRows = Constants.NUM_ROWS;
				double tmpLowerBound;
				double tmpUpperBound;
				String[][] outputData = new String[numRows][numCols];
				String[] pointLabels = new String[numRows];
				DataPoint tmpPoint = null;
				String prefix = null;
				String attName = null;
				int nominalIndex = 1;
				for (int dim = 0; dim < numCols; dim++)
				{
					prefix = Integer.toString(dim);
					writerData.write("@attribute dim" + prefix + " {");
					tmpNumMacroBins = outputBins[dim].size();
					for (int binIndex = 0; binIndex < tmpNumMacroBins; binIndex++)
					{
						tmpMacroBin = outputBins[dim].get(binIndex);
						//attName = prefix + "a" + Integer.toString(binIndex);
						attName = Integer.toString(nominalIndex);
						tmpMacroBin.name = attName;
						if (binIndex == 0)
							writerData.write(attName);
						else
							writerData.write("," + attName);
						nominalIndex++;
					} // end for
					
					writerData.write("}");
					writerData.newLine();
				} // end for
				
				// write categorical data
				for (int dim = 0; dim < numCatCols; dim++)
				{
					prefix = Integer.toString(numCols + dim);
					writerData.write("@attribute dim" + prefix + " {");
					tmpNumMacroBins = singleContexts[dim].size();
					for (int index = 0; index < tmpNumMacroBins; index++)
					{
						attName = singleContexts[dim].get(index).toString();
						if (index == 0)
							writerData.write(attName);
						else
							writerData.write("," + attName);
					} // end for
					
					writerData.write("}");
					writerData.newLine();
				} // end for
				
				// discretize data
				boolean attNameExists;
				for (int i = 0; i < numRows; i++)
				{
					tmpPoint = dataMatrix.data.get(i);
					pointLabels[i] = new String(tmpPoint.classID);
					for (int dim = 0; dim < numCols; dim++)
					{
						attNameExists = false;
						tmpNumMacroBins = outputBins[dim].size();
						for (int binIndex = 0; binIndex < tmpNumMacroBins; binIndex++)
						{
							//attName = dim + "a" + binIndex;
							tmpMacroBin = outputBins[dim].get(binIndex);
							tmpLowerBound = tmpMacroBin.lowerBound;
							tmpUpperBound = tmpMacroBin.upperBound;
							if (tmpPoint.measures[dim] <= tmpUpperBound && tmpPoint.measures[dim] >= tmpLowerBound)
							//if (tmpPoint.measures[dim] <= tmpUpperBound)
							{
								//outputData[i][dim] = attName;
								outputData[i][dim] = tmpMacroBin.name;
								attNameExists = true;
								break;
							}
						} // end for
						
						if (attNameExists == false)
						{
							//System.out.println(finalBins[dim].get(tmpNumMacroBins - 1).upperBound);
							System.out.println(i + " --- " + dim + " --- " + tmpPoint.measures[dim]);
							throw new Exception("Attribute nominal value not found");
						}
					} // end for
				} // end for
				
				// write class labels
				writerData.write("@attribute class {");
				int numLabels = Constants.CLASS_LABELS.size();
				for (int i = 0; i < numLabels; i++)
				{
					if (i == 0)
						writerData.write("\"" + Constants.CLASS_LABELS.get(i) + "\"");
					else
						writerData.write("," + "\"" +  Constants.CLASS_LABELS.get(i) + "\"");
				}
				writerData.write("}");
				writerData.newLine();
				writerData.newLine();
				
				// write the actual data
				writerData.write("@data");
				writerData.newLine();
				for (int i = 0; i < numRows; i++)
				{
					tmpPoint = dataMatrix.data.get(i);
					for (int j = 0; j < numCols; j++)
					{
						if (outputData[i][j] == null)
							throw new Exception("Null value");
							
						writerData.write(outputData[i][j] + ",");
					}
					for (int j = 0; j < numCatCols; j++)
						writerData.write(tmpPoint.catContexts[j] + ",");
					writerData.write("\"" + pointLabels[i] + "\"");
					
					if (i != numRows - 1)
						writerData.newLine();
				}
				
				writerData.flush();
				writerData.close();
				double t0 = 4 * Constants.MAX_VAL * (Constants.NUM_MEASURE_COLS - 1) * (2 + 1 / Math.log(2)) * Constants.DELTA;
				System.out.println("threshold " + t0);
			}
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}
	
	public static double getOriginalValue(double value, double min, double max)
	{
		double ret = min;
		if (max > min) 
			ret = (max - min) * (value + Constants.MAX_VAL) / (2 * Constants.MAX_VAL) + min;
		return ret;
	}
	
	public static void readInputString(String[] args) throws Exception
	{	
		int i;
		int total = args.length - 1;
		
		// take FILE_INPUT
		boolean found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-FILE_INPUT"))
			{
				Constants.FILE_INPUT = args[i + 1];
				found = true;
				break;
			}
		if (found == false)
			throw new Exception("Missing -FILE_INPUT");
		
		// take FILE_IDIST
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-FILE_IDIST"))
			{
				Constants.FILE_IDIST = args[i + 1];
				found = true;
				break;
			}
		
		// take FILE_CP_OUTPUT
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-FILE_CP_OUTPUT"))
			{
				Constants.FILE_CP_OUTPUT = args[i + 1];
				found = true;
				break;
			}
		if (found == false)
			throw new Exception("Missing -FILE_CP_OUTPUT");
		
		// take FILE_RUNTIME_OUTPUT
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-FILE_RUNTIME_OUTPUT"))
			{
				Constants.FILE_RUNTIME_OUTPUT = args[i + 1];
				found = true;
				break;
			}
		if (found == false)
			throw new Exception("Missing -FILE_RUNTIME_OUTPUT");
		
		// take FILE_DATA_OUTPUT
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-FILE_DATA_OUTPUT"))
			{
				Constants.FILE_DATA_OUTPUT = args[i + 1];
				found = true;
				break;
			}
		if (found == false)
			throw new Exception("Missing -FILE_DATA_OUTPUT");
		
		// take NUM_ROWS
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-NUM_ROWS"))
			{
				Constants.NUM_ROWS = Integer.parseInt(args[i + 1]);
				found = true;
				break;
			}
		if (found == false)
			throw new Exception("Missing -NUM_ROWS");
		
		// take NUM_MEASURE_COLS
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-NUM_MEASURE_COLS"))
			{
				Constants.NUM_MEASURE_COLS = Integer.parseInt(args[i + 1]);
				found = true;
				break;
			}
		if (found == false)
			throw new Exception("Missing -NUM_MEASURE_COLS");
		
		// take NUM_CAT_CONTEXT_COLS
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-NUM_CAT_CONTEXT_COLS"))
			{
				Constants.NUM_CAT_CONTEXT_COLS = Integer.parseInt(args[i + 1]);
				found = true;
				break;
			}
		if (found == false)
			throw new Exception("Missing -NUM_CAT_CONTEXT_COLS");
		
		// take FIELD_DELIMITER
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-FIELD_DELIMITER"))
			{
				Constants.FIELD_DELIMITER = args[i + 1];
				found = true;
				break;
			}
		//if (found == false)
		//	throw new Exception("Missing -FIELD_DELIMITER");
		
		// take MAX_VAL
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-MAX_VAL"))
			{
				Constants.MAX_VAL = Double.parseDouble(args[i + 1]);
				found = true;
				break;
			}
		if (found == false)
			throw new Exception("Missing -MAX_VAL");
		
		// take METHOD
		found = false;
		for (i = 0; i < total; i++)
			if (args[i].equals("-METHOD"))
			{
				Constants.METHOD = Integer.parseInt(args[i + 1]);
				found = true;
				break;
			}
		if (found == false)
			throw new Exception("Missing -METHOD");
		
		if (Constants.METHOD == Constants.MB_MDL_DP || Constants.METHOD == Constants.MB_MDL_GD || Constants.METHOD == Constants.MB_MDL_EGD || Constants.METHOD == Constants.MB_NM_GD || Constants.METHOD == Constants.DP_MEAN || Constants.METHOD == Constants.DP_MIN || Constants.METHOD == Constants.DP_MAX)
		{
			// take QUANTILE
//			found = false;
//			for (i = 0; i < total; i++)
//				if (args[i].equals("-QUANTILE"))
//				{
//					Constants.QUANTILE = Integer.parseInt(args[i + 1]);
//					found = true;
//					break;
//				}
//			if (found == false)
//				throw new Exception("Missing -QUANTILE");
//			
//			// take INIT_BIN_COUNT
//			found = false;
//			for (i = 0; i < total; i++)
//				if (args[i].equals("-INIT_BIN_COUNT"))
//				{
//					Constants.INIT_BIN_COUNT = Integer.parseInt(args[i + 1]);
//					found = true;
//					break;
//				}
//			
//			// take USE_CE
//			found = false;
//			for (i = 0; i < total; i++)
//				if (args[i].equals("-USE_CE"))
//				{
//					Constants.USE_CE = Boolean.parseBoolean(args[i + 1]);
//					found = true;
//					break;
//				}
			
			
		}
		else if (Constants.METHOD == Constants.J_MDL || Constants.METHOD == Constants.J_MDL_HM || Constants.METHOD == Constants.M_GD || Constants.METHOD == Constants.MMIC)
		{
			// take USE_CE
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-USE_CE"))
				{
					Constants.USE_CE = Boolean.parseBoolean(args[i + 1]);
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -USE_CE");
			
			// take NUM_CENTROIDS
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-NUM_CENTROIDS"))
				{
					Constants.NUM_CENTROIDS = Integer.parseInt(args[i + 1]);
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -NUM_CENTROIDS");
			
			// take INIT_BIN_COUNT
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-INIT_BIN_COUNT"))
				{
					Constants.INIT_BIN_COUNT = Integer.parseInt(args[i + 1]);
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -INIT_BIN_COUNT");
		}
		else if (Constants.METHOD == Constants.J_GD_HM || Constants.METHOD == Constants.M_GD_CAIM)
		{
			// take USE_CE
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-USE_CE"))
				{
					Constants.USE_CE = Boolean.parseBoolean(args[i + 1]);
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -USE_CE");
			
			// take NUM_CENTROIDS
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-NUM_CENTROIDS"))
				{
					Constants.NUM_CENTROIDS = Integer.parseInt(args[i + 1]);
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -NUM_CENTROIDS");
		}
		else if (Constants.METHOD == Constants.EF || Constants.METHOD == Constants.EW || Constants.METHOD == Constants.MVD)
		{
			// take INIT_BIN_COUNT
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-INIT_BIN_COUNT"))
				{
					Constants.INIT_BIN_COUNT = Integer.parseInt(args[i + 1]);
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -INIT_BIN_COUNT");
			
			// take GAMMA
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-GAMMA"))
				{
					Constants.GAMMA = Double.parseDouble(args[i + 1]);
					found = true;
					break;
				}
		}
		else if (Constants.METHOD == Constants.PCABINNING)
		{
			// take INIT_BIN_COUNT
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-INIT_BIN_COUNT"))
				{
					Constants.INIT_BIN_COUNT = Integer.parseInt(args[i + 1]);
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -INIT_BIN_COUNT");
			
			// take FILE_ORG_INPUT
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-FILE_ORG_INPUT"))
				{
					Constants.FILE_ORG_INPUT = args[i + 1];
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -FILE_ORG_INPUT");
			
			// take FILE_ORG_NORMAL_INPUT
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-FILE_ORG_NORMAL_INPUT"))
				{
					Constants.FILE_ORG_NORMAL_INPUT = args[i + 1];
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -FILE_ORG_NORMAL_INPUT");
			
			// take NUM_ORG_CAT_CONTEXT_COLS
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-NUM_ORG_CAT_CONTEXT_COLS"))
				{
					Constants.NUM_ORG_CAT_CONTEXT_COLS = Integer.parseInt(args[i + 1]);
					found = true;
					break;
				}
			if (found == false)
				throw new Exception("Missing -NUM_ORG_CAT_CONTEXT_COLS");
			
			// take GAMMA
			found = false;
			for (i = 0; i < total; i++)
				if (args[i].equals("-GAMMA"))
				{
					Constants.GAMMA = Double.parseDouble(args[i + 1]);
					found = true;
					break;
				}
		}
	} // end method
}
