package multibinning.business;

import java.io.*;
import java.util.*;

import weka.core.*;
import weka.filters.Filter;
import weka.filters.supervised.attribute.Discretize;

public class SupUniMDL 
{
	protected static Instances load(String filename) throws Exception 
	{
		Instances result;
	    BufferedReader reader;
	 
	    reader = new BufferedReader(new FileReader(filename));
	    result = new Instances(reader);
	    result.setClassIndex(result.numAttributes() - 1);
	    reader.close();
	 
	    return result;
	}
	
	protected static void save(Instances data, String filename) throws Exception 
	{
	    BufferedWriter writer;
	 
	    writer = new BufferedWriter(new FileWriter(filename));
	    writer.write(data.toString());
	    writer.newLine();
	    writer.flush();
	    writer.close();
	}
	
	public static void Discretize(String[] args) throws Exception 
	{
	    Instances inputTrain;
	    Instances outputTrain;
	    Discretize filter;
	    
	    // get number of categorical columns
	    int numCatCols = Integer.parseInt(args[1]);
	 
	    // load data (class attribute is assumed to be last attribute)
	    inputTrain = load(args[2]);
	    
	    // setup filter
	    filter = new Discretize();
	    filter.setInputFormat(inputTrain);
	    int method = Integer.parseInt(args[0]);
	    if (method == Constants.KoMDL)
	    	filter.setUseKononenko(true);
	    else
	    	filter.setUseKononenko(false);
	    
	    // apply filter
	    long start = System.currentTimeMillis();
	    outputTrain = Filter.useFilter(inputTrain, filter);
	    long end = System.currentTimeMillis();
	    BufferedWriter writerRuntime = new BufferedWriter(new FileWriter(new File(args[4])));
	    writerRuntime.write(Long.toString(end - start));
		writerRuntime.flush();
		writerRuntime.close();

		// open categorical data if any
		ArrayList<String> catHeaders = null;
		ArrayList<String> catData = null;
		if (numCatCols > 0)
		{
			catHeaders = new ArrayList<String>();
			catData = new ArrayList<String>();
			
			Scanner scanner = new Scanner(new FileInputStream(args[5]));
			String s = null;
			int count = 0;
			while (scanner.hasNextLine())
			{
				s = scanner.nextLine();
				count++;
				if (count <= numCatCols)
					catHeaders.add(new String(s));
				else
					catData.add(new String(s));
			}
		}
	    
	    // save output
		BufferedWriter writerData = new BufferedWriter(new FileWriter(new File(args[3])));
		writerData.write("@relation DB");
		writerData.newLine();
		writerData.newLine();
		
		int numRows = outputTrain.size();
		int numCols = outputTrain.get(0).toDoubleArray().length - 1;
		String prefix = null;
		Attribute curAtt = null;
		int tmpNumMacroBins;
		ArrayList<ArrayList<String>> attNames = new ArrayList<ArrayList<String>>();
		ArrayList<String> tmpAtt = null;
		String attName = null;
		int nominalIndex = 1;
		for (int dim = 0; dim < numCols; dim++)
		{
			prefix = Integer.toString(dim);
			writerData.write("@attribute dim" + prefix + " {");
			curAtt = outputTrain.attribute(dim);
			tmpNumMacroBins = curAtt.numValues();
			tmpAtt = new ArrayList<String>();
			for (int binIndex = 0; binIndex < tmpNumMacroBins; binIndex++)
			{
				//attName = curAtt.value(binIndex);
				attName = Integer.toString(nominalIndex);
				tmpAtt.add(attName);
				if (binIndex == 0)
					writerData.write(attName);
				else
					writerData.write("," + attName);
				nominalIndex++;
			} // end for
			attNames.add(tmpAtt);
			
			writerData.write("}");
			writerData.newLine();
		} // end for
		
		// write categorical headers
		for (int dim = 0; dim < numCatCols; dim++)
		{
			writerData.write(catHeaders.get(dim));
			writerData.newLine();
		} // end for
		
		// write class labels
		curAtt = outputTrain.attribute(numCols);
		writerData.write("@attribute class {");
		int numLabels = curAtt.numValues();
		for (int i = 0; i < numLabels; i++)
		{
			if (i == 0)
				writerData.write("\"" + curAtt.value(i) + "\"");
			else
				writerData.write(",\"" + curAtt.value(i) + "\"");
		}
		writerData.write("}");
		writerData.newLine();
		writerData.newLine();
		
		// write the actual data
		writerData.write("@data");
		writerData.newLine();
		Instance tmpPoint = null;
		String tmpValue = null;
		for (int i = 0; i < numRows; i++)
		{
			tmpPoint = outputTrain.get(i);
			for (int j = 0; j < numCols; j++)
			{
				tmpValue = tmpPoint.stringValue(j);
				curAtt = outputTrain.attribute(j);
				tmpNumMacroBins = curAtt.numValues();
				tmpAtt = attNames.get(j);
				for (int binIndex = 0; binIndex < tmpNumMacroBins; binIndex++)
				{
					if (tmpValue.equals(curAtt.value(binIndex)))
					{
						writerData.write(tmpAtt.get(binIndex) + ",");
						//writerData.write(tmpPoint.stringValue(j) + ",");
						break;
					}
				}
			}
			if (numCatCols > 0)
				writerData.write(catData.get(i) + ",");
			writerData.write("\"" + tmpPoint.stringValue(numCols) + "\"");
			
			if (i != numRows - 1)
				writerData.newLine();
		}
		
		writerData.flush();
		writerData.close();
	}
}
