package multibinning.business;

import java.io.*;
import java.util.ArrayList;

import multibinning.data.*;

import weka.filters.Filter;
import weka.attributeSelection.AttributeSelection;
import weka.attributeSelection.Ranker;
//import weka.core.Instance;
import weka.core.Instances;

public class PCABinning
{
	@SuppressWarnings("unchecked")
	public static ArrayList<MacroBin>[] discretizeData(DataMatrix dataMatrix, String orgInputFile) throws Exception
	{
		int numCols = dataMatrix.cols;
		//int numRows = dataMatrix.rows;
		
		ArrayList<MacroBin>[] pcaBins = MVD.discretizeData(dataMatrix);
		ArrayList<MacroBin>[] outputBins = new ArrayList[numCols];
		MacroBin tmpMacroBin = null;
		MacroBin tmpPrevMacroBin = null;
		int tmpNumMacroBins;
		for (int dim = 0; dim < numCols; dim++)
		{
			outputBins[dim] = new ArrayList<MacroBin>();
			tmpNumMacroBins = pcaBins[dim].size();
			tmpPrevMacroBin = pcaBins[dim].get(0);
			outputBins[dim].add(tmpPrevMacroBin);
			for (int binIndex = 1; binIndex < tmpNumMacroBins; binIndex++)
			{
				tmpMacroBin = pcaBins[dim].get(binIndex);
				if (tmpMacroBin.upperBound > tmpPrevMacroBin.upperBound)
					outputBins[dim].add(tmpMacroBin);
				tmpPrevMacroBin = tmpMacroBin;
			}
		}
		
		Instances inputData = getData(orgInputFile);
		int numOrgCols = inputData.get(0).toDoubleArray().length - 1;
		ArrayList<MacroBin>[] ret = new ArrayList[numOrgCols];
		/*String splitStr = "";
		for (int i = numOrgCols - 1; i >= 0; i--)
		{
			if (i == numOrgCols - 1)
				splitStr += "dim" + Integer.toString(i);
			else
				splitStr += "|dim" + Integer.toString(i);
		}*/
		
		// get coefficients of each attributes
		weka.filters.unsupervised.attribute.PrincipalComponents pca = new weka.filters.unsupervised.attribute.PrincipalComponents();
		pca.setVarianceCovered(0.95);
		pca.setMaximumAttributeNames(-1);
		pca.setInputFormat(inputData);
		Instances transformedData = Filter.useFilter(inputData, pca);
		//String[] vals = null;
		//String vals = null;
		String tmpVal = null;
		String str = null;
		//int pos;
		double[][] coefs = new double[numCols][numOrgCols];
		int curDim;
		int index;
		int[] dimIndices = new int[numOrgCols];
		int[] dims = new int[numOrgCols];
		boolean[] processedDims = new boolean[numOrgCols];
		boolean notExists;
		//String dimConst = "dim";
		int count;
		for (int i = 0; i < numCols; i++)
		{
			tmpVal = transformedData.attribute(i).name();
			System.out.println(tmpVal);
			
			for (int j = 0; j < numOrgCols; j++)
			{
				dimIndices[j] = -1;
				processedDims[j] = false;
			}
			
			count = 0;
			index = tmpVal.indexOf("dim");
			dimIndices[count++] = index;
			while (index >= 0) 
			{
				index = tmpVal.indexOf("dim", index + 1);
				dimIndices[count++] = index;
				if (count >= numOrgCols)
					break;
			}
			
			for (int j = 0; j < numOrgCols; j++)
			{
				notExists = false;
				if (j != numOrgCols - 1)
				{
					//System.out.println("dimIndices[j] = " + dimIndices[j] + " --- dimIndices[j + 1] = " + dimIndices[j + 1]);
					if (dimIndices[j] != -1)
					{
						if (dimIndices[j + 1] != -1)
							str = tmpVal.substring(dimIndices[j] + 3, dimIndices[j + 1]);
						else
							str = tmpVal.substring(dimIndices[j] + 3);
					}
					else
					{
						notExists = true;
						for (int k = 0; k < numOrgCols; k++)
							if (processedDims[k] == false)
							{
								str = Integer.toString(k);
								break;
							}
					}
					
					str = str.trim();
					if (str.length() != 0)
						str = str.split("\\+|-")[0];
					curDim = Integer.parseInt(str);
					if (curDim > numOrgCols)
						str = str.substring(0, str.length() - 1);
				}
				else
				{
					if (dimIndices[j] != -1)
						str = tmpVal.substring(dimIndices[j] + 3);
					else
					{
						notExists = true;
						for (int k = 0; k < numOrgCols; k++)
							if (processedDims[k] == false)
							{
								str = Integer.toString(k);
								break;
							}
					}
				}
				curDim = Integer.parseInt(str);
				processedDims[curDim] = true;
				dims[j] = curDim;
				System.out.println("curDim = " + curDim);
				
				if (notExists == false)
				{
					if (j == 0)
						str = tmpVal.substring(0, dimIndices[j]);
					else
					{
						str = tmpVal.substring(dimIndices[j - 1] + 3 + Integer.toString(dims[j - 1]).length(), dimIndices[j]);
						str = str.trim();
					}
					if (str.trim().length() == 0)
						coefs[i][curDim] = 0;
					else
						coefs[i][curDim] = Double.parseDouble(str.trim());
				}
				else
					coefs[i][curDim] = 0;
				System.out.println("coefs = " + coefs[i][curDim]);
				System.out.println("---");
			}
		}
		
		// get eigenvalues
		AttributeSelection selector = new AttributeSelection();
		weka.attributeSelection.PrincipalComponents pcaOther = new weka.attributeSelection.PrincipalComponents();
		pca.setVarianceCovered(0.95);
		pca.setMaximumAttributeNames(-1);
		Ranker ranker = new Ranker();
		selector.setEvaluator(pcaOther);
		selector.setSearch(ranker);
		selector.SelectAttributes(inputData);
		double[][] rankedAtts = selector.rankedAttributes();
		double[] eigenValues = new double[rankedAtts.length];
		for (int i = 0; i < rankedAtts.length; i++)
			eigenValues[i] = rankedAtts[i][1];
		
		// set up bins for each dimensions
		int numMacroBins;
		double scaledVal;
		double maxVal;
		int maxDim;
		double[] cps = null;
		for (int dim = 0; dim < numOrgCols; dim++)
		{
			maxVal = -Double.MAX_VALUE;
			maxDim = -1;
			for (int pcaDim = 0; pcaDim < numCols; pcaDim++)
			{
				scaledVal = Math.abs(coefs[pcaDim][dim] / eigenValues[pcaDim]);
				if (scaledVal > maxVal)
				{
					maxVal = scaledVal;
					maxDim = pcaDim;
				}
			}
			
			ret[dim] = new ArrayList<MacroBin>();
			numMacroBins = outputBins[maxDim].size();
			cps = new double[numMacroBins];
			for (int binIndex = 0; binIndex < numMacroBins; binIndex++)
			{
				tmpMacroBin = outputBins[maxDim].get(binIndex);
				cps[binIndex] = tmpMacroBin.upperBound * coefs[maxDim][dim];
			}

			ret[dim].add(new MacroBin(-Double.MAX_VALUE, cps[0]));
			for (int binIndex = 0; binIndex < numMacroBins - 1; binIndex++)
				ret[dim].add(new MacroBin(cps[binIndex], cps[binIndex + 1]));
			ret[dim].add(new MacroBin(cps[numMacroBins - 1], Double.MAX_VALUE));
		}
		
		Constants.NUM_MEASURE_COLS = numOrgCols;
		Constants.NUM_CAT_CONTEXT_COLS = Constants.NUM_ORG_CAT_CONTEXT_COLS;
		
		return ret;
	}
	
	public static Instances getData(String fileName) throws Exception
	{
		
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		Instances data = new Instances(reader);
		reader.close();
		data.setClassIndex(data.numAttributes() - 1);
		return data;
	} // end method
}
