package multibinning.data;

import java.util.*;

public class DataMatrix 
{
	public int rows;
	public int cols;
	public int catCols;
	public ArrayList<DataPoint> data;
	
	public DataMatrix()
	{
		// do nothing
	}
	
	public DataMatrix(int rows, int cols, int catCols)
	{
		this.rows = rows;
		this.cols = cols;
		this.catCols = catCols;
		this.data = new ArrayList<DataPoint>();
	}
	
	public DataMatrix(int rows, int cols, int catCols, ArrayList<DataPoint> data)
	{
		this.rows = rows;
		this.cols = cols;
		this.catCols = catCols;
		this.data = data;
	}
}