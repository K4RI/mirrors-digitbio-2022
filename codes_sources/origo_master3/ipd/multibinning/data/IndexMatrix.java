package multibinning.data;

public class IndexMatrix 
{
	public int rows;
	public int cols;
	public int[][] data;
	
	public IndexMatrix()
	{
		// do nothing
	}
			
	public IndexMatrix(int rows, int cols)
	{
		this.rows = rows;
		this.cols = cols;
		this.data = new int[rows][cols];
	}
	
	public IndexMatrix(int rows, int cols, int[][] data)
	{
		this.rows = rows;
		this.cols = cols;
		this.data = data;
	}
}