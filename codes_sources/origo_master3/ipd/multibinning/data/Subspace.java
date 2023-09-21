package multibinning.data;

import java.util.*;

public class Subspace 
{
	public int num_dims;
	public ArrayList<Integer> dims;
	
	public Subspace()
	{
		// do nothing
	}
	
	public Subspace(int num_dims)
	{
		this.num_dims = num_dims;
		dims = new ArrayList<Integer>();
	}
}
