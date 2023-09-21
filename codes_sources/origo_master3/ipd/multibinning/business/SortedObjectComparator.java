package multibinning.business;

import java.util.*;

public class SortedObjectComparator implements Comparator<SortedObject>
{
	public int compare(SortedObject o1, SortedObject o2)
	{
		if (o1.value > o2.value)
			return 1;
		
		if (o1.value == o2.value)
			return 0;
		
		return -1;
	}
}
