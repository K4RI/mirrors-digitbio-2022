from itemset import closure;
import sys;


sets = [frozenset(list(map(int, line.split()))) for line in open(sys.argv[1])];

res = closure(sets);

out = open(sys.argv[2], 'w');

for x in res:
	print(' '.join(map(str, x))); 
	print(' '.join(map(str, x)), file=out); 
