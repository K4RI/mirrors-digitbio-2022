from math import log;

def out_dag(dag, out, labels):
	print('digraph G {', file=out);
	
	for i in range(len(dag.edges)):
		print('\t%d [label = "%s"];' % (i, labels[i]), file=out);
		for p in dag.edges[i]:
			print('\t%d -> %d;' % (p, i), file=out);

	print('}', file=out);

def out_tree(tree, out, node, labels):
	if tree.leaf == False:
		print('\t%s [label = "%s:\\n%.1f (%.1f)"];' % (node, labels[tree.attr], tree.total / log(2), tree.saved / log(2)), file=out);
		nodepos = node + 'p' + str(tree.attr);
		nodeneg = node + 'n' + str(tree.attr);
		out_tree(tree.pos, out, nodepos, labels);
		out_tree(tree.neg, out, nodeneg, labels);
		print('\t%s -> %s [label = "%d"];' % (node, nodepos, 1), file=out);
		print('\t%s -> %s [label = "%d"];' % (node, nodeneg, 0), file=out);
	else:
		print('\t%s [label = "%.1f"];' % (node, tree.total / log(2)), file=out);


def out_scheme(schemes, out, orig_score, labels):
	print('digraph G {', file=out);
	print('\tlabel = "Needed bits: %.1f, originally: %.1f";' % (sum([s.tree.total for s in schemes]) / log(2), orig_score / log(2)), file=out);
	for i in range(len(schemes)):
		print('subgraph cluster%d {' % i, file=out);
		print('\tlabel = "%s:";' % labels[i], file=out);
		out_tree(schemes[i].tree, out, 't' + str(i), labels);
		print('}', file=out);
	print('}', file=out);

def set_comp(x, y):
	if len(x) < len(y):
		return -1;
	elif len(x) > len(y):
		return 1;
	elif x < y:
		return -1;
	elif x > y:
		return 1;
	else:
		return 0;
		
def compute_sets(tree, sets):
	if tree.leaf:
		return sets;
	else:
		n = frozenset([tree.attr]);
		sets2 = [n | x for x in sets];
		return compute_sets(tree.pos, sets2) + compute_sets(tree.neg, sets + sets2);

def out_sets(schemes, out):
	l = [];

	for s in schemes:
		l.extend(compute_sets(s.tree, [frozenset(), frozenset([s.target])]));

	l = [tuple(x) for x in set(l) if len(x) > 0];
	l.sort(set_comp);
	for x in l:
		print(' '.join(map(str, x)), file=out); 

def compute_leaves(tree, p):
	if tree.leaf == False:
		res = compute_leaves(tree.pos, [p[0] | set([tree.attr]), p[1]]);
		res.extend(compute_leaves(tree.neg, [p[0], p[1] | set([tree.attr])]));
		return res;
	else:
		return [[p[0], p[1], tree.pcount, tree.ncount]];

def out_model(schemes, out):
	print(len(schemes), file=out);
	for s in schemes:
		leaves = compute_leaves(s.tree, [set(), set()]);
		print(len(leaves), file=out);
		for l in leaves:
			print(' '.join(map(str, l[0])), file=out);
			print(' '.join(map(str, l[1])), file=out);
			print(l[2], l[3], file=out);
