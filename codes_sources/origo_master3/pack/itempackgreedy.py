import sys;
from math import log;
from mdl import binomial_mdl;
from graph import dmst, Dag;
from out import out_dag, out_scheme, out_sets, out_model;
import os, profile;
import itemset;
from functools import reduce

def binmdl(n):
	if bmdl[n] == -1:
		bmdl[n] = binomial_mdl(n);
	return bmdl[n];

def split(v, data):
	p = [];
	n = [];

	for t in data:
		if v in t:
			p.append(t);
		else:
			n.append(t);
	return [p, n];

class Tree:
	def __init__(s, sets, cands, path, fr, target, K):
		s.trivial = False;
		s.leaf = True;
		s.neg = None;
		s.pos = None;
		s.attr = None;
		s.cands = frozenset();

		s.sets = sets;
		s.path = path;
		s.fr = fr;
		s.target = target;
		s.K = K;

		if len(sets) > 0:
			s.cands = cands & reduce(lambda x, y: x | y, sets);

		s.saved = 0;
		s.score = s.bits();
		s.total = s.score;

	def split(s, neg, pos, attr):
		s.leaf = False;
		s.neg = neg;
		s.pos = pos;
		s.attr = attr;

	def update_total(s):
		if s.leaf == False:
			s.total = s.score + s.neg.update_total() + s.pos.update_total();
		else:
			s.total = s.score;
		return s.total;

	def deps(s):
		if s.leaf == False:
			return s.neg.deps() | s.pos.deps() | frozenset([s.attr]);
		else:
			return frozenset();

	def bits(s):
		entr = lambda x: -x*log(x);
		pen = lambda x: binmdl(int(x)) + log(2);

		x = itemset.freqexact(s.fr, s.path[0], s.path[1]);
		y = itemset.freqexact(s.fr, s.path[0] | frozenset([s.target]), s.path[1]);


		f = [y, x - y];
		s.pcount = y;
		s.ncount = x - y;

		n = sum(f);
		if min(f) == 0:
			s.trivial = True;
			return pen(n);
		if n == 0:
			return 0;

		return sum(map(entr, f)) - entr(n) + pen(n);

	def __str__(s):
		if not s.leaf:
			return str(s.attr) + ' (' + str(s.pos) + ', ' + str(s.neg) + ') ';
		else:
			return '';


	def best_tree(s):
		orig_score = s.total;
		pen = log(2) + log(s.K);

		if s.trivial or len(s.cands) == 0:
			return;

		trees = [];
		for c in s.cands:
			newcands = s.cands - set([c]);
			psets = split(c, s.sets)[0];

			pos = Tree(psets, newcands, [s.path[0] | set([c]), s.path[1]], s.fr, s.target, s.K);
			neg = Tree(psets, newcands, [s.path[0], s.path[1] | set([c])], s.fr, s.target, s.K);
			trees.append([pos.total + neg.total, pos, neg, c]);

		t = min(trees);
		pos = t[1];
		neg = t[2];

		pos.best_tree();
		neg.best_tree();

		if s.total > pos.total + neg.total + pen:
			s.split(neg, pos, t[3]);
			s.score = pen;
			s.total = pos.total + neg.total + s.score;

		s.saved = orig_score - s.total;



class Scheme:
	def __init__(s, sets, fr, target, K, cands):
		s.target = target;
		[psets, nsets] = split(target, sets);
		s.tree = Tree(psets, cands - frozenset([target]), [frozenset(), frozenset()], fr, target, K);

def find_dag(cands, fr, K):
	whites = [frozenset() for i in range(K)];
	cache = [{} for i in range(K)];

	splitsets = [split(t, cands)[0] for t in range(K)];

	deps = [frozenset() for i in range(K)];
	for i in range(K):
		if len(splitsets[i]) > 0:
			deps[i] = reduce(lambda x, y: x | y, splitsets[i]);
		deps[i] -= frozenset([i]);

	sinks = set();

	for i in range(K):
		s = Tree(splitsets[i], frozenset(), [frozenset(), frozenset()], fr, i, K);
		cache[i][frozenset()] = s.total;


	while True:
		g = [[]];

		print('Nodes left:', K - len(sinks));
		for i in range(K):
			print(i, '\r', end=' ')
			sys.stdout.flush();
			w = deps[i] & whites[i];
			edges = [[cache[i][w], 0]];
			if not i in sinks:
				t = deps[i] - whites[i];
				for j in deps[i] - whites[i]:
					src = frozenset([j]) | w;
					if src in cache[i]:
						total = cache[i][src];
					else:
						s = Tree(splitsets[i], src, [frozenset(), frozenset()], fr, i, K);
						#s = Tree(splitsets[i], src, splits[i], K);
						s.best_tree();
						total = cache[i][src] = s.total;
					if total < edges[0][0]:
						edges.append([total, j + 1]);
			g.append(edges);


		edges = dmst(g)

		sinks |= set([x for x in range(len(edges)) if edges[x] == 0]);

		if len(sinks) == K:
			return whites;

		dag = Dag(K);
		for i in range(K):
			if edges[i] > 0:
				dag.add(i, edges[i] - 1);

		for i in range(K):
			if not i in sinks:
				whites[i] |= dag.children[i] & sinks;

i = 1;
depfile = None;
srcfile = None;
treefile = None;
labels = None;
setfile = None;
mdlfile = None;
cands = frozenset();

while (i < len(sys.argv)):
	if (sys.argv[i] == '-out'):
		treefile = sys.argv[i + 1];
		i += 1;
	elif (sys.argv[i] == '-model'):
		mdlfile = sys.argv[i + 1];
		i += 1;
	elif (sys.argv[i] == '-cands'):
		cands = [frozenset(list(map(int, line.split()))) for line in open(sys.argv[i + 1])];
		i += 1;
	elif (sys.argv[i] == '-dep'):
		depfile = sys.argv[i + 1];
		i += 1;
	elif (sys.argv[i] == '-labels'):
		labels = [line.strip() for line in open(sys.argv[i + 1])];
		i += 1;
	elif (sys.argv[i] == '-set'):
		setfile = sys.argv[i + 1];
		i += 1;
	else:
		srcfile = sys.argv[i];
	i += 1;


data = [frozenset(list(map(int, line.split()))) for line in open(srcfile)];
K = 1 + max(list(map(max, list(filter(bool, data)))));

if labels == None:
	labels = list(map(str, list(range(K))));

starttime = os.times()[0];

print('Computing frequencies');
fr = itemset.freqs(itemset.transpose(data, K), cands, len(data));
print('Done');

bmdl = [-1] * (1 + len(data));


print('\nbuilding dag');
whites = find_dag(cands, fr, K);

# Recompute the tree with correct candidates - a little bit silly
# but the alternative is to save the trees.

schemes = [Scheme(cands, fr, i, K, whites[i]) for i in range(K)];
orig_score = sum([s.tree.total for s in schemes]);
for i in range(K):
	 schemes[i].tree.best_tree();

deps = Dag(K);
for i in range(K):
	p = schemes[i].tree.deps();
	for y in p:
		deps.add(i, y);

endtime = os.times()[0];

print('time:', endtime - starttime);

if depfile:
	out_dag(deps, open(depfile, 'w'), labels);

if treefile:
	out_scheme(schemes, open(treefile, 'w'), orig_score, labels);

if setfile:
	out_sets(schemes, open(setfile, 'w'));

if mdlfile:
	out_model(schemes, open(mdlfile, 'w'))


total_score = sum([s.tree.total for s in schemes]);
print('Ratio: %.1f / %.1f (%.2f)' % (total_score / log(2), orig_score / log(2), total_score / orig_score));
