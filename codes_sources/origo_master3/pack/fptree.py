import sys;

class TreeNode:
	def __init__(self, label, count):
		self.label = label;
		self.count = count;
		self.next = None;
		self.parent = None;
		self.children = [];

	def add(self, child):
		self.children.append(child);
		child.parent = self;

	def find(self, label):
		for child in self.children:
			if child.label == label:
				return child;
		return None;

class FPTree:
	def __init__(self, trans, labels, counts, support):
		[labs, ind] = self.sortlabels(trans, labels, counts, support);
		# labs contains the sorted labels of the nodes in fp tree.
		# ind contains the permutation s.t. ind[labs[i]] = i.
		# ind[j] = - 1 if there is no labs[i] such that labs[i] = j.

		# Create head node list.
		self.labels = labs;
		self.heads = [None] * len(labs);

		# Add root.
		self.root = TreeNode(-1, 0);

		# Go through transactions.
		for i in range(len(trans)):
			t = trans[i];
			n = self.root;
			n.count = n.count + counts[i];

			# prune and sort transaction.
			s = [];
			for v in t:
				if ind[v] == -1:
					continue;
				s.append(ind[v]);
			s.sort();
				
			for v in s:
				c = n.find(self.labels[v]);
				if (c == None):
					c = TreeNode(self.labels[v], 0);
					n.add(c);
					c.next = self.heads[v];
					self.heads[v] = c;
				n = c;
				n.count = n.count + counts[i];
				
	def sortlabels(self, trans, labels, counts, support):
		N = len(labels);
		if len(labels) == 0:
			M = 0;
		else:
			M = max(labels) + 1;

		# Compute means
		s = [0] * M;
		for i in range(len(trans)):
			for v in trans[i]:
				s[v] = s[v] + counts[i];

		# Sort margins.
		margins = [(s[labels[i]], labels[i]) for i in range(N)];
		margins.sort();
		margins.reverse(); # les plus hauts s[labels[i]] en premier

		# Get the reverse permutation.
		ind = [-1] * M;
		K = 0;
		for i in range(N):
			if (margins[i][0] < support): # si pas assez de counts
				break;
			ind[margins[i][1]] = i;
			K = K + 1;

		# Create head node list.
		labs = [m[1] for m in margins[0:K]];

		return [labs, ind];

def frequent(tree, support):
	sets = [[]];
	sups = [tree.root.count];

	for i in range(len(tree.labels)):
		# Create sub-database.
		D = [];
		counts = [];
		node = tree.heads[i];

		while (node != None):
			t = [];
			n = node.parent;
			while (n.label != -1):
				t.append(n.label);
				n = n.parent;
			D.append(t);
			counts.append(node.count);
			node = node.next;

		# Create new fp-tree from sub-database.
		subtree = FPTree(D, tree.labels[0:i], counts, support);

		# Call frequent with subtree.
		sets2, sup2 = frequent(subtree, support);
		for j in range(len(sets2)):
			sets2[j].append(tree.labels[i]);
			sets2[j].sort();

		# Add new items into known set.
		sets = sets + sets2;
		sups = sups + sup2;

	return sets, sups;

x = [list(map(int, line.split())) for line in open(sys.argv[1])];
thresh = int(sys.argv[2]);

tree = FPTree(x, list(range(1 + max(list(map(max, list(filter(bool, x))))))), [1] * len(x), thresh);
(sets, freqs) = frequent(tree, thresh);
out = open("output_pack.txt", 'w');
for s in sets:
	print(' '.join(map(str, s)));
	print(' '.join(map(str, s)), file=out);