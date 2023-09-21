from heapq import heappush, heappop;
import copy;

def scc(g, allowed):

	visited = [False] * len(g);

	stacked = [False] * len(g);
	comp = [];
	index = [-1] * len(g);
	lowlink = [-1] * len(g);
	res = [];

	while visited.count(False) > 0:
		res.extend(tarjan(g, visited.index(False), index, lowlink, comp, stacked, [0], visited, allowed));
		allowed = [allowed[x] or visited[x] for x in range(len(g))];

	return res;


def tarjan(g, v, index, lowlink, comp, stacked, cur, visited, allowed):
	res = [];
	index[v] = cur[0];
	lowlink[v] = cur[0];
	comp.append(v);
	stacked[v] = True;
	visited[v] = True;
	cur[0] += 1;
	
	for e in g[v]:
		if not allowed[e]:
			continue;
		if index[e] < 0:
			res.extend(tarjan(g, e, index, lowlink, comp, stacked, cur, visited, allowed));
			lowlink[v] = min(lowlink[v], lowlink[e]);
		elif stacked[e]:
			lowlink[v] = min(lowlink[v], index[e]);

	if (lowlink[v] == index[v]):
		cycle = [];
		while True:
			e = comp.pop();
			stacked[e] = False;
			cycle.append(e);
			if e == v:
				break;
		cycle.sort();
		res.append(cycle);
			

	return res;

def shortest_path(g, s, t, allowed):
	parent = [-1] * len(g);
	visit = [False] * len(g);
	stack = [s];
	visit[s] = True;

	while not visit[t]:
		n = stack.pop(0);
		for m in g[n]:
			if allowed[m] and not visit[m]:
				visit[m] = True;
				parent[m] = n;
				stack.append(m);

	path = [];
	n = t;
	while n != s:
		path.append(n);
		n = parent[n];
	path.append(s);
	return path;

def contractgraph(g, cycle, K, edges):
	mincost = min([edges[x][0] for x in cycle]);
	nodes = set(cycle);
	inedges = [];

	for x in cycle:
		inedges.extend([[y[0] - edges[x][0] + mincost] + y[1:] + [x] for y in g[x] if not y[1] in nodes]);

	for i in range(len(g)):
		if i in nodes:
			continue;
		for x in g[i]:
			if x[1] in nodes:
				x.append(x[1]);
				x[1] = K;

	g.append(inedges);

def dmst(g):
	N = len(g);
	g = copy.deepcopy(g);

	edges = [None] + [min(x) for x in g[1:]];
	visited = [False] * len(g);

	#print [x[1:] for x in edges[1:]];


	visited[0] = True;

	active = [False] * len(g);
	n = visited.index(False);

	while visited.count(False) > 0:

		while not visited[n]:
			visited[n] = True;
			active[n] = True;
			n = edges[n][1];

		if n == 0 or not active[n]:
			if visited.count(False) == 0:
				break;
			n = visited.index(False);
			active = [False] * len(g);
		elif active[n]: # Cycle found!
			x = edges[n][1];
			cycle = [n];
			while x != n:
				cycle.append(x);
				x = edges[x][1];

			#print cycle;
			K = len(g);
			contractgraph(g, cycle, K, edges);


			visited.append(False);
			active.append(False);
			edges.append(min(g[-1]));

			#print [x[1:] for x in edges[1:]];

			n = K;

	#print [x[1:] for x in edges[1:]];

	while len(edges) > N:
		k = len(edges) - 1;
		#print k;
		#print [x[1:] for x in edges[1:]];
		m = edges[-1][-1];
		edges[m] = [edges[m][0]] + edges[-1][1:-1];
		edges.pop();

		for i in range(1, len(edges)):
			if edges[i][1] == k:
				edges[i][1] = edges[i].pop();

	return [x[1] for x in edges[1:]];

class Dag:
	def __init__(s, K):
		s.edges = [frozenset([]) for i in range(K)];
		s.children = [frozenset([]) for i in range(K)];

	def add(s, par, child):
		s.edges[child] |= frozenset([par]);
		s.update(par, s.children[child] | frozenset([child]));

	def update(s, node, children):
		s.children[node] |= children;
		for c in s.edges[node]:
			s.update(c, children);


print (dmst([[], [[2, 0], [1, 3]], [[1, 1]], [[1, 2]],[[1, 2]]   ]));
