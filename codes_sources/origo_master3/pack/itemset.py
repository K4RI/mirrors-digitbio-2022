from heapq import heappush, heappop;
import sys;

class Track:
	def __init__(self, column):
		self.column = column;
		self.track = 0;

	def end(self):
		return self.track >= len(self.column);

	def value(self):
		return self.column[self.track];

	def advance(self):
		self.track += 1;

	def __cmp__(self, other):
		return self.value() - other.value();

def support(query, data):
	tracks = [];
	m = 0;
	sup = 0;

	for c in query:
		t = Track(data[c]);
		if t.end():
			return 0;
		m = max(m, t.value());
		heappush(tracks, t);

	while (1):
		t = heappop(tracks);
		if (t.value() == m):
			sup += 1;
		t.advance();
		if t.end():
			return sup;
		m = max(t.value(), m);
		heappush(tracks, t);

def transpose(data, K):
	res = [[] for i in range(K)];
	for i in range(len(data)):
		for x in data[i]:
			res[x].append(i);
	return res;

def freqs(data, sets, N):
	f = {};
	K = len(data);
	f[frozenset()] = N;

	for i in range(K):
		f[frozenset([i])] = support([i], data);

	for s in sets:
		if len(s) > 1:
			f[s] = support(list(s), data);
	return f;

def freqexact(f, pos, neg):
	if len(neg) == 0:
		return f[pos];
	for c in neg:
		n = neg - frozenset([c]);
		return freqexact(f, pos, n) - freqexact(f, pos | frozenset([c]), n);

def closeadd(fam, s):
	if not s in fam:
		fam[s] = 1;
		for x in s:
			closeadd(fam, s - set([x]));

def closure(sets):
	res = {frozenset(): 1};
	for s in sets:
		closeadd(res, s);
	return list(res.keys());
