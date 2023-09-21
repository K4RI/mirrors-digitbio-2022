from itemset import transpose;
from math import log;
import sys;

def entropy(x, M):
	N = float(len(x));
	p = [0] * M;
	for z in x:
		p[z] += 1;

	res = 0.0;
	for i in range(M):
		if p[i] > 0:
			res -= p[i]/N*log(p[i]/N);
	return res / log(2);


def mine(thresh, cond, prefix, data, attrs, out):
	ents = [];
	conds = [];
	cands = [];

	for x in attrs:
		t = [a << 1 for a in cond];
		for y in data[x]: t[y] += 1;
		e = entropy(t, 1 << (len(prefix) + 1));
		if e <= thresh:
			ents.append(e);
			conds.append(t);
			cands.append(x);

	for i in range(len(cands)):
		x = cands[i];
		leset = prefix | set([x]);
		out.append((ents[i], leset));
		mine(thresh, conds[i], leset, data, cands[(i + 1):], out);


def main():
	data = [list(map(int, line.split())) for line in open(sys.argv[1])];
	thresh = float(sys.argv[2]);

	K = 1 + max(list(map(max, list(filter(bool, data)))));
	N = len(data);

	lesets = [];
	mine(thresh, [0]*N, set(), transpose(data, K), list(range(K)), lesets);

	if len(sys.argv) >= 4:
		out = open(sys.argv[3], 'w');
		for s in lesets:
			print(' '.join(map(str, s[1])), file=out);

	if len(sys.argv) >= 5:
		out = open(sys.argv[4], 'w');
		for s in lesets:
			print(s[0], file=out);

if __name__ == "__main__":
    main()
