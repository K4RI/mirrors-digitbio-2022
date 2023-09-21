import sys;
import getopt;
from math import log;

def read_model(twist, name):
	f = open(name);
	K = int(f.readline());
	mdl = [None] * K;
	for i in range(K):
		L = int(f.readline());
		mdl[i] = [None] * L;
		for j in range(L):
			pos = set(map(int, f.readline().split()));
			neg = set(map(int, f.readline().split()));
			codes = list(map(float, f.readline().split()));
			n = sum(codes) + 2*twist;
			mdl[i][j] = [pos, neg, (twist + codes[0]) / n, (twist + codes[1]) / n];
	return mdl;

def score(model, t):
	K = len(model);
	sc = 0;
	for i in range(K):
		m = model[i];
		for s in m:
			if len(t & s[0]) == len(s[0]) and len(t & s[1]) == 0:
				if i in t:	
					sc -= log(s[2]);
				else:
					sc -= log(s[3]);
			break;
	return sc;

def main():
	twist = 0.0;
	opts, args = getopt.getopt(sys.argv[1:], 'i:t:');
	data=[];
	for o, a in opts:
		if o == '-i':
			data = [set(map(int, line.split())) for line in open(a)];
		elif o == '-t':
			twist = float(a);
	models = [None] * len(args);
	for i in range(len(args)):
		models[i] = read_model(twist, args[i]);

	for t in data:
		sc = [None] * len(models);
		for i in range(len(models)):
			sc[i] = (score(models[i], t), i);
		print(min(sc)[1]);


	#mdl = read_model(0, 'paleo.mdl');
	#print score(mdl, set());

if __name__ == "__main__":
    main()

