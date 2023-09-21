from math import log, exp;

def cumsum(i):
	s = 0;
	for x in i:
		s+=x;
		yield s;

def logfrac(k):
	return list(cumsum([0, 0] + list(map(log, list(range(2, k + 1))))));


def binomial_mdl(k):
	fr = logfrac(k);

	t = [fr[k] - fr[x] - fr[k - x] + x*(log(x) - log(k)) + (k - x)*(log(k - x) - log(k)) for x in range(1, k)];
	return log(sum(map(exp, t + [0, 0])));
