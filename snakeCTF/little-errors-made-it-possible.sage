from data import *

def norm2(v):
    return sum([x**2 for x in v])

def gaussian_solve_mod_p(Ab, p):
    ab2 = [a[:] for a in Ab] # keeping the original matrix intact
    m = len(ab2[0]) - 1
    n = len(ab2)
    for y in range(min(n, m)):
        t = ab2[y][y]
        denom = pow(t, -1, p)
        for y2 in range(y + 1, min(n,m)):
            num = ab2[y2][y]
            for x in range(y, m + 1):                         
                ab2[y2][x] = (ab2[y2][x] - num * denom * ab2[y][x]) % p

    solved = []
    if n < m:
        for i in range(m-n):
            solved.append(1)
    for y in range(min(n, m) - 1,-1,-1):
        b = ab2[y][m]
        for x, s in enumerate(solved):
            b = (b - ab2[y][m-1-x] * s) % p

        s = b * pow(ab2[y][y], -1, p) % p
        solved.append(s)

    solved.reverse()

    return solved

Ab = [A[i] + [1 if i == j else 0 for j in range(len(A))] + [b[i]] for i in range(len(A))]
X = gaussian_solve_mod_p(Ab, p)

E = X[5:] # the first 5 are potential (wrong) solutions to s

lattice = [A[i] + [E[i]] + [p if j == i else 0 for j in range(len(A))] for i in range(len(A))] # this also adds linearly depended vectors, but they won't cause harm

M = Matrix(ZZ, lattice)
reduced = M.transpose().LLL()

min_norm2 = None
min_index = None
for i, v in enumerate(reduced):
	v_length = norm2(v)
	if v_length > 0: # only non-zero vectors allowed
		if min_norm2 is None or min_norm2 > v_length:
			min_norm2 = v_length
			min_index = i

e = reduced[min_index]

Ab5 = [A[i][:] + [(b[i] - e[i]) % p] for i in range(5)]
s = gaussian_solve_mod_p(Ab5, p)

for x in s:
	print(int.to_bytes(int(x), 50, "big")) # this will output the flag
