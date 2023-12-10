## Scenario

We know the flag has been encoded into 5 integers, which form the secret vector s in the operation: $b = {(As + e)} \mod{p}$
Where p is prime, e is the errors vector of size n (50), A is a n x m (50 x 5) matrix and b is a vector of size n (50).
## Main Idea

To solve the challenge, we're gonna use the fact that we can rewrite such operation in two ways:
- $b = ([A|e][s,1]) \mod{p}$ 
(A concatenated with e as last vector, and s concatenated  with the scalar 1 to allow matrix multiplication).
- $b = ([A|I_{n}][s,e]) \mod{p}$
(A concatenated with an Identity matrix of size n x n, and s concatenated with the errors vector).

Please notice that when concatenating to A we're adding columns, while when we're concatenating to s we're adding rows.
## Solution

Starting from the second alternative representation, we can immediately see we've got a matrix which has more columns than rows, meaning we're going to face an infinite number of solutions. For the time being we will solve for $[s, e]$, fixing the free variables to 1. To achieve that we can use the classic Gaussian Reduction algorithm, keeping in mind we're working in modulo p (hence inverses are multiplicative inverses).

```python
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
```

Once we've found one potential solution X, we can split it into S = X[:5] and E = X[5:], which is one potential solution for s and e.
At this point, we need to realize that the error_bound given in the data.txt is several order of magnitudes smaller than the average value we find in E -> From this we can infere that the correct errors vector is probably the one with minimal length. How to achieve this? We need to go back and remember the first alternative representation of the formula. 

- $b = ([A|e][s,1]) \mod{p}$ 

That formula can be interpreted as following: "the result b is a linear combination with integer coefficients of the vectors forming A and the vector e (modulo p)". This sentence should ring a bell -> If we can build the <b>lattice</b> corresponding to this linear space, we can apply lattice reduction to minimize the basis and get the shortest possible vector for E.
(if you are not familiar with lattices, I suggest looking at the guided problems on cryptohack: https://cryptohack.org/challenges/maths/).<br>
Luckily for us, building such lattice is not too difficult, we just need the vectors in A, our potential solution E and to remember that in every dimension we can freely move by $kp$ thanks to the modulus. Hence, the set of vectors missing to complete the lattice basis is just the canonical basis element-wise multiplied by $p$.

```python
lattice = [A[i] + [E[i]] + [p if j == i else 0 for j in range(len(A))] for i in range(len(A))] # this also adds linearly depended vectors, but they won't cause harm
```

At this point, it's just a matter of running the LLL algorithm on our set of vectors, I chose to use the sagemath implementation because it was the most easily available for me.

```python
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
```

Once we have e, we can simply solve $A[0:5]s = (b[0:5] - e[0:5])$ for s reusing the Gaussian Reduction algorithm.

```python
Ab5 = [A[i][:] + [(b[i] - errors[i]) % p] for i in range(5)]
s = gaussian_solve_mod_p(Ab5, p)

for x in s:
	print(int.to_bytes(int(x), 50, "big")) # this will output the flag
```

## Full Code

