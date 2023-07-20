## Behavior

The program computes a faulty RSA encryption of the flag with only one prime number, <b>b</b>, which is provided to us with the ciphertext. By computing the totient (b-1) it's immediate to notice that its gcd with the public exponent (210) is the public exponent itself -> This means that recovering the plaintext is not as easy as computing pow(e, -1, b-1).

## Solution

To solve the challenge, we first must notice that the public exponent (210) is just 2 * 3 * 5 * 7. This tells us that the ciphertext is at the same time a quadratic, a cubic, a fifth power and seventh power residue. <br>To recover the roots modulo b, we must use a generalized version of Tonelli-Shank's algorithm that is able to find the rth root: such algorithm was presented by Adelman, Manders and Miller in 1977. <br>To get an idea of how the algorithm works please read this article by Scott Charles Lindhurst (Especially Chapter 3):<br> http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.144.2770 <br>
## Implementation

rth root finder for r prime: 

```
import math
from Crypto.Util.number import isPrime

def modular_roots(a: int, p: int, r: int) -> list[int] | None:
	
	if pow(a, (p-1)//r, p) != 1:
		print("n is not a rth-residue")
		return []
	
	if not isPrime(r):
		print("factorize r first")
		return
	
	q = p - 1
	n = 0
	while q % r == 0:
		n += 1
		q = q // r
	
	d = math.gcd(r, q)
	assert d == 1 # I'm implementing it only for prime r.
				  # If r is not prime just factorize it and run this more times
	
	m = 0
	for i in range(1, r): # for bigger r, this can be optimized with the euclidean algo
		if (i * q + 1) % r == 0:
			m = i
			break
	
	assert m != 0

	# find non-residue
	u = 2
	exp = (p-1) // r
	while True:
		z = pow(u, exp, p)
		if z != 1:
			break
			
		u += 1
	
	k = n
	z = pow(u, q, p)
	x = pow(a, (m*q + 1) // r, p)
	b = pow(a, m*q, p)
	
	# end of initialization-phase
	
	while b != 1:
		m = 0
		must_be_one = b
		for i in range(1, k):
			must_be_one = pow(must_be_one, r, p)
			if must_be_one == 1:
				m = i
				break
		
		assert m != 0
		
		t = pow(z, pow(r, k - m - 1, p-1), p)
		z = pow(t, r, p)
		
		# finding l such that bz^l is in S(m-1)
		l = 0
		ord = pow(r, m-1, p-1)
		bz = b
		
		for i in range(1, r):
			bz = (bz * z) % p 
			if pow(bz, ord, p) == 1:
				#then it's in the sub-group since the order is r^(m-1)
				l = i
				break
		
		assert l != 0
		
		b = (b * pow(z, l, p)) % p
		x = (x * pow(t, l, p)) % p
		k = m

	# when we arrive here we've got one solution (x), time to find the other r-1:
	g = p - 2
	while pow(g,(p-1) // r, p) == 1:
		g -= 1
		
	g = pow(g, (p-1) // r, p)
	return sorted([(x * pow(g, i, p) % p) for i in range(r)])
	
```

Decryption of given ciphertext:

```

b = 135954424160848276393136392848608760791498666756786983317146989739232222268153235587604168914827859099133726281621143020610041450200631778336472889038077986687446107427527703447531968569919642975653169056203851297117178187249653136191818357235077367060617558261023389453028554177668515375377299577050000000001

f210 = 13350651294793393689107684390908420972977381011191079381202728507002264420264784588373703945341668404762890725356808809021906408198983625375190500172144348596288910240548668158058030780501343680214713780242304547715977777103636873360269427453504233184515002477489763359569764117968027273137245802436961373256


f105 = modular_roots(f210, b, 2)
for root2 in f105:
	f35 = modular_roots(root2, b, 3)
	
	if len(f35) > 0:
		for root3 in f35:
			f7 = modular_roots(root3, b, 5)
			
			if len(f7) > 0:
				for root5 in f7:
					f = modular_roots(root5, b, 7)
					
					if len(f) > 0:
						for inverse_flag in f:
							
							# remember that the flag was inverted
							pflag = pow(inverse_flag, -1, b) 
							msg = pflag.to_bytes(150, "big")
							
							if b"amateursCTF" in msg:
								print(msg)
								print("")
```
