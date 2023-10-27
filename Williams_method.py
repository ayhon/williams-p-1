#!/usr/bin/env python
# coding: utf-8

# # Método $p+1$ de Williams

# # Contenidos
# ## Introducción: Factorización de números
# ## Método $p-1$ de Pollard
# ## Funciones de Lucas
# ## Método $p+1$ de William
# ## Adendum: Segundos pasos

# # Introducción: Factorización de números

# ### ¿Cómo podemos factorizar un número?

# In[81]:


def factor(N: int) -> dict[int,int]:
    factors = {}
    for number in range(2,sqrt(N)):
        multiplicity = 0
        while N % number == 0:
            multiplicity+=1
            N = N // number
        if multiplicity != 0:
            factors.update({number: multiplicity})
    return factors


# Complejidad: $O\left(\sqrt{N}\right)$

# In[82]:


def factor_once(N: int) -> int:
    for number in range(2,sqrt(N)):
        multiplicity = 0
        while N % number == 0:
            multiplicity+=1
            N = N // number
        if multiplicity != 0:
            return number**multiplicity
    return 1


# Si devolvemos `1` es que nos hemos rendido

# La complejidad sigue siendo $O\left(\sqrt{N}\right)$ en el peor caso

# ### ¿Hay una forma más rápida de sacar un factor de un número?

# In[14]:


def gcd(N: int, M: int) -> int:
    return N if M == 0 else gcd(M,N%M)


# Complejidad: $O(\log(N))$

# En verdad $O\left(\log^2{N}\right)$ si $N$ es de longitud arbitraria

# Esto es porque `rem` toma $O\left(\log^2(N)\right)$ en estos casos ([Source](https://en.wikipedia.org/wiki/Euclidean_algorithm#Algorithmic_efficiency))

# ### ¿Cómo conseguimos multiplos de un factor de N?

#  - Método $p-1$ de Pollard

#  - Método $p+1$ de William

# # Método $p-1$ de Pollard

# Teorema de Fermat:  
# Para $p$ primo y $a$ coprimo con $p$,
# $$
# a^{p-1} \equiv 1 \mod p
# $$

# $$
# a^{k\left(p-1\right)} \equiv 1 \mod p
# $$
# para $k\in \mathbb{Z}$

# Es decir
# $$a^{k(p-1)} - 1 \in p\mathbb{Z}$$

# $$
# \begin{align*}
# \varphi :\ &(p-1)\mathbb{Z} \to p\mathbb{Z} \\
#           & k(p-1) \mapsto a^{k(p-1)} - 1
# \end{align*}
# $$

# Lo que tenemos es una función que nos transforma múltiplos de $p-1$ en múltiplos de $p$

# ## ¿Como hallamos múltiplos de $p-1$?

# Supongamos que $p-1$ es $B$-suave

# > Un número $n$ se dice $B$-suave si todos sus factores potencia de primo son menores o iguales a $B$, es decir
# > $$ n = p_1^{\alpha_1} \cdot\dots\cdot p_r^{\alpha_r} \quad \Longrightarrow \quad p_i^{\alpha_i} \le B \quad \forall i \in \{1,\dots,r\} $$

# Múltiplos de $p-1$:

#  - $B! = \prod_{k=1}^Bk$  
#   si $\alpha_i = 1$ para $i\in\{1,\dots,r\}$

#  - $R = \prod_{p\in P} p^{\lfloor\log_pB\rfloor}$  
#   con $P = \{p \text{ primo } : p \le B\}$

# In[16]:


from math import gcd, log
from typing import Iterable, Callable
from itertools import count
from random import randint


# In[17]:


def metodo_pollard(N: int, B: int) -> int:
    a_R = randint(2,N)
    for prime in primes_until(B):
        r = pow(prime, int(log(B,prime)))
        a_R = pow(a_R,r,N)
        d = gcd(a_R - 1, N)
        if d not in (1,N):
            return d
    return 1


# Lo verdaderamente increible es que nunca necesitamos usar $R$, nos basta con construir $a^R$ de forma paulatina

# In[18]:


def primes_until(B: int) -> Iterable[int]:
    sequence = count(2)
    prime = next(sequence)
    while prime <= B:
        yield prime
        sequence = filter(lambda x, p=prime: x % p != 0, sequence)
        prime = next(sequence)


# In[19]:


import sympy
print(sympy.factorint(226799-1)) # Este puede servir para un ejemplo
assert sympy.isprime(226799)
print(sympy.factorint(347161-1))
assert sympy.isprime(347161)


# In[20]:


metodo_pollard(226799 * 347161, B=200)


# In[21]:


# metodo_pollard(, B=100) != 1
powersmooth_primes = [
    104723,
    226799,
    347161,
    799979,
    1952437,
    2501099,
]
Ns = [prime*70201 for prime in powersmooth_primes]
all(prime for N in Ns if metodo_pollard(N, B=10) != 1)


# In[34]:


def metodo_generico (
    N: int, 
    B: int, 
    initial_candidate: int,
    update_candidate: Callable[[int,int], int],
    extract: Callable[[int,int],int] = gcd,
) -> int:
    candidate = initial_candidate
    for prime in primes_until(B):
        r = pow(prime, int(log(B,prime)))
        candidate = update_candidate(candidate,r)
        d = extract(candidate, N)
        if d not in (1,N):
            return d
    return 1

def metodo_pollard(N: int, B: int) -> int:
    update = lambda a_R, r: pow(a_R,r,N)
    return metodo_generico ( N
                           , B
                           , initial_candidate = randint(2,N-1)
                           , update_candidate  = update
                           , extract = lambda a, b : gcd(a-1,b))


# # Funciones de Lucas

# Definimos las funciones de Lucas como
# $$
# U_n(P,Q) = \frac{a^n - b^n}{a-b} \quad\quad\quad 
# V_n(P,Q) = a^n + b^n
# $$
# donde $a$ y $b$ son las raices de $x^2 - Px + Q$

# Por lo tanto, tenemos que
# $$
# \begin{align*}
# P &= a + b \\
# Q &= ab
# \end{align*}
# $$

# Denotamos por $D$ al discriminante de $x^2 - Px + Q$
# $$
# \begin{align*}
# D &= P^2 - 4Q \\
#   &= (a - b)^2
# \end{align*}
# $$

# Las funciones de Lucas cumplen muchas propiedades

#  - $U_n = PU_{n-1} - QU_{n-2}$
#  - $V_n = PV_{n-1} - QV_{n-2}$

#  <video controls>
#   <source src="Move_a_formula.mp4" type="video/mp4">
#   Your browser does not support the video tag.
# </video> 

#  - $U_0 = \frac{1 - 1}{a-b} = 0$
#  - $U_1 = \frac{a-b}{a-b} = 1$
#  - $V_0 = 1 + 1 = 2$
#  - $V_1 = a + b = P$

# Con $P = -Q = 1$, $\{U_n\}_{n\ge1}$ es la serie de Fibonacci y $\{V_n\}_{n\ge1}$ es la serie de Lucas

# Las que nos interesan son:

# $$
# U_{k\left(p-\left(\frac{D}{p}\right)\right)}(P,Q) \equiv 0 \mod p
# $$

# $$
# V_{k\left(p-\left(\frac{D}{p}\right)\right)}(P,Q) \equiv 2Q^{\frac{1}{2}\left(1-\left(\frac{D}{p}\right)\right)} \mod p
# $$

# In[2]:


import sympy as sp
i,j,a,b = sp.symbols("i,j,a,b", reals=True, positive=True)
P = a + b
Q = a*b
D = (a-b)**2
U = lambda n : (a**n-b**n)/(a-b)
V = lambda n : a**n + b**n


# ## $$U_{p-\left(\frac{D}{p}\right)}(P,Q) \equiv 0 \mod p$$

# $$V_p \equiv V_1 = P \mod p$$

# $$U_p \equiv \left(\frac{D}{p}\right) \mod p$$

# $$2U_{i+j} = U_iV_j + U_jV_i$$

# In[4]:


sp.expand(U(i) * V(j) + U(j) * V(i)).simplify()


# In[5]:


2*U(i+j) 


# $$U_{p+1} \equiv 0\mod p\quad\text{ si }\left(\frac{D}{p}\right) = -1$$

# $$2Q^jU_{i-j} = U_iV_j - U_jV_i$$

# In[10]:


sp.expand(2*Q**j * U(i-j)).simplify()


# In[12]:


sp.expand(U(i) * V(j) - U(j) * V(i)).simplify()


# $$U_{p-1} \equiv 0\mod p\quad\text{ si }\left(\frac{D}{p}\right) = 1$$

# ## $$U_{k\left(p-\left(\frac{D}{p}\right)\right)}(P,Q) \equiv 0 \mod p$$

# $$ p \mid U_i \text{ y } p \mid U_j \Longrightarrow p \mid U_{i+j} $$

# ## $$ V_{p-\left(\frac{D}{p}\right)}(P,Q) \equiv 2Q^{\frac{1}{2}\left(1-\left(\frac{D}{p}\right)\right)} \mod p $$

# $$V_p \equiv V_1 = P \mod p$$
# $$U_p \equiv \left(\frac{D}{p}\right) \mod p$$

# $$2V_{i+j} = V_iV_j + DU_iU_j$$

# In[24]:


sp.expand(V(i)*V(j) + D*U(i)*U(j)).simplify()


# In[25]:


2*V(i+j)


# $$V_{p+1} \equiv 2Q \mod p\quad\text{ si } \left(\frac{D}{p}\right) = -1$$

# $$2Q^jV_{i-j} \equiv V_iV_j - DU_iU_j$$

# In[26]:


sp.expand(V(i)*V(j) - D*U(i)*U(j)).simplify()


# In[27]:


sp.expand(2*Q**j*V(i-j))


# $$V_{p-1} \equiv 2 \mod p\quad\text{ si }\left(\frac{D}{p}\right) = 1$$

# ## $$ V_{k\left(p-\left(\frac{D}{p}\right)\right)}(P,Q) \equiv 2Q^{\frac{1}{2}\left(1-\left(\frac{D}{p}\right)\right)} \mod p $$

# $$2U_{i+j} = U_iV_j + U_j V_i$$
# $$2V_{i+j} = V_iV_j + DU_iU_j$$

# $$2V_{i+2j} = V_iV_{2j} + DU_iU_jV_j$$

# $$p\mid V_i \text{ y } p \mid V_j \Longrightarrow p \mid V_{i+2j}$$

# $$V_{2i} = V_i^2 - 2Q^{i}$$

# In[29]:


sp.expand(V(i)**2 -2*Q**i).simplify()


# In[30]:


V(2*i)


# $$p \mid V_i \text{ y } p \mid 2Q^i \Longrightarrow p \mid V_{2i}$$

# # Método $p+1$ de William

# In[149]:


from typing import Callable


# In[31]:


def metodo_generico (
    N: int, 
    B: int, 
    initial_candidate: int,
    update_candidate: Callable[[int,int], int],
    extract: Callable[[int,int],int] = gcd,
) -> int:
    candidate = initial_candidate
    for prime in primes_until(B):
        r = pow(prime, int(log(B,prime)))
        candidate = update_candidate(candidate,r)
        d = extract(candidate, N)
        if d not in (1,N):
            return d
    return 1

def metodo_pollard(N: int, B: int) -> int:
    update = lambda a_R, r: pow(a_R,r,N)
    return metodo_generico ( N
                           , B
                           , initial_candidate = randint(2,N-1)
                           , update_candidate  = update
                           , extract = lambda a, b : gcd(a-1,b))


# ### ¿Cómo conseguimos algo parecido, aprovechándonos de las funciones de Lucas?

# _A partir de ahora asumimos que $Q = 1$. En el paper original, William prueba que esto se puede hacer sin pérdida de generalidad pues_

# $$U_{2m}(a+b,ab) \equiv PQ^{m-1}U_m\left(\frac{a}{b}+\frac{b}{a},1\right)\mod p$$ 

# _por lo que_

# $$p \mid U_m(P,Q) \Longrightarrow p \mid U_{2m}(P,Q) \Longrightarrow U_m(P',1)$$

# $$
# \begin{align*}
# U_{2n - 1} &= U_n^2 - QU_{n-1}^2 \\
# U_{2n} &= U_n\left(P U_n - 2QU_{n-1}\right) \\
# U_{2n + 1} &= PU_{2n} - QU_{2n-1} \\
# \end{align*}
# $$

# Funcionaría, pero tiene un problema. Necesitamos tener $R$ construido totalmente para hallar $U_R$

# $$V_n(V_k(P,Q),Q^k) = V_{nk}(P,Q)$$

# Si tenemos $R = r_1 \dots r_k$ y queremos construir $V_R(P) = V_R(P,1)$ entonces podemos considerar $v_0 = V_1 = P$ y $v_i = V_{r_1 \dots r_i}(P) = V_{r_i}(v_{i-1})$, donde $v_k = R$

# $$
# \begin{pmatrix}
# V_i \\ V_{i+1}
# \end{pmatrix}
# =
# \begin{pmatrix}
# 0 & 1 \\
# -Q & P
# \end{pmatrix}
# \begin{pmatrix}
# V_{i-1}\\ V_i
# \end{pmatrix}
#  = 
# \begin{pmatrix}
# 0 & 1 \\
# -Q & P
# \end{pmatrix} ^2
# \begin{pmatrix}
# V_{i-2}\\ V_{i-1}
# \end{pmatrix}
# =
# \cdots
# =
# \begin{pmatrix}
# 0 & 1 \\
# -Q & P
# \end{pmatrix} ^ i
# \begin{pmatrix}
# V_{0}\\ V_{1}
# \end{pmatrix}
# $$

# In[40]:


Matriz = list[list[int]]
def matmul(A: Matriz, B: Matriz,M:int) -> Matriz:
    return [[sum(A[i][k]*B[k][j]%M  for k in range(len(B)))%M
             for j in range(len(B[0]))]
             for i in range(len(A))]
def matpow(A: Matriz, n: int,M) -> Matriz:
    if n == 0:
        return [[int(i == j) for i in range(len(A[0]))]
                             for j in range(len(A))]
    if n == 1:
        return A
    if n%2 == 1:
        return matmul(A,matpow(A,n-1,M),M)
    aux = matpow(A,n//2,M)
    return matmul(aux,aux,M)

def V(P: int, Q: int, n: int, N: int) -> int:
    return matmul(
        [[2, P]],
        matpow([[0,-Q]
               ,[1, P] ], n, N),
        N
    )[0][0]


# In[41]:


def metodo_william(N: int, B: int) -> int:
    P_0 = randint(2,N)
    if (d := gcd(P_0**2-4,N)) == 1: return d
    update = lambda v_i,r : V(v_i,1,r,N)
    return metodo_generico( N
                          , B
                          , initial_candidate = P_0
                          , update_candidate  = update
                          , extract = lambda a, b : gcd(a-2,b))


# In[80]:


candidates = []
for prime in primes_until(1_000):
    factors = sp.factorint(prime+1)
    B = max(p**a for p, a in factors.items())
    candidates += [(prime,B)]
sorted(candidates, key=lambda c : c[1])


# In[156]:


next(d for _ in range(1000) if (d := metodo_william(769**10, 20)) != 1)


# In[157]:


next(d for _ in range(1000) if (d := metodo_william(839**10,8)) != 1)


# In[72]:


print(metodo_william(112729,23))
print(sympy.factorint(112729-1))
print(metodo_william(451889,10))
print(sympy.factorint(451889-1))
print(next(d for _ in range(100000) if (d := metodo_william(613**2,333)) != 1))


# In[1]:


import sympy as sp
i,j,a,b = sp.symbols("i,j,a,b", reals=True, positive=True)
P = a + b
Q = a*b
D = (a-b)**2
U = lambda n : (a**n-b**n)/(a-b)
V = lambda n : a**n + b**n


# # Fin de la presentación
# ### ¿Alguna pregunta?
