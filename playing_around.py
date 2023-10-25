from typing import Callable
from math import sqrt, sin, cos, pi


LucasFunction = Callable[[int], float]

def make_U(P: int, Q: int) -> LucasFunction:
    a = (P + sqrt(P**2-4*Q))/2
    b = (P - sqrt(P**2-4*Q))/2
    def U(n: int):
        if n <= 1: 
            return (a**n - b**n)/(a - b)
        return P*U(n-1) - Q*U(n-2)
    return U

def make_V(P: int, Q: int) -> LucasFunction:
    a = (P + sqrt(P**2-4*Q))/2
    b = (P - sqrt(P**2-4*Q))/2
    def V(n: int) -> float:
        if n <= 1: 
            return (a**n + b**n)
        return P*V(n-1) - Q*V(n-2)
    return V


def make_Ubar(P: int, Q: int) -> LucasFunction:
    a = P + sqrt(P**2-4*Q)
    b = P - sqrt(P**2-4*Q)
    def Ubar(n: int):
        if n <= 1: 
            return (a**n - b**n)/(a - b)
        return ( (1 + (P**2 - 1)*1j**(n-1) * sin(pi*n/2))*Ubar(n-1) - Q*Ubar(n-2) ).real
    return Ubar

def make_Vbar(P: int, Q: int) -> LucasFunction:
    a = P + sqrt(P**2-4*Q)
    b = P - sqrt(P**2-4*Q)
    def Vbar(n: int):
        if n <= 1: 
            return (a**n - b**n)/(a - b)
        return ( (1 + (P**2 - 1)*1j**n * cos(pi*n/2))*Vbar(n-1) - Q*Vbar(n-2) ).real
    return Vbar

def make_Lucas(P: int, Q: int) -> tuple[LucasFunction, LucasFunction, LucasFunction, LucasFunction]:
    return make_U(P,Q), make_V(P,Q), make_Ubar(P,Q),  make_Vbar(P,Q)

Matriz = list[list[int]]
def matmul(A: Matriz, B: Matriz) -> Matriz:
    return [[sum(A[i][k]*B[k][j]  for k in range(len(B))) 
             for j in range(len(B[0]))]
             for i in range(len(A))]
def matpow(A: Matriz, n: int) -> Matriz:
    if n == 0:
        return [[int(i == j) for i in range(len(A[0]))]
                             for j in range(len(A))]
    if n == 1:
        return A
    if n%2 == 1:
        return matmul(A,matpow(A,n-1))
    aux = matpow(A,n//2)
    return matmul(aux,aux)

def V(P: int, Q: int, n: int) -> int:
    return matmul(
        [[2, P]],
        matpow([[0,-Q]
               ,[1, P] ], n),
    )[0][0]
def U(P: int, Q: int, n: int) -> int:
    return matmul(
        [[0, 1]],
        matpow([[0,-Q]
               ,[1, P] ], n),
    )[0][0]

def runTests():
    import numpy as np
    from tqdm import tqdm
    from random import randint
    def random_mat(n):
        return [[randint(-10,10) for _ in range(n)] for _ in range(n)]
    for _ in tqdm(range(10), desc="Mul"):
        A, B = random_mat(100), random_mat(100)
        assert np.all( np.array(matmul(A,B)) == np.array(A)@np.array(B) )
    for _ in tqdm(range(10), desc="Pow"):
        A, n = random_mat(100), randint(1,10)
        assert np.all( (act := np.linalg.matrix_power(np.array(A), n)) == (exp := np.array(matpow(A,n)))), f"Failed for {n}, dif matrix = \n{act-exp}"

    lucas_numbers = [ 2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322, 521, 843, 1364, 2207, 3571, 5778, 9349, 15127, 24476, 39603, 64079, 103682, 167761, 271443, 439204, 710647, 1149851, 1860498, 3010349, 4870847, 7881196, 12752043, 20633239, 33385282, 54018521, 87403803 ]
    assert all(ln == V(1,-1,i) for i, ln in enumerate(lucas_numbers))
    for _ in range(10):
        P,Q = randint(-10,10), randint(-10,10)
        F, G = make_V(1,-1), lambda n : V(1,-1,n)
        assert all(int(F(x)) == G(x) for x in tqdm(range(20), desc=f"Checking for {P=} and {Q=}"))

F, G = make_V(1,-1), lambda n : V(1,-1,n)
