import numpy as np
import scipy.optimize as op

## CONSTANTS:
N = 10**6
C = 1   #Fix
c_dyn = 0.00129704
c_short = 0.0032426
c_leak = 0.002026625
lambda_mem = 36e-12     #Check units
lambda_net = lambda_mem
lambda_noc = .75e-12    #Check units
sigma_core = 0.0105e-6     #Check units
sigma_cache = 0.386e-6     #Check units
eps_node = 2

## CONSTRAINTS:
DIE_MAX = 141
POW_MAX = 173 #? Is this per processor?


## FUNCTION TO OPTIMIZE
'''
args: q, f, Z, beta_mem, beta_noc, beta_net, p
    q - Cores
    f - Frequency
    Z - Aggregate On-Chip Cache
    beta_mem - Mem Bandwidth
    beta_noc - On-chip Network Bandwidth
    beta_net - Network Bandwidth
    p - Number of nodes
'''
def matT(args):

    q, f, Z, beta_mem, beta_noc, beta_net, p = args

    T_comp = 2*(N**3)/(p*q*f)
    T_net = 2*(N**2)/(beta_net * (C*p)**0.5)
    T_mem = 2*(3**0.5)*(N**3) / (p*(Z**0.5)*beta_mem)
    T_noc = 2*(N**2)/(beta_noc*((C*p*q)**0.5))

    return max(T_comp,T_net,T_mem,T_noc)


## Power Constraint
def powCon(args):

    q, f, Z, beta_mem, beta_noc, beta_net, p = args

    P_comp = q*(c_dyn * f**3 + c_short * f**2 + c_leak * f)
    P_mem = beta_mem * lambda_mem
    P_net = beta_net * lambda_net * (3*p) #Check no of links
    P_noc = (4*q - 4*q**.5)* ((sigma_core + (Z/q)*sigma_cache)**.5 ) * beta_noc * lambda_noc

    print(P_comp,P_mem,P_net,P_noc,P_net)

    return POW_MAX*p - (p*(P_comp + P_mem + P_noc + eps_node) + P_net)


## Die Area Constraint
def dieCon(args):

    q, f, Z, beta_mem, beta_noc, beta_net, p = args

    return DIE_MAX - (q*sigma_core +Z*sigma_cache)


## Initializing array
arg0 = np.array([2048,      #q
                 10**6,     #f
                 2**20,     #Z
                 2**30,     #beta_mem
                 2**20,     #beta_noc
                 2**30,     #beta_net
                 10**6])    #p

print("Die Area: ",DIE_MAX - dieCon(arg0)," mm^2\nPower consumed: ",POW_MAX*arg0[6] - powCon(arg0))

print(matT(arg0))

#res = op.minimize(matT, arg0, method = 'COBYLA', constraints = [{"type": "ineq", "fun": powCon}, {"type": "ineq", "fun": dieCon}])
#print(res)

## C?

## T = max{T_comp, T_net, T_mem, T_noc}
## T_comp = 2(n^3)/(pqf)                        -> Vars: p,q,f
## T_net = 2(n^2)/(beta_net sqrt(Cp))           -> Vars: beta_net, p
## T_mem = 2root(3)(n^3)/(p beta_mem sqrt(Z))   -> Vars: Z, beta_mem, p
## T_noc = 2(n^2)/(beta_noc sqrt(Cpq))          -> Vars: beta_noc, p, q

## Pow = p(P_comp + P_mem + P_noc + eps_node) + P_net
## P_comp = q(c_dynamic f^3 + c_short f^2 + c_leakage f)    -> Vars: f
## P_mem = beta_mem lambda_mem                              -> Vars: beta_mem
## P_net = beta_net Links(p) lamda_net                      -> Vars: beta_net, Links(p)
## P_noc = (4q-4r(q)).....                                  -> Vars: q, Z, beta_noc
## Non-linear Constraint

## DieArea = f(q,Z)                                         -> q, Z
## Linear Constraint

## All p,q,f,... >= 0

## q, f, Z, beta_mem, beta_noc, beta_net, p