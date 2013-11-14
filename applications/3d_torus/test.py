#!/usr/bin/python
import math

R_outer=1.0
G = 1.0
M = 1.0
q = 0.499

def R_0(R_inner):
    A = R_inner**(2.0*(q-1.0)) - R_outer**(2.0*(q-1.0)) 
    B = -2.0*(q-1.0)*(1.0/R_inner - 1.0/R_outer)
    return (A/B)**(1.0/(2.0*q-1.0))
    

def j_0(R_inner):
    return math.sqrt(G*M*R_0(R_inner))

def Psi(R_inner,R):
    return -( (j_0(R_inner)/R_0(R_inner))**2.0) * ((R/R_0(R_inner))**(2.0*(q-1.0)))  /  (2.0*(q-1.0))

def C(R_inner):
#    return G*M*(( R_outer**(2.0*(q-1.0))/(R_0(R_inner))**(2.0*q-1.0) )/(2.0*(q-1.0)) - 1/R_outer)
    return Psi(R_inner,R_outer)-G*M/R_outer
   
def j_here(R_inner,R):
    return j_0(R_inner)*(R/R_0(R_inner))**q

def Psi_q0(R):
    return 0.5*(j_0(R_inner)/R)**2.0

def H(R_inner,R,Z):
    return G*M/math.sqrt(R**2.0+Z**2.0) - Psi(R_inner,R) + C(R_inner)

def z_max(R_inner,R):
#    temp1 = (G*M/(Psi(R_inner,R)+C(R_inner)))**2.0 
    A = G*M/(C(R_inner)-Psi(R_inner,R))
    B = A**2.0 - R**2.0
    return math.sqrt(B)
    

R_inner = 0.01



N = 0
start = R_inner
stop  = R_outer
while (N < 100):
    N=N+1
    R = start + (stop-start)*N/101.0
    
    print R, z_max(R_inner,R)
#    startz = 0.0
#    stopz = 0.1
#    NZ = 0
#    while (NZ < 100):
#        NZ = NZ + 1
#        Z = startz + (stopz-startz)*NZ/101.0
#        print R, Z, H(R_inner,R,Z)
    
        
        
#print "R_0 = ",R_0(R_inner)
#print "C = ",C(R_inner)
