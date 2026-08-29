import math                                                                                    
# --- YOUR MATERIAL PROPERTIES --- 
# Copy from properties.inp                                                 
K = 1000.0
C10 = 0.0
C01 = 0.0
PHINET = 1.0
a = 1.2
R0C = 0.014
ETAC = 0.5
MU0STRETCH = 38600000000000.0
BETA = 0.5
LP = 16.0
THETA = 25.0 + 273.0
DX = 0.1
BB = 0.000001
LAMBDA0 = 1.0
CACTIN = 0.0095
R = 0.1
RFMAX = 1.0
RBMAX = 0.25
CHI = 0.1
D = 0.1
MU0 = 0.0
VMOL = 0.15
KOFF0 = 0.05
KEQ = 0.25

# --- CAPACITIES ---
CR = CACTIN * R
F_MAX = CACTIN * RFMAX
C_MAX = CACTIN * RBMAX

RGAS = 8.31446261815324  # J/(mol*K)

# --- TRANSCENDENTAL EQUATION FOR cb0 ---
def evalh(cb0):
    cf0 = CR - cb0
    thetaf0 = cf0 / F_MAX
    # lhs - rhs = 0
    lhs = KEQ * cf0 * (1.0 - cb0 / C_MAX)
    rhs = cb0 * (1.0 - thetaf0) * math.exp(-CHI * (1.0 - 2.0 * thetaf0))            
    return lhs - rhs

# Bisection solver
def bisection(f, a, b, tol=1e-12):
    if f(a)*f(b) >= 0:
        return a
    while (b - a)/2.0 > tol:
        midpoint = (a + b)/2.0
        if f(midpoint) == 0:
            return midpoint
        elif f(a)*f(midpoint) < 0:
            b = midpoint
        else:
            a = midpoint
    return (a + b)/2.0

cb0_root = bisection(evalh, 0.0, min(CR, C_MAX))

# --- CALCULATE EXACT INITMU ---
cf0 = CR - cb0_root
thetaf0 = cf0 / F_MAX
JC = 1.0 + VMOL * CR
JE = 1.0 / JC

INITMU = MU0 + RGAS * THETA * (
    math.log(thetaf0 / (1.0 - thetaf0)) + 
    CHI * (1.0 - 2.0 * thetaf0) - 
    (K * VMOL / (RGAS * THETA)) * (math.log(JE) / JC)
)

print(f"THETAF_T for Abaqus = {thetaf0:.6f}")
print(f"INITMU for Abaqus = {INITMU:.6f}")
