import math                                                                                    
# --- YOUR MATERIAL PROPERTIES --- 
# Copy from properties.inp                                                 
density = 0.0005
K = 1000
C10 = 0.1
C01 = 0
PHINET = 0.9
a = 1.2
R0C = 0.014
ETAC = 0.333333333333
MU0STRETCH = 38600.0
BETA = 0.5
LP = 16
THETA = 298
DX = 0.005
BB = 1e-06
LAMBDA0 = 1
CACTIN = 0.0095
R = 0.1
RFMAX = 0.25
RBMAX = 0.25
CHI = 0.1
D = 0.01
MU0 = 0.0
VMOL = 0.15
KOFF0 = 0.05
KEQ = 0.25
KCATCH0 = 0.05
DXC = 0.005
# ** GAMMA = 0.2
UIND = -3.1


# Baseline
# koff0 = 0.05
# Keq = 0.25
# --- CAPACITIES ---
CR = CACTIN * R
F_MAX = CACTIN * RFMAX
C_MAX = CACTIN * RBMAX

RGAS = 8.31446261815324  # J/(mol*K)
KB = 1.380649e-5  # pN um / K

# --- REST FORCE OF THE PRESTRETCHED FILAMENT (lambda_bar = 1) ---
def rest_force(cb0):
    """f0 at F = I (as evalh.f90 / material.py::_init_equilibrium); 0 for LAMBDA0 = 1."""
    if cb0 <= 0.0:
        return 0.0
    r0f = 1.6 * (1.0e3 * cb0)**(-2.0 / 5.0)
    r0 = r0f + R0C
    L = a * r0f
    l0f = ETAC * (r0 / r0f) * (LAMBDA0 - 1.0) + 1.0 if 0.0 < ETAC <= 1.0 else LAMBDA0
    if l0f <= 1.0:
        return 0.0
    B0 = LP * THETA * KB
    def G(f):   # stretch-force relation of the extensible filament (evalg.f90), lambdaf = 1
        aux = f / MU0STRETCH
        rhs = 1.0 + aux - (1.0 - r0f / L) * (1.0 + aux)**BETA * (1.0 + 2.0 * aux) \
              / (1.0 + f * L**2 / (math.pi**2 * B0) * (1.0 + aux))**BETA
        return l0f * r0f / L - rhs
    return bisection(G, 0.0, MU0STRETCH * (1.0 + l0f * r0f / L))

def koff(f):
    """Catch-slip off-rate (koffcs.f90)."""
    return KCATCH0 * math.exp(-DXC / (KB * THETA) * f) + KOFF0 * math.exp(DX / (KB * THETA) * f)

# --- EQUATION FOR cb0 ---
def evalh(cb0):
    cf0 = CR - cb0
    thetaf0 = cf0 / F_MAX
    # prestress-aware: Keq = kon0/koff(0) becomes Keq*koff(0)/koff(f0)
    keq_eff = KEQ * koff(0.0) / koff(rest_force(cb0))
    # lhs - rhs = 0
    lhs = keq_eff * cf0 * (1.0 - cb0 / C_MAX)
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
JC = 1.0    # J^c = 1 + VMOL*(c - c0) is measured from the initial content c0 = CR (stress-free reference)
JE = 1.0 / JC

INITMU = MU0 + RGAS * THETA * (
    math.log(thetaf0 / (1.0 - thetaf0)) + 
    CHI * (1.0 - 2.0 * thetaf0) - 
    (K * VMOL / (RGAS * THETA)) * (math.log(JE) / JC)
)

print(f"cb0 for Abaqus = {cb0_root:.6f}")
print(f"THETAF_T for Abaqus = {thetaf0:.6f}")
print(f"INITMU for Abaqus = {INITMU:.6f}")
