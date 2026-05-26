from CDNozzleLib import *
import __pycache__ as self

# Air Ideal Gas Thermodynamic Constants:
GAMMA = 1.4  # Air Specific Heat
R = 287.0    # Air Gas Constant

# Nozzle  Geometric Parameters
Lc = 0.2     # Converging Lenght (m)
Ld = 0.7     # Diverginf Length (m)
Ri = 0.15    # Inlet Radius (m)
Rth = 0.05   # Throat Radius (m)
Re = 0.1     # Outlet Radius (m)

# Nozzle Condition Paramteres
P0 = 1.0     # Inlet Total Pressure

# Loading Parameters into Solver:
Constants.__init__(GAMMA, R, Lc, Ld, Ri, Rth, Re, P0)

print("="*60)
print("                 C-D NOZZLE FLOW ANALYSIS")
print("="*60)
L = Lc+Ld-0.001
# Analysis Cases:
cases = [("Throat Mach no. = 0.7", None, None, 0.7),
         ("Choked Subsonic", None, None, 1),
         ("Shock at x = 0.5m", None, 0.5, None),
         ("Fully Supersonic", 0.029787, None, None),
         ("Shock at Exhaust", None, L, None)]

# Solving for each Case
for case_name, Pe, shock_pos, M_throat in cases:
    print(f"\nANALYSING CASE: {case_name}:")
    print("-" * 60)
    results, flow_type, shock_x, P_total = FlowSolver.solve_flow(
        Pe, shock_position=shock_pos, M_throat=M_throat)
    PlotFunc.plot_nozzle_results(results, flow_type, shock_x, case_name)
    plt.savefig(f'./Plots/{case_name}.jpg')
    plt.show()
    print('Plots saved successfully.')
    print('\n')

Pe_start = 0.9999
Pe_mid = 0.985
Pe_end = 0


print("="*60)
print("CREATING ANIMATION: Back Pressure Effect")
print("="*60)
start = Pe_start
end = Pe_mid
n = 50
a = []
for i in np.linspace(1, n, n):
    a.append(start + ((end-start)/(n-1))*(i-1))

start = Pe_mid
end = Pe_end
n = 100
for i in np.linspace(1, n-1, n):
    a.append(start + ((end-start)/(n-1))*(i))

Pe_range = a
PlotFunc.create_flow_animation(Pe_range, 'nozzle_back_pressure_animation')
