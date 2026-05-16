import numpy as np
import matplotlib.pyplot as plt
import __pycache__ as self
from scipy.optimize import fsolve
import matplotlib.animation as animation
import os
import math


# =============================================================================
#    1. DEFINING THERMODYNAMIC CONSTANTS AND NOZZLE PARAMETERS
# =============================================================================
class Constants:
    def __init__(GAMMA, R, Lc, Ld, Ri, Rth, Re, P0, self=self):
        self.GAMMA = GAMMA  # Air Specific Heat
        self.R = R    # Air Gas Constant

        # Nozzle  Geometric Parameters
        self.Lc = Lc    # Converging Lenght (m)
        self.Ld = Ld    # Diverginf Length (m)
        self.Ri = Ri   # Inlet Radius (m)
        self.Rth = Rth  # Throat Radius (m)
        self.Re = Re    # Outlet Radius (m)

        # Nozzle Condition Paramteres
        self.P0 = P0    # Inlet Total Pressure

        # Areas
        self.A_inlet = NozzleGeometery.nozzle_area(0)
        self.A_throat = NozzleGeometery.nozzle_area(self.Lc)
        self.A_exit = NozzleGeometery.nozzle_area(self.Lc + self.Ld)
        pass


# =============================================================================
#    2. NOZZLE GEOMETRIC EQUATIONS
# =============================================================================
class NozzleGeometery:
    def nozzle_radius(x, nozzletype='parabola-cosine'):
        # Nozzle Profile:
        if x <= self.Lc:  # Converging Section
            return (self.Ri - self.Rth)/self.Lc**2 * (x - self.Lc)**2 + self.Rth
        else:  # Diverging Section
            return (self.Re - self.Rth)/2.0 * np.cos(np.pi/self.Ld * (x - (self.Lc - self.Ld))) + (self.Re + self.Rth)/2.0

    def nozzle_area(x, nozzletype='parabola-cosine'):
        """Nozzle Area at x"""
        R_x = NozzleGeometery.nozzle_radius(x, nozzletype=nozzletype)
        return np.pi * R_x**2


# =============================================================================
#    3. BASIC THERMODYNAMIC RELATIONS
# =============================================================================
class ThermoRelations:
    def area_ratio_from_mach(M):
        return (1.0/M) * ((2.0/(self.GAMMA+1.0))*(1.0 + 0.5*(self.GAMMA-1.0)*M**2))**(0.5*(self.GAMMA+1.0)/(self.GAMMA-1.0))

    def isentropic_pressure_ratio(M):
        return (1.0 + 0.5*(self.GAMMA-1.0)*M**2)**(-self.GAMMA/(self.GAMMA-1.0))

    def isentropic_temperature_ratio(M):
        return (1.0 + 0.5*(self.GAMMA-1.0)*M**2)**(-1.0)

    def isentropic_density_ratio(M):
        return (1.0 + 0.5*(self.GAMMA-1.0)*M**2)**(-1.0/(self.GAMMA-1.0))

    def normal_shock_relations(M1):
        if M1 <= 1.0:
            return 1.0, 1.0, 1.0, 1.0
        else:
            M2 = np.sqrt((1.0 + 0.5*(self.GAMMA-1.0)*M1**2) /
                         (self.GAMMA*M1**2 - 0.5*(self.GAMMA-1.0)))
            P2_P1 = 1.0 + (2.0*self.GAMMA/(self.GAMMA+1.0))*(M1**2 - 1.0)
            T2_T1 = (1.0 + (2.0*self.GAMMA/(self.GAMMA+1.0))*(M1**2 - 1.0)) * \
                ((2.0 + (self.GAMMA-1.0)*M1**2)/((self.GAMMA+1.0)*M1**2))
            P02_P01 = (((self.GAMMA+1.0)*M1**2)/((self.GAMMA-1.0)*M1**2 + 2.0))**(self.GAMMA/(self.GAMMA-1.0)) * \
                ((self.GAMMA+1.0)/(2.0*self.GAMMA*M1**2 -
                 (self.GAMMA-1.0)))**(1.0/(self.GAMMA-1.0))
            return M2, P2_P1, T2_T1, P02_P01

    def mach_from_area_ratio(A_ratio, subsonic=True):
        def equation(M):
            return ThermoRelations.area_ratio_from_mach(M) - A_ratio
        if subsonic:
            guess = 0.1
        else:
            guess = 2.0
        return fsolve(equation, guess)[0]


# =============================================================================
#    4. NOZZLE FLOW SOLVER
# =============================================================================
class FlowSolver:
    def throat_mach(P0, Pe):
        Pe_choked = ThermoRelations.isentropic_pressure_ratio(
            ThermoRelations.mach_from_area_ratio(self.A_exit/self.A_throat))*P0

        if P0 == Pe:
            M = 0
            Pi = P0

        elif Pe_choked >= Pe:
            M = 1
            Pi = ThermoRelations.isentropic_pressure_ratio(
                ThermoRelations.mach_from_area_ratio(self.A_inlet/self.A_throat))*P0

        else:
            def equation(M):
                A_ratio_exit = ThermoRelations.area_ratio_from_mach(
                    M)*self.A_exit/self.A_throat
                M_exit = ThermoRelations.mach_from_area_ratio(A_ratio_exit)
                P_ratio_exit = ThermoRelations.isentropic_pressure_ratio(
                    M_exit)
                return P_ratio_exit*P0 - Pe
            guess = 0.5
            M = fsolve(equation, guess)[0]
            Pi = ThermoRelations.isentropic_pressure_ratio(ThermoRelations.mach_from_area_ratio(
                ThermoRelations.area_ratio_from_mach(M)*self.A_inlet/self.A_throat))*P0
        return M, Pi

    def solve_flow(Pe=None, shock_position=None, M_throat=None, n=500):

        # CONDITION 1: Given Back Pressure:
        if Pe is not None:
            self.Pe = Pe
            self.M_throat, self.Pi = FlowSolver.throat_mach(self.P0, self.Pe)
            self.Pe_ratio = Pe/self.P0

            # CASE 1.1: Fully Subsonic:
            if self.M_throat < 1:
                flow_type = "Fully Subsonic"
                Output = FlowSolver.subsonic_solver(n)
                P_P0 = [r['P/P0'] for r in Output]
                self.P_total = Pe/P_P0[-1]
                return Output, flow_type, None, self.P_total

            # CASE 1.2: Chocked Subsonic:
            elif self.Pe_ratio == ThermoRelations.isentropic_pressure_ratio(ThermoRelations.mach_from_area_ratio(self.A_exit/self.A_throat, subsonic=True)):
                flow_type = "Choked Subsonic"
                Output = FlowSolver.subsonic_solver(n)
                P_P0 = [r['P/P0'] for r in Output]
                self.P_total = Pe/P_P0[-1]
                return Output, flow_type, None, self.P_total

            else:
                # Finding Shock Position
                self.x_shock, self.M1, self.M2, self.M_exit = FlowSolver.find_shock_location()

            # CASE 1.3: Normal Shock in Nozzle:
                if self.x_shock is not None:
                    flow_type = "Normal Shock in Nozzle"
                    Output = FlowSolver.shock_solver(n)
                    P_P0 = [r['P/P0'] for r in Output]
                    self.P_total = Pe/P_P0[-1]
                    return Output, flow_type, self.x_shock, self.P_total

            # CASE 1.4: Isentropic Supersonic:
                else:
                    Output = FlowSolver.supersonic_solver(n)
                    P_P0 = [r['P/P0'] for r in Output]
                    self.P_total = self.P0

                    if P_P0[-1]*self.P_total > Pe:
                        flow_type = 'Fully Supersonic (Under Expanded)'

                    elif P_P0[-1]*self.P_total == Pe:
                        flow_type = 'Fully Supersonic (Ideal Expansion)'

                    else:
                        flow_type = 'Fully Supersonic (Over Expanded)'

                    return Output, flow_type, None, self.P_total

        # CONDITION 2: Given Shock Position:
        elif shock_position is not None:
            self. x_shock = shock_position
            flow_type = "Shock in Nozzle"

            # Calculating Mach no. for Given Shock Position
            self.A_shock = NozzleGeometery.nozzle_area(self.x_shock)
            A_ratio_shock = self.A_shock/self.A_throat
            self.M1 = ThermoRelations.mach_from_area_ratio(
                A_ratio_shock, subsonic=False)
            self.M2, P2_P1, T2_T1, P02_P01 = ThermoRelations.normal_shock_relations(
                self.M1)
            M_exit = ThermoRelations.mach_from_area_ratio(
                self.A_exit/self.A_shock*ThermoRelations.area_ratio_from_mach(self.M2), subsonic=True)
            M_inlet = ThermoRelations.mach_from_area_ratio(
                self.A_inlet/self.A_throat, subsonic=True)

            Output = FlowSolver.shock_solver(n)
            P_P0 = [r['P/P0'] for r in Output]
            self.Pe = P_P0[-1]*P02_P01/P_P0[0]*self.P0
            self.P_total = self.Pe/P_P0[-1]
            return Output, flow_type, self.x_shock, self.P_total

        # CONDITION 3: Given Throat Mach No.:
        elif M_throat is not None:
            self.M_throat = M_throat
            Output = FlowSolver.subsonic_solver(n)
            P_P0 = [r['P/P0'] for r in Output]
            self.P_total = self.P0/P_P0[0]
            self.Pe = P_P0[-1]*self.P_total

            if M_throat == 1:
                flow_type = "Choked Subsonic"
            else:
                flow_type = "Fully Subsonic"
            return Output, flow_type, None, self.P_total

    def subsonic_solver(n):
        """Solving for Subsonic Flow"""
        x = np.linspace(0, self.Lc + self.Ld, n)
        results = []
        A_throat_ratio = ThermoRelations.area_ratio_from_mach(self.M_throat)
        for xi in x:
            A = NozzleGeometery.nozzle_area(xi)
            A_ratio = A_throat_ratio*A/self.A_throat

            # Subsonic Mach No.
            M = ThermoRelations.mach_from_area_ratio(A_ratio, subsonic=True)

            # Thermodynamic Parameters
            P_P0 = ThermoRelations.isentropic_pressure_ratio(M)
            T_T0 = ThermoRelations.isentropic_temperature_ratio(M)
            rho_rho0 = ThermoRelations.isentropic_density_ratio(M)
            P0_P0ref = 1.0

            results.append({'x': xi, 'M': M, 'P/P0': P_P0, 'T/T0': T_T0,
                           'rho/rho0': rho_rho0, 'P0/P0_ref': P0_P0ref})
        return results

    def supersonic_solver(n):
        """Solving for Supersonic Flow"""
        x = np.linspace(0, self.Lc + self.Ld, n)
        results = []
        A_throat_ratio = ThermoRelations.area_ratio_from_mach(self.M_throat)
        for xi in x:
            A = NozzleGeometery.nozzle_area(xi)
            A_ratio = A_throat_ratio*A/self.A_throat

            if xi <= self.Lc:
                # Subsonic Mach No.
                M = ThermoRelations.mach_from_area_ratio(
                    A_ratio, subsonic=True)
            else:
                # Supersonic Mach No.
                M = ThermoRelations.mach_from_area_ratio(
                    A_ratio, subsonic=False)

            # Thermodynamic Parameters
            P_P0 = ThermoRelations.isentropic_pressure_ratio(M)
            T_T0 = ThermoRelations.isentropic_temperature_ratio(M)
            rho_rho0 = ThermoRelations.isentropic_density_ratio(M)
            P0_P0ref = 1.0

            results.append({'x': xi, 'M': M, 'P/P0': P_P0, 'T/T0': T_T0,
                           'rho/rho0': rho_rho0, 'P0/P0_ref': P0_P0ref})
        return results

    def find_shock_location(num_points=1000):
        """Finding Shock Location with Trial and Error Method"""
        x_start = self.Lc
        x_end = self.Lc + self.Ld
        x_positions = np.linspace(x_start, x_end+x_end/num_points, num_points)

        for x_shock in x_positions:

            # Flow Before Shock (Isentropic)
            self.A_shock = NozzleGeometery.nozzle_area(x_shock)
            A_ratio_shock = self.A_shock / self.A_throat

            # Mach no. Before Shock
            M1 = ThermoRelations.mach_from_area_ratio(
                A_ratio_shock, subsonic=False)

            # Normal Shock Relations
            M2, P2_P1, T2_T1, P02_P01 = ThermoRelations.normal_shock_relations(
                M1)

            # Total Pressure After Shock
            P02 = self.P0 * P02_P01

            # Flow After Shock (Isentropic)
            A_exit_ratio = self.A_exit / self.A_shock
            M_exit = ThermoRelations.mach_from_area_ratio(
                A_exit_ratio*ThermoRelations.area_ratio_from_mach(M2), subsonic=True)

            # Exhaust Pressure Calculated
            P_exit_calc = P02*ThermoRelations.isentropic_pressure_ratio(M_exit)

            # Comparison with Back Pressure
            if abs(P_exit_calc - self.Pe) / P_exit_calc < 0.01:
                break

        if x_shock <= x_end:
            return x_shock, M1, M2, M_exit
        else:
            return None, None, None, None

    def shock_solver(n):
        """Solving for Flow with Specified Shock Position"""
        x = np.linspace(0, self.Lc + self.Ld, n)
        results = []
        # Shock Relations
        self.M2, P2_P1, T2_T1, P02_P01 = ThermoRelations.normal_shock_relations(
            self.M1)
        for xi in x:
            A = NozzleGeometery.nozzle_area(xi)
            A_ratio = A/self.A_throat

            # Before Throat: Subsonic
            if xi <= self.Lc:
                M = ThermoRelations.mach_from_area_ratio(
                    A_ratio, subsonic=True)
                P_P0 = ThermoRelations.isentropic_pressure_ratio(M)
                T_T0 = ThermoRelations.isentropic_temperature_ratio(M)
                rho_rho0 = ThermoRelations.isentropic_density_ratio(M)
                P0_P0ref = 1.0

            # After Throat, Before Shock: Supersonic
            elif xi <= self.x_shock:
                M = ThermoRelations.mach_from_area_ratio(
                    A_ratio, subsonic=False)
                P_P0 = ThermoRelations.isentropic_pressure_ratio(M)
                T_T0 = ThermoRelations.isentropic_temperature_ratio(M)
                rho_rho0 = ThermoRelations.isentropic_density_ratio(M)
                P0_P0ref = 1.0

            # After Shock: Subsonic
            else:
                # Effective Area Calculation After Shock
                A_ratio_effective = A/self.A_shock * \
                    ThermoRelations.area_ratio_from_mach(self.M2)
                M = ThermoRelations.mach_from_area_ratio(
                    A_ratio_effective, subsonic=True)

                # Total Pressure After Shock
                P0_after_shock = self.P0*P02_P01
                P_P0 = ThermoRelations.isentropic_pressure_ratio(
                    M)*(P0_after_shock/self.P0)
                T_T0 = ThermoRelations.isentropic_temperature_ratio(M)
                rho_rho0 = ThermoRelations.isentropic_density_ratio(M)
                P0_P0ref = P02_P01

            results.append({'x': xi, 'M': M, 'P/P0': P_P0, 'T/T0': T_T0,
                            'rho/rho0': rho_rho0, 'P0/P0_ref': P0_P0ref})

        return results


# =============================================================================
#    PLOTTING AND ANIMATION FUNCTIONS
# =============================================================================
class PlotFunc:
    def plot_nozzle_results(results, flow_type, shock_x=None, case_name=""):
        """Plotting Nozzle Results"""
        x = [r['x'] for r in results]
        M = [r['M'] for r in results]
        P_P0 = [r['P/P0'] for r in results]
        T_T0 = [r['T/T0'] for r in results]
        rho_rho0 = [r['rho/rho0'] for r in results]
        P0_P0ref = [r['P0/P0_ref'] for r in results]

        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(
            f'NOZZLE FLOW ANALYSIS - {case_name} ({flow_type})', fontsize=14, fontweight='bold')

        # Plot 1: Mach no.
        axes[0, 0].plot(x, M, 'k-', linewidth=2)
        axes[0, 0].set_xlabel('Axial Position x (m)')
        axes[0, 0].set_ylabel('Mach Number')
        axes[0, 0].set_title('Mach Number Distribution')
        axes[0, 0].grid(True, alpha=0.3)
        if shock_x is not None:
            axes[0, 0].axvline(x=shock_x, color='r', linestyle='--',
                               label=f'Shock at x={shock_x:.3f}m')
            axes[0, 0].legend()

        # Plot 2: Static Pressure Ratio
        axes[0, 1].plot(np.append(x, x[-1]+0.1), np.append(P_P0,
                        self.Pe/self.P_total), 'k-', linewidth=2)
        axes[0, 1].axvline(x=self.Lc+self.Ld, color='k', linestyle=':')
        axes[0, 1].set_xlabel('Axial Position x (m)')
        axes[0, 1].set_ylabel('Static Pressure Ratio P/P₀')
        axes[0, 1].set_title('Static Pressure Distribution')
        axes[0, 1].grid(True, alpha=0.3)
        if shock_x is not None:
            axes[0, 1].axvline(x=shock_x, color='r', linestyle='--')

        # Plor 3: Static Temperature Ratio
        axes[0, 2].plot(x, T_T0, 'k-', linewidth=2)
        axes[0, 2].set_xlabel('Axial Position x (m)')
        axes[0, 2].set_ylabel('Static Temperature Ratio T/T₀')
        axes[0, 2].set_title('Static Temperature Distribution')
        axes[0, 2].grid(True, alpha=0.3)
        if shock_x is not None:
            axes[0, 2].axvline(x=shock_x, color='r', linestyle='--')

        # Plot 4: Density Ratio
        axes[1, 0].plot(x, rho_rho0, 'k-', linewidth=2)
        axes[1, 0].set_xlabel('Axial Position x (m)')
        axes[1, 0].set_ylabel('Density Ratio ρ/ρ₀')
        axes[1, 0].set_title('Density Distribution')
        axes[1, 0].grid(True, alpha=0.3)
        if shock_x is not None:
            axes[1, 0].axvline(x=shock_x, color='r', linestyle='--')

        # Plot 5: Total Pressure Ratio
        axes[1, 1].plot(x, P0_P0ref, 'k-', linewidth=2)
        axes[1, 1].set_xlabel('Axial Position x (m)')
        axes[1, 1].set_ylabel('Total Pressure Ratio P₀/P₀,ref')
        axes[1, 1].set_title('Total Pressure Distribution')
        axes[1, 1].grid(True, alpha=0.3)
        if shock_x is not None:
            axes[1, 1].axvline(x=shock_x, color='r', linestyle='--')

        # Plot 6: Nozzle Profile
        radii = [NozzleGeometery.nozzle_radius(xi) for xi in x]
        axes[1, 2].plot(x, radii, 'k-', linewidth=2)
        axes[1, 2].plot(x, [-r for r in radii], 'k-', linewidth=2)
        axes[1, 2].set_xlabel('Axial Position x (m)')
        axes[1, 2].set_ylabel('Radius (m)')
        axes[1, 2].set_title('Nozzle Profile')
        axes[1, 2].grid(True, alpha=0.3)
        axes[1, 2].set_aspect('equal')
        if shock_x is not None:
            axes[1, 2].plot([shock_x, shock_x], [-NozzleGeometery.nozzle_radius(shock_x),
                                                 NozzleGeometery.nozzle_radius(shock_x)], color='r', linestyle='--')

        plt.tight_layout()

    def create_flow_animation(Pe_range, filename='nozzle_flow_animation'):
        """Creating Animation with Change in Back Pressure"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('C-D Nozzle Flow Animation - Effect of Back Pressure',
                     fontsize=16, fontweight='bold')
        plt.tight_layout(pad=3)

        # Axis and Style
        axes = [ax1, ax2, ax3, ax4]
        titles = ['Mach Number',
                  'Static Pressure (P/P₀)', 'Static Temperature (T/T₀)', 'Nozzle Profile']
        ylabels = ['Mach Number', 'P/P₀', 'T/T₀', 'Radius (m)']

        for ax, title, ylabel in zip(axes, titles, ylabels):
            ax.set_xlabel('Axial Position x (m)')
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0, self.Lc + self.Ld+0.1)

        # Y Axis Limits
        ax1.set_ylim(0, 3.2)
        ax2.set_ylim(0, 1.1)
        ax3.set_ylim(0.3, 1.1)
        ax4.set_ylim(-0.4, 0.2)
        ax4.set_aspect('equal')

        # Initial Lines
        x = np.linspace(0, self.Lc + self.Ld, 500)
        radii = [NozzleGeometery.nozzle_radius(xi) for xi in x]

        lines = []
        for ax in axes[:3]:
            line, = ax.plot([], [], 'k-', linewidth=2)
            lines.append(line)

        # Nozzle Profile
        nozzle_line1, = ax4.plot(x, radii, 'k-', linewidth=3)
        nozzle_line2, = ax4.plot(x, [-r for r in radii], 'k-', linewidth=3)
        exhaust_line = ax2.axvline(
            x=self.Lc+self.Ld, color='k', linestyle=':', linewidth=1)
        shock_line = ax4.axvline(
            x=0, ymin=0, ymax=1, color='r', linestyle='--', linewidth=2, alpha=0)

        # Status
        status_text = ax4.text(0.02, 0.2, '', transform=ax4.transAxes, fontsize=12,
                               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        pressure_text = ax4.text(0.02, 0.08, '', transform=ax4.transAxes, fontsize=10,
                                 bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

        # Animation DATA
        def animate(frame):
            Pe = Pe_range[frame]
            results, flow_type, shock_x, P_total = FlowSolver.solve_flow(Pe)
            PlotFunc.progress_bar(Pe-Pe_range[0], Pe_range[-1]-Pe_range[0], 60)

            x_vals = [r['x'] for r in results]
            M_vals = [r['M'] for r in results]
            P_vals = [r['P/P0'] for r in results]
            T_vals = [r['T/T0'] for r in results]

            # UPDATING LINES
            lines[0].set_data(x_vals, M_vals)
            lines[1].set_data(np.append(x_vals, x_vals[-1]+0.1),
                              np.append(P_vals, Pe/P_total))
            lines[2].set_data(x_vals, T_vals)

            # UPDATING SHOCK
            if shock_x is not None:
                shock_r = NozzleGeometery.nozzle_radius(shock_x)
                shock_line.set_xdata([shock_x, shock_x])
                shock_line.set_alpha(1.0)
            else:
                shock_line.set_alpha(0)

            # UPDATING STATUS
            status_text.set_text(
                f'Flow Type: {flow_type}\nBack Pressure: {Pe:.4f} bar')
            if shock_x is not None:
                pressure_text.set_text(f'Shock at x = {shock_x:.3f} m')
            else:
                pressure_text.set_text('No Shock')

            return lines + [nozzle_line1, nozzle_line2, shock_line, status_text, pressure_text]

        # Creating Animation
        anim = animation.FuncAnimation(
            fig, animate, frames=len(Pe_range),
            interval=50, blit=True, repeat=True
        )

        # Saving Animation
        # print("Creating animation...")
        anim.save(
            filename=f'./Plots/{filename}.mkv', writer='ffmpeg', dpi=300)
        print(f"\nAnimation saved as {filename}")
        return anim

    def progress_bar(progress, total, length):
        percent = (length) * (progress/total)
        bar = '█' * int(percent) + '░' * ((length) - int(percent))
        print(f'\r{bar} {percent*100/(length):.2f}%', end='\r')
