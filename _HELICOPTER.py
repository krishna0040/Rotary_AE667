from _ROTOR import *
from _MISSION_PLANNER import *
import numpy as np
from scipy.optimize import fsolve

class surface:
    def __init__(self, name, pos, area,
                 cl=0.0, cd=0.0,
                 cl_alpha=0.0, cl_delta=0.0,         ### cl vs alpha slope, cl vs delta slopt
                 alpha=0.0, delta=0.0,               ### alpha, delta [the deflection of the control surface]
                 surface_type="lifting",             ### "lifting", "control"
                 lift_direction=np.array([0.0, 0.0, -1.0])):
        
        self.name = name
        self.pos = np.array(pos)
        self.area = area
        self.cl = cl
        self.cd = cd
        self.cl_alpha = cl_alpha
        self.cl_delta = cl_delta
        self.alpha = alpha
        self.delta = delta
        self.surface_type = surface_type
        self.rho = 1.225  ### 1.225 kg/m3
        # Normalize the lift direction vector
        self.lift_direction = np.array(lift_direction, dtype=float)
        self.lift_direction /= np.linalg.norm(self.lift_direction)

    def compute_forces_and_moments(self, rho, V, ref_point):
        q = 0.5 * rho * V**2
        self.rho = rho
        # Compute total CL
        cl_total = self.cl + self.cl_alpha * self.alpha + self.cl_delta * self.delta
        
        # Aerodynamic forces
        L = q * self.area * cl_total
        D = q * self.area * self.cd

        # Drag direction (opposite to X-axis)
        drag_direction = np.array([-1.0, 0.0, 0.0])

        # Lift along the user-defined lift_direction
        F_lift = L * self.lift_direction
        F_drag = D * drag_direction
        
        # Total aerodynamic force
        F = F_lift + F_drag
        
        # Compute moment about reference point
        r = self.pos - np.array(ref_point)
        M = np.cross(r, F)
        
        return F, M

class helicopter:
    def __init__(self,
                 main_rotor,
                 tail_rotor, 
                 main_pos,      ### pass self.with correct sign conventions
                 tail_pos,
                 lifting_surfaces,
                 control_surfaces,
                 stabilizers,
                 weight,   # weight in kg
                 ### weight pos is considered the reference point
                 ref_point = np.array([0,0,0])):
        
        self.main_rotor = main_rotor
        self.tail_rotor = tail_rotor
        self.rho = main_rotor.rho
        self.main_pos = main_pos      # (x, y, z) of main rotor
        self.tail_pos = tail_pos      # (x, y, z) of tail rotor
        self.lifting_surfaces = lifting_surfaces
        self.control_surfaces = control_surfaces
        self.stabilizers = stabilizers
        self.ref_point = ref_point    # reference point for force/moment summation
        self.weight = weight*9.81     # self.weight in N 

    def calculate_total_forces(self,velocity, theta0, theta1c, theta1s, alpha_tpp, theta0_tail):
        main_rotor = self.main_rotor.forward(velocity, theta0, theta1c, theta1s, alpha_tpp)[0]   
        tail_rotor = self.tail_rotor.forward(velocity, theta0_tail, 0, 0, 0)[0]
        tail_rotor['FY'], tail_rotor['FZ'] = -tail_rotor['FZ'], tail_rotor['FY'] 
        tail_rotor['MY'], tail_rotor['MZ'] = -tail_rotor['MZ'], tail_rotor['MY']

        force = np.zeros(3)
        moment = np.zeros(3)

        # Main rotor
        force += np.array([main_rotor['FX'], main_rotor['FY'], main_rotor['FZ']])
        moment += np.array([main_rotor['MX'], main_rotor['MY'], main_rotor['MZ']])

        # Tail rotor
        force += np.array([tail_rotor['FX'], tail_rotor['FY'], tail_rotor['FZ']])
        moment += np.array([tail_rotor['MX'], tail_rotor['MY'], tail_rotor['MZ']])

        # Lifting surfaces
        for surf in self.lifting_surfaces:
            F, M = surf.compute_forces_and_moments(self.rho, velocity,self.ref_point)
            force += F
            moment += M

        # Control surfaces
        for surf in self.control_surfaces:
            F, M = surf.compute_forces_and_moments(self.rho, velocity,self.ref_point)
            force += F
            moment += M

        moment += np.cross(self.main_pos - self.ref_point,np.array([main_rotor['FX'], main_rotor['FY'], main_rotor['FZ']]))
        moment += np.cross(self.tail_pos - self.ref_point,np.array([tail_rotor['FX'], tail_rotor['FY'], tail_rotor['FZ']]))
        
        return {'FX': force[0], 'FY': force[1], 'FZ': force[2],
                'MX': moment[0], 'MY': moment[1], 'MZ': moment[2]}

    def find_trim(self, velocity, alpha_tpp,Drag):
        def trim_objective(x):
            theta0, theta1c, theta1s, theta0_tail = x
            fm = self.calculate_total_forces(
                velocity=velocity,
                theta0=theta0,
                theta1c=theta1c,
                theta1s=theta1s,
                alpha_tpp=alpha_tpp,
                theta0_tail=theta0_tail,
            )

            res = np.zeros(4)
            res[0] = fm['FZ'] + self.weight        # Thrust balances weight
            res[1] = fm['FX'] - Drag     # Rotor drag balances fuselage drag
            res[2] = fm['MY']            # Pitching moment = 0
            res[3] = fm['MZ']
            return res

        x0 = [0.1, 0.0, 0.0, 0.05]
        sol = fsolve(trim_objective, x0)
        theta0_trim, theta1c_trim, theta1s_trim, theta0_tail_trim = sol

        fm_trim = self.calculate_total_forces(
            velocity=velocity,
            theta0=theta0_trim,
            theta1c=theta1c_trim,
            theta1s=theta1s_trim,
            alpha_tpp=alpha_tpp,
            theta0_tail=theta0_tail_trim
        )

        print("\n=== Trim Results ===\n")
        print(f"Theta0 (main):  {np.degrees(theta0_trim):.2f}°")
        print(f"Theta1c:        {np.degrees(theta1c_trim):.2f}°")
        print(f"Theta1s:        {np.degrees(theta1s_trim):.2f}°")
        print(f"Theta0_tail:    {np.degrees(theta0_tail_trim):.2f}°")

        print("\nForces (N): FX, FY, FZ =", fm_trim['FX'], fm_trim['FY'], fm_trim['FZ'])
        print("Moments (N·m): MX, MY, MZ =", fm_trim['MX'], fm_trim['MY'], fm_trim['MZ'])

        return {
            'theta0': theta0_trim,
            'theta1c': theta1c_trim,
            'theta1s': theta1s_trim,
            'theta0_tail': theta0_tail_trim,
            'forces': fm_trim
        }







