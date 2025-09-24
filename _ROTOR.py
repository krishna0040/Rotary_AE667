import numpy as np
import matplotlib.pyplot as plt 
class rotor:
    def __init__(self):
        #change accordingly, initial set at 3500 m
        self.temp = 288.15    # ambient is 288.15 K
        self.rho = 1.225      # ambient is 1.006kg/m3
        self.pressure = 101325 # ambient is 79500 Pa
        self.sound_speed = 340.29 # ambient is 340.29 m/s

        # only the payload considered initially
        self.m = 8 * 70    # 8 people, 70kg each
        self.b = 3
        self.r = 0.762         # radius defined = 0.762m
        self.rc = 0.125        # inner radius = 0.125m

        self.twist = 0         # all the angles in radian
        self.twist_slope = 0

        self.chord = 0.0508
        self.chord_slope = 0

        self.omega = 30    # omega in rad/sec

        # cl = a0 * alpha, alpha (alpha in rad )
        self.cl = 5.75
        # cd = cdo + ep*alpha (alpha in rad )
        self.cd = 1.25
        self.cdo = 0.0113

        self.V = None           # climb velocity
        self.lambda_c = None

    def integrate(self,x, y):
        """
        @brief Compute the integral ydx from [a,b], take the input as np.array
        """
        return np.trapz(y,x)

    def set_rotor(self, radius = 0.762, root_cut_out = 0.125, blades = 3, chord0 = 0.0508,chord_slope = 0,
              twist0 = 0, twist_slope = 0,omega = 30):
        """
        @brief defines the parameters for the rotor, set the angles and everything in terms of rad-1 and rad
        @param radius: outer radius of the rotor
        @param root_cut_out: radius, where the rotors start
        @param blades: number of blades
        @param chord0 = chord at r = root_cut_out
        @param chord_slope: slope of the chord, chord at radius r = chord + slope*(r-root_cut_out), keep the slope negative
        @param twist0 : twist at r = root_cut_out
        @param twist_slope: slope of the twist, twist at radius r = twist0 + slope*(r-root_cut_out)
        """
        self.r = radius
        self.rc = root_cut_out
        self.b = blades
        self.chord = chord0
        self.chord_slope = chord_slope
        self.twist = twist0
        self.twist_slope = twist_slope
        self.omega = omega

    def set_airfoil(self, cl = 5.75, cdo = 0.0113,cd = 1.25):
        """
        @brief call this func to change the airfoil parameters
        @param cl: cl_val = cl*alph
        @param cd: cd_val = cdo + cd*alpha^2
        """
        self.cl = cl
        self.cdo = cdo
        self.cd = cd

    def atmosphere(self, h, dT=0):
        """
        @brief used to set the ambient conditions at the given altitude, and also returns them
        @return temp (K)  pressure (Pa), density (kgm/3), speed of sound(m/s)
        """
        T = 288.15 - 0.0065*h + dT
        p = 101325 * (1 - 0.00198*h/288.16)**5.2553
        rho = 1.225 * (1 - 0.00198*h/288.16)**4.2553
        a = (1.4*287.05*T)**0.5
        self.temp, self.pressure, self.rho, self.sound_speed = T, p, rho, a
        return T, p, rho, a


    def cal_cl(self,alpha):
        return self.cl * alpha
    def cal_cd(self,alpha):
        return self.cdo + self.cd*alpha*alpha
    def cal_chord(self,x):                          ### local chord 
        return self.chord + self.chord_slope*(x-self.rc)
    def cal_twist(self,x):                          ### local twist
        return  self.twist + self.twist_slope*(x-self.rc)
    def cal_solidity(self,r):                       ### local solidity
        return self.b  *self.cal_chord(r)/(np.pi * r)      
    def cal_F(self, val_lambda, r):
        x = np.exp(-self.b/2 * (1 - r/self.r) / max(val_lambda, 1e-8))
        x = np.clip(x, 0.0, 1.0)
        return 2 * np.arccos(x) / np.pi

    def find_lambda(self,sigma,r,theta, tol = 1e-8,beta = 0.6,iterations = 100):   ## this is the iterative step to calculate lambda for a particular r
        guess = self.lambda_c      ### start with the initial guess as lambda_c
        for i in range(iterations):
            F = self.cal_F(guess,r)
            nxt = np.sqrt((sigma*self.cl/(16*F) - self.lambda_c/2)**2  +  (sigma*self.cl/(8*F))*theta*(r/self.r)) + self.lambda_c/2 - sigma*self.cl/(16*F) 
            if abs(nxt-guess) < tol:
                guess = nxt
                break
            guess = nxt* beta + (1-beta)*guess

        v = guess*self.omega*self.r - self.V        # induced velocity, local induced
        return guess,v,self.cal_F(guess,r)

    def bemt(self,V,n=1000):       

        """
        @brief calculates everything - thrust, torque, lambda, induced velocity for hover or upward flight only
        @param V: climb velocity   [define the rotor parameters and all before]
        @return: lambdas,v,sigma,twist, thrust ,dT,Fs,d_torque, torque

        """
        r = np.linspace(self.rc,self.r - 1e-6,n)
        wt = self.omega*r
        self.V = V
        self.lambda_c = self.V/(self.omega * self.r)       ### global one

        sigma = np.array([self.cal_solidity(r[i]) for i in range(len(r))])
        twist = np.array([self.cal_twist(r[i]) for i in range(len(r))])

        ### calculating the parameters for thrust calculation
        lambdas = []
        Fs = []
        v = []
        for i in range(len(r)):
          temp_lambda , temp_v, temp_F = self.find_lambda (sigma[i],r[i],twist[i])
          lambdas.append(temp_lambda)
          v.append(temp_v)
          Fs.append(temp_F)

        ### thrust calculation
        v = np.array(v)
        lambdas = np.array(lambdas)
        Fs = np.array(Fs)
        dT = 4*np.pi*self.rho * r * Fs * v * (V+v)
        thrust = self.integrate(r,dT)

        ### torque calculation
        U2 = (self.omega * r)**2 + (self.V + v)**2
        chord = np.array([self.cal_chord(r[i]) for i in range(len(r))])

        phi = np.arctan((self.V + v) / (self.omega * r))
        alpha = twist - phi               ### net angle of attack
        cl = np.array([self.cal_cl(alpha[i]) for i in range(len(r))])
        cd = np.array([self.cal_cd(alpha[i]) for i in range(len(r))])
        d_torque = 0.5*self.b * r *chord* self.rho * U2*(cd * np.cos(phi) + cl*np.sin(phi))
        torque = self.integrate(r,d_torque)

        return lambdas,v,sigma,twist, thrust ,dT,Fs,d_torque, torque

    def integrate_2d(self, r, psi, f):
        """
        2D integral of f(r,psi) where:
         - r is 1D array of radial samples (len n_r)
         - psi is 1D array of azimuth samples (len n_psi)
         - f is shape (n_r, n_psi)
        Returns scalar integral using trapezoidal rule in psi then r.
        NOTE: this function already accounts for spacing (no extra dr*dpsi needed).
        """
        # integrate in psi (axis=1) -> returns array length n_r
        integral_over_psi = np.trapz(f, x=psi, axis=1)
        # integrate in r -> scalar
        total = np.trapz(integral_over_psi, x=r)
        return total

    def calculate_glauert_lambda(self, mu, alpha_tpp, CT, tol=1e-6, max_iter=100):
        """
        Keep your Glauert iteration - unchanged except return doc clarity.
        """
        lambda_G = mu * np.tan(alpha_tpp) / 2  # initial guess
        for i in range(max_iter):
            denominator = np.sqrt(mu**2 + lambda_G**2)
            lambda_G_new = mu * np.tan(alpha_tpp) / 2 + CT / (2 * denominator)
            if abs(lambda_G_new - lambda_G) < tol:
                lambda_G = lambda_G_new
                break
            lambda_G = lambda_G_new
        lambda_iG = lambda_G - mu * np.sin(alpha_tpp)
        return lambda_G, lambda_iG

    def calculate_non_uniform_inflow(self, r, psi, mu, lambda_G, lambda_iG):
        """
        Build full (r,psi) distribution of induced inflow.
        r is radial array (dimensional), psi is azimuth array.
        Note: use non-dimensional radius for R_grid in factor if desired.
        """
        R_nd, PSI_grid = np.meshgrid(r / self.r, psi, indexing='ij')  # nondim radius
        inflow_factor = (1 + ((4.0/3.0) * (mu / lambda_G))**1.2 + mu / lambda_G)
        # distribution has sign cos(psi) (longitudinal asymmetry)
        lambda_i = lambda_iG * inflow_factor * R_nd * np.cos(PSI_grid)
        # convert to dimensional induced velocity if you want: lambda_i * omega * r
        return lambda_i
    
    def calculate_flapping_angles(self, theta0, theta1c, theta1s, mu, alpha_tpp):
        """
        @brief Calculate flapping angles beta0, beta1c, beta1s
        """
        # Simplified flapping calculations (from blade dynamics)
        beta0 = theta0 / 8  # Coning angle approximation
        beta1c = (theta1s - mu * alpha_tpp) / (1 - mu**2/2)  # Longitudinal flapping
        beta1s = -(theta1c + mu * theta0/8) / (1 + mu**2/2)  # Lateral flapping
        
        return beta0, beta1c, beta1s

    def forward(self, forward_speed, theta0, theta1c, theta1s, alpha_tpp,
                CT_guess=0.008, n_r=100, n_psi=36):
        # discretization
        r = np.linspace(self.rc, self.r - 1e-6, n_r)
        # for periodic azimuth integrals: avoid duplicating 0 and 2pi
        psi = np.linspace(0.0, 2.0 * np.pi, n_psi, endpoint=False)

        mu = forward_speed / (self.omega * self.r)

        # Glauert inflow
        lambda_G, lambda_iG = self.calculate_glauert_lambda(mu, alpha_tpp, CT_guess)

        # flapping
        beta0, beta1c, beta1s = self.calculate_flapping_angles(theta0, theta1c, theta1s, mu, alpha_tpp)

        # lambda_i distribution (non-dimensional: lambda = w/(omega*r))
        lambda_i = self.calculate_non_uniform_inflow(r, psi, mu, lambda_G, lambda_iG)

        # meshes (dimensional radius used for velocities)
        R_grid, PSI_grid = np.meshgrid(r, psi, indexing='ij')  # shape (n_r, n_psi)

        # Tangential velocity UT:
        # baseline: omega * r
        # forward flight component projected into tangential direction: V * sin(psi)
        # (if TPP tilt included, incorporate cos/sin(alpha_tpp) per derivation)
        UT = self.omega * R_grid + forward_speed * np.sin(PSI_grid)

        # flapping rate -> tangential/axial contribution
        beta_dot = self.omega * (-beta1c * np.sin(PSI_grid) + beta1s * np.cos(PSI_grid))

        # Axial (perpendicular) velocity UP:
        # induced axial: lambda_i * omega * r  (lambda_i is nondim factor)
        # forward axial projection: V * cos(psi) * cos(alpha_tpp)  (approx)
        # flapping radial contribution: R_grid * beta_dot
        UP = (lambda_i * self.omega * self.r +
              forward_speed * np.cos(PSI_grid) * np.cos(alpha_tpp) +
              R_grid * beta_dot +
              forward_speed * np.sin(beta0) * np.cos(PSI_grid))

        # total local speed
        U_total = np.sqrt(UT**2 + UP**2)

        # inflow angle phi = atan2(UP, UT)   (signed and handles UT<0 correctly)
        phi = np.arctan2(UP, UT)

        # blade pitch distribution
        twist_grid = np.array([self.cal_twist(rv) for rv in r])[:, np.newaxis]  # (n_r,1)
        theta_total = (theta0 +
                       theta1c * np.cos(PSI_grid) +
                       theta1s * np.sin(PSI_grid) +
                       twist_grid)

        # angle of attack
        alpha = theta_total - phi

        # compute aero coefficients (ensure cal_cl/cal_cd accept arrays)
        cl = self.cal_cl(alpha)
        cd = self.cal_cd(alpha)

        chord_grid = np.array([self.cal_chord(rv) for rv in r])[:, np.newaxis]

        dL = 0.5 * self.rho * U_total**2 * chord_grid * cl
        dD = 0.5 * self.rho * U_total**2 * chord_grid * cd

        # sectional force components
        dT = dL * np.cos(phi) - dD * np.sin(phi)   # thrust/normal
        dH = dL * np.sin(phi) + dD * np.cos(phi)   # in-plane axial (x)
        dY = dD * np.sin(PSI_grid)                 # approx lateral

        # CORRECT torque contribution: torque about shaft mainly from drag (profile)
        # dQ (per unit span per blade) => r * dD  (plus other contributions if modelled)
        dQ = R_grid * dD

        # Moments due to thrust at lever arm r (components)
        dMx = R_grid * dT * np.sin(PSI_grid)
        dMy = R_grid * dT * np.cos(PSI_grid)

        # Integrate using integrate_2d (which already accounts for dx spacing)
        T_blade = self.integrate_2d(r, psi, dT)
        H_blade = self.integrate_2d(r, psi, dH)
        Y_blade = self.integrate_2d(r, psi, dY)
        Q_blade = self.integrate_2d(r, psi, dQ)
        Mx_blade = self.integrate_2d(r, psi, dMx)
        My_blade = self.integrate_2d(r, psi, dMy)

        # total rotor quantities (multiply by blade number)
        forces_moments = {
            'FZ': self.b * T_blade,             # thrust (positive up)
            'FX': -self.b * H_blade,            # forward x (sign convention as before)
            'FY': -self.b * Y_blade,
            'MX': self.b * Mx_blade,
            'MY': self.b * My_blade,
            'MZ': self.b * Q_blade,
            'power': self.b * Q_blade * self.omega,
            'mu': mu,
            'lambda_G': lambda_G,
            # supply both mean scalar and distribution
            'lambda_iG': lambda_iG,
            'lambda_i_distribution': lambda_i,  # nondimensional lambda (n_r x n_psi)
            'beta0': beta0,
            'beta1c': beta1c,
            'beta1s': beta1s
        }

        # also return grids & local arrays for diagnostics
        return forces_moments, R_grid, PSI_grid, alpha, dT, dH, dY



    # def forward(self, forward_speed , theta0, theta1c , theta1s , beta0, alpha_tpp , psi, head_wind = 0.0,n=1000):
    #     """ 
    #     @brief head_wind is taken positive if it is opposite to the direction of forward speed
    #     @param forward_speed: is the speed of the helicopter in the forward direction
    #     @param head_wind: is taken if it is opposite to the forward speed
    #     @return 
    #     """
    #     forward_speed = forward_speed + head_wind # net forward speed 
    #     r = np.linspace(self.rc,self.r - 1e-6,n)
    #     twist = np.array([self.cal_twist(r[i]) for i in range(len(r))])
    #     theta = theta0 + theta1c * np.cos(psi) + theta1s * np.sin(psi)   # collective + cyclic + geometric
    #     mew  = forward_speed/(self.omega*self.r)
    #     sigma = np.array([self.cal_solidity(r[i]) for i in range(len(r))])

        


    #     # v_induced = np.zeros_like(r)
    #     # tol = 1e-4
    #     # for i in range(iter):
    #     #     phi = np.arctan((forward_speed + v_induced) / (self.omega * self.r))
    #     #     alpha_sectional = 
    #     #     v_induced_temp = 

