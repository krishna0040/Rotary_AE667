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

    def bemt(self,V):       

        """
        @brief calculates everything - thrust, torque, lambda, induced velocity for hover or upward flight only
        @param V: climb velocity   [define the rotor parameters and all before]
        @return: lambdas,v,sigma,twist, thrust ,dT,Fs,d_torque, torque

        """
        r = np.linspace(self.rc,self.r - 1e-6,1000)
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
    
    