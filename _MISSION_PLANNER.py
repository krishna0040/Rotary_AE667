import numpy as np
class MissionPlanner:
    def __init__(self, rotor, engine_power_kw=1000, sfc=1.389e-4):
        """
        @param rotor: rotor object (defined by user)
        @param engine_power_kw: Maximum engine shaft power available [kW]
        @param sfc: Specific fuel consumption [kg/kW-s]
        """
        self.rotor = rotor
        self.engine_power_kw = engine_power_kw
        self.sfc = sfc
        self.mission_log = []
        
        # Mission state
        self.takeoff_weight = None
        self.fuel_weight = None
        self.altitude = None
        self.gross_weight = None

    def mission_inputs(self, takeoff_weight, fuel_weight, altitude):
        """
        Set initial mission conditions.
        """
        self.takeoff_weight = takeoff_weight
        self.fuel_weight = fuel_weight
        self.gross_weight = takeoff_weight
        self.altitude = altitude
        self.rotor.atmosphere(altitude)
        self.log_event(f"Mission initialized at {altitude} m with weight {takeoff_weight} kg")

    def log_event(self, message):
        self.mission_log.append(message)
        print(message)

    def power_required_hover(self, climb_rate=0.0):
        """
        Estimate hover/vertical climb power required using BEMT.
        """
        _, _, _, _, thrust, _, _, _, torque = self.rotor.bemt(climb_rate)
        weight_N = self.gross_weight * 9.81
        if thrust < weight_N:
            self.log_event("WARNING: Rotor thrust insufficient for hover at this weight!")
        power_req = torque * self.rotor.omega/1000
        return power_req

    def fuel_burn(self, power_kw, duration_s, eff=0.9):
        """
        Compute fuel burned using SFC in kg/kW-s.
        Accounts for mechanical/transmission efficiency.
        
        @param power_kw: shaft power required [kW]
        @param duration_s: segment duration [s]
        @param eff: efficiency (fraction of fuel power that reaches rotor), default = 0.9
        """
        # Effective fuel power must be higher due to losses
        fuel_used = (power_kw / eff) * self.sfc * duration_s
        self.fuel_weight -= fuel_used
        self.gross_weight -= fuel_used
        if self.fuel_weight < 0:
            self.log_event("WARNING: Fuel exhausted!")
        return fuel_used


    def segment_hover(self, duration_s=300, climb_rate=0.0):
        """
        Hover or vertical climb segment.
        @param duration_s: duration in seconds
        """
        P_req = self.power_required_hover(climb_rate)
        P_avail = self.engine_power_kw
        if P_req > P_avail:
            self.log_event("WARNING: Engine power insufficient for hover!")
        fuel_used = self.fuel_burn(P_req, duration_s)
        self.log_event(f"Hover/Climb {duration_s:.0f} s @ {climb_rate:.1f} m/s | "
                       f"Power req = {P_req:.1f} kW | Fuel used = {fuel_used:.4f} kg")

    def segment_payload(self, payload_change):
        """
        Pickup (+) or drop (-) payload.
        """
        self.gross_weight += payload_change
        self.log_event(f"Payload change: {payload_change:+.1f} kg | New gross weight = {self.gross_weight:.1f} kg")

    def report(self):
        """
        Print mission summary.
        """
        print("\n--- Mission Summary ---")
        for log in self.mission_log:
            print(log)

    def helicopter_performance(x, weight, fuel_weight, eta_p=0.8, P_avail=745.7*8000):
        """
        Compute stall-limited, power-limited, range-max, and endurance-max speeds and outputs.

        Parameters
        ----------
        x : rotor() instance
            Rotor class with thrust_power(V, alpha) defined.
        weight : float
            Gross weight of helicopter [N].
        fuel_weight : float
            Fuel weight [kg].
        eta_p : float
            Propulsive efficiency (default 0.8)
        P_avail : float
            Available engine power [W].

        Returns
        -------
        dict with:
        V_stall : m/s
        V_max_power : m/s
        range_max : km
        endurance_max : hr
        """
        rho = x.rho
        A = np.pi * x.radius**2

        alpha_stall = np.deg2rad(18)
        mu_values = np.linspace(0.01, 0.4, 100)
        V_values = mu_values * x.omega * x.radius

        alpha_values = []
        for mu in mu_values:
            _, _, _, _, _, _, _, alpha_tpp, _, _ = x.thrust_power(mu)
            alpha_values.append(alpha_tpp)
        alpha_values = np.array(alpha_values)

        # find where stall reached
        stall_idx = np.where(alpha_values >= alpha_stall)[0]
        if len(stall_idx) > 0:
            V_stall = V_values[stall_idx[0]]
        else:
            V_stall = V_values[-1]

        def power_required(V):
            T, P = x.thrust_power(V / (x.omega * x.radius))
            return P

        V_test = np.linspace(0, 120, 300)
        P_req = np.array([power_required(V) for V in V_test])
        idx = np.where(P_req < P_avail)[0]
        V_max_power = V_test[idx[-1]] if len(idx) > 0 else 0

        # Range ∝ η * (W_fuel/W) * (L/D)_max
        # Approx: (L/D)_max ≈ (T/W) at minimum power required
        i_min = np.argmin(P_req)
        T_min, P_min = x.thrust_power(V_test[i_min] / (x.omega * x.radius))
        LD_max = T_min / weight
        range_max = eta_p * (fuel_weight * 9.81 / weight) * LD_max * 3.6e3  # m → km

        # Endurance ∝ (η / P_min) * (W_fuel)
        endurance_max = eta_p * (fuel_weight * 9.81) / P_min / 3600  # sec → hr

        return {
            "V_stall": V_stall,
            "V_max_power": V_max_power,
            "range_max": range_max,
            "endurance_max": endurance_max
        }
