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
