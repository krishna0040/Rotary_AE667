import numpy as np
import matplotlib.pyplot as plt

class MissionPlanner:
    def __init__(self, helicopter, engine_power_kw=1800, sfc=1.389e-4):
        self.helicopter = helicopter
        self.main_rotor = helicopter.main_rotor
        self.tail_rotor = helicopter.tail_rotor
        self.engine_power_kw = engine_power_kw
        self.sfc = sfc
        
        self.takeoff_weight = None
        self.fuel_weight = None
        self.current_weight = None
        self.altitude = 0

        self.reset_logs()

    # ============================================================
    # Initialization and utilities
    # ============================================================
    def reset_logs(self):
        self.time_log = []
        self.fuel_log = []
        self.weight_log = []
        self.alt_log = []
        self.power_log = []
        self.speed_log = []
        self.climb_log = []
        self.dist_log = []
        self.seg_log = []
        self.fuel_burn_log = []
        self._time = 0.0
        self._dist = 0.0

    def initialize(self, takeoff_weight, fuel_weight, initial_altitude=0):
        self.takeoff_weight = takeoff_weight
        self.fuel_weight = fuel_weight
        self.current_weight = takeoff_weight
        self.altitude = initial_altitude
        
        self.main_rotor.atmosphere(initial_altitude)
        self.tail_rotor.atmosphere(initial_altitude)
        
        self.reset_logs()
        print(f"\nMission Start: Weight={takeoff_weight:.1f} kg | Fuel={fuel_weight:.1f} kg | Alt={initial_altitude:.0f} m")

    def burn_fuel(self, power_kw, duration_s):
        fuel_used = power_kw * self.sfc * duration_s
        self.fuel_weight -= fuel_used
        self.current_weight -= fuel_used
        return fuel_used

    def log_segment(self, duration_s, power_kw, velocity=0.0, climb_rate=0.0, segment=""):
        """Log only one entry per segment."""
        self._time += duration_s
        self._dist += velocity * duration_s
        self.time_log.append(self._time)
        self.fuel_log.append(self.fuel_weight)
        self.weight_log.append(self.current_weight)
        self.alt_log.append(self.altitude)
        self.power_log.append(power_kw)
        self.speed_log.append(velocity)
        self.climb_log.append(climb_rate)
        self.dist_log.append(self._dist / 1000)  # in km
        self.seg_log.append(segment)
        self.fuel_burn_log.append(power_kw * self.sfc)  


    # ============================================================
    # Hover
    # ============================================================
    def hover(self, duration_s):
        print(f"\n--- Hover: {duration_s:.0f}s ---")
        _, _, _, _, thrust, _, _, _, torque = self.main_rotor.bemt(0.0)
        power_main = torque * self.main_rotor.omega / 1000
        power_tail = 0.1 * power_main
        power_total = power_main + power_tail
        
        print(f"  Power: {power_total:.1f} kW")

        if power_total > self.engine_power_kw:
            print(f"  ❌ ABORT: Power {power_total:.1f} > {self.engine_power_kw:.1f} kW")
            return False
        
        fuel_used = self.burn_fuel(power_total, duration_s)
        print(f"  Fuel used: {fuel_used:.2f} kg | Remaining: {self.fuel_weight:.1f} kg")
        self.log_segment(duration_s, power_total, segment="Hover")
        return self.fuel_weight > 0

    # ============================================================
    # Climb
    # ============================================================
    def climb(self, height_m, climb_rate):
        duration_s = height_m / climb_rate
        print(f"\n--- Climb: {height_m:.0f} m @ {climb_rate:.1f} m/s ---")

        self.altitude += height_m
        self.main_rotor.atmosphere(self.altitude)
        self.tail_rotor.atmosphere(self.altitude)

        _, _, _, _, thrust, _, _, _, torque = self.main_rotor.bemt(climb_rate)
        power_main = torque * self.main_rotor.omega / 1000
        power_tail = 0.1 * power_main
        power_total = power_main + power_tail

        print(f"  Power: {power_total:.1f} kW | New Alt: {self.altitude:.0f} m")

        if power_total > self.engine_power_kw:
            print(f"  ❌ ABORT: Power {power_total:.1f} > {self.engine_power_kw:.1f} kW")
            return False

        fuel_used = self.burn_fuel(power_total, duration_s)
        print(f"  Fuel used: {fuel_used:.2f} kg | Remaining: {self.fuel_weight:.1f} kg")
        self.log_segment(duration_s, power_total, climb_rate=climb_rate, segment="Climb")
        return self.fuel_weight > 0

    # ============================================================
    # Forward Flight
    # ============================================================
    def forward_flight(self, velocity, distance_m, alpha_tpp=0.05, drag_coeff=2.14):
        duration_s = distance_m / velocity
        print(f"\n--- Forward Flight: {distance_m/1000:.1f} km @ {velocity:.1f} m/s ---")

        self.helicopter.weight = self.current_weight * 9.81
        drag = drag_coeff * velocity**2
        trim = self.helicopter.find_trim(velocity=velocity, alpha_tpp=alpha_tpp, Drag=drag)

        theta0, theta1c, theta1s, theta0_tail = trim['theta0'], trim['theta1c'], trim['theta1s'], trim['theta0_tail']
        print(f"  Trim: θ0={theta0:.3f}, θ1c={theta1c:.3f}, θ1s={theta1s:.3f}, θ_tail={theta0_tail:.3f}")

        forces_moments, *_ = self.main_rotor.forward(
            forward_speed=velocity,
            theta0=theta0,
            theta1c=theta1c,
            theta1s=theta1s,
            alpha_tpp=alpha_tpp
        )

        MZ = forces_moments['MZ']
        power_main = abs(MZ) * self.main_rotor.omega / 1000.0
        power_tail = 0.1 * power_main
        power_total = power_main + power_tail

        print(f"  Power: {power_total:.1f} kW | Duration: {duration_s/60:.1f} min")

        if power_total > self.engine_power_kw:
            print(f"  ❌ ABORT: Power {power_total:.1f} > {self.engine_power_kw:.1f} kW")
            return False

        fuel_used = self.burn_fuel(power_total, duration_s)
        print(f"  Fuel used: {fuel_used:.2f} kg | Remaining: {self.fuel_weight:.1f} kg")

        self.log_segment(duration_s, power_total, velocity=velocity, segment="Cruise")
        return self.fuel_weight > 0

    # ============================================================
    # Mission Execution
    # ============================================================
    def execute_mission(self, segments):
        print(f"\n{'='*60}\nMISSION EXECUTION\n{'='*60}")
        for seg in segments:
            t = seg['type']
            if t == 'hover':
                success = self.hover(seg['duration_s'])
            elif t == 'climb':
                success = self.climb(seg['height_m'], seg['climb_rate'])
            elif t == 'forward':
                success = self.forward_flight(seg['velocity'], seg['distance_m'])
            else:
                print(f"Unknown segment type: {t}")
                continue
            if not success:
                print("\n❌ MISSION ABORTED\n" + "="*60)
                return
        print("\n✅ MISSION COMPLETED\n" + "="*60)
        print(f"Final Fuel: {self.fuel_weight:.1f} kg")

    # ============================================================
    # Plotting
    # ============================================================
    def plot_results(self, mission_name="Mission", mission_segments=None):
        """
        Plot each performance metric in a separate figure window.
        Add labeled markers A, B, C, ... for each segment boundary,
        and a summary box listing segment descriptions (e.g., A–B: Hover 120 s).
        """
        time_min = np.array(self.time_log) / 60

        # Compute segment durations to get boundary times
        seg_durations_s = []
        for seg in mission_segments:
            if seg["type"] == "hover":
                seg_durations_s.append(seg["duration_s"])
            elif seg["type"] == "climb":
                seg_durations_s.append(seg["height_m"] / seg["climb_rate"])
            elif seg["type"] == "forward":
                seg_durations_s.append(seg["distance_m"] / seg["velocity"])

        seg_boundaries = [0]
        for d in seg_durations_s:
            seg_boundaries.append(seg_boundaries[-1] + d / 60.0)  # in minutes
        seg_labels = [chr(65 + i) for i in range(len(seg_boundaries))]  # A, B, C, ...

        # Prepare segment descriptions for legend box
        seg_texts = []
        for i, seg in enumerate(mission_segments):
            start_label, end_label = seg_labels[i], seg_labels[i + 1]
            if seg["type"] == "hover":
                desc = f"{start_label}–{end_label}: Hover {seg['duration_s']} s"
            elif seg["type"] == "climb":
                desc = f"{start_label}–{end_label}: Climb {seg['height_m']} m @ {seg['climb_rate']} m/s"
            elif seg["type"] == "forward":
                desc = f"{start_label}–{end_label}: Forward {seg['distance_m']/1000:.0f} km @ {seg['velocity']} m/s"
            seg_texts.append(desc)

        # Quantities to plot
        plots = [
            ("Gross Weight [kg]", self.weight_log),
            ("Fuel Weight [kg]", self.fuel_log),
            ("Total Power Required [kW]", self.power_log),
            ("Speed [m/s]", self.speed_log),
            ("Altitude [m]", self.alt_log),
            ("Distance [km]", self.dist_log),
            ("Climb Rate [m/s]", self.climb_log),
            ("Fuel Burn Rate [kg/s]", self.fuel_burn_log)
        ]

        # Plot each quantity separately
        for title, data in plots:
            plt.figure(figsize=(8, 5))
            plt.plot(time_min, data, 'o-', lw=2, color='steelblue', markerfacecolor='orange')
            plt.title(f"{title}", fontsize=13, fontweight="bold")
            plt.xlabel("Time [min]")
            plt.grid(True)

            # ✅ Add horizontal line at 1500 kW only for the total power plot
            if "Power" in title:
                plt.axhline(1500, color='red', linestyle='--', linewidth=2, label='Engine Limit (1500 kW)')
                plt.legend(loc='best')


            # Add markers A, B, C... at segment boundaries
            y_min, y_max = plt.ylim()
            for i, t in enumerate(seg_boundaries):
                y_pos = y_min + 0.02 * (y_max - y_min)
                plt.text(t, y_pos, seg_labels[i], fontsize=10, fontweight='bold',
                         ha='center', va='bottom', color='darkred',
                         bbox=dict(boxstyle="circle,pad=0.3", facecolor='white', edgecolor='darkred'))

            # Add mission summary box in top-right
            box_text = f"{mission_name}\n" + "\n".join(seg_texts)
            plt.text(0.98, 0.98, box_text,
                     transform=plt.gca().transAxes,
                     fontsize=9, ha='right', va='top', color='black',
                     bbox=dict(boxstyle="round,pad=0.4", facecolor="whitesmoke", edgecolor="gray"))

            plt.tight_layout()
            plt.show()
