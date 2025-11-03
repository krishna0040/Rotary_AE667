"""
Simple Mission Planner Examples
"""
import numpy as np
from _ROTOR import rotor
from _HELICOPTER import helicopter
from _MISSION_PLANNER import MissionPlanner

# ==================== SETUP HELICOPTER ====================

main_rotor = rotor()
main_rotor.set_rotor(radius=4,
                     root_cut_out=0.5,
                     blades =3,
                     chord0=0.5,
                     chord_slope=(0.3-0.5)/(5-0.5), 
                     twist0=np.deg2rad(5),
                     twist_slope=np.deg2rad(4-5) / (3.5),
                     omega= 50)   ### approx 500RPM

tail_rotor = rotor()
tail_rotor.set_rotor(
    radius=0.8,
    root_cut_out=0.2,
    blades=4,
    chord0=0.18,                                  
    chord_slope=(0.08 - 0.18) / (0.6),    
    twist0=np.deg2rad(3),
    twist_slope=np.deg2rad(-3/0.6),                      
    omega=250.0     ### approz  2388RPM                           
)

heli = helicopter(main_rotor,tail_rotor,main_pos = np.array([-1,0,-2]),
                  tail_pos=np.array([-8,0,0]),lifting_surfaces=[],control_surfaces=[],
                  stabilizers=[],weight=6000)
heli.find_trim(velocity=40,alpha_tpp=np.deg2rad(5),Drag = 8000)





# ==================== MISSION 1: Simple Patrol ====================
mission1 = MissionPlanner(heli, engine_power_kw=2000, sfc=1.389e-4)
mission1.initialize(takeoff_weight=6000, fuel_weight=500, initial_altitude=0)

segments = [
    {'type': 'hover', 'duration_s': 120},
    {'type': 'climb', 'height_m': 1500, 'climb_rate': 5},
    {'type': 'forward', 'velocity': 50, 'distance_m': 80000},
    {'type': 'hover', 'duration_s': 300},
    {'type': 'forward', 'velocity': 50, 'distance_m': 80000},
]

mission1.execute_mission(segments)
mission1.plot_results("Simple Patrol Mission",segments)
