Open the ipynb notebook in vscode, or upload on google colab.

Performance estimator tool using BEMT,          Rotor():-
- To give the inputs, first define an object of type rotor :- temp = rotor()
- Set the height of mission, and the ambient conditions are set according to the formulae for standard ISA temperature measurements at lower altitudes.
- Then define the parameters using the temp.set_airfoil() for setting the cl and cd for the airfoil. Cl = cl*alpha, alpha in rad, and cl in rad-1. cd = cdo + cd*alpha*alpha, here cdo is dimensionless, and cd is rad-2.
- Define the geometry of the airfoil using the temp.set_rotor(). Some default values are already set according to the paper. 
- Call the temp.bemt(), with the specified climb velocity V, which returns the values :- lambdas,v,sigma,twist, thrust ,dT,Fs,d_torque, torque. This discretizes the domain along the radius, and solves for lambda iteratively. 

Mission Planner tool , MissionPlanner():- 
- Define the class using mission = MisisonPlanner(), and give it the rotor, max power by the engine and the sfc of the helicopter. 
- The defined rotor is passed to the mission class, and the mission inputs such as takeoff_weight, fuel_weight, and altitude are defined. 
- Then give the inputs to the missionplanner using mission.segment_hover(), which takes inputs the duration and the upward velocity. This calculates if the thrust generated is enough and the power required can be provided or not. It also updates the remaining fuel in the tank. 
- The mission segments are stored in results, which can be obtained by calling mission.report() - this print the entire mission history. 


The rotor class is used for both the main rotor and the tail rotor, with the following specs:-

| Parameter                        | Rotor 1 (Main Rotor)             | Rotor 2 (Tail Rotor)           |
|----------------------------------|----------------------------------|--------------------------------|
| **Rotor Description (role)**     | Main Rotor                       | Tail Rotor                     |
| **Airfoil**                      | NACA 0012                        | NACA 0010                      |
| **Rotor Radius (m)**             | 5.5 m                            | 1 m                            |
| **Rotor Speed (rad/s)**          | 40 rad/s                         | 150 rad/s                      |
| **Number of Blades**             | 3                                | 2                              |
| **Chord Length Variation**       | 0.45 − 0.05 (r − 0.8)            | 0.18 − 0.09 (r − 0.2)          |
| **Twist Variation**              | 12 − 1.6 (r − 0.8)               | 10 − 6.4 (r − 0.2)             |
| **Root Cutout**                  | 0.8 m                            | 0.2 m                          
