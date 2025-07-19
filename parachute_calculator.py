# Model Rocket Parachute Calculator by David Presker

import math
import sys

# Custom Tee class to print and write to file
class Tee:
    def __init__(self, filename, mode="w"):
        self.file = open(filename, mode, encoding="utf-8")
        self.stdout = sys.stdout  # original terminal output

    def write(self, message):
        self.stdout.write(message)
        self.file.write(message)

    def flush(self):
        self.stdout.flush()
        self.file.flush()

def compute_diameter_from_area(A, d_spill_percent, shape):
    if shape == 'round':
        return 2 * math.sqrt(A / ((1-(d_spill_percent/100)**2) * math.pi))
    elif shape == 'hex':
        return math.sqrt(2*A/(3**(1/3) - (d_spill_percent/100)**2 * math.pi/2))
    elif shape == 'octo':
        return math.sqrt(A / ((2 * (math.sqrt(2)-1)) - (d_spill_percent/100)**2 * math.pi / 4))
    elif shape == 'square':
    		return math.sqrt((A) / (1 - ((d_spill_percent/100)**2 * math.pi / 4)))
    

#def compute_diameter_from_area(A, shape):
#    if shape == 'round':
#        return 2 * math.sqrt(A / math.pi)
#    elif shape == 'hex':
#        return math.sqrt(2*A/3**(1/3))
#    elif shape == 'octo':
#        return math.sqrt(A/2/(math.sqrt(2)-1))
#    elif shape == 'square':
#    		return math.sqrt(A)





def compute_parachute_area(mass, gravity, velocity, air_density, Cd):
		return 2 * mass / 1000 * gravity / (velocity**2 * air_density * Cd)

def model_rocket_parachute_calculator():
    print("=== Model Rocket Parachute Size Calculator ===")
    print("===            by David Presker            ===\n")
    try:
        
        project_name = str(input("Project name (optional): "))        
        
        mass_input = input("Rocket mass (grams) [Required]: ").strip()
        if not mass_input:
            print("❌ Error: Rocket mass is required.")
            return
        mass = float(mass_input)

        valid_shapes = ['round', 'hex', 'octo', 'square']
        shape_input = input("Parachute shape (round/hex/octo/square) [Default: octo]: ").strip().lower()
        shape = shape_input if shape_input in valid_shapes else 'octo'

        if shape_input and shape_input not in valid_shapes:
            print(f"⚠️ Unknown shape '{shape_input}', defaulting to 'octo'.")

        d_spill_percent = input("Spill hole size as % of diameter (flat to flat for polygons) [Default: 20]: ")
        d_spill_percent = float(d_spill_percent) if d_spill_percent.strip() else 20.0

        


        velocity = input("Descent rate v (m/s) [Default: 4.5]: ")
        velocity = float(velocity) if velocity.strip() else 4.5

        air_density = input("Air density ρ (kg/m³) [Default: 1.225]: ")
        air_density = float(air_density) if air_density.strip() else 1.225

        Cd = input("Drag coefficient Cd [Default: 0.75]: ")
        Cd = float(Cd) if Cd.strip() else 0.75

        gravity = input("Gravity g (m/s²) [Default: 9.8067]: ")
        gravity = float(gravity) if gravity.strip() else 9.8067

        # --- Calculations ---
        area_parachute = compute_parachute_area(mass, gravity, velocity, air_density, Cd)      
        d_no_spill_hole=compute_diameter_from_area(area_parachute, 0, shape) 
        d_final=compute_diameter_from_area(area_parachute, d_spill_percent, shape) 
        d_spill_hole=d_final*d_spill_percent/100        
        area_spill_hole=(d_final * (d_spill_percent / 100) /2)**2 * math.pi        
        area_ratio_spill_surface=area_spill_hole/area_parachute
        area_over_all=area_parachute+area_spill_hole
        area_ratio_spill_over_all=area_spill_hole/area_over_all
        
        # --- Print and Save Results ---
        print()
        tee = Tee("parachute.txt")        
        print(f"=== ✅ Results of Parachute Calculation, by David Presker  ===\n", file=tee)
        if project_name.strip():        
           print(f"Project name: {project_name}", file=tee)
        print(f"Parachute surface area: {(area_parachute*10000):.2f} cm²", file=tee)
        print(f"Parachute diameter without Spill hole: {(d_no_spill_hole*100):.1f} cm", file=tee)        
        print(f"Parachute diameter with Spill hole: {(d_final*100):.1f} cm", file=tee)
        print(f"Spill hole diameter: {(d_spill_hole*100):.1f} cm", file=tee)        
        print(f"Spill Hole Area: {(area_spill_hole*10000):.2f} cm²", file=tee)
        print(f"Parachute area over all: {(area_over_all*10000):.2f} cm²", file=tee)        
        print(f"Area ratio Spill hole/Surface area: {(area_ratio_spill_surface*100):.2f} %", file=tee)
        print(f"Area ratio Spill hole/Area over all: {(area_ratio_spill_over_all*100):.2f} %", file=tee)
        print(f"\nInputs → Rocket mass: {mass_input} g, Shape: {shape}, Spill Hole: {d_spill_percent} %, Cd: {Cd}, ρ: {air_density} kg/m³, g: {gravity} m/s², v: {velocity} m/s", file=tee)
        tee.flush()

        tee.file.close()

        print("✅ Results also saved to 'parachute.txt'")

    except ValueError:
        print("❌ Invalid input. Please enter valid numbers.")

# Run the calculator
model_rocket_parachute_calculator()

