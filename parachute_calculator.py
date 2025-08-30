# Model Rocket Parachute Calculator by David Presker

# Model Rocket Parachute Calculator + PDF Template
import math, sys, re
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ------------------------
# Tee class for TXT output
# ------------------------
class Tee:
    def __init__(self, filename, mode="w"):
        self.file = open(filename, mode, encoding="utf-8")
        self.stdout = sys.stdout
    def write(self, message):
        self.stdout.write(message)
        self.file.write(message)
    def flush(self):
        self.stdout.flush()
        self.file.flush()

# ------------------------
# Calculations
# ------------------------
def compute_diameter_from_area(A, d_spill_percent, shape):
    if shape == 'round':
        return 2 * math.sqrt(A / ((1-(d_spill_percent/100)**2) * math.pi))
    elif shape == 'hex':
        return math.sqrt(2*A/(3**(1/3) - (d_spill_percent/100)**2 * math.pi/2))
    elif shape == 'octo':
        return math.sqrt(A / ((2 * (math.sqrt(2)-1)) - (d_spill_percent/100)**2 * math.pi / 4))
    elif shape == 'square':
        return math.sqrt(A / (1 - ((d_spill_percent/100)**2 * math.pi / 4)))

def compute_parachute_area(mass, gravity, velocity, air_density, Cd):
    return 2 * mass / 1000 * gravity / (velocity**2 * air_density * Cd)

# ------------------------
# Page sizes (mm)
# ------------------------
A_SERIES = {
    "A0": (841, 1189),
    "A1": (594, 841),
    "A2": (420, 594),
    "A3": (297, 420),
    "A4": (210, 297),
}

# ------------------------
# Shape vertices
# ------------------------
def get_vertices(shape, D_mm):
    if shape=="round":
        n_sides = 8
        R = D_mm/2
        verts = [(R*math.cos(2*math.pi*i/n_sides), R*math.sin(2*math.pi*i/n_sides)) for i in range(n_sides)]
    else:
        if shape=="square": 
            n_sides = 4
            rot = math.pi/4  # rotate 45° so flat-to-flat matches page edges
        elif shape=="hex": 
            n_sides = 6
            rot = math.pi/6  # rotate 30° so flats are horizontal
        elif shape=="octo": 
            n_sides = 8
            rot = math.pi/8
        R = D_mm/2/math.cos(math.pi/n_sides)
        verts = [(R*math.cos(2*math.pi*i/n_sides + rot), R*math.sin(2*math.pi*i/n_sides + rot)) for i in range(n_sides)]
    return verts, n_sides

# ------------------------
# Draw parachute
# ------------------------
def draw_parachute(ax, shape, D_mm, spill_mm, project_name, rocket_mass, page_w, page_h, page_size_name, k_factor=0.01, draw_cutout=True, author_name="David Presker", border_mm=5):
    verts, n_sides = get_vertices(shape, D_mm)
    
    # Main shape
    if shape=="round":
        R = D_mm/2
        ax.add_patch(patches.Circle((0,0), R, fill=False, lw=1))
    else:
        ax.add_patch(patches.Polygon(verts, closed=True, fill=False, lw=1))
    
    # Border rectangle around parachute
    max_x = max(x for x,_ in verts)
    min_x = min(x for x,_ in verts)
    max_y = max(y for _,y in verts)
    min_y = min(y for _,y in verts)
    ax.add_patch(patches.Rectangle((min_x-border_mm, min_y-border_mm),
                                   (max_x-min_x)+2*border_mm,
                                   (max_y-min_y)+2*border_mm,
                                   fill=False, lw=0.7, linestyle='dotted', edgecolor='gray'))
    
    # Helping lines
    margin = max(page_w,page_h)*2
    for i,(x,y) in enumerate(verts):
        length = math.hypot(x,y)
        if length>0:
            ux,uy = x/length, y/length
            ax.plot([-ux*margin, ux*margin], [-uy*margin, uy*margin], 'k--', lw=0.3)
        x2,y2 = verts[(i+1)%n_sides]
        xm,ym = (x+x2)/2,(y+y2)/2
        length = math.hypot(xm,ym)
        if length>0:
            ux,uy = xm/length, ym/length
            ax.plot([-ux*margin, ux*margin], [-uy*margin, uy*margin], 'k--', lw=0.3)

    # Cut-out reinforcement and spill hole
    if draw_cutout:
        small_diameter = 6.5
        outer_diameter = 13
        R_out = max(math.hypot(x,y) for x,y in verts)
        for x, y in verts:
            length = math.hypot(x, y)
            if length == 0: continue
            ux, uy = x/length, y/length
            R_center = R_out - outer_diameter/2 - k_factor*D_mm
            cx, cy = ux*R_center, uy*R_center
            ax.add_patch(patches.Circle((cx,cy), small_diameter/2, fill=False, color='blue', lw=0.6))
            ax.add_patch(patches.Circle((cx,cy), outer_diameter/2, fill=False, color='red', lw=0.6))
        ax.add_patch(patches.Circle((0,0), spill_mm/2, fill=False, color='gray', lw=1))

    # Labels – move left with safe margin
    left_margin = -page_w/2 + 15  # 15 mm inside left page
    top_y = page_h/2 - 15
    line_spacing = 15

    if project_name.strip():
        ax.text(left_margin, top_y, f"Project: {project_name}", ha='left', va='top', fontsize=10, fontweight='bold')
    if rocket_mass:
        ax.text(left_margin, top_y - line_spacing, f"Rocket mass: {rocket_mass} g", ha='left', va='top', fontsize=10, fontweight='normal')
    ax.text(left_margin, -page_h/2 + 30, f"Diameter (/flat-to-flat): {D_mm:.1f} mm", ha='left', va='bottom', fontsize=8)
    ax.text(left_margin, -page_h/2 + 15, f"Spill hole diameter: {spill_mm:.1f} mm", ha='left', va='bottom', fontsize=8)

    # Format in top right corner
    right_margin = page_w/2 - 15
    ax.text(right_margin, page_h/2 - 15, f"Format: {page_size_name}", ha='right', va='top', fontsize=10, fontweight='bold')

    # Author just above bottom-right
    ax.text(right_margin, -page_h/2 + 20, author_name, ha='right', va='bottom', fontsize=8, fontweight='normal')

# ------------------------
# Split PDF vector
# ------------------------
def split_pdf_vector_full(original_pdf_name, page_w, page_h):
    print(f"📄 Original page is ready: {original_pdf_name}")
    print(f"📌 For large prints, use pdfposter or similar for vector tile splitting.")

# ------------------------
# Default k per shape
# ------------------------
DEFAULT_K = {"round":0.007,"square":0.018,"hex":0.012,"octo":0.01}

# ------------------------
# Main calculator
# ------------------------
def model_rocket_parachute_calculator():
    print("=== Model Rocket Parachute Calculator ===")
    print("===            by David Presker            ===\n")
    try:
        # Inputs
        project_name = input("Project name (optional): ").strip()
        mass_input = input("Rocket mass (grams) [Required]: ").strip()
        if not mass_input: print("❌ Rocket mass required"); return
        mass = float(mass_input)

        # Dual input for shapes
        SHAPE_MAP = {"0": "round", "4": "square", "6": "hex", "8": "octo"}
        VALID_SHAPES = ["round", "square", "hex", "octo"]
        shape_input = input("Parachute shape: 0=round, 4=square, 6=hex, 8=octo [Default:8]: ").strip().lower()
        if shape_input in VALID_SHAPES:
            shape = shape_input
        else:
            shape = SHAPE_MAP.get(shape_input, "octo")  # default octo

        d_spill_percent = float(input("Spill hole size % of diameter [Default:20]: ").strip() or 20.0)
        velocity = float(input("Descent rate v (m/s) [Default:4.5]: ") or 4.5)
        air_density = float(input("Air density ρ (kg/m³) [Default:1.225]: ") or 1.225)
        Cd = float(input("Drag coefficient Cd [Default:0.75]: ") or 0.75)
        gravity = float(input("Gravity g (m/s²) [Default:9.8067]: ") or 9.8067)
        draw_cutout = input("Draw cut-out? (y/n) [Default:y]: ").strip().lower() != 'n'
        default_k = DEFAULT_K.get(shape,0.01)
        k_input = input(f"Reinforcement offset factor k (0–0.05) [Default:{default_k:.3f}]: ")
        k_factor = float(k_input) if k_input.strip() else default_k

        # Calculations
        area = compute_parachute_area(mass, gravity, velocity, air_density, Cd)
        d_no_spill = compute_diameter_from_area(area, 0, shape)
        d_final = compute_diameter_from_area(area, d_spill_percent, shape)
        d_spill = d_final*d_spill_percent/100
        area_spill = math.pi*(d_spill/2)**2
        area_total = area + area_spill
        ratio1 = area_spill/area*100
        ratio2 = area_spill/area_total*100

        # TXT filename
        fname_safe = re.sub(r'[^\w\d-]', '_', project_name.strip()) if project_name.strip() else "parachute"
        txt_filename = f"parachute_{fname_safe}.txt"
        tee = Tee(txt_filename)
        print("=== ✅ Results of Parachute Calculation, by David Presker ===\n", file=tee)
        if project_name: print(f"Project name: {project_name}", file=tee)
        print(f"Parachute surface area: {area*10000:.2f} cm²", file=tee)
        print(f"Parachute diameter without Spill hole: {d_no_spill*100:.1f} cm", file=tee)
        print(f"Parachute diameter with Spill hole: {d_final*100:.1f} cm", file=tee)
        print(f"Spill hole diameter: {d_spill*100:.1f} cm", file=tee)
        print(f"Spill Hole Area: {area_spill*10000:.2f} cm²", file=tee)
        print(f"Parachute area over all: {area_total*10000:.2f} cm²", file=tee)
        print(f"Area ratio Spill hole/Surface area: {ratio1:.2f} %", file=tee)
        print(f"Area ratio Spill hole/Area over all: {ratio2:.2f} %\n", file=tee)
        print(f"Inputs → Rocket mass: {mass_input} g, Shape: {shape}, Spill Hole: {d_spill_percent} %, Cd: {Cd}, ρ: {air_density} kg/m³, g: {gravity} m/s², v: {velocity} m/s", file=tee)
        tee.flush(); tee.file.close()
        print(f"✅ Results saved to '{txt_filename}'")

        # Original page size
        D_mm = d_final*1000
        margin=10
        for size_name,(w,h) in sorted(A_SERIES.items(), key=lambda x:x[1][0]*x[1][1]):
            if D_mm+2*margin <= w and D_mm+2*margin <= h:
                page_w,page_h = w,h
                page_size_name=size_name
                break
        else: page_w,page_h=D_mm+2*margin,D_mm+2*margin; page_size_name=f"{page_w:.0f}x{page_h:.0f}mm"

        # Draw PDF
        fig, ax = plt.subplots(figsize=(page_w/25.4,page_h/25.4))
        ax.set_xlim(-page_w/2,page_w/2); ax.set_ylim(-page_h/2,page_h/2)
        ax.set_aspect('equal'); ax.axis('off')
        draw_parachute(ax, shape, d_final*1000, d_spill*1000, project_name, mass_input, page_w, page_h, page_size_name, k_factor, draw_cutout)
        pdf_filename = f"parachute_{fname_safe}_{page_size_name}.pdf"
        fig.savefig(pdf_filename, bbox_inches='tight', pad_inches=0)
        plt.close(fig)
        print(f"✅ Original page saved as {pdf_filename}")
        split_pdf_vector_full(pdf_filename, page_w, page_h)

    except ValueError:
        print("❌ Invalid input, please enter numbers correctly.")

if __name__=="__main__":
    model_rocket_parachute_calculator()
