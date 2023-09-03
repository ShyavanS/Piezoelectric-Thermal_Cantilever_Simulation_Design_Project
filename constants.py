"""Global constants and variables to be held for other programs to access.

   Inputs:
   		 auto_sweep.py: Stores output of data in global variables

   Outputs:
   		 auto_sweep.py: Provides TEMPLATE_SCRIPT and FLEX_VERSION for use in running FlexPDE
		 solver.py: Provides data output from global variables and inital sweep parameters
"""

# Template FlexPDE script to be modified for sweep
TEMPLATE_SCRIPT = """
TITLE 'Design Project - Steady State'    												! The Problem Identification

COORDINATES cartesian2 ("x", "z") 														! Coordinate System, 1D,2D,3D, etc

VARIABLES
V(threshold=1e-6)																		! Voltage (V)
rho_free(threshold=1e-6)																! Free Charge Density (C/m^3)
Temp(threshold=1e-6)																	! Temperature (C)
u(threshold=1e-6)																		! Displacement in x (m)
w(threshold=1e-6)																		! Displacement in z (m)

SELECT
ngrid = 20     																			! Method Controls

DEFINITIONS
! Dynamic
alpha_T																					! Coefficient of Thermal Expansion (K^-1)
k																						! Thermal Conductivity (W/m K^-1)

epsilon_r																				! Relative Permittivity (unitless ratio)
sigma_e																					! Electrical Conductivity (S/m)

! Piezoelectric Coupling Coefficients (m/V)
d_11 d_13 d_15
d_31 d_33 d_35

! Stiffness Matrix (Pa^-1)
C_11 C_13 C_15
C_31 C_33 C_35
C_51 C_53 C_55

! Static
mag = 0.1*globalmax(magnitude(x, z))/globalmax(magnitude(u, w))							! Magnification factor for plot (unitless ratio)

epsilon_0 = 8.85e-12																	! Permittivity of Free Space (F/m)

h_air = 200																				! Convection Coefficient for Air (W/m^2 K^-1)

E_metal = {0}																			! Young's Modulus for Metal (Pa)
nu_metal = {1}																			! Poisson's Ratio for Metal (unitless ratio)

V_max = 1																				! Electric Potential over Piezoelectric (V)

l_x = 2e-3																				! Length of Piezoelectric (m)
t_piezo = 2e-4																			! Thickness of Piezoelectric (m)

Delta_T = Temp - 25																		! Temperature Difference (K)

T_inf = 25																				! Temperature of Surroundings (C)

! Sweep
t_metal_top = {2}																		! Thickness of Metal on Top (m)
t_metal_bottom = {3}																	! Thickness of Metal on Bottom (m)
t_metal_side = {4}																		! Thickness of Metal on Free Side (m)

! Equations
E_field = -grad(V)																		! Electric Field (V/m)
J = sigma_e*E_field																		! Current Density (A/m^2)
q_dot = -k*grad(Temp)																	! Heat Flux for Conduction (W/m^2)
q_dotc = h_air*(T_inf - Temp)															! Heat Flux for Convection (W/m^2)
q_dotvol = dot(J, E_field)																! Volumetric Heat Generation (W/m^3)

E_x = xcomp(E_field)																	! Electric Field in x (V/m)
E_z = zcomp(E_field)																	! Electric Field in z (V/m)

! Strain Definitions from Displacements (unitless ratio)
epsilon_x = dx(u)
epsilon_z = dz(w)
gamma_xz = dx(w) + dz(u)

! Mechanical Strain (unitless ratio)
epsilon_xm = epsilon_x - alpha_T*Delta_T - (d_11*E_x + d_31*E_z)
epsilon_zm = epsilon_z - alpha_T*Delta_T - (d_13*E_x + d_33*E_z)
gamma_xzm = gamma_xz - (d_15*E_x + d_35*E_z)

! Hookes Law for Stresses (Pa)
s_x = C_11*epsilon_xm + C_13*epsilon_zm + C_15*gamma_xzm
s_z = C_31*epsilon_xm + C_33*epsilon_zm + C_35*gamma_xzm
s_xz = C_51*epsilon_xm + C_53*epsilon_zm + C_55*gamma_xzm

! Electric Displacement Field due to Piezoelectric Effect (C/m^2)
D_piezox = d_11*s_x + d_13*s_z + d_15*s_xz
D_piezoz = d_31*s_x + d_33*s_z + d_35*s_xz

D_piezo = vector(D_piezox, 0, D_piezoz)

! Total Electric Displacement Field (C/m^2)
D_field = epsilon_0*epsilon_r*E_field + D_piezo

INITIAL VALUES
V = V_max/2																				! Voltage begins at half V_max
Temp = T_inf																			! Cantilever begins at same temperature as surroundings
rho_free = 0																			! No free charge on cantilever at start

EQUATIONS        																		! PDEs; one for each variable
V: div(J) = 0
rho_free: div(D_field) = rho_free
u: dx(s_x) + dz(s_xz) = 0
w: dx(s_xz) + dz(s_z) = 0
Temp: q_dotvol = div(q_dot)

BOUNDARIES
	REGION "Metal"
			alpha_T = {5}																! Coefficient of Thermal Expansion (K^-1)
			k = {6}																		! Thermal Conductivity (W/m K^-1)

			epsilon_r = 1e6																! Relative Permittivity (unitless ratio)
			sigma_e = 1/{7}																! Electrical Conductivity (S/m)

			! Piezoelectric Coupling Coefficients (m/V)
			d_11 = 0 d_13 = 0 d_15 = 0
			d_31 = 0 d_33 = 0 d_35 = 0

			! Stiffness Matrix (Pa^-1)
			C_11 = E_metal*(1 - nu_metal)/((1 + nu_metal)*(1 - 2*nu_metal)) C_13 =  E_metal*nu_metal/((1 + nu_metal)*(1 - 2*nu_metal)) C_15 = 0
			C_31 = E_metal*nu_metal/((1 + nu_metal)*(1 - 2*nu_metal))  		C_33 = C_11												   C_35 = 0
			C_51 = 0 														C_53 = 0												   C_55 = E_metal/(2*(1 + nu_metal))

			START (0, 0)
				load(Temp) = q_dotc														! Set up convection on all free sides
			LINE TO (l_x + t_metal_side, 0)
			LINE TO (l_x + t_metal_side, t_piezo + t_metal_top + t_metal_bottom)
			LINE TO (0, t_piezo + t_metal_top + t_metal_bottom)							! Hold one side fixed for cantilever
				value(V) = V_max														! Set V to V_max on top contact
				value(rho_free) = 0														! Set free charge density to zero on contact
				load(Temp) = 0															! Thermally insulate fixed side
				value(u) = 0
				value(w) = 0
			LINE TO (0, t_piezo/2 + t_metal_bottom)										! Hold one side fixed for cantilever
				value(V) = 0															! Set V to 0 on bottom contact
				value(rho_free) = 0														! Set free charge density to zero on contact
				load(Temp) = 0															! Thermally insulate fixed side
				value(u) = 0
				value(w) = 0
			LINE TO CLOSE

	REGION "Piezoelectric"
			alpha_T = 1.15e-6															! Coefficient of Thermal Expansion (K^-1)
			k = 7.7																		! Thermal Conductivity (W/m K^-1)

			epsilon_r = 4.68															! Relative Permittivity (unitless ratio)
			sigma_e = 1/10e5															! Electrical Conductivity (S/m)

			! Piezoelectric Coupling Coefficients (m/V)
			d_11 = 2e-11   d_13 = 0 	d_15 = 0
			d_31 = 4.5e-11 d_33 = 7e-11 d_35 = 0

			! Stiffness Matrix (Pa^-1)
			C_11 = 8.66766243465273e10 C_13 = 1.12023898431665e10 C_15 = 0
			C_31 = C_13 			   C_33 = 1.02688573562360e11 C_35 = 0
			C_51 = 0 				   C_53 = 0 				  C_55 = 5.78034682080925e10

			START (0, t_metal_bottom)
			LINE TO (l_x, t_metal_bottom)
			LINE TO (l_x, t_metal_bottom + t_piezo)
			LINE TO (0, t_metal_bottom + t_piezo)										! Hold one side fixed for cantilever
				load(Temp) = 0															! Thermally insulate fixed side
				value(u) = 0
				value(w) = 0
			LINE TO CLOSE

PLOTS          	  																		! Save result displays
	grid(x + mag*u, z + mag*w)
SUMMARY export file="output_{8}.txt"
	report val(w, l_x + t_metal_side, t_metal_bottom + t_piezo/2) as "Delta_tip "		! Export Transverse Tip Displacement (m)
	report val(Temp, l_x, t_metal_bottom + t_piezo/2) as "T_junc "						! Export Tip Junction Temperature (m)
END
"""

# Second Template FlexPDE script to be modified for sweep
SECONDARY_TEMPLATE_SCRIPT = """
TITLE 'Design Project - Resonance'    								! The Problem Identification

COORDINATES cartesian2 ("x", "z") 									! Coordinate System, 1D,2D,3D, etc

VARIABLES
u(threshold=1e-6)													! Displacement in x (m)
w(threshold=1e-6)													! Displacement in z (m)

SELECT    															! Method Controls
ngrid = 20
modes = 1

DEFINITIONS
! Dynamic
alpha_T																! Coefficient of Thermal Expansion (K^-1)

rho																	! Density (kg/m^3)

E_x																	! Electric Field in x (V/m)
E_z																	! Electric Field in z (V/m)

! Piezoelectric Coupling Coefficients (m/V)
d_11 d_13 d_15
d_31 d_33 d_35

! Stiffness Matrix (Pa^-1)
C_11 C_13 C_15
C_31 C_33 C_35
C_51 C_53 C_55

! Static
mag = 0.1*globalmax(magnitude(x, z))/globalmax(magnitude(u, w))		! Magnification factor for plot (unitless ratio)

epsilon_0 = 8.85e-12												! Permittivity of Free Space (F/m)

E_metal = {0}														! Young's Modulus for Metal (Pa)
nu_metal = {1}														! Poisson's Ratio for Metal (unitless ratio)

V_max = 1															! Electric Potential over Piezoelectric (V)

l_x = 2e-3															! Length of Piezoelectric (m)
t_piezo = 2e-4														! Thickness of Piezoelectric (m)

Delta_T = 25														! Temperature Difference (K)

! Sweep
t_metal_top = {2}													! Thickness of Metal on Top (m)
t_metal_bottom = {3}												! Thickness of Metal on Bottom (m)
t_metal_side = {4}													! Thickness of Metal on Free Side (m)

! Equations
! Strain Definitions from Displacements (unitless ratio)
epsilon_x = dx(u)
epsilon_z = dz(w)
gamma_xz = dx(w) + dz(u)

! Mechanical Strain (unitless ratio)
epsilon_xm = epsilon_x - alpha_T*Delta_T - (d_11*E_x + d_31*E_z)
epsilon_zm = epsilon_z - alpha_T*Delta_T - (d_13*E_x + d_33*E_z)
gamma_xzm = gamma_xz - (d_15*E_x + d_35*E_z)

! Hookes Law for Stresses (Pa)
s_x = C_11*epsilon_xm + C_13*epsilon_zm + C_15*gamma_xzm
s_z = C_31*epsilon_xm + C_33*epsilon_zm + C_35*gamma_xzm
s_xz = C_51*epsilon_xm + C_53*epsilon_zm + C_55*gamma_xzm

EQUATIONS        													! PDEs; one for each variable
u: dx(s_x) + dz(s_xz) = -rho*lambda*u
w: dx(s_xz) + dz(s_z) = -rho*lambda*w

BOUNDARIES
	REGION "Metal"
			alpha_T = {5}											! Coefficient of Thermal Expansion (K^-1)

			rho = {6}												! Density (kg/m^3)

			E_x = 0													! Electric Field in x (V/m)
			E_z = 0													! Electric Field in z (V/m)

			! Piezoelectric Coupling Coefficients (m/V)
			d_11 = 0 d_13 = 0 d_15 = 0
			d_31 = 0 d_33 = 0 d_35 = 0

			! Stiffness Matrix (Pa^-1)
			C_11 = E_metal*(1 - nu_metal)/((1 + nu_metal)*(1 - 2*nu_metal)) C_13 =  E_metal*nu_metal/((1 + nu_metal)*(1 - 2*nu_metal)) C_15 = 0
			C_31 = E_metal*nu_metal/((1 + nu_metal)*(1 - 2*nu_metal))  		C_33 = C_11												   C_35 = 0
			C_51 = 0 														C_53 = 0												   C_55 = E_metal/(2*(1 + nu_metal))

			START (0, 0)
			LINE TO (l_x + t_metal_side, 0)
			LINE TO (l_x + t_metal_side, t_piezo + t_metal_top + t_metal_bottom)
			LINE TO (0, t_piezo + t_metal_top + t_metal_bottom)		! Hold one side fixed for cantilever
				value(u) = 0
				value(w) = 0
			LINE TO CLOSE

	REGION "Piezoelectric"
			alpha_T = 1.15e-6										! Coefficient of Thermal Expansion (K^-1)

			rho = 2650												! Density (kg/m^3)

			E_x = 0													! Electric Field in x (V/m)
			E_z = V_max/t_piezo										! Electric Field in z (V/m)

			! Piezoelectric Coupling Coefficients (m/V)
			d_11 = 2e-11   d_13 = 0 	d_15 = 0
			d_31 = 4.5e-11 d_33 = 7e-11 d_35 = 0

			! Stiffness Matrix (Pa^-1)
			C_11 = 8.66766243465273e10 C_13 = 1.12023898431665e10 C_15 = 0
			C_31 = C_13 			   C_33 = 1.02688573562360e11 C_35 = 0
			C_51 = 0 				   C_53 = 0 				  C_55 = 5.78034682080925e10

			START (0, t_metal_bottom)
			LINE TO (l_x, t_metal_bottom)
			LINE TO (l_x, t_metal_bottom + t_piezo)
			LINE TO (0, t_metal_bottom + t_piezo)					! Hold one side fixed for cantilever
				value(u) = 0
				value(w) = 0
			LINE TO CLOSE

PLOTS          	  													! Save result displays
	grid(x + mag*u, z + mag*w)
SUMMARY export file="frequency.txt"
	report sqrt(lambda)/(2*pi) as "f "								! Export Resonant Frequency (Hz)
END
"""

# Minimum and maximum thickness to begin sweep range
MIN_THICK = 5e-5
MAX_THICK = 5e-4

# Initialize empty lists to store width range to sweep and results of sweep
sweep_range = []
displacement, temp_junc, t_metal_top, t_metal_bottom, t_metal_side = [], [], [], [], []

# Metal Layer Properties
E_LS = [120e9, 79e9, 100e9, 300e9]
NU_LS = [0.34, 0.42, 0.35, 0.49]
ALPHA_T_LS = [8.6e-6, 4.2e-5, 1.7e-5, 9.5e-6]
K_LS = [17, 310, 150, 54]
SIGMA_E_LS = [4.2e-7, 2.2e-8, 1.1e-7, 1.3e-7]
RHO_LS = [4500, 19300, 8500, 497]

# List of all metals
MATERIAL_LS = ["Ti", "Au", "Brass", "Vb"]

# Path where your FlexPDE executable is installed
# (or just name of executable if added to PATH correctly)
FLEX_VERSION = "FlexPDE6s"

# Ensuring the script is only imported into other programs as intended
if __name__ == "__main__":
    print("This script is only meant to be implemented as a module, not run independently.")
