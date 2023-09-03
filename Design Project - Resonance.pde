TITLE 'Design Project - Resonance'    																! The Problem Identification

COORDINATES cartesian2 ("x", "z") 																! Coordinate System, 1D,2D,3D, etc

VARIABLES
u(threshold=1e-6)																									! Displacement in x (m)
w(threshold=1e-6)																									! Displacement in z (m)

SELECT    																													! Method Controls
ngrid = 20
modes = 1

DEFINITIONS
! Dynamic
alpha_T																														! Coefficient of Thermal Expansion (K^-1)

rho																																! Density (kg/m^3)

E_x																																! Electric Field in x (V/m)
E_z																																! Electric Field in z (V/m)

! Piezoelectric Coupling Coefficients (m/V)
d_11 d_13 d_15
d_31 d_33 d_35

! Stiffness Matrix (Pa^-1)
C_11 C_13 C_15
C_31 C_33 C_35
C_51 C_53 C_55

! Static
mag = 0.1*globalmax(magnitude(x, z))/globalmax(magnitude(u, w))		! Magnification factor for plot (unitless ratio)

epsilon_0 = 8.85e-12																								! Permittivity of Free Space (F/m)

E_metal = 120e9																										! Young's Modulus for Metal (Pa)
nu_metal = 0.34																										! Poisson's Ratio for Metal (unitless ratio)

V_max = 1																												! Electric Potential over Piezoelectric (V)

l_x = 2e-3																													! Length of Piezoelectric (m)
t_piezo = 2e-4																											! Thickness of Piezoelectric (m)

Delta_T = 25																												! Temperature Difference (K)

! Sweep
t_metal_top = 1e-4																									! Thickness of Metal on Top (m)
t_metal_bottom = 2e-4																							! Thickness of Metal on Bottom (m)
t_metal_side = 1e-4																								! Thickness of Metal on Free Side (m)

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

EQUATIONS        																										! PDEs; one for each variable
u: dx(s_x) + dz(s_xz) = -rho*lambda*u
w: dx(s_xz) + dz(s_z) = -rho*lambda*w

BOUNDARIES
	REGION "Metal"
			alpha_T = 8.6e-6																							! Coefficient of Thermal Expansion (K^-1)

			rho = 4500																									! Density (kg/m^3)

			E_x = 0																											! Electric Field in x (V/m)
			E_z = 0																											! Electric Field in z (V/m)

			! Piezoelectric Coupling Coefficients (m/V)
			d_11 = 0 d_13 = 0 d_15 = 0
			d_31 = 0 d_33 = 0 d_35 = 0

			! Stiffness Matrix (Pa^-1)
			C_11 = E_metal*(1 - nu_metal)/((1 + nu_metal)*(1 - 2*nu_metal)) C_13 =  E_metal*nu_metal/((1 + nu_metal)*(1 - 2*nu_metal)) C_15 = 0
			C_31 = E_metal*nu_metal/((1 + nu_metal)*(1 - 2*nu_metal))  		C_33 = C_11																							 C_35 = 0
			C_51 = 0 																										C_53 = 0																								 C_55 = E_metal/(2*(1 + nu_metal))

			START (0, 0)
			LINE TO (l_x + t_metal_side, 0)
			LINE TO (l_x + t_metal_side, t_piezo + t_metal_top + t_metal_bottom)
			LINE TO (0, t_piezo + t_metal_top + t_metal_bottom)						! Hold one side fixed for cantilever
				value(u) = 0
				value(w) = 0
			LINE TO CLOSE

	REGION "Piezoelectric"
			alpha_T = 1.15e-6																						! Coefficient of Thermal Expansion (K^-1)

			rho = 2650																									! Density (kg/m^3)

			E_x = 0																											! Electric Field in x (V/m)
			E_z = V_max/t_piezo																					! Electric Field in z (V/m)

			! Piezoelectric Coupling Coefficients (m/V)
			d_11 = 2e-11 	d_13 = 0 		 d_15 = 0
			d_31 = 4.5e-11 d_33 = 7e-11 d_35 = 0

			! Stiffness Matrix (Pa^-1)
			C_11 = 8.66766243465273e10 C_13 = 1.12023898431665e10 C_15 = 0
			C_31 = C_13 				    			  C_33 = 1.02688573562360e11 C_35 = 0
			C_51 = 0 										  C_53 = 0 						   	 				C_55 = 5.78034682080925e10

			START (0, t_metal_bottom)
			LINE TO (l_x, t_metal_bottom)
			LINE TO (l_x, t_metal_bottom + t_piezo)
			LINE TO (0, t_metal_bottom + t_piezo)													! Hold one side fixed for cantilever
				value(u) = 0
				value(w) = 0
			LINE TO CLOSE

PLOTS          	  																											! Save result displays
	grid(x + mag*u, z + mag*w)
SUMMARY export file="stdout.txt"
	report sqrt(lambda)/(2*pi) as "f "																	! Export Resonant Frequency (Hz)
END
