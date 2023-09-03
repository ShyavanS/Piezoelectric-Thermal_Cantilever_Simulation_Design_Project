"""Main program to be run. Uses functions from other programs and
   solves for optimum in given problem.

   Inputs:
   		 constants.py: Retrieves output data from global variables
         auto_sweep.py: Uses functions to run FlexPDE scripts in a multithreaded manner
         advanced_plotter.py: Uses functions to plot data, fit curves, and identify optima

   Outputs:
   		 None
"""

# Imports
import constants as ct
from advanced_plotter import plotter_3d_contour
from auto_sweep import *

# Function Definitions


def main():
    """Main function of main program. Uses functions from other programs and
       solves for optimum in given problem.

       Args:
           None

       Returns:
           None
    """

    for (ind, material) in enumerate(ct.MATERIAL_LS):
        # Begin range at minimum possible and maximum possible thicknesses
        initialize_range([ct.MIN_THICK, ct.MAX_THICK], [
                         ct.MIN_THICK, ct.MAX_THICK], [ct.MIN_THICK, ct.MAX_THICK])

        # Sweep and organize data
        sweep(material)
        collate_data()

        # Initialize variables for figure plots
        fig_counter = ind + 1
        title = f"Displacement of a Piezoelectric Cantilever Varying {material} Thickness"
        x_label = f"Top {material} Thickness (\u03BCm)"
        y_label = f"Bottom {material} Thickness (\u03BCm)"
        z_label = f"Side {material} Thickness (\u03BCm)"
        c_label = "Tip Displacement (\u03BCm)"

        # Plot data and output maximum
        maximum = plotter_3d_contour(ct.t_metal_top, ct.t_metal_bottom, ct.t_metal_side,
                                     ct.displacement, fig_counter, title, x_label, y_label, z_label, c_label)

        # Open FlexPDE script and write formatted template
        with open("frequency.pde", "w") as stdin:
            print(ct.SECONDARY_TEMPLATE_SCRIPT.format(
                ct.E_LS[ind], ct.NU_LS[ind], maximum[0]*1e-6, maximum[1]*1e-6, maximum[2]*1e-6, ct.ALPHA_T_LS[ind], ct.RHO_LS[ind]), file=stdin)

        # Adjust timeout according to computing power
        # (if scripts return with timeout, increase timeout)
        completed = run([ct.FLEX_VERSION, "-S", "frequency.pde"], timeout=30)

        # Output runcode of file to confirm successful completion
        print(f"File frequency returned {completed.returncode}")

        with open(f"frequency_0.txt", "r") as stdin:
            freq = float(stdin.readlines()[5].split(" = ")[1])

        print(f"--{material} Maximum--")
        print(
            f"t_metal_top: {maximum[0]} (\u03BCm), t_metal_bottom: {maximum[1]} (\u03BCm), t_metal_side: {maximum[2]} (\u03BCm)")
        print(
            f"displacement: {maximum[3]} (\u03BCm), temp_junc: {ct.temp_junc[ct.displacement.index(maximum[3])]} (\u00b0C), res_freq: {freq} (Hz)")


# Run main() when run as main program; Not meant to be imported as a module
if __name__ == "__main__":
    main()
