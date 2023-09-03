"""Program contaning all necessary functions to plot figures and
   somewhat automate that process. Held for other programs to access.

   Inputs:
         solver.py: Plotting functions are called from within solver.py with any parameters necessary passed to them
         
   Outputs:
   		 solver.py: Outputs optimum values from plots to solver.py to be printed to user
         figure_*.png: Saves plots of output data as figures in png format
"""

# Imports
import numpy as np
from matplotlib import pyplot as plt
from typing import List
from auto_sweep import *

# Function Definitions


def adj_r(x: List[float], y: List[float], degree: int) -> float:
    """Takes a set of data and tries to fit a polynomial curve to it.
       Outputs the adjusted R Squared value of the fit to see how well
       the curve fits the data.

       Args:
           x (List[float]): The x-axis data of the set to be plotted
           y (List[float]): The y-axis data of the set to be plotted
           degree (int): Degree of polynomial to try and fit to data

       Returns:
           r_squared (float): The adjusted R Squared value of the fit
                              representing how well the curve fits the data
    """

    # Fit polynomial of provided degree to data
    coeffs = np.polyfit(x, y, degree)
    p = np.poly1d(coeffs)

    # Get y hat and y bar from curve fit to x-axis data and given y-axis
    # data resepectively to see what the curve predicts should be the
    # y-values in comparison to the actual y-values
    y_hat = p(x)
    y_bar = np.sum(y)/len(y)

    # Calculate Sum of Squares Regression and Total Sum of Squares
    ssreg = np.sum((y_hat - y_bar)**2)
    sstot = np.sum((y - y_bar)**2)

    # Calculate adjusted R Squared value from previous computations
    r_squared = 1 - (((1 - (ssreg/sstot))*(len(y) - 1))/(len(y) - degree - 1))

    return r_squared


def plotter_2d(x: List[float], y: List[float], count: int, title: str, x_label: str, y_label: str, maximize=True, max_degree=5) -> List[float]:
    """Takes a set of data and plots a 2D figure with it given chart titles,
       a figure number, and other optional parameters. Returns the
       optimum point in the data after fitting to a curve. Saves the
       resulting figure as a png image in the working directory.

       Args:
           x (List[float]): The x-axis data of the set to be plotted
           y (List[float]): The y-axis data of the set to be plotted
           count (int): The figure number to output in the filename
           title (str): The figure title
           x_label (str): The x-axis data label for the figure
           y_label (str): The y-axis data label for the figure
           maximize (bool): Default value is True; Determines whether
                            optimum value to be searched for is a minimum
                            or maximum of the given data; True yields
                            a maximum, False yields a minimum
           max_degree (int): Default value is 5; Maximum degree of
                             polynomial to try and fit to the data

       Returns:
           optimum (List[float]): The optimum datapoint that matches
                                  the requested criteria based on the curve
                                  that best fits the data
    """

    # Calculate adjusted R Squared value for data using degrees
    # of polynomials until the maximum degree specified
    r_squared = [adj_r(x, y, i)
                 for i in range(1, max_degree + 1)]

    # Use the degree that best fit the data going forward
    degree = r_squared.index(max(r_squared)) + 1

    # Fit a model to the data with the best degree polynomial
    model = np.poly1d(np.polyfit(x, y, degree))
    polyline = np.linspace(min(x), max(x), 100)

    # Find critical points of curve fit using first derivative test
    bounds = [min(x), max(x)]
    crit_x = bounds + [i for i in model.deriv().r if i.imag ==
                       0 and bounds[0] < i.real < bounds[1]]
    crit_pts = zip(model(crit_x), crit_x)

    # Depending on the requested optimum, find
    # the maimum or minimum of the critical points
    if maximize == True:
        optimum = max(crit_pts)
    else:
        optimum = min(crit_pts)

    # Re-organize the optimal point so it is
    # in the format (x, y) and is entirely real
    optimum = [i.real for i in optimum][::-1]

    # Initialize figure
    fig = plt.figure()

    # Plot FlexPDE data as a scatter in default blue
    plt.scatter(x, y, marker="o")

    # Plot polynomial as a curve in red
    plt.plot(polyline, model(polyline), color="red")

    # Mark optimum point with green triangle
    plt.scatter(optimum[0], optimum[1], color="green", marker="^")

    # Add legend and chart titles
    plt.legend(["FlexPDE Data", "Polynomial Fit", "Maximum"], loc='best')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Set up chart grid
    plt.grid(visible=True, which="major")
    plt.minorticks_on()
    plt.grid(visible=True, which="minor", linewidth=0.25)

    # Save figure and clear plotter for next figure
    plt.savefig(f"figure_{count}.png", format="png", dpi=500)
    fig.clear()
    plt.close(fig)

    return optimum


def plotter_3d(x: List[float], y: List[float], z: List[float], count: int, title: str, x_label: str, y_label: str, z_label: str, maximize=True) -> List[float]:
    """Takes a set of data and plots a 3D figure with it given chart titles,
       a figure number, and other optional parameters. Returns the
       optimum point in the data. Saves the resulting figure
       as a png image in the working directory.

       Args:
           x (List[float]): The x-axis data of the set to be plotted
           y (List[float]): The y-axis data of the set to be plotted
           z (List[float]): The z-axis data of the set to be plotted
           count (int): The figure number to output in the filename
           title (str): The figure title
           x_label (str): The x-axis data label for the figure
           y_label (str): The y-axis data label for the figure
           z_label (str): The z-axis data label for the figure
           maximize (bool): Default value is True; Determines whether
                            optimum value to be searched for is a minimum
                            or maximum of the given data; True yields
                            a maximum, False yields a minimum

       Returns:
           optimum (List[float]): The optimum datapoint that matches
                                  the requested criteria
    """

    # Initialize figure
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    # Depending on the requested optimum, find
    # the maimum or minimum of the points
    if maximize == True:
        optimum = [x.pop(z.index(max(z))), y.pop(
            z.index(max(z))), z.pop(z.index(max(z)))]
    else:
        optimum = [x.pop(z.index(min(z))), y.pop(
            z.index(min(z))), z.pop(z.index(min(z)))]

    # Plot FlexPDE data as a scatter in default blue
    ax.scatter(x, y, z, marker="o")

    # Mark optimum point with green triangle
    ax.scatter(optimum[0], optimum[1], optimum[2], color="green", marker="^")

    # Add legend and chart titles
    ax.legend(["FlexPDE Data", "Maximum"], loc='best')
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)

    # Set up chart grid
    ax.grid(visible=True, which="major")
    ax.minorticks_on()
    ax.grid(visible=True, which="minor", linewidth=0.25)

    # Save figure and clear plotter for next figure
    plt.savefig(f"figure_{count}.png", format="png", dpi=500)
    ax.clear()
    plt.close(fig)

    return optimum


def plotter_3d_contour(x: List[float], y: List[float], z: List[float], c: List[float], count: int, title: str, x_label: str, y_label: str, z_label: str, c_label: str, maximize=True) -> List[float]:
    """Takes a set of data and plots a 3D figure with it with colour
       representing a fourth dataset given chart titles, a figure number,
       and other optional parameters. Returns the optimum point in the data.
       Saves the resulting figure as a png image in the working directory.

       Args:
           x (List[float]): The x-axis data of the set to be plotted
           y (List[float]): The y-axis data of the set to be plotted
           z (List[float]): The z-axis data of the set to be plotted
           c (List[float]): The data to be plotted over the colour gradient
           count (int): The figure number to output in the filename
           title (str): The figure title
           x_label (str): The x-axis data label for the figure
           y_label (str): The y-axis data label for the figure
           z_label (str): The z-axis data label for the figure
           c_label (str): The colour gradient data label for the figure
           maximize (bool): Default value is True; Determines whether
                            optimum value to be searched for is a minimum
                            or maximum of the given data; True yields
                            a maximum, False yields a minimum

       Returns:
           optimum (List[float]): The optimum datapoint that matches
                                  the requested criteria
    """

    # Initialize figure
    fig, (ax, cax) = plt.subplots(ncols=2, figsize=(
        7.2, 4.2), gridspec_kw={"width_ratios": [1, 0.035]})
    fig.tight_layout(pad=3)
    fig.subplots_adjust(left=0, right=0.89)

    ax.remove()
    ax = fig.add_subplot(1, 1, 1, projection="3d")

    # Depending on the requested optimum, find
    # the maimum or minimum of the points
    if maximize == True:
        ind_c = c.index(max(c))
    else:
        ind_c = c.index(min(c))

    optimum = [x[ind_c], y[ind_c], z[ind_c], c[ind_c]]

    # Plot FlexPDE data as a scatter with colour varying with capacitance
    plot = ax.scatter(x, y, z, c=c, marker="o", cmap="gist_rainbow_r")

    # Add colour bar and chart titles
    fig.colorbar(plot, cax=cax, orientation="vertical",
                 extend="neither", label=c_label)

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)

    # Set up chart grid
    ax.grid(visible=True, which="major")
    ax.minorticks_on()
    ax.grid(visible=True, which="minor", linewidth=0.25)

    # Save figure and clear plotter for next figure
    plt.savefig(f"figure_{count}.png", format="png", dpi=500)
    ax.clear()
    plt.close(fig)

    return optimum


# Ensuring the script is only imported into other programs as intended
if __name__ == "__main__":
    print("This script is only meant to be implemented as a module, not run independently.")
