"""Program containing all necessary functions to run flex,
   sweep through data, and organize results. Held for other programs to access.

   Inputs:
         constants.py: Loads TEMPLATE_SCRIPT and FLEX_VERSION for use in running FlexPDE
         solver.py: Global functions are called from within solver.py with any parameters necessary passed to them
         output*.txt: Takes output data from FlexPDE and collates it
    
   Outputs:
         constants.py: Saves output data in constants.py to facilitate sharing across multiple scripts
         output*.txt: Directs FlexPDE to store output data in text files
"""

# Imports
import constants as ct
from threading import Thread
from itertools import groupby
from numpy import linspace
from subprocess import run
from typing import List

# Initialize empty list to keep track of how many processes are alive at once
threads = list()

# Function Definitions


def initialize_range(*args: List[int], points=7):
    """Initializes given a list of upper and lower bounds to sweep from for
       each value to be varied. Creates a specified number of
       sweep points (default is 7).

       Args:
           *args (List[List[int]]): A list of lists, each list within
           containing an upper and lower bound to be unpacked and swept
           through as a range for its respective variable
           points (int): Defualt value is 7; Number of points to make for
                         sweeping through range for each variable

       Returns:
           None
    """

    # Set up range for each variable provided with bounds
    var_ranges = [linspace(i[0], i[1], points) for i in args]

    # Place each valid value in range in order and grouped for running in a FlexPDE script
    # Round all values to 6 decimal places and remove any duplicates after rounding
    rounded_vars = [[round(i, 6), round(j, 6), round(
        k, 6)] for i in var_ranges[0] for j in var_ranges[1] for k in var_ranges[2]]

    ct.sweep_range = [i for i, _ in groupby(rounded_vars)]


def run_flex(name: str, *args: List[str]):
    """Runs FlexPDE given a name for the script and any arguments to pass to
       the template script located in constants. Outputs the file name and
       returncode. Saves and runs the FlexPDE script in the working directory.

       Args:
           name (str): A name for the FlexPDE script excluding the file extension
           *args (List[str]): Parameters to be passed onto the script in constants
                        through the format method and in the order they were given

       Returns:
           None     
    """

    # Open FlexPDE script and write formatted template
    with open(f"{name}.pde", "w") as stdin:
        print(ct.TEMPLATE_SCRIPT.format(*args), file=stdin)

    # Adjust timeout according to computing power
    # (if scripts return with timeout, increase timeout)
    completed = run([ct.FLEX_VERSION, "-S", f"{name}.pde"], timeout=1000)

    # Output runcode of file to confirm successful completion
    print(f"File {name} returned {completed.returncode}")


def sweep(material: str):
    """Sweeps through all values present in the pre-defined sweep range by
       instantiating a new thread for each datapoint and running the
       corresponding number of FlexPDE instances to produce the datapoints.
       Outputs when a new thread has started and when all threads have completed.

       Args:
           None

       Returns:
           None     
    """

    # Set maximum number of threads to run at once, change depending on CPU core count
    chunk_size = 50

    # Get length of range to break into chunks
    sweep_length = len(ct.sweep_range)

    # Break range into manageable chunks
    sweep_blocks = [ct.sweep_range[i:i + chunk_size]
                    for i in range(0, sweep_length, chunk_size)]

    # Taking all items in the sweep range and indexes to be each
    # passed into a thread with the target set to run_flex()
    for (i, chunk) in enumerate(sweep_blocks):
        for (j, t_list) in enumerate(chunk):
            identifier = j + i*chunk_size

            if material == "Ti":
                t = Thread(target=run_flex, args=(f"script_{identifier}", ct.E_LS[0], ct.NU_LS[0], t_list[
                           0], t_list[1], t_list[2],  ct.ALPHA_T_LS[0], ct.K_LS[0], ct.SIGMA_E_LS[0], identifier,))
            elif material == "Au":
                t = Thread(target=run_flex, args=(f"script_{identifier}", ct.E_LS[1], ct.NU_LS[1], t_list[
                           0], t_list[1], t_list[2],  ct.ALPHA_T_LS[1], ct.K_LS[1], ct.SIGMA_E_LS[1], identifier,))
            elif material == "Brass":
                t = Thread(target=run_flex, args=(f"script_{identifier}", ct.E_LS[2], ct.NU_LS[2], t_list[
                           0], t_list[1], t_list[2],  ct.ALPHA_T_LS[2], ct.K_LS[2], ct.SIGMA_E_LS[2], identifier,))
            elif material == "Vb":
                t = Thread(target=run_flex, args=(f"script_{identifier}", ct.E_LS[3], ct.NU_LS[3], t_list[
                           0], t_list[1], t_list[2],  ct.ALPHA_T_LS[3], ct.K_LS[3], ct.SIGMA_E_LS[3], identifier,))

            threads.append(t)
            t.start()
            print(f"Thread script_{identifier} Started")

        # Wait for all threads to finish before moving on to next step
        # to prevent FileNotFound errors
        while True:
            alive = [thread.is_alive() for thread in threads]
            if not any(alive):
                print("All threads complete!")
                break


def collate_data():
    """Collates all data from the output files produced by FlexPDE
       and saves them in blank lists stored in constants.

       Args:
           None

       Returns:
           None     
    """

    # Re-initalize output variables as empty lists to ensure
    # previous data from other iterations is erased
    ct.displacement, ct.temp_junc, ct.t_metal_top, ct.t_metal_bottom, ct.t_metal_side = [], [], [], [], []

    # Taking all items in the sweep range and indexes to view output files
    # produced by FlexPDE
    for (i, t_list) in enumerate(ct.sweep_range):
        with open(f"output_{i}.txt", "r") as stdin:
            data = stdin.readlines()

            d = abs(float(data[6].split(" = ")[1]))
            t = float(data[7].split(" = ")[1])

            ct.displacement.append(d*1e6)
            ct.temp_junc.append(t)
            ct.t_metal_top.append(t_list[0]*1e6)
            ct.t_metal_bottom.append(t_list[1]*1e6)
            ct.t_metal_side.append(t_list[2]*1e6)


# Ensuring the script is only imported into other programs as intended
if __name__ == "__main__":
    print("This script is only meant to be implemented as a module, not run independently.")
