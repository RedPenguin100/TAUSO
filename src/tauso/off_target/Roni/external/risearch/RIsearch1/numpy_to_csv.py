import pandas as pd
import numpy as np
import os

def numpy_to_csv(matrix = np.random.rand(36, 36) , filename = "./modified_dsm.csv"):
    # Get the directory from the provided path
    module_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(module_dir, filename)
    assert matrix.shape == (36, 36), f"Expected shape (36, 36), but got {matrix.shape}"
    df = pd.DataFrame(matrix)
    df.to_csv(file_path, index=False, header=False)

#func and func_extend: function that modified np_array with shape (6,6,6,6)
def dsm_variable_to_csv(dsm_name = "dsm_su95_rev_woGU_pos",func=None, kwargs=None, func_extend=None, kwargs_extend=None):
    import ctypes
    import numpy as np

    module_dir = os.path.dirname(os.path.abspath(__file__))
    dsm_so_path = os.path.join(module_dir, "./dsm.so")
    # Load the shared library
    lib = ctypes.CDLL(dsm_so_path)  # Use "dsm.dll" on Windows

    # Define a function to load and reshape the 4D array
    def load_and_reshape_array(var_name,func=None,kwargs=None):
        """Load a 6x6x6x6 C array and reshape it into 36x36 NumPy array."""
        # Define the array type (4D short array)
        ArrayType = ctypes.c_short * 6 * 6 * 6 * 6

        # Get a reference to the array from the shared library
        c_array = ArrayType.in_dll(lib, var_name)

        # Convert to a NumPy array
        np_array = np.ctypeslib.as_array(c_array).reshape((6, 6, 6, 6))
        if func:
            if kwargs:
                np_array = func(np_array,**kwargs)
            else:
                np_array = func(np_array)
        # Reshape it to 36x36
        np_reshaped = np_array.reshape((36, 36))

        return np_reshaped

    # Call the function for both variables
    dsm_name_array = load_and_reshape_array(dsm_name,func,kwargs)
    
    dsm_extend = load_and_reshape_array("dsm_extend",func_extend,kwargs_extend)

    numpy_to_csv(matrix = dsm_name_array)
    numpy_to_csv(matrix = dsm_extend ,filename ="./modified_dsm_extend.csv")
