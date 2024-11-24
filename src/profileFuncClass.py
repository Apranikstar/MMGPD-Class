import numpy as np

class ProfileFunction:
    def __init__(self, parameters, x):
        # Calculate the function at initialization
        self.funcAtPoint = parameters[0] * np.power(1 - x, 3) * np.log(1 / x) \
                         + parameters[1] * np.power(1 - x, 3) \
                         + parameters[2] * x * np.power(1 - x, 2)

    def __call__(self):
        # Recompute the function with new parameters and x
        return self.funcAtPoint
    

