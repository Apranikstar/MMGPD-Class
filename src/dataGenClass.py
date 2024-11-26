import numpy as np

class SkewedDataGenerator:
    def __init__(self, initial, final, num_points, scattering_rate):
        """
        Initialize the generator with the given parameters and generate data immediately.
        
        :param initial: The start of the range.
        :param final: The end of the range.
        :param num_points: Number of points to generate.
        :param scattering_rate: Rate at which values are skewed towards the initial.
        """
        if initial >= final:
            raise ValueError("Initial value must be less than final value")
        
        self.initial = initial
        self.final = final
        self.num_points = num_points
        self.scattering_rate = scattering_rate
        
        # Generate the skewed data immediately
        self.data = self.__generate_data__()
        
    def __call__(self):
        return self.data
    def __generate_data__(self):
        """
        Internal method to generate skewed data.
        
        :return: A sorted array of skewed data points in ascending order.
        """
        # Generate values in the range [initial, final]
        uniform_random_values = np.random.uniform(0, 1, self.num_points)
        
        # Apply the transformation to skew values towards initial
        skewed_values = np.power(uniform_random_values, self.scattering_rate)
        
        # Transform values back to the range [initial, final]
        scaled_values = self.initial + (self.final - self.initial) * skewed_values
        
        # Remove zero values
        non_zero_values = scaled_values[scaled_values != self.initial]
        
        # Sort the values in ascending order
        sorted_values = np.sort(non_zero_values)
        
        return sorted_values

    def __repr__(self):
        """
        Return a string representation of the generated data for easy visualization.
        """
        return f"{self.data}"

    def __getitem__(self, index):
        """
        Allow indexing into the generated data.
        
        :param index: The index or slice to access.
        :return: The corresponding data point(s).
        """
        return self.data[index]
