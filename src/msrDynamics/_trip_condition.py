TRIP_TYPES = ['state', 'diff']

class TripCondition:
    """
    Class to store information about trip conditions for state variables.

    Attributes:
        trip_type (str): The type of trip condition, either 'state' or 'diff'.
        bounds (tuple): The bounds for the trip condition, default is (-inf, inf).
        idx (int or None): The index of the state variable the trip condition applies to.
        check_after (float or None): Time after which the trip condition should be checked.
        delay (float or None): Delay time for the trip condition.
    """

    def __init__(self, 
                 trip_type=None, 
                 bounds=(-float('inf'), float('inf')),
                 idx=None,
                 check_after=None,
                 delay=None) -> None:
        """
        Initialize a TripCondition instance.

        Parameters:
            trip_type (str, optional): The type of trip condition, either 'state' or 'diff'. Default is 'state'.
            bounds (tuple, optional): The bounds for the trip condition. Default is (-inf, inf).
            idx (int or None, optional): The index of the state variable the trip condition applies to. Default is None.
            check_after (float or None, optional): Time after which the trip condition should be checked. Default is None.
            delay (float or None, optional): Delay time for the trip condition. Default is None.
        """
        # Set trip type 
        if trip_type is None:
            self.trip_type = 'state'
        else:
            self._check_type(trip_type)
            self.trip_type = trip_type

        self.bounds = bounds
        self.idx = idx
        self.check_after = check_after
        self.delay = delay

    def _check_type(self, trip_type):
        """
        Check if the provided trip type is valid.

        Parameters:
            trip_type (str): The type of trip condition to check.

        Raises:
            ValueError: If the trip_type is not 'state' or 'diff'.
        """
        if trip_type not in TRIP_TYPES:
            raise ValueError('Invalid trip type.')
