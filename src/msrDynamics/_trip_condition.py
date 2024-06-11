
TRIP_TYPES = ['state', 'diff']

class TripCondition:
    '''
    Object to store information about trip conditions for state variables
    '''

    def __init__(self, 
                 trip_type = None, 
                 bounds = (-float('inf'),float('inf')),
                 idx = None,
                 check_after = None
                 ) -> None:
        
        # set trip type 
        if trip_type is None:
            self.trip_type = 'state'
        else:
            self._check_type(trip_type)
            self.trip_type = trip_type

        self.bounds = bounds
        self.idx = idx
        self.check_after = check_after

    def _check_type(self,trip_type):
        if trip_type not in TRIP_TYPES:
            raise ValueError('Invalid trip type.')