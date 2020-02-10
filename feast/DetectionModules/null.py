from .abstract_detection_method import DetectionMethod
from ..GeneralClassesFunctions.simulation_functions import set_kwargs_attrs


class Null(DetectionMethod):
    """
    This class specifies a null detection method. It includes a detection method and several parameters.
    """
    def __init__(self, time, gas_field, **kwargs):
        """
        Inputs:
              gas_field    a gas_field object (Defined in feast_classes)
              time         a time object (Defined in feast_classes)
              kwargs       optional input dictionary that will override default parameters
        """
        DetectionMethod.__init__(self, time, gas_field)
        # The null repair rate defaults to the same setting as the gas field, but can be set to a different value.
        # self.repair_rate = gas_field.null_repair_rate

        # -------------- Process details: None --------------
        # -------------- Financial Properties --------------
        # Capital costs are zero in this module.
        self.capital = [0]*time.n_timesteps  # dollars

        # maintenance costs are zero in the null module
        self.maintenance = [0]*time.n_timesteps  # $

        # Find cost is zero in the null module
        self.find_cost = [0]*time.n_timesteps  # $

        # Set attributes defined by kwargs. Only set attributes that already exist.
        set_kwargs_attrs(self, kwargs, only_existing=True)

    def detection(self, time, gas_field):
        """
        The null detection method is simply the null detection method defined in the super class DetectionMethod
        Inputs:
            time        an object of type Time (defined in feast_classes)
            gas_field   an object of type GasField (defined in feast_classes)
        """
        DetectionMethod.null_detection(self, time, gas_field)
