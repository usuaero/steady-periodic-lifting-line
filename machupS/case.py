import json

import numpy as np

from machupS.wing import Wing


class Case:
    """Defines a steady-periodic wing simulation."""

    def __init__(self, **kwargs):

        # Get kwargs
        wing_kwargs = kwargs.get("wing")
        state_kwargs = kwargs.get("state")
        self._V_inf = state_kwargs.get("V_inf")
        self._w = state_kwargs.get("frequency")
        self._y_A = state_kwargs.get("plunge_amplitude")
        self._phi_P = state_kwargs.get("plunge_phase")
        self._a_hat_A = state_kwargs.get("pitch_amplitude")
        self._phi_a_hat = state_kwargs.get("pitch_phase")
        self._dt = kwargs.get("case_options").get("dt")

        # Initialize wing
        self.wing = Wing(**wing_kwargs)