import machupS as ms
import numpy as np


if __name__=="__main__":

    # Initialize input
    input_dict = {
        "wing" : {
            "span" : 10.0,
            "chord" : 1.0,
            "section" : {
                "CLa" : 6.28,
                "aL0" : 0.0
            },
            "N" : 20
        },
        "state" : {
            "V_inf" : 10.0,
            "frequency" : 1.0,
            "plunge_amplitude" : 1.0,
            "plunge_phase" : 0.0,
            "pitch_amplitude" : 5.0,
            "pitch_phase" : 0.0
        },
        "case_options" : {
            "dt" : 0.01
        }
    }

    # Intiialize case
    my_case = ms.Case(**input_dict)