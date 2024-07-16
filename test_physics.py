import unittest
import numpy as np
import physics as p


class TestHello(unittest.TestCase):
    def test_calculate_bouyancy(self):
        with self.assertRaises(TypeError):
            p.calculate_bouyancy("hi", 5)

    def test_plot(self):
        p.plot_auv2_motion(
            np.array([0, 0, 300, 300]),
            np.pi / 4,
            1,
            1,
            t_final=1,
            dt=0.05,
            inertia=1,
            mass=1,
        )
