import unittest
from physics import *


class TestHello(unittest.TestCase):
    def test_calculate_bouyancy(self):
        self.assertRaises(ValueError, calculate_bouyancy("hi", 5))
