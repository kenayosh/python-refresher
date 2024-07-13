import unittest
import psetp1


class TestHello(unittest.TestCase):
    def setUp(self):
        self.account1 = psetp1.Bank("Ken", "1")

    def test_name(self):
        account2 = psetp1.Bank("Bob", "2")
        self.assertEqual(self.account1.name, "Ken")
        self.assertEqual(account2.name, "Bob")
        self.assertNotEqual(self.account1.name, "Bob")
        self.assertNotEqual(account2.name, "Ken")

    def test_number(self):
        account2 = psetp1.Bank("Bob", "2")
        self.assertEqual(self.account1.number, "1")
        self.assertEqual(account2.number, "2")
        self.assertNotEqual(self.account1.number, "2")
        self.assertNotEqual(account2.number, "1")

    # def test_sin(self):
    #     self.assertEqual(hello.sin(0), 0)
    #     self.assertEqual(hello.sin(1), 0.8414709848078965)

    # def test_cos(self):
    #     self.assertEqual(hello.cos(0), 1)
    #     self.assertEqual(hello.cos(1), 0.5403023058681398)

    # def test_tan(self):
    #     self.assertEqual(hello.tan(0), 0)
    #     self.assertEqual(hello.tan(1), 1.5574077246549023)

    # def test_cot(self):
    #     self.assertEqual(hello.cot(0), float("inf"))
    #     self.assertEqual(hello.cot(1), 0.6420926159343306)


if __name__ == "__main__":
    unittest.main()
