class Bank:
    def __init__(self, name, number):
        # Give warning if name isn't string and set as name
        if not isinstance(name, str):
            print("Warning: Name is not a string")
            self.name = name.toString()
        else:
            self.name = name

        # If its not an instance of a string
        if not isinstance(number, str):
            print("Error: Please ensure number is an string")
            self.name = name.toString()
            return
        else:
            self.number = number

        # Set balance to 0
        self.balance = 0

    def withdraw(self, amount):
        if amount < 0:
            print("Please withdrawl a negative ")
            return
        if amount > self.balance:
            print("Not enough money in balance")
        try:
            self.balance = self.balance - amount
        except TypeError:
            print("Please enter a number")

    def deposit(self, amount):
        if amount < 0:
            print("Please withdrawl a negative ")
            return
        try:
            self.balance = self.balance + amount
        except TypeError:
            print("Please enter a number")

    def printBalance(self):
        return self.balance
