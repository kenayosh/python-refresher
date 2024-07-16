# import numpy as np


def hello():
    return "Hello, world!"


def add(a, b):
    return a + b


def sub(a, b):
    return a - b


def mul(a, b):
    return a * b


def div(a, b):
    if b == 0:
        raise ValueError("Can't divide by zero!")
    return a / b


# def sqrt(a):
#     return np.sqrt(a)


# def power(a, b):
#     return np.power(a, b)


# def log(a):
#     return np.log(a)


# def exp(a):
#     return np.exp(a)


# def sin(a):
#     return np.sin(a)


# def cos(a):
#     return np.cos(a)


# def tan(a):
#     return np.tan(a)


# def cot(a):
#     return 1 / np.tan(a)


def __main__():
    hello()


if __name__ == "__main__":
    __main__()
    
    #ros2 launch mavros apm.launch fcu_url:=udp://127.0.0.1:14550@14555 gcs_url:=udp://:14550@169.254.66.80:14550 tgt_system:=1 tgt_component:=1 system_id:=255 component_id:=240