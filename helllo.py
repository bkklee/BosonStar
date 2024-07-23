from ctypes import *
file = "./testing.so"
my_func = CDLL(file)

print(my_func.square(10))