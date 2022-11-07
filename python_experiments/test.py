from OHFSS.Kernel import PyKernel, PyTestCPP
from OHFSS.Hello import call_hello
from math import pi

w01 = 3 * 2 * pi * 1e9
w12 = w01 - 0.25 * 2 * pi * 1e9
wt = 25 * 2 * pi * 1e9
w = 4e-12
T = pi / wt * 2
tstep = 5e-14
h = 1.054e-34

print("1")
cl = PyTestCPP(h, w01, w12)
print(cl.get())
cl.eig(3)
print(cl.get())


# print(tstep, w01, w12, wt, w, T, 3)

# obj = PyKernel(tstep, w01, w12, wt, w, T, 3)

# print(tstep, w01, w12, wt, w, T, 3)
# print(obj.CreateStartSCALLOP(120, 0.001))
