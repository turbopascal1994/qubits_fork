with open('output_c++.txt', 'r') as f:
    c = f.read().split()

with open('C:/Users/User/Documents/MATLAB/output_matlab.txt', 'r') as f:
    matlab = f.read().split()

print(len(c), len(matlab))

if len(c) != len(matlab):
    print('beda')
    exit(0)

cnt = 0

for i in range(len(c)):
    if c[i] != matlab[i]:
        print(f'index = {i};c++ number = {c[i]};matlab number = {matlab[i]}')
        print(f'diff = {abs(float(c[i]) - float(matlab[i]))}')
        cnt += 1
        # exit(0)

print(cnt)
print('everything is ok')