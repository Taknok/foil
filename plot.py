import sys

with open(sys.argv[1]) as f1:
    lines = f1.readlines()
    name0 = lines[0].strip()
    x1 = []
    y1 = []
    for i in lines[1:]:
        x, y = i.split()
        x1.append(float(x))
        y1.append(float(y))
        
    for i in range(len(x1)-1):
        print(f"{x1[i+1]:.5f} {y1[i+1]:.5f} ", end="")
              
        