import numpy as np

filename = "out.log"

x = []
y = []

seen_x = False

with open(filename, 'r') as f:
  for line in f:
    if line.strip() == "" and seen_x:
        break
    if line.strip() == "A" or line.strip() == "B":
        continue
    if line.strip() == "":
        seen_x = True
        continue
    if not seen_x:
        row = line.strip().split("  ")
        x.append([float(i) for i in row])
    else:
        row = line.strip().split("  ")
        y.append([float(i) for i in row])

res = np.dot(x, y)

print(res)