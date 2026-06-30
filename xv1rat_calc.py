x1min = 1e8
x1max = 6e13
nx1 = 1000

xv1rat = (x1max / x1min) ** (1.0 / nx1)
print(f"xv1rat = {xv1rat:.6f}")