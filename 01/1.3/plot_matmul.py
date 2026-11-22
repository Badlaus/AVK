import csv
import statistics

import matplotlib.pyplot as plt

data = {}
with open("matmul-O2.csv") as f:
    for row in csv.DictReader(f):
        key = (row["method"], int(row["n"]))
        data.setdefault(key, []).append(float(row["seconds"]))

methods = ["ijk", "ikj", "jik", "jki", "kij", "kji"]
colors  = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]
markers = ["o", "s", "^", "D", "v", "P"]

for method, color, marker in zip(methods, colors, markers):
    ns = sorted(n for (m, n) in data if m == method)
    med = [statistics.median(data[(method, n)]) for n in ns]
    plt.plot(ns, med, marker=marker, color=color, label=method)

plt.xscale("log")
plt.yscale("log")
plt.xlabel("размер матрицы n")
plt.ylabel("медиана времени, с")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("matmul-orders.png", dpi=150)
