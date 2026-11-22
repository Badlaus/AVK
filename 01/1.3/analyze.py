
import csv
import statistics

data = {}

with open("matmul-O2.csv") as f:
    for row in csv.DictReader(f):
        key = (row["method"], int(row["n"]))
        data.setdefault(key, []).append(float(row["seconds"]))

methods = ["ijk", "ikj", "jik", "jki", "kij", "kji"]
sizes = [128, 256, 384, 512]

print(f'{"Порядок":<10} {"128":>12} {"256":>12} {"384":>12} {"512":>12}')

for method in methods:
    values = [
        statistics.median(data[(method, n)])
        for n in sizes
    ]
    print(f'{method:<10} ' + ' '.join(f'{v:>12.6f}' for v in values))
