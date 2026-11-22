
import csv
import statistics

data = {}

with open("strassen.csv") as f:
    for row in csv.DictReader(f):
        key = (int(row["n"]), row["method"])
        data.setdefault(key, []).append(float(row["seconds"]))

print(f'{"Размер":<10} {"Naive, с":>12} {"Strassen, с":>14} {"Ускорение":>12}')

for n in [128, 256, 512, 1024]:
    naive = statistics.median(data[(n, "naive")])
    strassen = statistics.median(data[(n, "strassen")])

    print(f'{n:<10} {naive:>12.6f} {strassen:>14.6f} {naive / strassen:>11.2f}x')
