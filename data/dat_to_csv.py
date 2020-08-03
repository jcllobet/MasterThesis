import csv

with open("airports.dat") as f:
    with open("airports.csv", "w") as f1:
        for line in f:
            f1.write(line)