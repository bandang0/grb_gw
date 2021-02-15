from os import system
from sys import argv

exe = "./kn_pop.exe"
gwh_d = [157, 229, 472, 181, 715]
runs_d = ['O3', 'O4', 'O5', '190425', 'max']

if 'test' in argv:
    print("Testing mode...")
    N = int(1.e6)
    filenames_d = [f"data/kn_{r}_test.data" for r in runs_d]
else:
    print("Prod. mode...")
    N = int(3e7)
    filenames_d = [f"data/kn_{r}.data" for r in runs_d]


for i in range(len(gwh_d)):
    com = f"{exe} {gwh_d[i]} {N:d} > {filenames_d[i]}"
    print(f"Running '{com}'")
    system(com)
