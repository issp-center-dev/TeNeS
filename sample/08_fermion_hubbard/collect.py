import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--Dmax", type=int, default=5, help="Maximum bond dimension")
parser.add_argument("--mu", type=float, default=2.0, help="Chemical potential")

args = parser.parse_args()

file_res = open("result.dat", "w")
print("# D E n doublon SS", file=file_res)

for D in range(2, args.Dmax+1):
    E = 0.0
    n = 0.0
    doublon = 0.0
    SS = 0.0
    with open(f"output_D{D}_measure/density.dat") as f:
        for line in f:
            words = line.strip().split()
            if len(words) == 0:
                continue
            if words[0].strip() == "Energy":
                E += float(words[2])
            elif words[0].strip() == "n" :
                E += args.mu * float(words[2])
                n += float(words[2])
            elif words[0].strip() == "doublon":
                doublon += float(words[2])
            elif words[0].strip() == "SzSz":
                SS += float(words[2])
            elif words[0].strip() == "SxSx":
                SS += float(words[2])
            elif words[0].strip() == "SySy":
                SS += float(words[2])
    print(f"{D} {E} {n} {doublon} {SS/2}", file=file_res)

