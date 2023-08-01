import pyalps

params = []

for L in [32]:
    for g in [0.0, 0.5, 0.8, 2.0]:
        for t in [
            0.05,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0,
            1.25,
            1.5,
            1.75,
            2.0,
            5.0,
            10.0,
        ]:
            params.append(
                {
                    "LATTICE": "square lattice",
                    "MODEL": "spin",
                    "local_S": 0.5,
                    "T": t,
                    "Jz": -1,
                    "Jxy": 0,
                    "THERMALIZATION": 5000,
                    "SWEEPS": 50000,
                    "Gamma": g,
                    "L": L,
                    "ALGORITHM": "loop",
                }
            )

pyalps.writeInputFiles("params", params)
