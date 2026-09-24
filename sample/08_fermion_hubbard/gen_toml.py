import subprocess
import toml


def update(p, D, chi, stage):
    if stage == "opt":
        opt(p, D, chi)
    else:
        measure(p, D, chi)

def opt(p, D, chi):
    p["lattice"]["virtual_dim"] = D
    p["parameter"]["ctm"]["dimension"] = chi
    p["parameter"]["general"]["output"] = f"output_D{D}_opt"

    p["parameter"]["general"]["measure"] = False
    if D == 2:
        p["parameter"]["general"]["tensor_load"] = ""
    else:
        p["parameter"]["general"]["tensor_load"] = f"tensor_D{D-1}"

    p["parameter"]["general"]["tensor_save"] = f"tensor_D{D}"

    p["parameter"]["simple_update"]["num_step"] = 3000
    p["parameter"]["full_update"]["num_step"] = 0

    return p

def measure(p, D, chi):
    p["lattice"]["virtual_dim"] = D
    p["parameter"]["ctm"]["dimension"] = chi
    p["parameter"]["general"]["output"] = f"output_D{D}_measure"

    p["parameter"]["general"]["measure"] = True
    p["parameter"]["general"]["tensor_load"] = f"tensor_D{D}"
    p["parameter"]["general"]["tensor_save"] = ""

    p["parameter"]["simple_update"]["num_step"] = 0
    p["parameter"]["full_update"]["num_step"] = 0

with open("common.toml") as f:
    p = toml.load(f)

def main(D, chi):
    for stage in ("opt", "measure"):
        simple_toml = f"simple-{stage}-D{D}.toml"
        std_toml = f"std-{stage}-D{D}.toml"
        input_toml = f"input-{stage}-D{D}.toml"
        update(p, D, chi, stage)
        with open(simple_toml, "w") as f:
            toml.dump(p, f)
        cmd = f"tenes_simple {simple_toml} -o {std_toml}"
        subprocess.call(cmd.split())
        cmd = f"tenes_std {std_toml} -o {input_toml}"
        subprocess.call(cmd.split())

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--Dmax", type=int, default=5, help="Maximum virtual dimension")
    parser.add_argument("--chi-coeff", type=int, default=1, help="Coefficient of envionment virtual dimension chi = coeff*D*D")

    args = parser.parse_args()
    Dmax = args.Dmax
    chi_coeff = args.chi_coeff

    for D in range(2, Dmax+1):
        chi = D*D*chi_coeff
        main(D, chi)

