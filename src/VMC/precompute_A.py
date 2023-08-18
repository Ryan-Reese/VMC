from precompute.create_A import create_A
from precompute.run_stab_A import run_stab_A

def main():

    # take PDB code from user input
    pdb_code = input("PDB Code: ")

    # create A matrix
    create_A(pdb_code)

    # run Markov Stability
    run_stab_A(pdb_code)

main()
