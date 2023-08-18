import pickle

def get_results():


    # take inputs
    pdb_code = input("PDB code: ")
    a_m = input("A/M/M_noncovalent: ")
    Constructor = input("Constructor: ")
    datetime = input("Datetime: ")

    # if matrix is an adjacency matrix
    try:
        with open(f"./pygenstability/results_{pdb_code}_{a_m}_{Constructor}_{datetime}.pkl", "rb") as results_file:
            return pickle.load(results_file)
    except:
        raise FileNotFoundError(f"File './PyGenStability_results/results_{pdb_code}_{a_m}_{Constructor}_{datetime}.pkl' does not exist")

def header(words):
    print("\n" + words)
    print('-' * len(words))

def print_results(results):

    # print run_params
    header("run_params")
    run_params = results["run_params"]
    for param in run_params.keys():
        print(f"{param}: {run_params[param]}")

    # print time_scales
    header("scales")
    scales = results["scales"]
    print(scales)
    print(f"Number of time scales = {len(scales)}")

    # print number_of_communities
    header("number_of_communities")
    number_of_communities = results["number_of_communities"]
    print(number_of_communities)

    max_communities = max(number_of_communities)
    max_locations = []
    for i in range(number_of_communities.count(max_communities)):
        if len(max_locations) == 0:
            max_locations.append(number_of_communities.index(max_communities, 0))
        else:
            max_locations.append(number_of_communities.index(max_communities, max_locations[-1]+1))
    print(f"Maximum number of communities = {max_communities} at time scales: {max_locations}")

    min_communities = min(number_of_communities)
    min_locations = []
    for i in range(number_of_communities.count(min_communities)):
        if len(min_locations) == 0:
            min_locations.append(number_of_communities.index(min_communities, 0))
        else:
            min_locations.append(number_of_communities.index(min_communities, min_locations[-1]+1))
    print(f"Minimum number of communities = {min_communities} at time scales: {min_locations}")

    # print stabilities
    header("stability")
    stability = results["stability"]
    print(stability)

    # print NVI
    header("NVI")
    nvi = results["NVI"]
    print(nvi)

    # print selected_partitions
    header("selected_partitions")
    selected_partitions = results["selected_partitions"]
    print(selected_partitions)

    # print community_id at specified time scale
    scale = int(input("\nTime scale for community_id: "))
    header(f"community_id (time scale = {scale})")
    community_id = results["community_id"][scale]
    print(community_id.tolist())
    print(f"Number of communities = {number_of_communities[scale]}")
    print(f"Number of bonds = {len(community_id)}")

def main():
    results = get_results()
    print_results(results)

main()
