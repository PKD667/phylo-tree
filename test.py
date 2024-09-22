from tree import PhyloTree, build_tree_dist, print_dists

from compare import random_string, taxon_distances,fasta_distances,parralel_fasta_distances

import numpy as np


    

def test_tree():
    distances = [
        [0,1,1,5], # Species A
        [1,0,2,7], # Species B
        [3,2,0,9], # Species C
        [2,7,9,0]  # Species D
    ]

    taxons = ["A","B","C","D"]

    # Theoretically we should have:
    # A
    # |---B
    # |   |---C
    # D

    tree = build_tree_dist(taxons,distances)

    tree.visualize()

def mutate(genome, p):

    new_genome = ""

    c = 0

    # mutate n characters of the genome
    for g in genome:
        # DELETE
        if np.random.random() < (p/3):
            c += 1
            continue
        # REPLACE
        elif np.random.random() < (2*p/3):
            c += 1
            new_genome += np.random.choice(['A','C','G','T'])
            continue
        # INSERT
        elif np.random.random() < p:
            c += 1
            new_genome += np.random.choice(['A','C','G','T'])
        
        new_genome += g

    print(f"Mutated {c} characters")

    return new_genome   

def get_eval_data():
    # here we are gonna define species with a custom genome

    base = random_string(1000)

    species = {
        "A" : base,
        "B" : mutate(base,0.4),
        "C" : mutate(base,0.3),
        "D" : mutate(base,0.25)
    }

    # create sub species
    taxons = list(species.keys())
    
    species_b = {
        "E" : mutate(species["B"],0.20),
        "F" : mutate(species["B"],0.15)
    }

    species_c = {
        "G" : mutate(species["C"],0.21),
        "H" : mutate(species["C"],0.25)
    }

    species_d = {
        "I" : mutate(species["D"],0.2),
        "J" : mutate(species["D"],0.4)
    }

    full_taxons = list(species.keys()) + list(species_b.keys()) + list(species_c.keys()) + list(species_d.keys())
    full_genomes = list(species.values()) + list(species_b.values()) + list(species_c.values()) + list(species_d.values())

    print(full_taxons)
    print(full_genomes)

    return full_taxons,full_genomes

def test_distances():

    full_taxons,full_genomes = get_eval_data()

    # calculate distances
    distances_matrix = parralel_fasta_distances(full_genomes)

    print(distances_matrix)

    return full_taxons,distances_matrix



if __name__ == "__main__":
    t,d = test_distances()
    print_dists(t,d)
    # un-numpy the distances
    d = [[d[i][j] for j in range(0,len(d[i]))] for i in range(0,len(d))]
    tree = build_tree_dist(t,d)
    tree.visualize()