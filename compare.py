
import numpy as np
import random
import time
from get import get_fasta

import matplotlib.pyplot as plt

def matrix_render(matrix):
    plt.imshow(matrix, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.show()


def levenshtein_matrix(s, t):
    m = len(s)
    n = len(t)
    d = np.zeros((m+1, n+1))
    for i in range(1, m+1):
        d[i,0] = i
    for j in range(1, n+1):
        d[0,j] = j
    for j in range(1, n+1):
        for i in range(1, m+1):
            if s[i-1] == t[j-1]:
                d[i,j] = d[i-1,j-1]
            else:
                d[i,j] = min(d[i-1,j]+1, d[i,j-1]+1, d[i-1,j-1]+1)
    return d

# implementation the levenshtein distance using the matrix method
def mcompare(s,t):
    
    m = len(s)
    n = len(t)

    d = levenshtein_matrix(s,t)

    return int(d[m,n])

# compare two set of genomic data
# Input: a - the first set of genomic data
#        b - the second set of genomic data
# Output: Levenshtein distance between the two sets of genomic data
# Genome strings are huge (millions of characters) so we need to be memory efficient.
# The matrices methods uses O(mn) space, where m and n are the lengths of the two strings.
# basicaly Exabytes of memory for a genome of 10^9 characters.
# Algo from here : https://en.wikipedia.org/wiki/Levenshtein_distance#Iterative_with_two_matrix_rows
def compare(s, t) -> int:
    
    m = len(s)
    n = len(t)

    v0 = np.asarray(range(0,n+1))
    v1 = np.zeros(n+1)

    for i in range(0, m):
        v1[0] = i+1

        for j in range(0, n):
            v1[j+1] = min(
                v0[j+1] + 1, # deletion cost
                v1[j] + 1, # insertion cost
                v0[j] + (s[i] != t[j]) # substitution cost
            )

        v0 = np.copy(v1)
        print(f"{i}/{m}", end="\r")
    
    return int(v1[n])

# approximation of the levenshtein distance
def approx(s, t, k=1024) -> int:

    # get small length
    m = min(len(s), len(t))
    # get big length
    n = max(len(s), len(t))

    print(f"Approx with k = {k}")

    f = m // k

    d = 0

    for i in range(0, f):
        j = k * i  
        d += compare(s[j:j+k], t[j:j+k])
        print(f"{i}/{f} {d}", end="\r")
    d+= compare(s[f*k:], t[f*k:])
    print()

    #d *= 0.9

    #d += n-m

    return int(d)

def rapprox(s, t, k=1024,n=1000) -> int:
    # get small length
    m = min(len(s), len(t))
    # get big length
    n = max(len(s), len(t))



    d = 0

    chunk_size = m // n

    for i in range(0, n):
        j = chunk_size * i

        d += compare(s[j:j+chunk_size], t[j:j+chunk_size])
        print(f"{i}/{n} {d}", end="\r")
    
    # convert the distance to get an approximation of the levenshtein distance
    # we need to multiply the distance by the ratio of the two strings
    d *= m / n

    # Add the difference in length between the two strings
    #d += n - m

    return int(d)

# chosen compare function
cmp=compare

# the goal of this function is to provide dimensional Edit-distance between two strings
# We need to get a norm (edit distance) between two strings
# and this norm should have dimensionality in order to distinguish between different mutations
# We can use the Levenshtein distance as a norm, but it is not dimensional.
def analyse_mutations(s, t,d=100):

    # divide the string into d chunks
    
    if len(s) > len(t):
        s, t = t, s
    
    m = len(s)
    n = len(t)

    chunk_size = m // d

    d_vec = np.zeros(d)

    for i in range(0, d):
        j = chunk_size * i

        d_vec[i] = cmp(s[j:j+chunk_size], t[j:j+chunk_size])
    
    return d_vec


def visualize_analisys(ld):
    # ld is a list of d-dimensional levenshtein distances
    # plot all of them in a bar chart with different colors
    # the x axis is the dimension, the y axis is the distance

    if not ld:
        print("The input list is empty.")
        return

    l = len(ld[0])
    print(f"l = {l}")

    bar_width = 0.5  # Width of the bars
    indices = range(l)

    for i in range(len(ld)):
        plt.bar([x + i * bar_width for x in indices], ld[i], bar_width, label=f'Series {i+1}')
    
    plt.xlabel('Dimension')
    plt.ylabel('Levenshtein Distance')
    plt.title('Levenshtein Distances Across Dimensions')
    plt.legend()
    plt.show()

# compare two set of genomic data
def fasta_distances(fastas):
    n = len(fastas)
    d = np.zeros((n,n))
    for i in range(0, n):
        for j in range(i+1, n):
            d[i,j] = cmp(fastas[i], fastas[j])
            d[j,i] = d[i,j]
    
    return d


from concurrent.futures import ThreadPoolExecutor
def parralel_fasta_distances(fastas):
    def compute_distance(fastas, i, j):
        return i, j, cmp(fastas[i], fastas[j])
    n = len(fastas)
    d = np.zeros((n, n))
    
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(compute_distance, fastas, i, j) for i in range(n) for j in range(i+1, n)]
        
        for future in futures:
            i, j, distance = future.result()
            d[i, j] = distance
            d[j, i] = distance
    
    return d

# compute a set of distances between taxons.
# Input: taxons - a list of taxons
# Output: a matrix of distances between the taxons
def taxon_distances(taxons):

    codes = []
    for i in range(0, len(taxons)):
        codes.append(get_fasta(taxons[i]))

    print(f"Got {len(codes)} codes")

    d = fasta_distances(codes)

    return d




def random_string(n):
    return ''.join(np.random.choice(['A','C','G','T'], n))


def test(n=10,l=100,ap_k=16):

    if not test.eval_data:
        test.eval_data = get_fasta("soybean")


    s = []
    # cut n sequences of 1000 characters
    for i in range(0, n):
        s.append(test.eval_data[i*l:(i+1)*l])
    

    t = {
        "c" : [],
        "a" : [],
        "r": [],
        "m": []
    }

    o = {
        "c" : [],
        "a" : [],
        "r": [],
        "m": []
    }

    # run compare
    for i in range(0, n):

        start = time.time()
        c = compare(s[i], s[(i+1)%n])
        t["c"].append(time.time() - start)

        start = time.time()
        a = approx(s[i], s[(i+1)%n], k=ap_k)
        t["a"].append(time.time() - start)

        start = time.time()
        r = rapprox(s[i], s[(i+1)%n], k=ap_k)
        t["r"].append(time.time() - start)

        start = time.time()
        m = mcompare(s[i], s[(i+1)%n])
        t["m"].append(time.time() - start)

        o["c"].append(c)
        o["a"].append(a)
        o["r"].append(r)
        o["m"].append(m)


    print(f'{"Method":<10} {"Time":<10} {"Result":<10}')
    print(f'{"-"*10} {"-"*10} {"-"*10}')

    for i in range(0, n):
        print(f'{"Compare":<10} {t["c"][i]:<10.5f} {o["c"][i]:<10}')
        print(f'{"MCompare":<10} {t["m"][i]:<10.5f} {o["m"][i]:<10}')
        print(f'{"Approx":<10} {t["a"][i]:<10.5f} {o["a"][i]:<10}')
        print(f'{"RApprox":<10} {t["r"][i]:<10.5f} {o["r"][i]:<10}')

        print()
    
    # compare the output from approx, papprox and compare to get the average error
    approx_error = (sum([abs(o["a"][i] - o["c"][i]) for i in range(0, n)]) / n)
    rapprox_error = (sum([abs(o["r"][i] - o["c"][i]) for i in range(0, n)]) / n)

    print(f'{"Method":<10} {"Error":<10}')
    print(f'{"-"*10} {"-"*10}')
    print(f'{"Approx":<10} {approx_error:<10.5f}')
    print(f'{"RApprox":<10} {rapprox_error:<10.5f}')

test.eval_data = None