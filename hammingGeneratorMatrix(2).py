##Helper

def decimalToVector(i, l):
    bits = []
    
    for k in range(l, 0, -1):
        v = 2 ** (k - 1)
        if i - v >= 0:
            i -= v
            bits.append(1)
        else:
            bits.append(0)
            
    return bits

##Repetition Codes

def repetitionEncoder(m, n):
    return m * n

def repetitionDecoder(v):
    num_of_ones = sum(v)
    num_of_zeroes = len(v) - num_of_ones

    if num_of_ones == num_of_zeroes:
        return []

    if num_of_ones > num_of_zeroes:
        return [1]

    return [0]

##Hamming

def message(a):
    l = len(a)
    r = 2
    r_sum = 2 ** r - 2 * r - 1

    while r_sum < l:
        r += 1
        r_sum = 2 ** r - 2 * r - 1

    k = r_sum + r
    bin_l = decimalToVector(l, r)
    result = bin_l + a

    while len(result) != k:
        result.append(0)

    return result

def hammingEncoder(m):
    #derive r
    total_len = len(m)
    r = 2
    r_calc = 2 ** r - r - 1
    
    while r_calc < total_len:
        r += 1
        r_calc = 2 ** r - r - 1

    #if r isn't exact this cant be encoded
    if r_calc != total_len:
        return []

    #multiply the matrices to get the encoded message
    generator_matrix = hammingGeneratorMatrix(r)
    encoded_message = multiply_matrices([m], generator_matrix)
    
    return encoded_message[0]

def multiply_matrices(a, b):
    dim_a = (len(a), len(a[0]))
    dim_b = (len(b), len(b[0]))

    assert dim_a[1] == dim_b[0]
    dim_result = (dim_a[0], dim_b[1])

    b_cols = []

    for i in range(dim_b[1]):
        b_col = []
        for row in b:
            b_col.append(row[i])
        b_cols.append(b_col)

    result_matrix = []

    for row in a:
        result_row = []
        for col in b_cols:
            element = 0
            for i in range(0, len(row)):
                element += row[i] * col[i]
            result_row.append(element)
        result_matrix.append(result_row)

    return result_matrix

#function HammingG
#input: a number r
#output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code
def hammingGeneratorMatrix(r):
    n = 2**r-1
    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    return G



print(hammingEncoder([1, 0, 0, 0]))
