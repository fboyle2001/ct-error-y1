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
    encoded_message = multiply_matrices_m2([m], generator_matrix)
    
    return encoded_message[0]

def multiply_matrices_m2(a, b):
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
            element %= 2
            result_row.append(element)
        result_matrix.append(result_row)

    return result_matrix

def hammingDecoder(v):
    r = 2

    #check we have a valid vector
    
    while 2 ** r - 1 < len(v):
        r += 1
        
    if 2 ** r - 1 != len(v):
        return []

    #find the parity check transpose
    #multiply it with v to find any errors
    
    pc_transpose = parity_check_transpose(r)
    result = multiply_matrices_m2([v], pc_transpose)[0]

    #if the vector's components are all 0 then there are no errors

    error = False

    for element in result:
        if element != 0:
            error = True
            break

    if not error:
        return v

    #otherwise we have an error
    #correct this by finding the position of the error and flipping the bit

    error_position = 0

    for i in range(0, len(result)):
        error_position += result[i] * (2 ** (len(result) - i - 1))

    #correct the error and fix the binary

    v[error_position - 1] += 1
    v[error_position - 1] %= 2
    
    return v

def parity_check_transpose(r):
    rows = []

    #calculate the binary for every number 1...2**r - 1
    #these are the rows of the transpose
    
    for i in range(1, 2 ** r):
        bin_vector = decimalToVector(i, r)
        rows.append(bin_vector)
        
    return rows

def messageFromCodeword(c):
    r = 2

    #check we have a valid vector
    
    while 2 ** r - 1 < len(c):
        r += 1
        
    if 2 ** r - 1 != len(c):
        return []

    k = 1
    message = []

    #skip each index if it is 2 ** k - 1

    for i in range(0, len(c)):
        if i == k - 1:
            k *= 2
            continue
        message.append(c[i])
    
    return message

def dataFromMessage(m):
    r = 2

    #check we have a valid vector
    
    while 2 ** r - r - 1 < len(m):
        r += 1
        
    if 2 ** r - r - 1 != len(m):
        return []

    #first r bits represent the length of the data

    len_bin = m[0:r]
    data_pos = 0

    #calculate the length from the binary

    for i in range(0, len(len_bin)):
        data_pos += len_bin[i] * (2 ** (len(len_bin) - i - 1))

    #ensure that the length of the data is actually possible

    if data_pos > 2 ** r - 2 * r - 1:
        return []

    #select and return the data

    data = m[r:r+data_pos]
    
    return data

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
