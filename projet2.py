from utils import *
import matplotlib.pyplot as plt

# Question 1
class Tree:
    def __init__(self, value=None, left=None, right=None):
        self.value = value
        self.left = left
        self.right = right

    def __get__(self, instance, owner):
        return self.value

    def __lt__(self, other):
        return self.value < other.value


    def has_children(self):
        return self.left != None or self.right != None

    def get_children(self):
        return (self.left, self.right)

# Traverse an huffman tree from the root to each leaf assigning to each child either the value '0' or '1'
def huffman_rec(node):
    if not node.has_children():
        return [""]

    i = 0
    codes = []
    for child in node.get_children():
        for code in huffman_rec(child):
            codes.append(str(i) + code) # Progressively concatenate the codes with the additional values '0' and '1'
        i += 1

    return codes

# Assign a Huffman code to each probability
def huffman(probabilities):
    trees = np.array([])
    for proba in probabilities:
        trees = np.append(trees, Tree(proba)) # For each given probability, we create a new leaf node

    while len(trees) > 1:
        sum_indices = np.argsort(trees)[:2] # Find the 2 nodes (by their indices) having the smallest probability values

        sum_val = 0
        for node in trees[sum_indices]:
            sum_val += node.value  # Sum their value up ...
        node = Tree(sum_val, trees[sum_indices[0]], trees[sum_indices[1]]) # ... and instantiate a new internal node with the result.
                                                                           # The 2 nodes are assigned to this new node left and right child.
        trees = np.append(trees, node) # add the newly created parent node
        trees = np.delete(trees, sum_indices) # remove children nodes


    root = trees[0]
    if not root.has_children():
        return {probabilities[0], "0"}
    else:
        codes = huffman_rec(root)
        sorted_codes = np.array(sorted(codes, key=len))
        sorted_indices = np.argsort(-np.fromiter(probabilities, dtype=float))
        ranked_indices = np.argsort(sorted_indices)
        return sorted_codes[ranked_indices]


# Question 2
def lempelZiv(input):
    code = input[0]
    decimal_address = 1
    read_head = 0
    lempel_ziv_dict = {input[read_head]:decimal_address} # dictionary of (input_subsequence, decimal_address)
                                                         # initialized with the first character of the input sequence
    decimal_address += 1
    read_head += 1

    while read_head < len(input): # Read the whole input sequence
        read_tail = read_head + 1
        read_block = input[read_head:read_tail]
        while read_block in lempel_ziv_dict:
            # Progressively increase the current read block of the input sequence until
            # it can not be found in the dictionary of previously ennb_bits_flippedered blocks of
            # input sequences
            read_tail += 1
            if read_tail > len(input): # border effect for the last input subsequence read
                binary_address_length = int(np.ceil(np.log2(decimal_address - 1))) # compute the last address length (for on-line encoding)
                code += binarize(lempel_ziv_dict[read_block], binary_address_length) # do not add the last bit of the input subsequence as it is not a new unique one
                return code
            read_block = input[read_head:read_tail]

        lempel_ziv_dict[read_block] = decimal_address # Add the address of the new (until-now-)unique input sequence in the dictionary
        binary_address_length = int(np.ceil(np.log2(decimal_address))) # compute the current address length (for on-line encoding)
        if len(read_block) > 1:
            code += binarize(lempel_ziv_dict[read_block[:-1]], binary_address_length) + input[read_tail - 1] # add address and last bit of source words to current code
        else: # border effect for the second unique length-one sequence appearing in the code
            code += '0' + input[read_tail - 1]

        read_head = read_tail
        decimal_address += 1 # yes ! 1 more encoded entries

    return code, lempel_ziv_dict


# Question 4
# def Lz77(): 
# Result: Write here the result
# A sliding window size l;
# An input string;
# while input is not empty do
# prefix := longest prefix of input that begins in window;
# if prefix exists in window then
#   d := distance to the start of the prefix;
#   l := length of prefix;
#   c := char following the prefix in input;
# else
#   d := 0;
#   l := 0;
#   c := first symbol of input;
# end
# append (d, l, c) to encoded input;
# shift the sliding window by l + 1 symbols (i.e., discard l + 1 symbols from the beginning of
# window and add the l + 1 first symbols of the input at the end of the window). end
def LZ77(input, sliding_window_size, decimal=False):
    binary_number_length = int(np.ceil(np.log2(sliding_window_size + 1)))
    input_index = 0
    sliding_window = "".ljust(sliding_window_size)
    code = ""
    while input_index < len(input): # read the whole input sequence
        look_ahead_buffer = input[input_index:input_index + sliding_window_size] # shift the look ahead buffer after the last encoded sequence
        sliding_window_search_index = 0
        look_ahead_buffer_index = 0
        longest_prefix = ""
        prefix = ""
        d = 0
        l = 0
        c = look_ahead_buffer[look_ahead_buffer_index]
        while sliding_window_search_index <= sliding_window_size \
                and look_ahead_buffer_index <= len(look_ahead_buffer) \
                and len(longest_prefix) <= (sliding_window_size - sliding_window_search_index): # read the whole sliding window (+ border effects)
            if sliding_window_search_index < sliding_window_size \
                    and look_ahead_buffer_index < len(look_ahead_buffer) \
                    and sliding_window[sliding_window_search_index] == look_ahead_buffer[look_ahead_buffer_index]:
                # If the read heads are still within the sliding window and in the look ahead buffer,
                # check for similar sequence blocks in the sliding window
                # and at the beginning of the look ahead buffer
                prefix += look_ahead_buffer[look_ahead_buffer_index]
                look_ahead_buffer_index += 1
            else: # if sequences are not matching (anymore) or if we are out of the sliding window or look ahead buffer
                if len(prefix) > 0 and len(longest_prefix) <= len(prefix): # if sequences are not matching anymore
                    longest_prefix = prefix
                    l = look_ahead_buffer_index # length of the prefix
                    d = sliding_window_size - sliding_window_search_index + l # distance to the start of the prefix;
                    if l < len(look_ahead_buffer):
                        c = look_ahead_buffer[l] # char following the prefix in input
                    else:
                        c = '' # border effect: when arriving at the end of the sequence, no more character can be read further

                prefix = ""
                look_ahead_buffer_index = 0

            sliding_window_search_index += 1


        # Append new coded sequence
        if decimal:
            code += '(' + str(d) + ',' + str(l) + ',' + c + ')' # Readable (decimal_d, decimal_l, char) tuples
        else:
            code += str(binarize(d, binary_number_length)) + str(binarize(l, binary_number_length)) + c # <binary_d><binary_l><c> sequence

        sliding_window = sliding_window[l + 1:] + input[input_index:input_index + l + 1] # shift sliding window after the newly coded sequence
        input_index += l + 1

    return code


# Question 12
# Combination of LZ77 and Huffman encoding
def LZ77_Huffman(input, sliding_window_size, is_LZ77_input=False):
    LZ77_res = input
    if not is_LZ77_input: # if the sequence is already encoded in LZ77
        LZ77_res = LZ77(input, sliding_window_size = sliding_window_size)

    # Compute the frequency of each l value from LZ77 output
    binary_number_length = int(np.ceil(np.log2(sliding_window_size + 1)))
    l_frequencies = dict()
    LZ77_res_length = len(LZ77_res)
    i = binary_number_length
    while i < LZ77_res_length:
        binary_l = LZ77_res[i:i + binary_number_length]
        if binary_l not in l_frequencies:
            l_frequencies[binary_l] = 1
        else:
            l_frequencies[binary_l] += 1
        i += 1 + 2 * binary_number_length
    # Assign Huffman codes to the different l values
    huffman_codes = huffman(l_frequencies.values())
    huffman_codes = dict(zip(l_frequencies.keys(), huffman_codes))

    # Encode the l values of the LZ77 sequence with Huffman
    LZ77_Huffman_l = ""
    LZ77_Huffman_res = ""
    i = binary_number_length
    while i <= LZ77_res_length - binary_number_length:
        binary_l = LZ77_res[i:i + binary_number_length] # <l> in binary
        LZ77_Huffman_l += huffman_codes[binary_l] # encode the current <l> value with Huffman
        LZ77_Huffman_res += LZ77_res[i - binary_number_length:i] # concatenate the current <d> value form LZ77 output
        if i + binary_number_length < LZ77_res_length:
            LZ77_Huffman_res += LZ77_res[i + binary_number_length] # concatenat the current <c> value from LZ77 output
        i += 1 + 2 * binary_number_length

    LZ77_Huffman_res += LZ77_Huffman_l # concatenate all <l> values from LZ77 output after being encoded with Huffman
    return LZ77_Huffman_res


# Question 6
def length_bin(arr):
    length = np.ones(len(arr))
    for i in range(len(arr)):
        length[i] = len(arr[i])
    return length

def expected_huffman_length(probabilities, arr = None):
    if not arr:
        arr = huffman(probabilities)
    length = length_bin(arr)
    return np.sum(probabilities * length)

def entropy(probabilities):
    tmp = 0
    for i in probabilities:
        tmp -= i * np.log2(i)
    return tmp


# Question 18
def add_redundancy(arr):
    s1 = int(arr[0])
    s2 = int(arr[1])
    s3 = int(arr[2])
    s4 = int(arr[3])
    p1 = (s1 + s2 + s3) % 2
    p2 = (s2 + s3 + s4) % 2
    p3 = (s1 + s3 + s4) % 2
    arr = arr + str(p1) + str(p2) + str(p3)
    return arr

# Question 19
def decode_redundancy(arr):
    signal = arr[0:4]
    s1 = int(signal[0])
    s2 = int(signal[1])
    s3 = int(signal[2])
    s4 = int(signal[3])
    p1 = (s1 + s2 + s3) % 2
    p2 = (s2 + s3 + s4) % 2
    p3 = (s1 + s3 + s4) % 2
    #codeword = signal + str(p1) + str(p2) + str(p3)
    syndrom = str((int(arr[4]) + p1) % 2) + str((int(arr[5]) + p2) % 2) + str((int(arr[6]) + p3) % 2)
    p_false = np.zeros(3)
    for i in range(3):
        p_false[i] = int(syndrom[i])
    if np.sum(p_false) == 2:
        if p_false[0]:
            if p_false[1]:
                s2 = (s2 + 1) % 2
            else:
                s1 = (s1 + 1) % 2
        else:
            s4 = (s4 + 1) % 2
    elif np.sum(p_false) == 3:
        s3 = (s3 + 1) % 2
    decode = str(s1) + str(s2) + str(s3) + str(s4)
    return decode




if __name__ == "__main__":
    # --- Source Coding and Reversible (Lossless) Data Compression ---

    print("\nQuestion 1")
    probabilities = [0.05,0.10,0.15,0.15,0.2,0.35] #[0.05,0.15,0.35,0.15,0.2,0.10]
    print(probabilities)
    huffman_codes = huffman(probabilities)
    print(huffman_codes)


    print("\nQuestion 2")
    input = "1011010100010"
    res, lempel_ziv_dict = lempelZiv(input)
    print(res)
    print(lempel_ziv_dict)


    print("\nQuesion 4")
    input = "abracadabrad" #"azertya"#"abroabro"#"azertyy"#"abraabra"
    res1 = LZ77(input, sliding_window_size = 7, decimal = True)
    res2 = LZ77(input, sliding_window_size=7)
    print(res1)
    print(res2)


    print("\nLoad Genome")
    codon_length = 3
    genome = load_text_sample()
    genome_length = len(genome)
    if genome_length > 0:
        print('Genome successfully loaded (starts with {})'.format(genome[:10]))

    binary_genome = load_binary_text_sample(spaces = False)
    binary_genome_length = len(binary_genome)
    if binary_genome_length > 0:
        print('Binary genome successfully loaded (starts with {})'.format(binary_genome[:20]))
    print(f"Length of binary genome: {binary_genome_length}")


    print("\nQuestion 5")
    # Compute the probability of occurence of each codon
    n_codons = int(np.ceil(len(genome) / codon_length))
    marginal_probabilities = dict()
    for i in range(0, len(genome), codon_length):
        codon = genome[i:i+codon_length]
        if codon not in marginal_probabilities:
            marginal_probabilities[codon] = 1
        else:
            marginal_probabilities[codon] += 1
    marginal_probabilities = {k: v / n_codons for k, v in marginal_probabilities.items()}
    #print(marginal_probabilities)
    
    # Assign Huffman codes
    huffman_codes = huffman(marginal_probabilities.values())
    huffman_codes = dict(zip(marginal_probabilities.keys(), huffman_codes))
    #print(codes)
    
    # Translate the genome with the Huffman codes
    huffman_genome = ""
    for i in range(0, len(genome), codon_length):
        codon = genome[i:i + codon_length]
        huffman_genome += huffman_codes[codon]
    
    # avg_coded_codon_length = sum([len(code) * proba for code, proba in zip(codes.values(), marginal_probabilities.values())])
    # print(f"Average coded codon length: {avg_coded_codon_length}")
    # gene_alphabet_size = 4 # 'A', 'T', 'G', 'C'
    # code_alphabet_size = 2 # '0', '1'
    # Length and compression rate
    huffman_genome_length = len(huffman_genome)
    print("Length of coded genome: {}".format(huffman_genome_length))
    huffman_compression_rate = binary_genome_length / huffman_genome_length #* np.log(gene_alphabet_size) / np.log(code_alphabet_size)
    print(f"Compression rate: {huffman_compression_rate}")


    print("\nQuestion 6")
    probabilities = np.fromiter(marginal_probabilities.values(), dtype=float)
    huffman_expected_length = sum([len(code) * proba for code, proba in zip(huffman_codes.values(), marginal_probabilities.values())]) #expected_huffman_length(probabilities)
    print(f"Expected huffman code length: {huffman_expected_length}")
    huffman_empirical_avg_length = huffman_genome_length / n_codons
    print(f"Empirical huffman code average length: {huffman_empirical_avg_length}")
    lower_bound_huffman_length = entropy(probabilities)
    upper_bound_huffman_length = lower_bound_huffman_length + 1
    print(f"Bounds on the huffman code length: [{lower_bound_huffman_length}, {upper_bound_huffman_length}[")


    print("\nQuestion 7")
    huffman_genome_empirical_avg_length = {}
    huffman_genome = ""
    for i in range(0, len(genome), codon_length):
        genome_length = i + codon_length
        codon = genome[i:genome_length]
        huffman_genome += huffman_codes[codon]
        n_codons = int(np.ceil(genome_length / codon_length))
        huffman_genome_empirical_avg_length[i + codon_length] = len(huffman_genome) / n_codons
    
    plt.figure()
    plt.plot(huffman_genome_empirical_avg_length.keys(), huffman_genome_empirical_avg_length.values())
    plt.ylabel("Huffman Code Empirical Average Length")
    plt.xlabel("Input Genome Length")
    #plt.savefig("Huffman Code Empirical Average Length.pdf")
    plt.show()


    print("\nQuestion 9")
    lempel_ziv_genome = lempelZiv(genome)
    binary_lempel_ziv_genome = ''.join(lempel_ziv_genome[i] if lempel_ziv_genome[i] == '0' or
                                                               lempel_ziv_genome[i] == '1'
                                       else binarize(ord(lempel_ziv_genome[i]), nb=8)
                                       for i in range(len(lempel_ziv_genome)))
    binary_lempel_ziv_genome_length = len(binary_lempel_ziv_genome)
    print(f"Total length of the lempel ziv encoded genome: {binary_lempel_ziv_genome_length}")
    binary_lempel_ziv_compression_rate = binary_genome_length / binary_lempel_ziv_genome_length
    print(f"Compression rate with lempel ziv: {binary_lempel_ziv_compression_rate}")


    # print("\nQuestion 10")
    # sliding_window_size = 4000
    # #readable_LZ77_genome = LZ77(genome, sliding_window_size=sliding_window_size, decimal=True)
    # LZ77_genome = LZ77(genome, sliding_window_size=sliding_window_size)
    # binary_LZ77_genome = ''.join(LZ77_genome[i] if LZ77_genome[i] == '0' or
    #                                                LZ77_genome[i] == '1'
    #                                    else binarize(ord(LZ77_genome[i]), nb=8)
    #                                    for i in range(len(LZ77_genome)))
    # binary_LZ77_genome_length = len(binary_LZ77_genome)
    # print(f"Total length of the LZ77 encoded genome: {binary_LZ77_genome_length}")
    # binary_LZ77_compression_rate = binary_genome_length / binary_LZ77_genome_length
    # print(f"Compression rate with LZ77: {binary_LZ77_compression_rate}")


    # print("\nQuestion 12")
    # # Takes a long time ...
    # sliding_window_size = 4000
    # LZ77_Huffman_genome = LZ77_Huffman(genome, sliding_window_size=sliding_window_size)
    # binary_LZ77_Huffman_genome = ''.join(LZ77_Huffman_genome[i] if LZ77_Huffman_genome[i] == '0' or
    #                                                                LZ77_Huffman_genome[i] == '1'
    #                              else binarize(ord(LZ77_Huffman_genome[i]), nb=8)
    #                              for i in range(len(LZ77_Huffman_genome)))
    # binary_LZ77_Huffman_genome_length = len(binary_LZ77_Huffman_genome)
    # print(f"Total length of the LZ77-Huffman encoded genome: {binary_LZ77_Huffman_genome_length}")
    # binary_LZ77_Huffman_compression_rate = binary_genome_length / binary_LZ77_Huffman_genome_length
    # print(f"Compression rate with LZ77-Huffman: {binary_LZ77_Huffman_compression_rate}")


    # print("\nQuestion 13")
    # # Takes a very very long time...
    # sliding_window_size_max = 6000
    # sliding_window_size_step = 200
    # sliding_window_sizes = np.arange(sliding_window_size_step, sliding_window_size_max, sliding_window_size_step)
    # binary_LZ77_genome_lengths = []
    # binary_LZ77_compression_rates = []
    # binary_LZ77_Huffman_genome_lengths = []
    # binary_LZ77_Huffman_compression_rates = []
    #
    # # readable_LZ77_genome = LZ77(genome, sliding_window_size = sliding_window_size, decimal = True)
    # for sliding_window_size in sliding_window_sizes:
    #     LZ77_genome = LZ77(genome, sliding_window_size=sliding_window_size)
    #     binary_LZ77_genome = ''.join(LZ77_genome[i] if LZ77_genome[i] == '0' or
    #                                                    LZ77_genome[i] == '1'
    #                                  else binarize(ord(LZ77_genome[i]), nb=8)
    #                                  for i in range(len(LZ77_genome)))
    #     binary_LZ77_genome_length = len(binary_LZ77_genome)
    #     binary_LZ77_genome_lengths.append(binary_LZ77_genome_length)
    #     binary_LZ77_compression_rate = binary_genome_length / binary_LZ77_genome_length
    #     binary_LZ77_compression_rates.append(binary_LZ77_compression_rate)
    #
    #     LZ77_Huffman_genome = LZ77_Huffman(LZ77_genome, sliding_window_size=sliding_window_size, is_LZ77_input=True)
    #     binary_LZ77_Huffman_genome = ''.join(LZ77_Huffman_genome[i] if LZ77_Huffman_genome[i] == '0' or
    #                                                                    LZ77_Huffman_genome[i] == '1'
    #                                          else binarize(ord(LZ77_Huffman_genome[i]), nb=8)
    #                                          for i in range(len(LZ77_Huffman_genome)))
    #     binary_LZ77_Huffman_genome_length = len(binary_LZ77_Huffman_genome)
    #     binary_LZ77_Huffman_genome_lengths.append(binary_LZ77_Huffman_genome_length)
    #     binary_LZ77_Huffman_compression_rate = binary_genome_length / binary_LZ77_Huffman_genome_length
    #     binary_LZ77_Huffman_compression_rates.append(binary_LZ77_Huffman_compression_rate)
    #
    # print("Sliding Window Size \t LZ77 Length \t LZ77-Huffman Length \t LZ77 Compression Rate \t LZ77-Huffman Compression Rate \n")
    # for i in range(len(binary_LZ77_genome_lengths)):
    #     print(f"{str(sliding_window_sizes[i]):23s}"
    #           f"{str(binary_LZ77_genome_lengths[i]):23s}"
    #           f"{str(binary_LZ77_Huffman_genome_lengths[i]):23s}"
    #           f"{str(binary_LZ77_compression_rates[i]):23s}"
    #           f"{str(binary_LZ77_Huffman_compression_rates[i]):23s}")
    #
    # plt.figure()
    # plt.plot(sliding_window_sizes, binary_LZ77_genome_lengths)
    # plt.ylabel("LZ77 Genome Length")
    # plt.xlabel("Sliding Window Size")
    # plt.savefig("LZ77 Genome Length.pdf")
    #
    # plt.figure()
    # plt.plot(sliding_window_sizes, binary_LZ77_compression_rates)
    # plt.ylabel("LZ77 Compression Rate")
    # plt.xlabel("Sliding Window Size")
    # plt.savefig("LZ77 Compression Rate.pdf")
    #
    # plt.figure()
    # plt.plot(sliding_window_sizes, binary_LZ77_Huffman_genome_lengths)
    # plt.ylabel("LZ77-Huffman Genome Length")
    # plt.xlabel("Sliding Window Size")
    # plt.savefig("LZ77-Huffman Genome Length.pdf")
    #
    # plt.figure()
    # plt.plot(sliding_window_sizes, binary_LZ77_Huffman_compression_rates)
    # plt.ylabel("LZ77-Huffman Compression Rate")
    # plt.xlabel("Sliding Window Size")
    # plt.savefig("LZ77-Huffman Compression Rate.pdf")


    # --- Channel Coding ---

    print("\nChannel Coding")

    # Load Sound
    r, s = load_wav()
    
    # Question 15
    # Plot the Sound
    print("\nQuestion 15")
    plt.figure()
    plt.plot(s)
    plt.title("Sound without noise")
    plt.show()
    
    # Question 16
    # Encode in Binary
    # Number of bits = 8 since the values range from 0 to 255, so it fits on 8 bits (2^8)
    print("Question 16")
    word1 = []
    for i in s:
        word1.extend(binarize(i, 8))
    # print(word1[0:16])
    
    # Question 17
    # Add Noise
    print("\nQuestion 17")
    nb_bits_flipped = 0
    for i in range(len(word1)):
        if np.random.binomial(n = 1, p = 0.01, size = 1):
            #print(i)
            word1[i] = str((int(word1[i]) + 1) % 2) # flip a bit
            nb_bits_flipped += 1
    # Noise added
    print(f"Number of value flipped due to Noise: {nb_bits_flipped}")
    #print(word1)
    
    # Decode the Noisy Sound
    res = []
    for i in range(len(s)):
        tmp = word1[8 * i]
        for j in np.arange(1, 8):
            tmp += word1[8 * i + j]
        res.append(bin_to_dec(tmp))
    #print(res)
    
    # Plot the Decoded Noisy Sound
    plt.figure()
    plt.plot(res)
    plt.title("Sound with noise without Hamming")
    plt.show()
    save_wav("Q17.wav", r, np.array(res, dtype='uint8'))
    
    # Question 18
    print("\nQuestion 18")
    word2 = []
    for i in s:
        word2.extend(binarize(i, 8))
    
    add_red = []
    for i in range(2 * len(s)):
        tmp = word2[4 * i]
        for j in np.arange(1, 4):
            tmp += word2[4 * i + j]
        add_red.extend(add_redundancy(tmp))
    #print(add_red)
    
    # Question 19
    print("\nQuestion 19")
    # Noise
    nb_bits_flipped = 0
    for i in range(len(word2)):
        if np.random.binomial(n = 1, p = 0.01, size = 1):
            #print(i)
            nb_bits_flipped += 1
            add_red[i] = str((int(add_red[i]) + 1) % 2)
    # Noise added
    print(f"Number of value flipped due to Noise: {nb_bits_flipped}")
    
    # Decode the Noisy Sound
    dec_redundancy = []
    for i in range(2 * len(s)):
        tmp = add_red[7 * i]
        for j in np.arange(1, 7):
            tmp += add_red[7 * i + j]
        dec_redundancy.extend(decode_redundancy(tmp))
    #print(dec_redundancy)
    
    res2 = []
    for i in range(len(s)):
        tmp = dec_redundancy[8 * i]
        for j in np.arange(1, 8):
            tmp += dec_redundancy[8 * i + j]
        res2.append(bin_to_dec(tmp))
    #print(res2)
    print(f"Number of value flipped due to Noise, but with Haming Code: {np.sum(res2!=s)}")
    
    plt.figure()
    plt.plot(res2)
    plt.title("Sound with noise with Hamming")
    plt.show()
    save_wav("Q19.wav", r, np.array(res2, dtype='uint8'))

    # # --- Channel Coding Annex on 7bits ---

    tpm_array = []
    tmp = np.array(s)
    print("\nOn 7 bits")
    for i in range(256):
        tpm_array.append(len(np.where(tmp==i)[0]))
    # print(tpm_array) 
    tmp = np.asarray(tpm_array)
    decodage_lenght = np.where(tmp!=0)[0] # in order to switch to 7bits
    # print(decodage_lenght)
    s_bis = np.zeros(len(s),dtype=int)
    for i in range(len(s)):
        s_bis[i] = np.where(decodage_lenght==s[i])[0]
    # print(s_bis)
    # for the graphic of value used
    # tmp = []
    # for i in range(len(tpm_array)):
    #     tmp.append(1 if tpm_array[i] != 0 else 0)
    
    # plt.figure()
    # plt.plot(tmp)
    # plt.show()
    # plt.title("Occurancies of values")

    word1 = []
    for i in s_bis:
        word1.extend(binarize(i,7))
    # print(word1[0:16])

    print("\nQuestion 17 bis")
    nb_bits_flipped = 0
    for i in range(len(word1)):
        if np.random.binomial(n = 1, p = 0.01, size = 1):
            nb_bits_flipped += 1
            word1[i] = str((int(word1[i]) + 1) % 2)

    print(f"Number of value flipped due to Noise: {nb_bits_flipped}")
    # print(word1)
    res = []
    for i in range(len(s)):
        tmp = word1[7 * i]
        for j in np.arange(1, 7):
            tmp += word1[7 * i + j]
        res.append(bin_to_dec(tmp))
    # print(res)

    res_bis = []
    for i in range(len(res)):
        if ( res[i] >= len(decodage_lenght)):
            res_bis.append(decodage_lenght[res[i]-64])
        else :
            res_bis.append(decodage_lenght[res[i]])
    # print(res_bis)
    
    plt.figure()
    plt.plot(res_bis)
    # print(np.sum(res_bis!=s))
    # plt.show()
    save_wav('Q17bis.wav', r, np.array(res_bis, dtype='uint8'))

    # Question 18
    print("Question 18 bis")
    word2 = []
    for i in s_bis:
        word2.extend(binarize(i,7))
    # print(word2[0:16])
    
    add_red = []
    for i in range(2 * len(s)):
        if(len(word2) < 4 * i+4):
            # it means that the number is not on 4 bits, we choose to not decode it
            tmp = "0000"
        else :
            tmp = word2[4 * i]
            for j in np.arange(1, 4):
                tmp += word2[4 * i + j]
        add_red.extend(add_redundancy(tmp))
    # print(add_red)
    
    # Question 19
    print("\nQuestion 19 bis")
    # Noise
    for i in range(len(word2)):
        if np.random.binomial(n = 1, p = 0.01, size = 1):
            nb_bits_flipped += 1
            add_red[i] = str((int(add_red[i]) + 1) % 2)
    # Noise added
    dec_redundancy = []
    for i in range(2 * len(s)):
        tmp = add_red[7 * i]
        for j in np.arange(1, 7):
            tmp += add_red[7 * i + j]
        dec_redundancy.extend(decode_redundancy(tmp))
    # print(dec_redundancy)
    
    res2 = []
    for i in range(len(s)):
        tmp = dec_redundancy[7 * i]
        for j in np.arange(1, 7):
            tmp += dec_redundancy[7 * i + j]
        res2.append(bin_to_dec(tmp))
    # print(res2)
    
    res_bis2 = []
    for i in range(len(res2)):
        if ( res2[i] >= len(decodage_lenght)):
            res_bis2.append(decodage_lenght[res2[i]-64])
        else :
            res_bis2.append(decodage_lenght[res2[i]])
    # print(res_bis2)
    print(f"Number of value flipped due to Noise, but with Haming Code: {np.sum(res_bis2!=s)}")

    plt.figure()
    plt.plot(res_bis2)
    # print(np.sum(res_bis2!=s))
    save_wav('Q19bis.wav', r, np.array(res_bis2, dtype='uint8'))
    plt.show()