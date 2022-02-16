 #!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read, write
from scipy.io import loadmat

""" 
This function gives the binary representation of value v (dec)
	using nb bits (0 by default corresponds to the minimal number of 
	bits required to represent v).

	examples: 
	binarize(2)
	>>> '10'

	binarize(2,nb=8)
	>>> '00000010'

"""
def binarize(v, nb=0):
    if nb == 0:
        return bin(v)[2:]
    else:
        return np.binary_repr(v,width=nb)

# This functions returns the decimal representation of a sequence (str) of bits b
def bin_to_dec(b):
    return int(b,2)

# This function loads the text sample
def load_text_sample(name="genome.txt"):
	f = np.genfromtxt("genome.txt",dtype='str')
	f = "".join(f)
	return f

# This function loads the text sample and outputs the binary version (8-bits representation).
# If spaces=True, then each byte is separated by a space 
def load_binary_text_sample(name="genome.txt",spaces=True):
	f = load_text_sample(name)
	binary_text = ''
	for c in f:
		if spaces:
			binary_text = binary_text + " " + binarize(ord(c),nb=8)
		else:
			binary_text = binary_text + binarize(ord(c),nb=8)
	if spaces:
		return binary_text[1:]
	else:
		return binary_text

# This function loads the sound signal (.wav)
def load_wav():
	rate, data = read('sound.wav')
	return rate, data

# This function save the sound signal (.wav)
def save_wav(filename,rate,data):
    write(filename, rate, data)

if __name__ == "__main__":

	f = load_text_sample()
	if len(f) > 0:
		print('Text successfully loaded (starts with {})'.format(f[:10]))

	# Takes some time.
	#f = load_binary_text_sample()
	#if len(f) > 0:
	#	print('Binary text successfully loaded (starts with {})'.format(f[:10]))

	r,s = load_wav()
	print(s)
	save_wav("test.wav",r,s)
