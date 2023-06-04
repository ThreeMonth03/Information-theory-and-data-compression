import contextlib,sys
import argparse
import six
class ArithmeticCoderBase:
	
	# Constructs an arithmetic coder, which initializes the code range.
	def __init__(self, numbits):
		if numbits < 1:
			raise ValueError("State size out of range")
		
		# -- Configuration fields --
		# Number of bits for the 'low' and 'high' state variables. Must be at least 1.
		# - Larger values are generally better - they allow a larger maximum frequency total (maximum_total),
		#   and they reduce the approximation error inherent in adapting fractions to integers;
		#   both effects reduce the data encoding loss and asymptotically approach the efficiency
		#   of arithmetic coding using exact fractions.
		# - But larger state sizes increase the computation time for integer arithmetic,
		#   and compression gains beyond ~30 bits essentially zero in real-world applications.
		# - Python has native bigint arithmetic, so there is no upper limit to the state size.
		#   For Java and C++ where using native machine-sized integers makes the most sense,
		#   they have a recommended value of num_state_bits=32 as the most versatile setting.
		self.num_state_bits = numbits
		# Maximum range (high+1-low) during coding (trivial), which is 2^num_state_bits = 1000...000.
		self.full_range = 1 << self.num_state_bits
		# The top bit at width num_state_bits, which is 0100...000.
		self.half_range = self.full_range >> 1  # Non-zero
		# The second highest bit at width num_state_bits, which is 0010...000. This is zero when num_state_bits=1.
		self.quarter_range = self.half_range >> 1  # Can be zero
		# Minimum range (high+1-low) during coding (non-trivial), which is 0010...010.
		self.minimum_range = self.quarter_range + 2  # At least 2
		# Maximum allowed total from a frequency table at all times during coding. This differs from Java
		# and C++ because Python's native bigint avoids constraining the size of intermediate computations.
		self.maximum_total = self.minimum_range
		# Bit mask of num_state_bits ones, which is 0111...111.
		self.state_mask = self.full_range - 1
		
		# -- State fields --
		# Low end of this arithmetic coder's current range. Conceptually has an infinite number of trailing 0s.
		self.low = 0
		# High end of this arithmetic coder's current range. Conceptually has an infinite number of trailing 1s.
		self.high = self.state_mask
	
	
	# Updates the code range (low and high) of this arithmetic coder as a result
	# of processing the given symbol with the given frequency table.
	# Invariants that are true before and after encoding/decoding each symbol
	# (letting full_range = 2^num_state_bits):
	# - 0 <= low <= code <= high < full_range. ('code' exists only in the decoder.)
	#   Therefore these variables are unsigned integers of num_state_bits bits.
	# - low < 1/2 * full_range <= high.
	#   In other words, they are in different halves of the full range.
	# - (low < 1/4 * full_range) || (high >= 3/4 * full_range).
	#   In other words, they are not both in the middle two quarters.
	# - Let range = high - low + 1, then full_range/4 < minimum_range
	#   <= range <= full_range. These invariants for 'range' essentially
	#   dictate the maximum total that the incoming frequency table can have.
	def update(self, freqs, symbol):
		# State check
		low = self.low
		high = self.high
		if low >= high or (low & self.state_mask) != low or (high & self.state_mask) != high:
			raise AssertionError("Low or high out of range")
		range = high - low + 1
		if not (self.minimum_range <= range <= self.full_range):
			raise AssertionError("Range out of range")
		
		# Frequency table values check
		total = freqs.get_total()
		symlow = freqs.get_low(symbol)
		symhigh = freqs.get_high(symbol)
		if symlow == symhigh:
			raise ValueError("Symbol has zero frequency")
		if total > self.maximum_total:
			raise ValueError("Cannot code symbol because total is too large")
		
		# Update range
		newlow  = low + symlow  * range // total
		newhigh = low + symhigh * range // total - 1
		self.low = newlow
		self.high = newhigh
		
		# While low and high have the same top bit value, shift them out
		while ((self.low ^ self.high) & self.half_range) == 0:
			self.shift()
			self.low  = ((self.low  << 1) & self.state_mask)
			self.high = ((self.high << 1) & self.state_mask) | 1
		# Now low's top bit must be 0 and high's top bit must be 1
		
		# While low's top two bits are 01 and high's are 10, delete the second highest bit of both
		while (self.low & ~self.high & self.quarter_range) != 0:
			self.underflow()
			self.low = (self.low << 1) ^ self.half_range
			self.high = ((self.high ^ self.half_range) << 1) | self.half_range | 1
	
	
	# Called to handle the situation when the top bit of 'low' and 'high' are equal.
	def shift(self):
		raise NotImplementedError()
	
	
	# Called to handle the situation when low=01(...) and high=10(...).
	def underflow(self):
		raise NotImplementedError()



# Encodes symbols and writes to an arithmetic-coded bit stream.
class ArithmeticEncoder(ArithmeticCoderBase):
	
	# Constructs an arithmetic coding encoder based on the given bit output stream.
	def __init__(self, numbits, bitout):
		super(ArithmeticEncoder, self).__init__(numbits)
		# The underlying bit output stream.
		self.output = bitout
		# Number of saved underflow bits. This value can grow without bound.
		self.num_underflow = 0
	
	
	# Encodes the given symbol based on the given frequency table.
	# This updates this arithmetic coder's state and may write out some bits.
	def write(self, freqs, symbol):
		if not isinstance(freqs, CheckedFrequencyTable):
			freqs = CheckedFrequencyTable(freqs)
		self.update(freqs, symbol)
	
	
	# Terminates the arithmetic coding by flushing any buffered bits, so that the output can be decoded properly.
	# It is important that this method must be called at the end of the each encoding process.
	# Note that this method merely writes data to the underlying output stream but does not close it.
	def finish(self):
		self.output.write(1)
	
	
	def shift(self):
		bit = self.low >> (self.num_state_bits - 1)
		self.output.write(bit)
		
		# Write out the saved underflow bits
		for _ in range(self.num_underflow):
			self.output.write(bit ^ 1)
		self.num_underflow = 0
	
	
	def underflow(self):
		self.num_underflow += 1



# ---- Frequency table classes ----

# A table of symbol frequencies. The table holds data for symbols numbered from 0
# to get_symbol_limit()-1. Each symbol has a frequency, which is a non-negative integer.
# Frequency table objects are primarily used for getting cumulative symbol
# frequencies. These objects can be mutable depending on the implementation.
class FrequencyTable:
	
	# Returns the number of symbols in this frequency table, which is a positive number.
	def get_symbol_limit(self):
		raise NotImplementedError()
	
	# Returns the frequency of the given symbol. The returned value is at least 0.
	def get(self, symbol):
		raise NotImplementedError()
	
	# Sets the frequency of the given symbol to the given value.
	# The frequency value must be at least 0.
	def set(self, symbol, freq):
		raise NotImplementedError()
	
	# Increments the frequency of the given symbol.
	def increment(self, symbol):
		raise NotImplementedError()
	
	# Returns the total of all symbol frequencies. The returned value is at
	# least 0 and is always equal to get_high(get_symbol_limit() - 1).
	def get_total(self):
		raise NotImplementedError()
	
	# Returns the sum of the frequencies of all the symbols strictly
	# below the given symbol value. The returned value is at least 0.
	def get_low(self, symbol):
		raise NotImplementedError()
	
	# Returns the sum of the frequencies of the given symbol
	# and all the symbols below. The returned value is at least 0.
	def get_high(self, symbol):
		raise NotImplementedError()



# An immutable frequency table where every symbol has the same frequency of 1.
# Useful as a fallback model when no statistics are available.
class FlatFrequencyTable(FrequencyTable):
	
	# Constructs a flat frequency table with the given number of symbols.
	def __init__(self, numsyms):
		if numsyms < 1:
			raise ValueError("Number of symbols must be positive")
		self.numsymbols = numsyms  # Total number of symbols, which is at least 1
	
	# Returns the number of symbols in this table, which is at least 1.
	def get_symbol_limit(self):
		return self.numsymbols
	
	# Returns the frequency of the given symbol, which is always 1.
	def get(self, symbol):
		self._check_symbol(symbol)
		return 1
	
	# Returns the total of all symbol frequencies, which is
	# always equal to the number of symbols in this table.
	def get_total(self):
		return self.numsymbols
	
	# Returns the sum of the frequencies of all the symbols strictly below
	# the given symbol value. The returned value is equal to 'symbol'.
	def get_low(self, symbol):
		self._check_symbol(symbol)
		return symbol
	
	
	# Returns the sum of the frequencies of the given symbol and all
	# the symbols below. The returned value is equal to 'symbol' + 1.
	def get_high(self, symbol):
		self._check_symbol(symbol)
		return symbol + 1
	
	
	# Returns silently if 0 <= symbol < numsymbols, otherwise raises an exception.
	def _check_symbol(self, symbol):
		if not (0 <= symbol < self.numsymbols):
			raise ValueError("Symbol out of range")
	
	# Returns a string representation of this frequency table. The format is subject to change.
	def __str__(self):
		return "FlatFrequencyTable={}".format(self.numsymbols)
	
	# Unsupported operation, because this frequency table is immutable.
	def set(self, symbol, freq):
		raise NotImplementedError()
	
	# Unsupported operation, because this frequency table is immutable.
	def increment(self, symbol):
		raise NotImplementedError()



# A mutable table of symbol frequencies. The number of symbols cannot be changed
# after construction. The current algorithm for calculating cumulative frequencies
# takes linear time, but there exist faster algorithms such as Fenwick trees.
class SimpleFrequencyTable(FrequencyTable):
	
	# Constructs a simple frequency table in one of two ways:
	# - SimpleFrequencyTable(sequence):
	#   Builds a frequency table from the given sequence of symbol frequencies.
	#   There must be at least 1 symbol, and no symbol has a negative frequency.
	# - SimpleFrequencyTable(freqtable):
	#   Builds a frequency table by copying the given frequency table.
	def __init__(self, freqs):
		if isinstance(freqs, FrequencyTable):
			numsym = freqs.get_symbol_limit()
			self.frequencies = [freqs.get(i) for i in range(numsym)]
		else:  # Assume it is a sequence type
			self.frequencies = list(freqs)  # Make copy
		
		# 'frequencies' is a list of the frequency for each symbol.
		# Its length is at least 1, and each element is non-negative.
		if len(self.frequencies) < 1:
			raise ValueError("At least 1 symbol needed")
		for freq in self.frequencies:
			if freq < 0:
				raise ValueError("Negative frequency")
		
		# Always equal to the sum of 'frequencies'
		self.total = sum(self.frequencies)
		
		# cumulative[i] is the sum of 'frequencies' from 0 (inclusive) to i (exclusive).
		# Initialized lazily. When it is not None, the data is valid.
		self.cumulative = None
	
	
	# Returns the number of symbols in this frequency table, which is at least 1.
	def get_symbol_limit(self):
		return len(self.frequencies)
	
	
	# Returns the frequency of the given symbol. The returned value is at least 0.
	def get(self, symbol):
		self._check_symbol(symbol)
		return self.frequencies[symbol]
	
	
	# Sets the frequency of the given symbol to the given value. The frequency value
	# must be at least 0. If an exception is raised, then the state is left unchanged.
	def set(self, symbol, freq):
		self._check_symbol(symbol)
		if freq < 0:
			raise ValueError("Negative frequency")
		temp = self.total - self.frequencies[symbol]
		assert temp >= 0
		self.total = temp + freq
		self.frequencies[symbol] = freq
		self.cumulative = None
	
	
	# Increments the frequency of the given symbol.
	def increment(self, symbol):
		self._check_symbol(symbol)
		self.total += 1
		self.frequencies[symbol] += 1
		self.cumulative = None
	
	
	# Returns the total of all symbol frequencies. The returned value is at
	# least 0 and is always equal to get_high(get_symbol_limit() - 1).
	def get_total(self):
		return self.total
	
	
	# Returns the sum of the frequencies of all the symbols strictly
	# below the given symbol value. The returned value is at least 0.
	def get_low(self, symbol):
		self._check_symbol(symbol)
		if self.cumulative is None:
			self._init_cumulative()
		return self.cumulative[symbol]
	
	
	# Returns the sum of the frequencies of the given symbol
	# and all the symbols below. The returned value is at least 0.
	def get_high(self, symbol):
		self._check_symbol(symbol)
		if self.cumulative is None:
			self._init_cumulative()
		return self.cumulative[symbol + 1]
	
	
	# Recomputes the array of cumulative symbol frequencies.
	def _init_cumulative(self):
		cumul = [0]
		sum = 0
		for freq in self.frequencies:
			sum += freq
			cumul.append(sum)
		assert sum == self.total
		self.cumulative = cumul
	
	
	# Returns silently if 0 <= symbol < len(frequencies), otherwise raises an exception.
	def _check_symbol(self, symbol):
		if not (0 <= symbol < len(self.frequencies)):
			raise ValueError("Symbol out of range")
	
	
	# Returns a string representation of this frequency table,
	# useful for debugging only, and the format is subject to change.
	def __str__(self):
		result = ""
		for (i, freq) in enumerate(self.frequencies):
			result += "{}\t{}\n".format(i, freq)
		return result



# A wrapper that checks the preconditions (arguments) and postconditions (return value) of all
# the frequency table methods. Useful for finding faults in a frequency table implementation.
class CheckedFrequencyTable(FrequencyTable):
	
	def __init__(self, freqtab):
		# The underlying frequency table that holds the data
		self.freqtable = freqtab
	
	
	def get_symbol_limit(self):
		result = self.freqtable.get_symbol_limit()
		if result <= 0:
			raise AssertionError("Non-positive symbol limit")
		return result
	
	
	def get(self, symbol):
		result = self.freqtable.get(symbol)
		if not self._is_symbol_in_range(symbol):
			raise AssertionError("ValueError expected")
		if result < 0:
			raise AssertionError("Negative symbol frequency")
		return result
	
	
	def get_total(self):
		result = self.freqtable.get_total()
		if result < 0:
			raise AssertionError("Negative total frequency")
		return result
	
	
	def get_low(self, symbol):
		if self._is_symbol_in_range(symbol):
			low   = self.freqtable.get_low (symbol)
			high  = self.freqtable.get_high(symbol)
			if not (0 <= low <= high <= self.freqtable.get_total()):
				raise AssertionError("Symbol low cumulative frequency out of range")
			return low
		else:
			self.freqtable.get_low(symbol)
			raise AssertionError("ValueError expected")
	
	
	def get_high(self, symbol):
		if self._is_symbol_in_range(symbol):
			low   = self.freqtable.get_low (symbol)
			high  = self.freqtable.get_high(symbol)
			if not (0 <= low <= high <= self.freqtable.get_total()):
				raise AssertionError("Symbol high cumulative frequency out of range")
			return high
		else:
			self.freqtable.get_high(symbol)
			raise AssertionError("ValueError expected")
	
	
	def __str__(self):
		return "CheckedFrequencyTable (" + str(self.freqtable) + ")"
	
	
	def set(self, symbol, freq):
		self.freqtable.set(symbol, freq)
		if not self._is_symbol_in_range(symbol) or freq < 0:
			raise AssertionError("ValueError expected")
	
	
	def increment(self, symbol):
		self.freqtable.increment(symbol)
		if not self._is_symbol_in_range(symbol):
			raise AssertionError("ValueError expected")
	
	
	def _is_symbol_in_range(self, symbol):
		return 0 <= symbol < self.get_symbol_limit()



# ---- Bit-oriented I/O streams ----

# A stream of bits that can be read. Because they come from an underlying byte stream,
# the total number of bits is always a multiple of 8. The bits are read in big endian.
class BitInputStream:
	
	# Constructs a bit input stream based on the given byte input stream.
	def __init__(self, inp):
		# The underlying byte stream to read from
		self.input = inp
		# Either in the range [0x00, 0xFF] if bits are available, or -1 if end of stream is reached
		self.currentbyte = 0
		# Number of remaining bits in the current byte, always between 0 and 7 (inclusive)
		self.numbitsremaining = 0
	
	
	# Reads a bit from this stream. Returns 0 or 1 if a bit is available, or -1 if
	# the end of stream is reached. The end of stream always occurs on a byte boundary.
	def read(self):
		if self.currentbyte == -1:
			return -1
		if self.numbitsremaining == 0:
			temp = self.input.read(1)
			if len(temp) == 0:
				self.currentbyte = -1
				return -1
			self.currentbyte = temp[0]
			self.numbitsremaining = 8
		assert self.numbitsremaining > 0
		self.numbitsremaining -= 1
		return (self.currentbyte >> self.numbitsremaining) & 1
	
	
	# Reads a bit from this stream. Returns 0 or 1 if a bit is available, or raises an EOFError
	# if the end of stream is reached. The end of stream always occurs on a byte boundary.
	def read_no_eof(self):
		result = self.read()
		if result != -1:
			return result
		else:
			raise EOFError()
	
	
	# Closes this stream and the underlying input stream.
	def close(self):
		self.input.close()
		self.currentbyte = -1
		self.numbitsremaining = 0



# A stream where bits can be written to. Because they are written to an underlying
# byte stream, the end of the stream is padded with 0's up to a multiple of 8 bits.
# The bits are written in big endian.
class BitOutputStream:
	
	# Constructs a bit output stream based on the given byte output stream.
	def __init__(self, out):
		self.output = out  # The underlying byte stream to write to
		self.currentbyte = 0  # The accumulated bits for the current byte, always in the range [0x00, 0xFF]
		self.numbitsfilled = 0  # Number of accumulated bits in the current byte, always between 0 and 7 (inclusive)
	
	
	# Writes a bit to the stream. The given bit must be 0 or 1.
	def write(self, b):
		if b not in (0, 1):
			raise ValueError("Argument must be 0 or 1")
		self.currentbyte = (self.currentbyte << 1) | b
		self.numbitsfilled += 1
		if self.numbitsfilled == 8:
			towrite = bytes((self.currentbyte,))
			self.output.write(towrite)
			self.currentbyte = 0
			self.numbitsfilled = 0
	
	
	# Closes this stream and the underlying output stream. If called when this
	# bit stream is not at a byte boundary, then the minimum number of "0" bits
	# (between 0 and 7 of them) are written as padding to reach the next byte boundary.
	def close(self):
		while self.numbitsfilled != 0:
			self.write(0)
		self.output.close()

# Returns a frequency table based on the bytes in the given file.
# Also contains an extra entry for symbol 256, whose frequency is set to 0.
def get_frequencies(filepath):
	freqs = SimpleFrequencyTable([0] * 257)
	with open(filepath, "rb") as input:
		while True:
			b = input.read(1)
			if len(b) == 0:
				break
			freqs.increment(b[0])
	return freqs


def write_frequencies(bitout, freqs):
	for i in range(256):
		write_int(bitout, 32, freqs.get(i))


def compress(freqs, inp, bitout):
	enc = ArithmeticEncoder(32, bitout)
	while True:
		symbol = inp.read(1)
		if len(symbol) == 0:
			break
		enc.write(freqs, symbol[0])
	enc.write(freqs, 256)  # EOF
	enc.finish()  # Flush remaining code bits


# Writes an unsigned integer of the given bit width to the given stream.
def write_int(bitout, numbits, value):
	for i in reversed(range(numbits)):
		bitout.write((value >> i) & 1)  # Big endian

class ArithmeticCoding():
    def __init__(self, bits=8):
        self.read_bits = bits
        if(self.read_bits == 8):
            self.interval_bits = 32
        if(self.read_bits == 1):
            self.interval_bits = 64
        
        #Interval control
        self.total_num = 1 << self.interval_bits #2^interval_bits
        self.msb = self.total_num >> 1 #MSB
        self.msb2 = self.msb >> 1 #2nd-MSB
        self.mask = self.total_num - 1 #Get the valid number
        self.low = 0 #lower bound
        self.high = self.mask #upper bound
        self.E3 = 0
        
    def compress(self,infilename, outfilename, print_info):
        freq = {}
        cum_freq = {}
        total_cum = 0
        if self.read_bits == 8:
            #Calculate the frequency and cumulative frequency
            freq, cum_freq, total_cum = self.get_frequencies_256(infilename, print_info)
            #Compress
            self.encode_256(infilename, outfilename, cum_freq, total_cum, print_info)
        #        if self.readbits == 1:
        #    freq = self.get_frequencies_2(infilename, print_info)
        
        #Write the header
        #First byte is read bit
        #Second byte is interval_bits
        #Third byte is 
        
        #Operate the algorithm
        #self.encode()

    def get_frequencies_256(self,infilename, print_info):
        frequencies = {}
        cum_frequencies = {}
        cum = 0
        #update frequency table
        for i in range(1<<self.read_bits):
            frequencies[i] = 1
            
        with open(infilename, "rb") as f:
            infile_data = f.read()
            infile_size = f.tell()
            for i in range(0, infile_size):
                symbol = infile_data[i]
                frequencies[symbol] += 1
                    
        #update cumulative frequency table
        cum_frequencies[-1] = cum
        for symbol in sorted(frequencies.keys()):
            cum += frequencies[symbol]
            cum_frequencies[symbol] = cum
        total_cum = cum
        
        if print_info:
            print("Frequency table")
            for symbol in sorted(frequencies.keys()):
                print(symbol,' : ',frequencies[symbol])

            print("Cumulative frequency table")
            for symbol in sorted(cum_frequencies.keys()):
                print(symbol,' : ',cum_frequencies[symbol])
                            
        return frequencies, cum_frequencies ,total_cum

    def encode_256(self, infilename, outfilename, cum_freq, total_cum, print_info):
        with open(infilename, "rb") as f:
            with open(outfilename, 'wb') as output:
                #Write the header
                output.write(six.int2byte(self.read_bits))#First byte is read bit
                output.write(six.int2byte(self.interval_bits))#Second byte is interval_bits
                #Then we store cum_freq table
                for i in range(1<<self.read_bits):
                    freq = cum_freq[i]
                    o4 = freq&255
                    freq = freq>>8
                    o3 = freq&255
                    freq = freq>>8
                    o2 = freq&255
                    freq = freq>>8
                    o1 = freq&255
                    output.write(six.int2byte(o1))
                    output.write(six.int2byte(o2))
                    output.write(six.int2byte(o3))
                    output.write(six.int2byte(o4))                    
                #read and compress    
                infile_data = f.read()
                infile_size = f.tell()
                coding = ""
                for i in range(0, infile_size):                
                    symbol = infile_data[i]
                    out_int = 0
                    coding += self.update(symbol,cum_freq,total_cum)
                    while len(coding)>8:
                        for j in range(8):
                            out_int = out_int<<1
                            if coding[j] == '1':
                                out_int = out_int|1
                        coding = coding[8:]
                        output.write(six.int2byte(out_int))
                        out_int = 0
                #deal with remaining code
                remaining_length = len(coding)
                out_int = 0
                for j in range(remaining_length):
                    out_int = out_int<<1
                    if coding[j] == '1':
                        out_int = out_int|1        
                out_int = out_int<<(8 - remaining_length)
                output.write(six.int2byte(out_int))
                #write the remaining length
                output.write(six.int2byte(remaining_length))
                #write the lower bound
                lb = []
                temp = self.low
                for i in range(self.interval_bits//8):
                    lb.append(temp&255)
                    temp = temp>>8
                lb.reverse()
                for bytes in lb:
                    output.write(six.int2byte(bytes))
                output_size = output.tell()
                if print_info:
                    header_bit_size = 8 + 8 + 256*32 + 8 + self.interval_bits
                    encoding_bit_size = output_size*8 - header_bit_size + remaining_length - 8
                    expected_coding_length = encoding_bit_size / infile_size
                    print("The number of bits that one symbol have: ",self.read_bits)
                    print("The number of bits that bound have: ",self.interval_bits)
                    print("Header size: ",header_bit_size/8,"bytes , ",header_bit_size,"bits")
                    print("Padding size: ",(8 - remaining_length)/8,"bytes , ",8 - remaining_length,"bits")
                    print("Encoding size: ",encoding_bit_size/8,"bytes , ",encoding_bit_size,"bits")
                    print("Expected coding length: ",expected_coding_length/8,"bytes , ",expected_coding_length,"bits")
                    print("Total size: ",output_size,"bytes , ",output_size*8,"bits")
                    print("Lower bound: {:032b}".format(self.low))
                    print("Upper bound: {:032b}".format(self.high))
    def update(self, symbol, cum_freq, total_cum):
        #update range
        cum_low = cum_freq[symbol - 1]
        cum_high = cum_freq[symbol]
        range = self.high - self.low + 1
        new_low = self.low + cum_low * range // total_cum
        new_high = self.low + cum_high * range // total_cum - 1
        self.low = new_low
        self.high = new_high
        #E1/E2/E3
        coding = ""
        #E1/E2
        while ((self.low ^ self.high) & self.msb) == 0:
            #update coding
            bit = self.low >> (self.interval_bits - 1)            
            coding += str(bit)
            coding += (str(bit ^ 1) * self.E3)
            self.E3 = 0 
            #update interval
            self.low  = ((self.low  << 1) & self.mask)
            self.high = ((self.high << 1) & self.mask) | 1
        #E3 
        while (self.low & ~self.high & self.msb2) != 0:
            self.E3 += 1
            self.low = (self.low << 1) ^ self.msb
            self.high = ((self.high ^ self.msb) << 1) | self.msb | 1        
        return coding
          
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('infilename' , help='the file we want to deal with')
    parser.add_argument('outfilename' , help='the output file name')
    parser.add_argument('--algorithm1', default='256_AC', help='Choose the first compress algorithm.')
    parser.add_argument('--algorithm2', default=None, help='Choose the second compress algorithm.')
    parser.add_argument('--print_info' , default=True, help='Print the compress information.')
    args = parser.parse_args()
    return args

def choose_algorithm(algorithm):
    if algorithm == '256_AC':
        return ArithmeticCoding(bits=8)
        #print("ArithmeticCoding(bits=8)")
    #if algorithm == 'binary_AC':
    #    return ArithmeticCoding(bits=1)
    #if algorithm == 'HC':
    #    return HuffmanCoding()
    else:
        raise NameError("The algorithm should be 256_AC, binary_AC, HC.")
    
def main():
    #parse args
    args = parse_args()
    
    ##Operate algorithm1
    algo1 = choose_algorithm(args.algorithm1)
    print("Compress ", args.infilename)
    ###Compress
    if (args.algorithm2!=None):
        algo1.compress(args.infilename, "temp.cps" , args.print_info)
        #print("algo1.compress(",args.infilename, "temp.cps" , args.print_info)
    else:
        algo1.compress(args.infilename, args.outfilename , args.print_info)
        #print("algo1.compress(",args.infilename, args.outfilename , args.print_info)
        
    if (args.algorithm2!=None):
        pass
    '''
	if len(args) != 2:
		sys.exit("Usage: python arithmetic-compress.py InputFile OutputFile")
	inputfile, outputfile = args
	
	# Read input file once to compute symbol frequencies
	freqs = get_frequencies(inputfile)
	freqs.increment(256)  # EOF symbol gets a frequency of 1
	
	# Read input file again, compress with arithmetic coding, and write output file
	with open(inputfile, "rb") as inp, \
			contextlib.closing(BitOutputStream(open(outputfile, "wb"))) as bitout:
		write_frequencies(bitout, freqs)
		compress(freqs, inp, bitout)
    '''
    
if __name__ == "__main__":
	main()
