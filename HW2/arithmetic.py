import argparse
import six

####################################################
##                                                ##
##    Arithmetic coding with fixed probability    ##
##                                                ##
#################################################### 
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
        
        if self.read_bits == 1:
            #Calculate the frequency and cumulative frequency
            freq, cum_freq, total_cum = self.get_frequencies_2(infilename, print_info)
            #Compress
            self.encode_2(infilename, outfilename, cum_freq, total_cum, print_info)
        

    def get_frequencies_256(self,infilename, print_info):
        frequencies = {}
        cum_frequencies = {}
        cum = 0
        #update frequency table
        for i in range(1<<self.read_bits):
            frequencies[i] = 0
            
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

    def get_frequencies_2(self,infilename, print_info):
        frequencies = {}
        cum_frequencies = {}
        cum = 0
        #update frequency table
        for i in range(1<<self.read_bits):
            frequencies[i] = 0
            
        with open(infilename, "rb") as f:
            infile_data = f.read()
            infile_size = f.tell()
            for i in range(0, infile_size):
                symbol = infile_data[i]
                #for j in range(8):
                #    bit = (symbol >> (7-j)) & 1
                #    frequencies[bit] += 1
                frequencies[1] += symbol.bit_count()
            frequencies[0] = 8 * infile_size - frequencies[1]


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
                for symbol in range(1<<self.read_bits):
                    freq = cum_freq[symbol]
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
                    #header_bit_size = read bit + interval_bits + table size + padding size
                    header_bit_size = 8 + 8 + 256*32 + 8
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
                
    def encode_2(self, infilename, outfilename, cum_freq, total_cum, print_info):
        with open(infilename, "rb") as f:
            with open(outfilename, 'wb') as output:
                #Write the header
                output.write(six.int2byte(self.read_bits))#First byte is read bit
                output.write(six.int2byte(self.interval_bits))#Second byte is interval_bits
                #Then we store cum_freq table
                for symbol in range(1<<self.read_bits):
                    freq = cum_freq[symbol]
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
                
                #buffer = memoryview(infile_data)
                #np_bit_data = np.reshape(np.frombuffer(buffer, dtype=np.uint8),(1,-1))
                #buffer = None
                #print(np_data.shape)
                #print(np_data)   
                #np_bit_data = np.unpackbits(np_bit_data, axis=1)
                #print(np_bit_data.shape)
                #print(np_bit_data)
                #np_bit_data = np.reshape(np_bit_data, -1)
                #bit_size = np_bit_data.shape[0]
                #print(np_bit_data.shape)
                #print(np_bit_data)
                
                #for symbol in np_bit_data:  
                for i in range(0, infile_size):
                    symbol = infile_data[i]
                    for j in range(7,-1,-1):
                        bit = (symbol >> j) & 1
                        out_int = 0
                        coding += self.update(bit,cum_freq,total_cum)
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
                    #header_bit_size = read bit + interval_bits + table size + padding size
                    header_bit_size = 8 + 8 + 2*32 + 8
                    encoding_bit_size = output_size*8 - header_bit_size + remaining_length - 8
                    expected_coding_length = encoding_bit_size / infile_size/8
                    print("The number of bits that one symbol have: ",self.read_bits)
                    print("The number of bits that bound have: ",self.interval_bits)
                    print("Header size: ",header_bit_size/8,"bytes , ",header_bit_size,"bits")
                    print("Padding size: ",(8 - remaining_length)/8,"bytes , ",8 - remaining_length,"bits")
                    print("Encoding size: ",encoding_bit_size/8,"bytes , ",encoding_bit_size,"bits")
                    print("Expected coding length: ",expected_coding_length/8,"bytes , ",expected_coding_length,"bits")
                    print("Total size: ",output_size,"bytes , ",output_size*8,"bits")
                    print("Lower bound: {:064b}".format(self.low))
                    print("Upper bound: {:064b}".format(self.high))
                
                    
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

######################################
##                                  ##
##    Arithmetic coding with PPM    ##
##                                  ##
######################################
class PPM():
    def __init__(self, bits=8, order= -1):
        self.read_bits = bits
        self.escape_symbol = 1<<self.read_bits
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
        #Table
        self.order = order
        self.table = Table(self.order, self.escape_symbol)

        
    def compress(self,infilename, outfilename, print_info): 
        if self.read_bits == 8:
            self.encode_256(infilename, outfilename, print_info)   

        if self.read_bits == 1:
            self.encode_2(infilename, outfilename, print_info) 

    def encode_256(self, infilename, outfilename, print_info):
        with open(infilename, "rb") as f:
            with open(outfilename, 'wb') as output:
                #Write the header
                output.write(six.int2byte(self.read_bits))#First byte is read bit
                output.write(six.int2byte(self.interval_bits))#Second byte is interval_bits
                infile_data = f.read()
                infile_size = f.tell()
                history = []
                coding = ""
                for i in range(0, infile_size):
                    symbol = infile_data[i]
                    out_int = 0
                    coding += self.encoding_table(symbol, self.table, history)
                    self.table.increment_table(history, symbol)
                    while len(coding)>8:
                        for j in range(8):
                            out_int = out_int<<1
                            if coding[j] == '1':
                                out_int = out_int|1
                        coding = coding[8:]
                        output.write(six.int2byte(out_int))
                        out_int = 0
                    #update history
                    if self.order >= 1:
                        if len(history) == self.order:
                            history.pop()
                        history.insert(0, symbol)
                
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
                    #header_bit_size = read bit + interval_bits + padding size
                    header_bit_size = 8 + 8 + 8
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
    
    def encode_2(self, infilename, outfilename, print_info):
        with open(infilename, "rb") as f:
            with open(outfilename, 'wb') as output:
                #Write the header
                output.write(six.int2byte(self.read_bits))#First byte is read bit
                output.write(six.int2byte(self.interval_bits))#Second byte is interval_bits
                infile_data = f.read()
                infile_size = f.tell()
                history = []
                coding = ""
                for i in range(0, infile_size):
                    symbol = infile_data[i]
                    for j in range(7,-1,-1):
                        bit = (symbol >> j) & 1
                        out_int = 0
                        coding += self.encoding_table(bit, self.table, history)
                        self.table.increment_table(history, bit)
                        while len(coding)>8:
                            for j in range(8):
                                out_int = out_int<<1
                                if coding[j] == '1':
                                    out_int = out_int|1
                            coding = coding[8:]
                            output.write(six.int2byte(out_int))
                            out_int = 0
                        #update history
                        if self.order >= 1:
                            if len(history) == self.order:
                                history.pop()
                            history.insert(0, bit)
                
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
                    #header_bit_size = read bit + interval_bits + padding size
                    header_bit_size = 8 + 8 + 8
                    encoding_bit_size = output_size*8 - header_bit_size + remaining_length - 8
                    expected_coding_length = encoding_bit_size / infile_size/8
                    print("The number of bits that one symbol have: ",self.read_bits)
                    print("The number of bits that bound have: ",self.interval_bits)
                    print("Header size: ",header_bit_size/8,"bytes , ",header_bit_size,"bits")
                    print("Padding size: ",(8 - remaining_length)/8,"bytes , ",8 - remaining_length,"bits")
                    print("Encoding size: ",encoding_bit_size/8,"bytes , ",encoding_bit_size,"bits")
                    print("Expected coding length: ",expected_coding_length/8,"bytes , ",expected_coding_length,"bits")
                    print("Total size: ",output_size,"bytes , ",output_size*8,"bits")
                    print("Lower bound: {:064b}".format(self.low))
                    print("Upper bound: {:064b}".format(self.high))
    
    def encoding_table(self, symbol, table, history):
        coding = ""
        for order in reversed(range(len(history) + 1)):
            now_context = table.root_context
            for cond_symbol in history[ : order]:
                now_context = now_context.childern[cond_symbol]
                if now_context is None:
                    break
            else:
                if symbol != self.escape_symbol and now_context.get_count(symbol) > 0:
                    coding += self.encoding_context(symbol, now_context)
                    return coding
                
                coding += self.encoding_context(self.escape_symbol, now_context)
        coding += self.encoding_context(symbol, table.order_minus1_context)
        return coding
                
    def encoding_context(self, symbol, context):
        #update range
        cum_low = context.get_low(symbol)
        cum_high = context.get_high(symbol)
        range = self.high - self.low + 1
        new_low = self.low + cum_low * range // context.get_cum()
        new_high = self.low + cum_high * range // context.get_cum() - 1
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
    
class Table():
    def __init__(self, order, symbol_num):
        #parameter
        self.order = order
        self.context_size = symbol_num + 1
        self.escape_symbol = symbol_num

        #context_node
        if order >= 0:
            self.root_context = Context_node(self.context_size, order >= 1)
            self.root_context.increment(self.escape_symbol)
        else:
            self.root_context = None
        #-1 order context
        self.order_minus1_context = Context_node(self.context_size, 0)
        self.order_minus1_context.set_order_minus1_context()
        
    def increment_table(self, condtion, symbol):
        if self.order == -1:
            return

        now_context = self.root_context
        now_context.increment(symbol)
        for (i, cond_symbol) in enumerate(condtion):
            children = now_context.childern
            
            if children[cond_symbol] is None:
                children[cond_symbol] = Context_node(self.context_size, (i+1) < self.order)
                children[cond_symbol].increment(self.escape_symbol)
            now_context = children[cond_symbol]
            now_context.increment(symbol)
                
class Context_node():
    def __init__(self , context_size, is_not_leaf):
        #Frequency table
        self.frequencies = {}
        self.cum_frequencies = {}
        self.cum = 0
        self.context_size = context_size
        
        self.cum_frequencies[-1] = 0
        for i in range(context_size):
            self.frequencies[i] = 0
            self.cum_frequencies[i] = 0
        
        #Children
        self.childern = ([None] * context_size) if is_not_leaf else None
        
    def increment(self, symbol):
        self.frequencies[symbol] += 1
        for i in range(symbol,self.context_size):
            self.cum_frequencies[i] += 1
        self.cum += 1
    
    def get_cum(self):
        return self.cum
    
    def get_high(self, symbol):
        return self.cum_frequencies[symbol]
    
    def get_low(self, symbol):
        return self.cum_frequencies[symbol - 1]
    
    def get_count(self, symbol):
        return self.frequencies[symbol]
    
    def set_order_minus1_context(self):
        for i in range(self.context_size):
            self.frequencies[i] = 1
            self.cum_frequencies[i] = i+1
        
        self.cum = self.context_size
             
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('infilename' , help='the file we want to deal with')
    parser.add_argument('outfilename' , help='the output file name')
    parser.add_argument('--algorithm', default='256_AC', help='Choose the first compress algorithm.')
    parser.add_argument('--order' , default=2, type=int, help='The order of PPM, the value should be 0~2.')
    parser.add_argument('--print_info' , default=True, help='Print the compress information.')
    args = parser.parse_args()
    return args

def choose_algorithm(algorithm, order):
    if algorithm == '256_AC':
        return ArithmeticCoding(bits = 8)
    
    if algorithm == 'binary_AC':
        return ArithmeticCoding(bits = 1)
    
    if algorithm == '256_PPM':
        if order < 0:
            raise NameError("The order should be an integer and >= 0 .")
        return PPM(bits = 8, order = order)
    
    if algorithm == 'binary_PPM':
        if order < 0:
            raise NameError("The order should be an integer and >= 0 .")
        return PPM(bits = 1, order = order)    
    
    else:
        raise NameError("The algorithm should be 256_AC, binary_AC, 256_PPM, binary_PPM.")
    
def main():
    #Parse args
    args = parse_args()
    #Operate algorithm
    algorithm = choose_algorithm(args.algorithm, args.order)
    print("Compress ", args.infilename)
    #Compress
    algorithm.compress(args.infilename, args.outfilename, args.print_info)

    
if __name__ == "__main__":
	main()
