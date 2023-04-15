import six
import argparse
import os

class BasicNode:
    def isleaf(self):
        raise NotImplementedError("")


class Leaf(BasicNode):
    def __init__(self, symbol=0, frequency=0):
        super(Leaf, self).__init__()
        self.symbol = symbol
        self.weight = frequency

    def isleaf(self):
        return True
    


class NotLeaf(BasicNode):
    def __init__(self, left_child=None, right_child=None):
        super(NotLeaf, self).__init__()
        self.weight = left_child.weight + right_child.weight
        self.left_child = left_child
        self.right_child = right_child

    def isleaf(self):
        return False

class Tree(object):
    def __init__(self, leaf, symbol =0, freq=0, left_tree=None, right_tree=None):
        super(Tree, self).__init__()
        if leaf == 1:
            self.root = Leaf(symbol, freq)
        else:
            self.root = NotLeaf(left_tree.root, right_tree.root)

    def traverse_tree(self, root, code, symbol_frequency):
        if root.isleaf():
            symbol_frequency[root.symbol] = code
            print("it = ",root.symbol,"  and  freq = ",root.weight,"  code = ", code)
        else:
            self.traverse_tree(root.left_child, code+'0', symbol_frequency)
            self.traverse_tree(root.right_child, code+'1', symbol_frequency)

class HuffmanCoding:
    def buildHuffmanTree(self,list_hufftrees):
        while len(list_hufftrees) >1 :
            list_hufftrees.sort(key=lambda x: x.root.weight) 
            
            node1 = list_hufftrees[0]
            node2 = list_hufftrees[1]
            list_hufftrees = list_hufftrees[2:]

            newed_hufftree = Tree(0, 0, 0, node1, node2)

            list_hufftrees.append(newed_hufftree)
        return list_hufftrees[0]
    
    def compress(self,infile_name, outfile_name, symbol_length):
        #get the information of input file
        f = open(infile_name,'rb')
        infile_data = f.read()
        infile_size = f.tell()
        
        #statistic the symbol
        symbol_frequency = {}
        infile_size_wo_redundant = infile_size//symbol_length * symbol_length
        redundent_symbol = ''
        
        ##calculate most of the symbol except the last symbol
        for i in range(0, infile_size_wo_redundant, symbol_length):
            symbol = infile_data[i:i+symbol_length]
            if symbol in symbol_frequency:
                symbol_frequency[symbol] += 1
            else:
                symbol_frequency[symbol] = 1
                
        ##read redundant symbol
        if(infile_size != infile_size_wo_redundant):
            redundent_symbol = infile_data[infile_size_wo_redundant:]
            #print(redundent_symbol)
            #print(len(redundent_symbol))
            
        
        ##output statistic
        #for symbol in symbol_frequency:
        #    print(symbol,' : ',symbol_frequency[symbol])
        print(len(symbol_frequency))
            
        #store the num wo the redundant
        length = len(symbol_frequency)
        output = open(outfile_name, 'wb')
        #int = 4byte
        o4 = length&255
        length = length>>8
        o3 = length&255
        length = length>>8
        o2 = length&255
        length = length>>8
        o1 = length&255
        output.write(six.int2byte(o1))
        output.write(six.int2byte(o2))
        output.write(six.int2byte(o3))
        output.write(six.int2byte(o4))
        
        #store the redundant length
        length = len(redundent_symbol)
        #int = 4byte
        o4 = length&255
        length = length>>8
        o3 = length&255
        length = length>>8
        o2 = length&255
        length = length>>8
        o1 = length&255
        output.write(six.int2byte(o1))
        output.write(six.int2byte(o2))
        output.write(six.int2byte(o3))
        output.write(six.int2byte(o4))
        
        #store the redundant byte
        for res in redundent_symbol:
            #print(res)
            output.write(six.int2byte(res))
        
        #store the table
        for symbol in symbol_frequency:
            for sb in symbol:
                output.write(six.int2byte(sb))
            freq = symbol_frequency[symbol]
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
            
        #create the huffman node list
        node_list = []
        for symbol in symbol_frequency:
            node = Tree(1, symbol, symbol_frequency[symbol])
            node_list.append(node)
        #establish huffman tree
        tree = self.buildHuffmanTree(node_list)
        #get the coding
        tree.traverse_tree(tree.root,'',symbol_frequency)
        #compress
        coding=''
        for i in range(0, infile_size_wo_redundant, symbol_length):
            symbol = infile_data[i:i+symbol_length]
            coding += symbol_frequency[symbol]
            out_int = 0
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
        #print(six.int2byte(remaining_length))
        #close the output
        output.close()
        
    def decompress(self,infile_name, outfile_name, symbol_length):
        f = open(infile_name,'rb')
        infile_data = f.read()
        infile_size = f.tell()
        
        #Get the symbol number
        symbol_frequency = {}
        symbol_number = 0
        for i in range(4):
            symbol_number = symbol_number<<8
            a = infile_data[i]
            symbol_number = symbol_number|a
        #print(symbol_number)
        
        #Get the redundant symbol and length
        redundant_len = 0
        for i in range(4,8):
            redundant_len = redundant_len<<8
            a = infile_data[i]
            redundant_len = redundant_len|a
        redundant_symbol = infile_data[8:8+redundant_len]
        #print(redundant_symbol)
        #print(len(redundant_symbol))
        
        #get the remaining code length
        last_byte = infile_data[infile_size-1:infile_size]
        remaining_length = six.byte2int(last_byte)
        #print(remaining_length)
        
        #read the dictionary
        for i in range(symbol_number):
            symbol = infile_data[8+redundant_len+i*(4+symbol_length):8+redundant_len+i*(4+symbol_length)+symbol_length]
            freq = 0
            for j in range(4):   
                freq = freq<<8        
                a = infile_data[8+redundant_len+i*(4+symbol_length)+symbol_length+j]
                freq = freq|a
            #print(symbol)
            #print(freq)
            symbol_frequency[symbol] = freq

        #create the huffman node list
        node_list = []
        for symbol in symbol_frequency:
            node = Tree(1, symbol, symbol_frequency[symbol])
            node_list.append(node)
        #establish huffman tree
        tree = self.buildHuffmanTree(node_list)
        #get the coding
        tree.traverse_tree(tree.root,'',symbol_frequency)
        
        #decoding
        output = open(outfile_name, 'wb')
        now_node = tree.root
        coding = ''
        for cd in range((8+redundant_len + symbol_number*(4+symbol_length)),(infile_size-1)):
            byte_code = infile_data[cd]
            for i in range(8):
                if (byte_code&128):
                    coding += '1'
                else:
                    coding += '0'
                byte_code = byte_code << 1

            while len(coding) > 16:
                if now_node.isleaf():
                    output.write(now_node.symbol)
                    #print(now_node.symbol)
                    now_node = tree.root

                if coding[0] == '1':
                    now_node = now_node.right_child
                else:
                    now_node = now_node.left_child
                coding = coding[1:]
        
        #deal with the last 16 bit
        coding = coding[:-8 + remaining_length]
        #print(coding)
        while len(coding) > 0:
            if now_node.isleaf():
                output.write(now_node.symbol)
                #print(now_node.symbol)
                now_node = tree.root
            if coding[0] == '1':
                now_node = now_node.right_child
            else:
                now_node = now_node.left_child
            coding = coding[1:]
        if now_node.isleaf():
            output.write(now_node.symbol)
            now_node = tree.root

        
        output.write(redundant_symbol)
        output.close()

class AdaptiveHuffmanCoding:
    def compress(self,infile_name, outfile_name, symbol_length):
        pass
    def decompress(self,infile_name, outfile_name, symbol_length):
        pass
            
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode' , help='compress or decompress')
    parser.add_argument('infilename' , help='the file we want to deal with')
    parser.add_argument('--algorithm', default='huffman_coding', help='choose algorithm')
    parser.add_argument('--outfilename' , help='the output file name')
    parser.add_argument('--symbol_length', default=1, type=int, help='determine the length of source coding every time')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    #choose algorithm    
    if(args.algorithm == 'huffman_coding'):
        algorithm = HuffmanCoding()
    elif(args.algorithm == 'adaptive_huffman_coding'):
        algorithm = AdaptiveHuffmanCoding()
    else:
        raise NameError("The algorithm should be huffman_coding or adaptive_huffman_coding.")
    
    #compress or decompress
    if(args.mode=='compress'):
        #setting output filename
        if(args.outfilename==None):
            #add new file extension
            args.outfilename = args.infilename
        if(args.algorithm == 'huffman_coding'):
            args.outfilename += '.'  + 'hcps' + str(args.symbol_length)
        if(args.algorithm == 'adaptive_huffman_coding'):            
            args.outfilename += '.'  + 'ahcps' + str(args.symbol_length)
            
        print('Compress '+args.infilename)
        print("outputfile: ", args.outfilename)
        algorithm.compress(args.infilename,args.outfilename, args.symbol_length)
    
    elif(args.mode=='decompress'):
        if(args.outfilename==None):
            #remove new file extension
            args.outfilename = os.path.splitext(args.infilename)[0]
        file_extension = os.path.splitext(args.infilename)[1]
        #renew source_coding_length
        if(file_extension[1] == 'a'):
            args.symbol_length = int(file_extension[6:])
        else:
            args.symbol_length = int(file_extension[5:])     
        print(args.symbol_length)
        print('decompress file '+args.infilename)
        print("outputfile: ", args.outfilename)
        algorithm.decompress(args.infilename,args.outfilename, args.symbol_length)
    else:
        raise NameError("The mode should be compress or decompress.")
        
        

if __name__ == '__main__':
    main()
