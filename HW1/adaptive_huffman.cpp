#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cstdint>
#include <vector>
#include <algorithm>

constexpr size_t leaf_number = 256;
constexpr size_t node_number = 512;
constexpr size_t buffer_size = 1024;
constexpr size_t symbol_bit = 8;
int symbol_num = 0;

struct Node {
    Node* parent;
	Node* left;
	Node* right;
    int      symbol;
	int      order;
	uint64_t weight;

	Node(int Symbol, int Order, uint64_t Weight = 0, Node* Parent = nullptr, Node* Left = nullptr, Node* Right = nullptr)
		: symbol(Symbol), order(Order), weight(Weight), parent(Parent), left(Left), right(Right){ 
        }
};

class AdaptiveHuffmanCoding{
    private:
		Node*    node_nyt;
		Node*    node_root;
		Node*    array_leaves[leaf_number];
		Node*    array_nodes[node_number + 1];
		Node*    decoder_node_curr;
		std::vector<bool> bit_buffer;

		void Swap(Node* a, Node* b) {
			if (!((a == b) || (a == node_root) || (a->parent == b) || (b == node_root) || (b->parent == a))){
			    if (a->parent->left == a && b->parent->left == b)
			    	std::swap(a->parent->left, b->parent->left);
			    else if (a->parent->right == a && b->parent->left == b)
			    	std::swap(a->parent->right, b->parent->left);
			    else if (a->parent->left == a && b->parent->right == b)
			    	std::swap(a->parent->left, b->parent->right);
			    else if (a->parent->right == a && b->parent->right == b)
			    	std::swap(a->parent->right, b->parent->right);
			    std::swap(a->parent, b->parent);
			    std::swap(a->order, b->order);
			    std::swap(array_nodes[a->order], array_nodes[b->order]);
            }
		}

		void Delete(const Node* node_root) {
			if (node_root != nullptr){
			    Delete(node_root->left);
			    Delete(node_root->right);
			    delete node_root;
            }
		}

		void Update(Node* node) {
			while (node){
			    Node* swap_node = node;
			    int num = node->order;
			    for (size_t i = num + 1; i <= node_number ; i++){
                    if(array_nodes[i]->weight != node->weight)
                        break;
				    swap_node = array_nodes[i];
				}
                Swap(swap_node, node);
				node->weight += 1;
				node = node->parent;
			}
		}

		std::vector<bool> GetLeafCode(const Node* leaf) const {
			std::vector<bool> leaf_code;
			leaf_code.reserve(16);
			const Node* node_curr = leaf;
			while (node_curr != node_root) {
                leaf_code.push_back(node_curr->parent->right == node_curr);
				node_curr = node_curr->parent;
			}
			std::reverse(leaf_code.begin(), leaf_code.end());
			return leaf_code;
		}

		std::vector<bool> Encoder(uint8_t byte) {
			if (array_leaves[byte]) {
				std::vector<bool> encoding_code = GetLeafCode(array_leaves[byte]);

			    Update(array_leaves[byte]);

				return encoding_code;
			}
			else {
				symbol_num ++ ;
				std::vector<bool> encoding_code = GetLeafCode(node_nyt);
				
				for (size_t i = 0; i < symbol_bit; i++)
					encoding_code.push_back(byte & (1 << (7 - i)));

			    node_nyt->symbol = -2;
			    node_nyt->right  = new Node(byte, node_nyt->order - 1, 0, node_nyt);
			    array_leaves[byte] = node_nyt->right;
			    array_nodes[node_nyt->right->order] = node_nyt->right;

			    node_nyt->left   = new Node(-1, node_nyt->order - 2, 0, node_nyt);
			    array_nodes[node_nyt->left->order] = node_nyt->left;
			    node_nyt = node_nyt->left;

			    Update(array_leaves[byte]);

				return encoding_code;
			}
		}

		std::vector<uint8_t> Decoder(const std::vector<bool>& bit_temp) {
			std::vector<uint8_t> symbol_temp;

			std::copy(std::begin(bit_temp), std::end(bit_temp), std::back_inserter(bit_buffer));

			size_t i = 0;
			while (i < bit_buffer.size() || (!decoder_node_curr->left && !decoder_node_curr->right))
				if (decoder_node_curr->symbol == -1) {
					if (bit_buffer.size() - i >= symbol_bit) {
						uint8_t ascii_char = 0;

						for (size_t j = 0; j < symbol_bit; j++)
							ascii_char |= bit_buffer[i + j] << (7 - j);

						symbol_temp.push_back(ascii_char);

			            node_nyt->symbol = -2;
			            node_nyt->right  = new Node(ascii_char, node_nyt->order - 1, 0, node_nyt);
			            array_leaves[ascii_char] = node_nyt->right;
			            array_nodes[node_nyt->right->order] = node_nyt->right;
			            
                        node_nyt->left   = new Node(-1, node_nyt->order - 2, 0, node_nyt);
			            array_nodes[node_nyt->left->order] = node_nyt->left;
			            node_nyt = node_nyt->left;

                        Update(array_leaves[ascii_char]);
						
                        decoder_node_curr = node_root;
						i += symbol_bit;
					}
					else 
                        break; 
				}
				else if ((decoder_node_curr->left == nullptr) && (decoder_node_curr->right == nullptr)){
					symbol_temp.push_back(decoder_node_curr->symbol);

			        Update(array_leaves[decoder_node_curr->symbol]);

					decoder_node_curr = node_root;
				}
				else {
					if (bit_buffer[i])
                        decoder_node_curr = decoder_node_curr->right;
					else
                        decoder_node_curr = decoder_node_curr->left;
					i++;
				}

			bit_buffer.erase(std::begin(bit_buffer), std::begin(bit_buffer) + i);

			return symbol_temp;
		}

		void WriteOutputfile(std::ostream& outputfile, std::vector<bool>& outputbuffer) {
			size_t pos = 0;
			std::vector<uint8_t> out_bytes;
			out_bytes.reserve((buffer_size * symbol_bit) / 2);

			while (pos + symbol_bit <= outputbuffer.size()){
				uint8_t byte = 0;
				for (size_t i = 0; i < symbol_bit; i++)
					if (outputbuffer[i + pos])
						byte |= 1 << (7 - i);
				pos += symbol_bit;
				out_bytes.push_back(byte);
			}
			outputfile.write(reinterpret_cast<const char*>(&out_bytes.front()), out_bytes.size());
			outputbuffer.erase(std::begin(outputbuffer), std::begin(outputbuffer) + pos);

		}
		void statistic(size_t extra_bits,long long inputfile_size,long long outputfile_size){
			long long encoding_bit_size = outputfile_size*8-extra_bits-symbol_num*8;
			float expected_coding_length = float(encoding_bit_size)/float(inputfile_size);
			std::cout<<"The detail of compressed file: "<<std::endl;
			std::cout<<"Header size: "<<symbol_num<<"bytes, "<<symbol_num*8<<"bits"<<std::endl;	
			std::cout<<"Encoding size: "<<double(encoding_bit_size)/8<<"bytes, "<<encoding_bit_size<<"bits"<<std::endl;		
			std::cout<<"Padding size: "<<double(extra_bits)/8<<"bytes, "<<extra_bits<<"bits"<<std::endl;
			std::cout<<"Total size: "<<outputfile_size<<"bytes, "<<outputfile_size*8<<"bits"<<std::endl;
			std::cout<<"Expected coding length: "<<expected_coding_length/8<<"bytes, "<<expected_coding_length<<"bits"<<std::endl;	
			std::cout<<"Compression ratio: "<<float(outputfile_size)*100/float(inputfile_size)<<"%"<<std::endl;		
		}
    public:
		AdaptiveHuffmanCoding(){
        	node_nyt = new Node(-1, node_number);
            node_root = node_nyt;

            for (auto&& node: array_nodes)
				node = nullptr;
			array_nodes[node_number] = node_nyt;

			for (auto&& leaf : array_leaves)
				leaf = nullptr;

            decoder_node_curr = node_root;
		}

		~AdaptiveHuffmanCoding(){ 
            Delete(node_root); 
        }

        void Compress(std::ifstream& inputfile, std::ofstream& outputfile){
            uint8_t inputbuffer[buffer_size];
			inputfile.seekg(0,inputfile.end);
			long long inputfile_size = inputfile.tellg();
			inputfile.seekg(0,inputfile.beg);			
            std::vector<bool> outputbuffer;
			outputbuffer.reserve(buffer_size * symbol_bit + 64);
			size_t extra_bits;
            while(!inputfile.eof()){
				inputfile.read(reinterpret_cast<char*>(inputbuffer), buffer_size);
				size_t read_bytes_num = inputfile.gcount();
				if (outputbuffer.size() >= buffer_size * symbol_bit){
					WriteOutputfile(outputfile, outputbuffer);
                }
				for (size_t i = 0; i < read_bytes_num ; i++) {
					std::vector<bool> encoding_code = Encoder(inputbuffer[i]);
					std::copy(std::begin(encoding_code), std::end(encoding_code), std::back_inserter(outputbuffer));
				}
            }

			if (!outputbuffer.empty()) {
				if (outputbuffer.size() % symbol_bit != 0) {
					std::vector<bool> nyt_code = GetLeafCode(node_nyt);
					extra_bits = symbol_bit - outputbuffer.size() % symbol_bit;
					for (size_t i = 0; i < extra_bits; i++)
						outputbuffer.push_back(nyt_code[i % nyt_code.size()]);
				}
				WriteOutputfile(outputfile, outputbuffer);
			}
			outputfile.seekp(0,outputfile.end);
			long long outputfile_size = outputfile.tellp();			
			statistic(extra_bits,inputfile_size,outputfile_size);
        }

        void Decompress(std::ifstream& inputfile, std::ofstream& outputfile){
			uint8_t inputbuffer[buffer_size];

			while (!inputfile.eof()) {
				inputfile.read(reinterpret_cast<char*>(inputbuffer), buffer_size);
				std::vector<bool> bit_temp;
                bit_temp.reserve(buffer_size * symbol_bit);
				size_t read_bytes_num = inputfile.gcount();

				for (size_t i = 0; i < read_bytes_num; i++) {
					std::vector<bool> encoding_code;
                    encoding_code.resize(symbol_bit);
					for (size_t j = 0; j < symbol_bit; j++){
                        encoding_code[j] = inputbuffer[i] & (1 << (7 - j));
                    }
					std::copy(std::begin(encoding_code), std::end(encoding_code), std::back_inserter(bit_temp));
				}

				std::vector<uint8_t> outputbuffer = Decoder(bit_temp);
				outputfile.write(reinterpret_cast<const char*>(&outputbuffer.front()), outputbuffer.size());
			}        
        }

};

int main(int argc, char* argv[]) {
    //check the number of argv
    if (argc<3){
        std::cout<<"The argv should be at least 3 argument, including mode and infilename."<<std::endl;        
        return -1;
    }

    //read the mode is compress and decompress
    std::string mode = argv[1];
    if((mode != "compress")&&(mode != "decompress")){
        std::cout<<"Mode should be compress or decompress"<<std::endl;
        return -1;
    }

    //read the input file name
    std::string infilename = argv[2];

    //read the output file name
    std::string outfilename;
    if(argc>3){
        outfilename = argv[3];
        if(mode=="compress"){
            outfilename += ".ahcps";
        }
    }
    else{
        outfilename = infilename;
        if(mode=="compress"){
            outfilename += ".ahcps";
        }
        else{
            outfilename = outfilename.substr(0,outfilename.length()-6);
        }     
    }

    //std::cout<<mode<<std::endl;
    //std::cout<<infilename<<std::endl;
    //std::cout<<outfilename<<std::endl;

	std::ifstream inputfile;
	inputfile.open(infilename.c_str(), std::ios::binary);
    if (!inputfile.is_open()) {
        std::cerr<<"Fail to read file"<<std::endl;
        return -1;
	}

	std::ofstream outputfile;
	outputfile.open(outfilename.c_str(), std::ios::binary);
    if (!outputfile.is_open()) {
        std::cerr<<"Fail to write file"<<std::endl;
        return -1;
	}

    AdaptiveHuffmanCoding algorithm;
    
    if(mode=="compress"){
		std::cout<<"Compress "<<infilename<<std::endl;
    	std::cout<<"Outputfile: "<<outfilename<<std::endl;
        algorithm.Compress(inputfile,outputfile);
    }

    if(mode=="decompress"){
		std::cout<<"Decompress "<<infilename<<std::endl;
    	std::cout<<"Outputfile: "<<outfilename<<std::endl;
        algorithm.Decompress(inputfile,outputfile);
    }
    //close the file
	inputfile.close();
	outputfile.close();    
    return 0;
}