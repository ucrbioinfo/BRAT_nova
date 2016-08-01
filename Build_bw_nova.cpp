//  BRAT is a whole genome bisulfite sequence mapping program
//  Copyright (C) 2009  Elena Yavorska Harris
//	Mapping is done to Concatenation of two genome versions:
//	First is Complement of genome ( reverse of reverse-complement of genome is complement of genome)
//	Second is Reverse of genome
//	Reverse is taken to be able to map 5' end of reads (since BWT alignment is done from right-to-left)
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
/*
Same seeds as in Brat-1.1.4, but less space since it treats chroms as one continuous genome
  */

#include<iomanip>
#include<map>
#include<bitset>
#include<queue>
#include<ctime>
#include<list>
#include<cmath>
#include<vector>
#include<string.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<assert.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

using namespace std;


typedef unsigned long int _uli;
typedef unsigned int _ui;
typedef unsigned short int _usi;
typedef vector<_uli>::iterator _uli_iter;



_ui min(_ui x, _ui y){
	if(x < y)
		return x;
	else
		return y;
}//min

void usage(){

	string mes("\nUSAGE: build_bw -r <references> -P <prefix> \nOptions:\n");
	string gen_ref_opt("  -r <references-file>       file with file names containing references (one\n                          reference per file)");
	string prefix(     "  -P <prefix>                name of absolute path name for output files");
	cout << mes << endl;
	cout << gen_ref_opt << endl << prefix << endl;
 
	string mes2(      "Limitations:\n");
	string ref_size(  "  a single reference size <= 2^45");
	string chrom_ids( "  number of references <= 2^18");
	string space(     "  memory <= 2GB; space on hard disk <= 3GB");
	
	cout << mes2 << endl;
	cout << ref_size << endl << chrom_ids << endl << space << endl;

}//usage()

class Anchor{

public:
	Anchor(int argc, char* argv[]);
	~Anchor();
//TOP
	vector<string> names;//for references

	vector<unsigned long int> mask;	

	vector<_uli> C;//#occ of As and Ts in ta-genome

	unsigned long int* ht_indices;//ht_idices
	
	vector<_uli> ht_sizes;

	vector< pair<_ui, _ui> > ranges;
	vector<_ui> ranges_total;
	vector< pair<_ui, _ui> > sub_bucket;//ranges<st, en> for sub-buckets
	vector<_ui> sub_bucket_sum;//sum of indices in each sub-bucket range

	vector<_uli> bucket_components;
	vector<_uli> bc_indices;//pointers to indices of bucket_components
	vector<_uli> end_of_suffix;//if > 0, then this component contains suffix after *value* bits, the less bits, the earlier in sorting order
	vector<_ui> prev_same_comp_ids;

	vector<_uli> bwt;//cat-representation, C->T

	vector<_ui> bwt_marked;//mark those 1-bit bwt positions per char whose pos indices are stored
	vector<_ui> pos_strand;//nova 1bit per base: 0 for 1st half of genome (- strand), 1 for 2nd half of genome with Reverse(original strand) (+ strand)
	//pos_strand should be as bwt_marked
	vector<_ui> marked_ns;

	vector<_ui> num_occ_sqr;//to keep # occ of As in a current bucket of size 4096 bits, cumulative sum from 0 to current - 1 block
	//need to make sure that neither of C[i] for any character is greater than 2^32 (otherwise, num_occ_sqr must be _uli)
	vector< vector<_usi> > num_occ; // # occ of As in 64 sub-buckets of current bucket
	vector<_uli> num_occ_sqr_sum;//  num_occ_sqr_sum[k] = sum of entries of num_occ_sqr from 0 to k-1


	vector< _ui > pos_indices;//key = pos inside a 2^16-size bucket, value = genome pos; 
	vector<_usi> pos_num_occ;//as num_occ for char count in bwt, this one is for count of 1s in bwt_marked
	_ui pos_num_sqr;//as num_occ_sqr for char count in bwt, this one is for count of 1s in bwt_marked
	//pos_num_sqr will not work for STEP=1, for STEP = 1, must be _uli
	_uli pos_num_occ_sqr_sum;//similar to num_occ_sqr_sum
	_ui pos_indices_ind;//index to reference pos_indices from 0 to 65,536 (ok _ui)

	vector<_ui> size_chrom;//contains size of each chromosome, ok _ui
	vector<_uli> starts_size_chrom;
	map<_ui, string> chroms_names;//don't need map: can do in vector, _ui is ok since is in range 0...num_chroms
	vector<_uli> chroms_starts;//for ta/cg/ns-genomes only
	vector<_uli> starts_chroms_starts;

	//auxiliary to prevent huge stacking during bucket sort
	queue< pair< pair<_uli_iter, _uli_iter>, pair<_ui, _ui> > > not_sorted;

	unsigned long int* ta_genome;
	unsigned long int* cg_genome;
	unsigned long int* ns_genome;
	unsigned long int* full_genome;
	unsigned long int* bwt_genome;
	

	int parse_options(int argc, char* argv[]);

	void read_genome();

	void init_num_occ();//add here pos_num_occ too
	void init_pos_indices();//pos_marked
	void init_marked();//for bwt_marked, marked_ns, pos_strand

	_uli get_ta_lmer(_uli cur_ind, _uli asize);
	_uli get_ns_lmer(_uli cur_ind, _uli asize);
	void make_mask();
	void collect_ranges(_ui st, _ui en, _ui two_in_bucket, vector< pair<_ui, _ui> > &buckets, vector<_ui> &bucket_sum, vector<_uli> &ht_sizes_dif);
	void collect_bucket_indices(_ui st, _ui en, _ui range_ind);
	void init_comp_indices(_ui st, _ui en, _ui sub_bucket_size);
	void update_component(_uli_iter st, _uli_iter en, _ui &component, const _ui &sub_bucket_st);
	void update_component_first(_uli_iter st, _uli_iter en, _ui component, const _ui &sub_bucket_st);
	_uli get_component(_uli cur_ind, _uli asize, _uli comp_id, _uli &bits_before_suffix);
	void get_pivot_value(_uli_iter st, _uli_iter en, _uli &value, _ui &total_left);
	void quick_multikeys(_ui component, _uli_iter start, _uli_iter fin, const _ui &sub_bucket_st, _ui total_comp);
	void quick_multikeys(vector< _uli > &lengths, _uli_iter start, _uli_iter fin);
	_uli_iter partition_multikeys(vector< _uli > &lengths, _uli_iter st, _uli_iter en, _ui &total_left);
	void sort_small_range(_uli_iter st, _uli_iter en);
	_uli_iter partition_multikeys(_uli_iter st, _uli_iter en, _ui &total_left);
	void update_end_of_suffix(_uli_iter &st, _uli_iter en);
	_uli get_char( _uli index );
	void update_bwt(_ui sub_bucket_st);
	void finish_update_bwt(_ui sub_bucket_st);
	void collect_subbuckets(vector<_ui> &bucket_sum);
	_uli_iter find_first_equal_range(_uli_iter st, _uli_iter en, _ui &total_left);

	void flash_num_occ(_uli upto);
	void flash_bwt();
	void flash_pos_indices(_uli upto);
	void flash_marked_ns();
	void flash_marked();
	void flash_genome(ofstream &out, _uli *genome, _uli asize);
	void output_info();
	void usi_to_char(_usi, char abuffer[]);
	void uli_to_char(_uli, char abuffer[]);
	void ui_to_char(_ui, char abuffer[]);

	_uli Occ(_uli ch, _uli sp, vector< vector<_ui> > &num_occ_sqr);	
	void exact_match(_uli ta_read, _uli &sp, _uli &ep, vector< vector<_ui> > &num_occ_sqr);//sp=1, ep = 0, remainer = 0: init
	_uli get_bwt_lmer(_uli sp, _uli asize, _uli ch);
	void read_num_occ_sqr(vector< vector<_ui> > &num_occ_sqr);
	void read_num_occ(vector< vector<_usi> > &num_occ);
	void read_bwt_genome();
	void precalculate_seeds();
	_uli count_mism(unsigned long int bits);

	ofstream out_bwt;
	ofstream out_num_occ;
	ofstream out_num_sqr;	
	ofstream out_ta_genome;
	ofstream out_cg_genome;
	ofstream out_ns_genome;
	ofstream out_pos_strand;//nova
	ofstream out_pos_indices;
	ofstream out_bwt_marked;
	ofstream out_pos_num_sqr;
	ofstream out_pos_num_occ;

	ofstream out_chrom;

	ofstream out_seeds;
	ifstream in_bwt;

	void build_bwt();

	unsigned long int len_mer;
	unsigned long int read_len;
	_uli orig_num_sqr;
	_uli orig_num;

	string genome_names;
	string prefix;

	//all of the below refers to size of concatenation of genome to itself
	_uli binary_total;//size of binary representation of the real genome 
	_uli ta_size;//size of ta_genome
	_uli full_size;
	_uli orig_genome_size;//real size of original genome in bases

	unsigned long int byte_size;
	unsigned long int SIZE ;//genome char per an unsigned long int
	unsigned long int two_in_Size; //log_2(SIZE)
	unsigned long int MASK_SIZE;
	_uli MASK_LMER;
	_uli MASK_BYTE;
	_uli block_64;//64*64=4096
	_uli block_32;//32*32 = 1024
	_uli block_256; //256*256
	_uli two_in_block_64;//2^12=4096
	_uli two_in_block_32;//2^10
	_uli two_in_block_256;//2^16 = 256^2
	_uli two_in_32;
	_uli curL_bwt;//current index of BWT full size
	_uli cur_ind_bwt;

	_uli SEED;//24
	_uli MASK_SEED;
	_uli STEP;//16
	_uli two_in_Step;
	_uli orig_gen_inL;//index of original genome in bwt
	_ui abyte;
	_ui two_bytes;
	_ui three_bytes;
	_ui four_bytes;
	_ui five_bytes;
	_ui six_bytes;
	_ui seven_bytes;
	_ui two_in_byte;
	_ui num_chroms;
	_uli hash_table_size;
	_uli total_components;
	_uli genome_byte_size;
	void print();
	vector< vector<_usi> > num_occ_bytes;//hash table for the number of 1s in a byte

	//first entry of hash table might be hash table huge, sort separately
	vector< _uli > ends_ranges_of_As;
	vector< _uli > lengths;//similar to bucket components, holds lengths of regions of As at the beg of suffixes
	vector< _uli > ht_sizes_lengths;//holds number of suffixes with As by lengths subdivided into 2^16 entries
	void sort_first_entry();
	void collect_bucket_indices_first(_ui st, _ui en, _ui asize, vector<_uli> &count_lengths);
	void collect_indices_first_eos(_uli &count_eos);//end of suffix char
	_uli small_ht_size;
	_uli shift_ht_lengths;


	_uli two_in_half_word;//half_word = 32, two in it = 5
	_uli full_binary_total;

	void count_char(_uli kmer, vector<_usi> & characters, int asize);
	_uli pos_indices_total;

	void space_release();
	_uli half_Size;
	_uli half_SEED;

	_uli max_range_As;

};	


Anchor::~Anchor(){


	delete [] bwt_genome;

	out_bwt.close();
	out_num_occ.close();
	out_num_sqr.close();

	out_pos_strand.close();
	out_ta_genome.close();
	out_cg_genome.close();
	out_ns_genome.close();
	out_pos_indices.close();
	out_chrom.close();
	out_pos_num_occ.close();
	out_pos_num_sqr.close();
	out_bwt_marked.close();

	if(ht_indices != 0)
		delete [] ht_indices;
	if(ns_genome != 0)
		delete [] ns_genome;
	if(ta_genome != 0)
		delete [] ta_genome;
	if(cg_genome != 0)
		delete [] cg_genome;
	if(full_genome != 0)
		delete [] full_genome;

}
void Anchor::space_release()
{

	out_bwt.close();
	out_num_occ.close();
	out_num_sqr.close();

	out_pos_strand.close();
	out_ta_genome.close();
	out_cg_genome.close();
	out_ns_genome.close();
	out_pos_indices.close();
	out_chrom.close();
	out_pos_num_occ.close();
	out_pos_num_sqr.close();
	out_bwt_marked.close();

	if(ht_indices != 0)
		delete [] ht_indices;
	ht_indices = 0;
	if(ns_genome != 0)
		delete [] ns_genome;
	ns_genome = 0;
	if(ta_genome != 0)
		delete [] ta_genome;
	ta_genome = 0;
	if(cg_genome != 0)
		delete [] cg_genome;
	cg_genome = 0;
	if(full_genome != 0)
		delete [] full_genome;
	full_genome = 0;
	if(bwt_genome != 0)
		delete [] bwt_genome;
	bwt_genome = 0;
}//space_release

Anchor::Anchor(int argc, char* argv[])
{
	half_SEED = 12;
	pos_indices_total = 0;
	STEP = 16;
	SIZE = 64;//genome char per an unsigned long int
	half_Size = SIZE >> 1;

	make_mask();

	int res = parse_options(argc, argv);

	if((STEP > 1) && ((STEP & mask[1]) != 0))
		STEP &= ~(mask[1]);
	if(STEP > 32)
		STEP = 32;
	if(STEP == 0 || STEP == 1)
		STEP = 2;//does not support STEP = 1: some count will be off (number of occurrences or positions)
cout << "STEP = " << STEP << endl;
	if(res == -1){
		usage();
		space_release();
		exit(0);
	}

	if(prefix.empty()){
		usage();
		space_release();
		exit(0);
	}
	string genome_repres("CT");
	
	string rem = (prefix[prefix.length() - 1] == '/' ? "" : "/");
	string output_ta_genome = prefix + rem + "ta.ta";
	string output_cg_genome = prefix + rem +  "cg.cg";
	string output_pos_strand = prefix + rem + "pos_strand.txt";//nova 1bit per base: 0 for 1st half of genome (- strand), 1 for 2nd half of genome with Reverse(original strand) (+ strand)
	string output_chrom = prefix + rem + genome_repres + ".chr";

	string output_pos_indices = prefix + rem + genome_repres + ".pos_ind";
	string output_pos_num_sqr = prefix + rem + genome_repres + ".pos_num_sqr";
	string output_pos_num_occ = prefix + rem + genome_repres + ".pos_num_occ";
	string output_bwt_marked = prefix + rem + genome_repres + ".bwt_marked";//marks 1-bit bwt entries that have genome pos stored

	string output_ns_genome = prefix + rem + genome_repres + ".ns";//for each 1-bit bwt entry it marks 1 if 24-SEED prefix of this suffix contains N
	string output_num_occ = prefix + rem + genome_repres + ".no";//charactor count in each 64-bit word
	string output_num_sqr = prefix + rem + genome_repres + ".nosqr";//character count in square blocks: sum of counts in 0...i-1 (64*64)blocks
	string output_bwt = prefix + rem + genome_repres + ".bwt";

	string output_seeds = prefix + rem + genome_repres + ".seeds";

	out_bwt.open(output_bwt.c_str(), ios::out);
	out_num_occ.open(output_num_occ.c_str(), ios::out);
	out_num_sqr.open(output_num_sqr.c_str(), ios::out);
	out_seeds.open(output_seeds.c_str(), ios::out);
	out_ns_genome.open(output_ns_genome.c_str(), ios::out);
	out_pos_indices.open(output_pos_indices.c_str(), ios::out);

	out_pos_num_sqr.open(output_pos_num_sqr.c_str(), ios::out);
	out_pos_num_occ.open(output_pos_num_occ.c_str(), ios::out);
	out_bwt_marked.open(output_bwt_marked.c_str(), ios::out);
	out_chrom.open(output_chrom.c_str(), ios::out);	

	
		out_ta_genome.open(output_ta_genome.c_str(), ios::out);
		out_cg_genome.open(output_cg_genome.c_str(), ios::out);
		out_pos_strand.open(output_pos_strand.c_str(), ios::out);//nova
	

	total_components = 10000;//stack has limitation on sgi-2 it is 42000, total stack depth = total_comp * log(binary_total);


	unsigned long int one = 1;
	orig_gen_inL = 0;
	abyte = 8;


	two_in_Size = log(static_cast<float>(SIZE))/log(static_cast<float>(2));
	MASK_SIZE = (one << two_in_Size) - 1;
	MASK_BYTE = (one << abyte) - 1;
	two_in_32 = 5;


	two_bytes = 16;
	three_bytes = 24;
	four_bytes = 32;
	five_bytes = 40;
	six_bytes = 48;
	seven_bytes = 56;

	byte_size = 8;//byte size is characters(genome) per a byte
	two_in_byte = 3;//log_2(byte_size)
	two_in_half_word = 5;

	SEED = 24;
	len_mer = SEED;
	MASK_SEED = (one << 24) -1;
	MASK_LMER = MASK_SEED;
	block_64 = 4096;
	block_32 = 1024;
	two_in_block_64 = 12;
	two_in_block_32 = 10;	
	
	two_in_Step = log(static_cast<float>(STEP))/log(static_cast<float>(2));//2^4 = 16

	unsigned long int ta_ind = 0;

	float base = 2;
	float apow = static_cast<float>(SEED);

	hash_table_size = static_cast<_ui>(pow(base, apow));

	ifstream in;
	in.open(genome_names.c_str(), ios::in);
	if(!in){
		cout << "ERROR: cannot open " << genome_names << endl;
		space_release();
		exit(0);
	}

	string aname;
	num_chroms = 0;
	in >> aname;
	while(!in.eof()){
	
		if(aname.length() == 0){
			cout << "\nERROR: No input in the file " << genome_names << endl;
			space_release();
			exit(0);
		}
		if(aname[0] == '>' || aname[0] == '@'){
			cout << "\nERROR: " << genome_names << " must contain the names of the files with references." << endl;
			space_release();
			exit(0);
		}

		if(aname[aname.size() - 1] == '\n' || aname[aname.size() - 1] == '\r')//nova AWS treats end of line char as a valid character
			aname.erase(aname.size() - 1);

		names.push_back(aname);
		num_chroms++;
		in >> aname;
	}//while
	in.close(); in.clear();

	string line;

	//find out sizes of chromosomes (or references)
//	ta_genome.resize(1);
	vector<string> chrom_refs_names;

	_uli cur_size;

	for(_ui u = 0; u < num_chroms; u++){
		in.open(names[u].c_str(), ios::in);
		if(!in){
			cout << "ERROR: cannot open " << names[u] << endl;
		space_release();
			exit(0);
		}


		getline(in, line);

		if(line.length() == 0){
			cout << "\nERROR: No input in the file " << names[u] << endl;
					space_release();
			exit(0);
		}
		if(line[0] != '>' ){
			cout << "\nERROR: " << names[u] << " must be in FASTA format." << endl;
					space_release();
			exit(0);
		}

		string achrom;
		istringstream istr(line.substr(1, line.length() -1));
		istr >> achrom;
		chroms_names.insert(make_pair(u, achrom));
		chrom_refs_names.push_back(achrom);
		cur_size = 0;
		getline(in, line);

		while(!in.eof()){
			long int line_len = line.length();
			char last_char = line[line_len - 1] ;
			if(last_char == '\n' || last_char == '\r')//nova changed because AWS reads in EOL character
				line_len--;
			cur_size += line_len;
			getline(in, line);
		}//while

		size_chrom.push_back(cur_size);//in char/bits one char per bit
		_uli new_size = (cur_size >> two_in_Size) ;// =/(8*64) 1pb takes 1 bit, and 1 long takes 64 bits
		if((cur_size & mask[two_in_Size]) > 0)
			new_size += 1;
		chroms_starts.push_back(new_size);


		in.close(); in.clear();
	}//for u

	_ui total_chroms = num_chroms;
	num_chroms <<= 1;// multiply by 2
	for(int j = chrom_refs_names.size() - 1; j >= 0; j--){
			string achrom = chrom_refs_names[j];
			chroms_names.insert(make_pair(total_chroms, achrom));
			total_chroms++;
	}//for j since we have concatenation of genome to itself in this order

	//nova: since we concatenate two strands S_complement and S_reverse, where S is original string
	//chroms will be in order 1,2,3,3,2,1
	_uli last_chr = size_chrom.size() - 1;//index
	_uli size_of_orig_genome = 0;
	for(int j = last_chr; j >= 0; j--){
		size_of_orig_genome += size_chrom[j];
		size_chrom.push_back(size_chrom[j]);
		chroms_starts.push_back(chroms_starts[j]);
	}//for j
	orig_genome_size = size_of_orig_genome;//nova in bases: we use this to convert genomic positions within 
	//a single genome version (ConcatenationOf(genome) or Reverse(genome)) to genomic position within entire Concatenated genome

	_uli cur_total = 0;
	_uli starts_total = 0;
	starts_size_chrom.push_back(cur_total);
	starts_chroms_starts.push_back(starts_total);
	_ui i = 1;
	for(i = 1; i < size_chrom.size(); i++){//nova size_chrom contains 1,2,3,...,n,n,...,3,2,1 chroms
		cur_total += size_chrom[i - 1];
		starts_size_chrom.push_back(cur_total);
		
		starts_chroms_starts.push_back(starts_size_chrom[i] >> two_in_Size);

	}
	binary_total = cur_total + size_chrom[i - 1];//binary size of genome
	
	starts_chroms_starts.push_back(binary_total >> two_in_Size);//the end of the last chrom

	cur_total = (binary_total >> two_in_Size) + 1 ;//index of bucket plus 1 to get size, + 1 extra

	ta_size = cur_total;
	full_binary_total = binary_total << 1;
	full_size = ((binary_total << 1) >> two_in_Size) + 1;
cout << "References sizes are calculated" << endl;

	genome_byte_size = ((binary_total + 1)  >> two_in_byte) + 1;
	ta_genome = new unsigned long int[cur_total];
	assert(ta_genome != 0);
	cg_genome = new unsigned long int[cur_total];
	assert(cg_genome != 0);
	ns_genome = new unsigned long int[cur_total];
	assert(ns_genome != 0);
	full_genome = new unsigned long int[full_size];
	assert(full_genome != 0);
	
cout << "Space is allocated for genome representation" << endl;

	ht_sizes.resize(hash_table_size);
	_usi alphabet_size = 4;
	vector<_usi> dum(alphabet_size, 0);
	num_occ.resize(SIZE, dum);
	num_occ_sqr.resize(alphabet_size, 0);
	num_occ_sqr_sum.resize(alphabet_size, 0);

	C.resize(alphabet_size, 0);
	
	int fast_size = static_cast<int>(pow(static_cast<double>(2), static_cast<double>(STEP)));

	block_256 = pow(2.0, 16.0);
	two_in_block_256 = 16;
	pos_indices.resize(block_256, 0);
	pos_num_sqr = 0;
	pos_num_occ_sqr_sum = 0;
	pos_num_occ.resize(SIZE, 0);
	bwt_marked.resize(block_32, 0);
	pos_strand.resize(block_32, 0);
	pos_indices_ind = 0;

	ht_indices = 0;
	
	bwt.resize(block_32);//64*64

	marked_ns.resize(block_32);
	for(i = 0; i < block_32; i++)
	{	bwt[i] = 0;
		marked_ns[i] = 0;
	}//for i
	curL_bwt = 0;
	cur_ind_bwt = 0;

	int apow_seed = 16;
	_uli hash_small_size = pow(2.0, apow_seed + 0.0);
	num_occ_bytes.resize(hash_small_size);//2^SEED	
	for(_uli i = 0; i < hash_small_size; i++){
		vector< _usi > characters(alphabet_size, 0);
		count_char(i, characters, apow_seed);
        num_occ_bytes[i] = characters;

	}//for

	small_ht_size = pow(2.0, 24.0);
	ht_sizes_lengths.resize(small_ht_size, 0);

	cout << "Initialization is completed" << endl;
}//Anchor()

void Anchor::count_char(_uli kmer, vector<_usi> & characters, int asize){

	for(int i = 0; i < asize; i += 2){
		_uli cur_char = kmer & mask[2];
		characters[cur_char]++;
		kmer >>= 2;
	}
}//count_char

void Anchor::print()
{
	for(_uli i = 0; i < ta_size; i++){
		_uli cur = cg_genome[i];
//		bitset<64> z(cur);
//		cout << z << endl;
	}
}//for test only

void Anchor::build_bwt()
{
	//construct ta-, cg- and ns-genomes
	read_genome();//counts ta, cg, ns and ht_sizes (Hash Table holds # of 24-bit mers in full genome), full_genome

cout << "Genome is read " << endl;

	_uli cur_bwt_char = get_char(binary_total - 1); 

	bwt[0] = cur_bwt_char;
	_uli no_ind = (curL_bwt & mask[two_in_block_64]) >> two_in_Size;//index of one of 64 sub-buckets inside 4096 bucket
	num_occ_sqr[cur_bwt_char]++;
	num_occ[no_ind][cur_bwt_char]++;

	curL_bwt = 1;
	//write to files cg- and ns-genomes

		flash_genome(out_cg_genome, cg_genome, ta_size);
		flash_genome(out_ta_genome, ta_genome, ta_size);
		//WHY flash_genome on ns_genome is not done here? Answer: because we use it
		//If a seed at genome pos has an N char, it is not used (in mapping)

	//free space for cg- and ns-genomes
	delete [] cg_genome;
	cg_genome = 0;
	delete [] ta_genome;
	ta_genome = 0;

	//collect ranges of indices of HashTable representing buckets of size N/32, where N is genome size
	
	//first entry of hash table has to be processed separately
	//NOVA_DEBUG: using total_components = ta_size might cause problem of stack overflow
//	total_components = ta_size ;//sort_first_entry uses quick_multikeys conditioning on total_components
//	sort_first_entry();
	total_components = 5000;

//cout << "first entry is sorted" << endl;

	//now start from pos 1 not 0, because we have processed entry 0 already
	collect_ranges(0, hash_table_size - 1, two_in_32, ranges, ranges_total, ht_sizes);

	_ui ranges_size = ranges.size();
	for(_ui i = 0; i < ranges_size; i++){
		
		collect_bucket_indices(ranges[i].first, ranges[i].second, ranges_total[i]);

		//subdivide indices into sub buckets of size N/(32*32)
		vector< pair<_ui, _ui> > temp_buckets;
		vector<_ui> temp_sum;
		collect_ranges(ranges[i].first, ranges[i].second, two_in_block_64, temp_buckets, temp_sum, ht_sizes);
		collect_subbuckets(temp_sum);

		//for each sub-bucket 
		_ui sub_size = sub_bucket.size();
		for(_ui j = 0; j < sub_size; j++){
			//initialize bucket components and bc-indices
			init_comp_indices(sub_bucket[j].first, sub_bucket[j].second, sub_bucket_sum[j]);

			//sort indices in each subrange
			quick_multikeys(0, bc_indices.begin(), (bc_indices.end() -1) , sub_bucket[j].first, total_components);
			//sort not sorted until everything is sorted
			_uli not_sorted_size = not_sorted.size();
			_uli not_sorted_count = 0;
				total_components += total_components;//*2
			_uli count_total_ta_size = 0;
//DEBUG
cout << "quick_multikeys for j= " << j << " of " << sub_size << ", i = " << i << " of "
	<< ranges_size << ", not_sorted size = " << not_sorted_size << endl;
			while(!not_sorted.empty()){

				pair< pair<_uli_iter, _uli_iter>, pair<_ui, _ui> > cur_pair = not_sorted.front();
				not_sorted.pop();

				quick_multikeys(cur_pair.second.first, cur_pair.first.first, cur_pair.first.second, cur_pair.second.second, total_components); 
				not_sorted_count++;
				//each of ranges in not_sorted has been processed with next value of total_components
				//so increment total_components
				if(not_sorted_count == not_sorted_size){
					total_components += min(total_components + 5000, full_size );
					not_sorted_size = not_sorted.size();
					not_sorted_count = 0;
				}//
				if(total_components == full_size)
					count_total_ta_size++;
				if(count_total_ta_size > 1)
					break;
			}//while
			total_components = 5000;
			//update bwt, num_occ, pos_indices
			update_bwt(sub_bucket[j].first);

		}//for j

	}//for i
cout << "Sorting is complete" << endl;
	//check if the there is something needed to flash in bwt, num_occ and pos_indices
	finish_update_bwt(sub_bucket[sub_bucket.size() - 1].first);

	//output last genome: we don't need to output full genome


	//finish additional infor output
	output_info();
//DEBUG
cout << "output_info is done " << endl;
	precalculate_seeds();
//DEBUG
cout << "precalculated_seeds are done" << endl;
}//build_bwt

void Anchor::update_bwt(_ui sub_bucket_st)
{
	_uli azero = 0;
	_uli asize = bc_indices.size();

	for(_uli i = 0; i < asize; i++){
		_uli ht_ind = bc_indices[i] + sub_bucket_st;
		_uli gen_ind = ht_indices[ht_ind];//real concatenated genome index
		_uli cur_bwt_char = get_char(gen_ind - 1);//prev char
		marked_ns[cur_ind_bwt] <<= 1;
		bwt_marked[cur_ind_bwt] <<= 1;
		pos_strand[cur_ind_bwt] <<= 1;

		_uli start_lmer = get_ns_lmer(gen_ind, SEED);//if 'N' in genome at this pos within 12char
		if(start_lmer > 0)
			marked_ns[cur_ind_bwt] |= 1; //if 'N' in genome at this pos within 12char
		if(gen_ind > 0){
			bwt[cur_ind_bwt] = (bwt[cur_ind_bwt] << 2) | cur_bwt_char;
			//update no_occ
			_uli no_ind = (curL_bwt & mask[two_in_block_64]) >> two_in_Size;//index of one of 64 sub-buckets inside 4096 bucket
			num_occ_sqr[cur_bwt_char]++;
			num_occ[no_ind][cur_bwt_char]++;
		}//gen_ind > 0
		else{
			bwt[cur_ind_bwt] = (bwt[cur_ind_bwt] << 2) | 0;// $ character which is 0 here
			orig_gen_inL = curL_bwt; 
		}//else
	
		//update pos indices
		if((gen_ind & mask[two_in_Step]) == 0){//each 16-th position is saved
			//check if gen_ind is in the first half of the concatenation (strand "-", corresponds to ComplementOf(OriginalGenome))
			//If gen_ind is in the second half, we start pos indices from start, because
			//for human genome, genomic positions can go as high as 2^32 only
			//to fit genomic pos into _ui, we need to count genomic positions for the second half of genome starting with 0 
			//If pos_strand[i] = 1, that means that correspondin genomic position comes from second half
			//We mark only positions pos_strand[i] for which we actually store genomic positions
			_uli converted_gen_ind = gen_ind;
			if(gen_ind >= orig_genome_size){
				converted_gen_ind = gen_ind - orig_genome_size;
				pos_strand[cur_ind_bwt] |= 1; //genomic position in the second half of concatenated genome
			}
			pos_indices[pos_indices_ind] = converted_gen_ind;
			pos_indices_ind++;
			pos_indices_total++;
			bwt_marked[cur_ind_bwt] |= 1;
		
			_uli no_ind = (curL_bwt & mask[two_in_block_64]) >> two_in_Size;//index of one of 64 sub-buckets inside 4096 bucket
			pos_num_occ[no_ind]++;
			pos_num_sqr++;
		}

		curL_bwt++;
		if((curL_bwt & mask[two_in_block_64]) == 0){ //NOVA_NEXT
			num_occ_sqr_sum[0] += num_occ_sqr[0];
			num_occ_sqr_sum[1] += num_occ_sqr[1];
			num_occ_sqr_sum[2] += num_occ_sqr[2];
			num_occ_sqr_sum[3] += num_occ_sqr[3];

			num_occ_sqr[0] = num_occ_sqr_sum[0];
			num_occ_sqr[1] = num_occ_sqr_sum[1];
			num_occ_sqr[2] = num_occ_sqr_sum[2];
			num_occ_sqr[3] = num_occ_sqr_sum[3];
			for(_ui i = 1; i < SIZE; i++){
				num_occ[i][0] += num_occ[i - 1][0];
				num_occ[i][1] += num_occ[i - 1][1];
				num_occ[i][2] += num_occ[i - 1][2];
				num_occ[i][3] += num_occ[i - 1][3];
			}//for
			
			pos_num_occ_sqr_sum += pos_num_sqr;
			pos_num_sqr = pos_num_occ_sqr_sum;
			for(_ui i = 1; i < SIZE; i++)
				pos_num_occ[i] += pos_num_occ[i-1];

			flash_num_occ(SIZE);//pos_num_sqr and pos_num_occ flashed here too
			init_num_occ();
		}//if

		if((curL_bwt & mask[two_in_half_word]) == 0){
			cur_ind_bwt++;
			if(cur_ind_bwt == block_32){
				flash_bwt();
				flash_marked_ns();
				init_marked();
				cur_ind_bwt = 0;
			}
		}
		if(pos_indices_ind == block_256){
			flash_pos_indices(block_256);
			init_pos_indices();
		}//if
	}//for

}//update_bwt

void Anchor::finish_update_bwt(_ui sub_bucket_st)
{
		if((curL_bwt & mask[two_in_block_64]) > 0){
			num_occ_sqr_sum[0] += num_occ_sqr[0];
			num_occ_sqr_sum[1] += num_occ_sqr[1];
			num_occ_sqr_sum[2] += num_occ_sqr[2];
			num_occ_sqr_sum[3] += num_occ_sqr[3];

			num_occ_sqr[0] = num_occ_sqr_sum[0];
			num_occ_sqr[1] = num_occ_sqr_sum[1];
			num_occ_sqr[2] = num_occ_sqr_sum[2];
			num_occ_sqr[3] = num_occ_sqr_sum[3];
			for(_ui i = 1; i < SIZE; i++){
				num_occ[i][0] += num_occ[i - 1][0];
				num_occ[i][1] += num_occ[i - 1][1];
				num_occ[i][2] += num_occ[i - 1][2];
				num_occ[i][3] += num_occ[i - 1][3];
			}//for
			pos_num_occ_sqr_sum += pos_num_sqr;
			pos_num_sqr = pos_num_occ_sqr_sum;
			for(_ui i = 1; i < SIZE; i++)
				pos_num_occ[i] += pos_num_occ[i-1];

			_uli cur_ind = ((curL_bwt - 1) & mask[two_in_block_64]) >> two_in_Size;//index of one of 64 sub-buckets inside 4096 bucket
			flash_num_occ(cur_ind + 1);//up to not included this index
		}//if
			
		if((curL_bwt & mask[two_in_half_word]) > 0){
				//curLbwt points to unused base or bit
				bwt[cur_ind_bwt] <<= (SIZE - ((curL_bwt << 1 ) & mask[two_in_Size]) );//DEBUG POSSIBLE
				marked_ns[cur_ind_bwt] <<= ((SIZE >> 1) - (curL_bwt & mask[two_in_half_word]));
				bwt_marked[cur_ind_bwt] <<= ((SIZE >> 1) - (curL_bwt & mask[two_in_half_word]) );
				pos_strand[cur_ind_bwt] <<= ((SIZE >> 1) - (curL_bwt & mask[two_in_half_word]) );
				cur_ind_bwt++;

			flash_bwt(); //upto and including cur_ind_bwt
			flash_marked_ns();
			cur_ind_bwt = 0;
		}//if bwt is not full
		else if((curL_bwt & mask[two_in_half_word]) == 0){
			if(cur_ind_bwt > 0){
				//cur_ind_bwt points to an unused 64-bit word, so just flush upto cur ind that is unused
				flash_bwt(); //upto and including cur_ind_bwt
				flash_marked_ns();
				cur_ind_bwt = 0;
			}
		}//else

		if(pos_indices_ind < block_256){
			flash_pos_indices(pos_indices_ind);
			init_pos_indices();
		}//if
        //output one additional SIZE-word
        bwt[0] = 0;
        marked_ns[0] = 0;
		bwt_marked[0] = 0;
		pos_strand[0] = 0;
        cur_ind_bwt = 1;
        flash_bwt(); //upto not including cur_ind_bwt
		flash_marked_ns();
		// Not efficient: flash_genome(out_pos_strand, pos_strand, ta_size);

}//finish_update_bwt

_uli Anchor::get_char(_uli ind)//
{
	//index is given in 1-bit representation
	if((ind << 1) > full_binary_total)
		return 0;
	else if((ind << 1) == full_binary_total)
		return ((full_genome[0] >> (SIZE - 2)) & mask[2]);
	_uli ta_ind = (ind << 1) >> two_in_Size;
	_uli inside_ind = (ind << 1) & mask[two_in_Size];
	_uli res = full_genome[ta_ind];
	res = (res >> (SIZE - inside_ind - 2)) & mask[2];//2 for 2-bit full-genome representation
	return res;
}//get_char

void Anchor::init_num_occ()
{
	for(_ui i = 0; i < SIZE; i++){
		num_occ[i][0] = 0;
		num_occ[i][1] = 0;
		num_occ[i][2] = 0;
		num_occ[i][3] = 0;
	}
	num_occ_sqr[0] = 0;
	num_occ_sqr[1] = 0;
	num_occ_sqr[2] = 0;
	num_occ_sqr[3] = 0;

	pos_num_sqr = 0;
	for(_uli i = 0; i < SIZE; i++)
		pos_num_occ[i] = 0;

}//init_num_occ
//lookup8.c, by Bob Jenkins, January 4 1997, Public Domain.

void Anchor::init_pos_indices()
{
	for(_ui i = 0; i < block_256; i++)
	{ 
		pos_indices[i] = 0;
	}//for i
	pos_indices_ind = 0;
}//init_pos_indices

void Anchor::init_marked()
{
	for(_ui i = 0; i < block_32; i++){
		marked_ns[i] = 0;
		bwt_marked[i] = 0;
		pos_strand[i] = 0;
	}
}

void Anchor::make_mask(){
	_uli one = 1;
	mask.push_back(0);
	for(_uli i = 1; i <= SIZE; i++){
		_uli res = (one << i) - 1;
		mask.push_back(res);
	}//for
	mask[SIZE] = 0xFFFFFFFFFFFFFFFF;
}//make mask


//read genome: First read original genome S, and take its complement
//Since we built ReverseOf(ReverseComplementOf(S))*ReverseOf(S) with Cs converted to Ts
//ReverseOf(ReverseComplementOf(S)) = ComplementOf(S) with Cs converted to Ts.
//So, read S, take its complement and convert Cs to Ts
// G's complement is C, convert it to T, so if "G" is read, convert it to T, the rest of characters have usual complement
//Another conversion is Ns to As
//Next function read_genome2() will read S and take its reverse and convert Cs to Ts, and Ns to As,
//and concatenate this version of genome to the first version, ComplementOf(S)
void Anchor::read_genome()
{

	unsigned long int index = 0;
	unsigned long int full_ind = 0;
	ifstream in;

	const long int buffer_size = 10000000;
	char *buffer = new char[buffer_size];
	assert( buffer != 0);

	unsigned long int real_size = 0;
	unsigned long int count_bits = 0;

	unsigned long int lmer_1 = len_mer - 1;
	unsigned long int ns = 0;
	unsigned long int cg = 0;
	unsigned long int ta = 0;
	unsigned long int ta_seed = 0;

	unsigned long int full = 0;
	unsigned long int full_seed = 0;
	unsigned long int full_count_bits = 0;


	//to calculate consecutive ranges of As whose length is greater than or equal to SEED
	_uli range_As_len = 0;
	bool consecutive_a = false;
	max_range_As = 0;

	for(_ui u = 0; u < names.size(); u++){ 

		in.open(names[u].c_str(), ios::in);
		if(!in){
			cout << "can't open " << names[u] << " names[u] reading genome" << endl;
					space_release();
			exit(0);
		}
		string line;
		getline(in, line); //the first line has >

		in.read(buffer, buffer_size);
		real_size = in.gcount();
		if(real_size == 0){
			cout << "\nERROR: " << names[u] << " has no input " << endl;
					space_release();
			exit(0);
		}

		_uli last_lmer_char = 0;


		if(u == 0){
			string temp_seed("");
			_ui i = 0;
			while(i < ((SEED >> 1) - 1)) {// i < lmer_1){
				if(buffer[last_lmer_char] != '\n' && buffer[last_lmer_char] != '\r'){
					temp_seed = temp_seed + buffer[last_lmer_char];
					i++;
				}
				last_lmer_char++;
			}//while

			for(i = 0; i < ((SEED >> 1) - 1); i++){//i = 0; i < lmer_1; i++){	
				char next_char = temp_seed[i];
				cg <<= 1;
				ta <<= 1;
				ns <<= 1;
				full <<=2; 

				if(next_char == 'T' || next_char == 't'){//T's complement is A
					/*
					//OLD for char 'T'
					ta |= 1;
					C[3]++;
					full |= 3;

					consecutive_a = false;
					if(range_As_len >= half_SEED){
						if(max_range_As < range_As_len)
							max_range_As = range_As_len;
						ends_ranges_of_As.push_back(count_bits);
					}
					range_As_len = 0;
					// END OF OLD for char 'T'
					*/
					//OLD for 'A' 
					//Now we take complement of T, it becomes character A, so use "Old for 'A'"
					range_As_len++;//length is 
					consecutive_a = true;
					C[0]++;
					//END of OLD for 'A'

				}
				else if(next_char == 'C' || next_char == 'c' ){
					/*
					//OLD for 'C'
					cg |= 1;
					ta |= 1;//for c In BRAT-bw ta-representation has T and C represented as 1, and (A and G) represented with bit 0
					if(Genome == 1){//C->T conversion
						full |= 3;
						C[3]++;					
					}
					else{//normal or G->A conversion
						full |= 1; 
						C[1]++;						
					}
					consecutive_a = false;
					if(range_As_len >= half_SEED){
						if(max_range_As < range_As_len)
							max_range_As = range_As_len;
						ends_ranges_of_As.push_back(count_bits);
					}
					range_As_len = 0;
					//END of OLD for 'C'
					*/

					//OLD for 'G';
					//Complement of C is G
					cg |=1;
					full |= 2; 
					C[2]++;		

					consecutive_a = false;
					if(range_As_len >= half_SEED){
						if(max_range_As < range_As_len)
							max_range_As = range_As_len;
						ends_ranges_of_As.push_back(count_bits);
					}//if
					range_As_len = 0;						
					//END of OLD for 'G'

				}
				else if(next_char == 'G' || next_char == 'g'){
					/*
					//OLD for 'G'
					cg |=1;
						if(Genome == 2){
							range_As_len++;//length is in bits
							consecutive_a = true;
						}
						else{
							consecutive_a = false;
							if(range_As_len >= half_SEED){
								if(max_range_As < range_As_len)
									max_range_As = range_As_len;
								ends_ranges_of_As.push_back(count_bits);
							}
							range_As_len = 0;						
						}//else if normal or C->T
					if(Genome == 2){//G->A conversion
						full |= 0;
						C[0]++;					
					}
					else{
						full |= 2; 
						C[2]++;						
					}
					//END of OLD for 'G'
					*/

					//OLD for 'C'	Now G's complement is C (here we need to take complement)
					cg |= 1;
					ta |= 1;//for c In BRAT-bw ta-representation has T and C represented as 1, and (A and G) represented with bit 0
					//C->T conversion
					full |= 3;
					C[3]++;					

					consecutive_a = false;
					if(range_As_len >= half_SEED){
						if(max_range_As < range_As_len)
							max_range_As = range_As_len;
						ends_ranges_of_As.push_back(count_bits);
					}
					range_As_len = 0;
					//END of OLD for 'C'

				}//if char 'G', take complement(G) = 'C', convert 'C' to 'T'
				else if(next_char == 'A' || next_char == 'a'){
					/*
					//OLD for 'A'
					range_As_len++;//length is 
							consecutive_a = true;
						C[0]++;
					//END of OLD for 'A'
					*/
					//OLD for char 'T',
					//since complement('A') = 'T'
					ta |= 1;
					C[3]++;
					full |= 3;

					consecutive_a = false;
					if(range_As_len >= half_SEED){
						if(max_range_As < range_As_len)
							max_range_As = range_As_len;
						ends_ranges_of_As.push_back(count_bits);
					}
					range_As_len = 0;
					// END OF OLD for char 'T'
				}
				else{
					/*
					//OLD for 'N'
					ns |= 1;
					C[0]++;
							range_As_len++;//length is in bases
							consecutive_a = true;
					//END of OLD for 'N'
					*/
					ns |= 1; //this is N's position 
					//OLD for char 'T', 
					//'N' is converted to 'A' in original genome
					//But now we take complement of original genome,
					//Hence, complement('A') = 'T' 
					ta |= 1;
					C[3]++;
					full |= 3;

					consecutive_a = false;
					if(range_As_len >= half_SEED){
						if(max_range_As < range_As_len)
							max_range_As = range_As_len;
						ends_ranges_of_As.push_back(count_bits);
					}
					range_As_len = 0;
					// END OF OLD for char 'T'					
				}

				count_bits++;
				full_count_bits += 2;
				if((count_bits & mask[two_in_Size]) == 0){//BUG possible: change byte_size = 4
					ta_genome[index] =ta ;
					cg_genome[index] = cg ;
					ns_genome[index] = ns;
					index++;
				}//if
				if((count_bits & mask[two_in_half_word]) == 0){
					full_genome[full_ind] = full;
					full_ind++;
				}
			}//for i
		}//if chrom_id = 0

		while(true){
			for(; last_lmer_char < real_size; last_lmer_char++){
				char next_char = buffer[last_lmer_char];
				if(next_char != '\n' && next_char != '\r'){
					cg <<= 1;
					ta <<= 1;
					ns <<= 1;
					full <<= 2;
					if(next_char == 'T' || next_char == 't'){
						/*
						//OLD for 'T'
						ta |= 1;
						C[3]++;
						full |= 3;					
						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}
						range_As_len = 0;
						//END of OLD for 'T'
						*/
					
						//OLD for 'A' 
						//Now we take complement of T, it becomes character A, so use "Old for 'A'"
						range_As_len++;//length is 
						consecutive_a = true;
						C[0]++;
						//END of OLD for 'A'
					}//if 'T'
					else if(next_char == 'C' || next_char == 'c' ){
						/*
						//OLD for 'C'
						cg |= 1;
						ta |= 1;//for c
						if(Genome == 1){
							full |= 3;
							C[3]++;					
						}
						else{
							full |= 1; 
							C[1]++;						
						}
						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}
						range_As_len = 0;
						//END of OLD for 'C'
						*/

						//OLD for 'G';
						//Complement of C is G
						cg |=1;
						full |= 2; 
						C[2]++;		

						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}//if
						range_As_len = 0;						
						//END of OLD for 'G'
					}//if 'C' in genome
					else if(next_char == 'G' || next_char == 'g'){
						/*
						//OLD for 'G'
						cg |=1;
						if(Genome == 2){
							full |= 0;
							C[0]++;					
						}
						else{
							full |= 2; 
							C[2]++;						
						}
						if(Genome == 2){
							range_As_len++;//length is in bits
							consecutive_a = true;
						}
						else{
						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}
						range_As_len = 0;	
						}//else if normal or C->T
						
						//END of OLD for 'G'
						*/

						//OLD for 'C'	Now G's complement is C (here we need to take complement)
						cg |= 1;
						ta |= 1;//for c In BRAT-bw ta-representation has T and C represented as 1, and (A and G) represented with bit 0
						//C->T conversion
						full |= 3;//full is used to sort suffixes, so need converted C to T
						C[3]++;	//for BWT Cs converted to Ts				

						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}
						range_As_len = 0;
						//END of OLD for 'C'
					}//if 'G' in genome
					else if(next_char == 'A' || next_char == 'a'){
						/*
						//OLD for 'A'
						C[0]++;
							range_As_len++;//length is in 
							consecutive_a = true;	
						//END OLD for 'A'
						*/

						//OLD for char 'T',
						//since complement('A') = 'T'
						ta |= 1;
						C[3]++;
						full |= 3;

						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}
						range_As_len = 0;
						// END OF OLD for char 'T'
					}//if 'A' in genome
					else{
						/* 
						//OLD for 'N'
						ns |= 1;//
						C[0]++;
							range_As_len++;//length is in 
							consecutive_a = true;
						//END of OLD for 'N'
						*/

						ns |= 1; //Since this is 'N' is genome

						//OLD for char 'T',
						//since complement('A') = 'T'
						ta |= 1;
						C[3]++;
						full |= 3;

						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}
						range_As_len = 0;
						// END OF OLD for char 'T'
					}
		
					full_seed = full & MASK_LMER;
					ht_sizes[full_seed]++;

					count_bits++;
					full_count_bits += 2;

					if((count_bits & mask[two_in_Size]) == 0){//BUG possible: change byte_size = 4
						ta_genome[index] = ta;
						cg_genome[index] = cg;
						ns_genome[index] = ns;
						index++;
					}//if
					if((count_bits & mask[two_in_half_word]) == 0){
						full_genome[full_ind] = full;
						full_ind++;
					}
				}//if not end of line char
			}//for i buffer size
			last_lmer_char = 0;
			if(in.eof())
				break;
			in.read(buffer, buffer_size);
			real_size = in.gcount();
			if(real_size == 0)
				break;
		}//while
		in.close(); in.clear();
		cout << names[u] << " is preprocessed " << endl;
	}//for u all chroms

	cout << "Complement of genome is taken with Cs converted to Ts" << endl;

/**********NOVA: Now Read in genome, take its reverse and convert Cs to Ts
			And concatenate the result to the previous result (complement of original genome with Cs changed to Ts)
**************************************************************************************************/

	for(int u = names.size() - 1; u >= 0;  u--){ //in reverse order 

		in.open(names[u].c_str(), ios::in);
		if(!in){
			cout << "can't open " << names[u] << " names[u] reading genome" << endl;
					space_release();
			exit(0);
		}
		string line;
		getline(in, line); //the first line has >

		in.read(buffer, buffer_size);
		real_size = in.gcount();
		if(real_size == 0){
			cout << "\nERROR: " << names[u] << " has no input " << endl;
					space_release();
			exit(0);
		}

		//read in entire chromosome
		const long int current_size = size_chrom[u];
		char *cur_chromosome = new char[current_size];
		assert( cur_chromosome != 0);
		_uli cur_chrom_ind = 0;

		_uli last_lmer_char = 0;

		while(true){
			for(; last_lmer_char < real_size; last_lmer_char++){
				char next_char = buffer[last_lmer_char];
				if(next_char != '\n' && next_char != '\r'){
					cur_chromosome[cur_chrom_ind] = next_char;
					cur_chrom_ind++;
				}//if not end of line char
			}//for i buffer size
			last_lmer_char = 0;
			if(in.eof())
				break;
			in.read(buffer, buffer_size);
			real_size = in.gcount();
			if(real_size == 0)
				break;
		}//while
		in.close(); in.clear();

		//check that chromosome size is correct
		if(cur_chrom_ind != current_size){
			cout << "ERROR: Reading genome second time, total number of characters read differs\n"
				<< "from previously calculated size of chromosome " << names[u] << endl;
			space_release();
			exit(0);
		}//if
		
		//Process current chromosome in reverse order, converting Cs to Ts at the same time
		//continue building ta-, cg-, ns-, full- genome representatives
		for(long int ii = cur_chrom_ind - 1; ii >= 0; ii--){
			
			char next_char = cur_chromosome[ii];
			cg <<= 1;
			ta <<= 1;
			ns <<= 1;
			full <<= 2;
			if(next_char == 'T' || next_char == 't'){
						//OLD for 'T'
						ta |= 1;
						C[3]++;
						full |= 3;					
						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}
						range_As_len = 0;
						//END of OLD for 'T'

					}//if 'T'
					else if(next_char == 'C' || next_char == 'c' ){
						
						//OLD for 'C'
						cg |= 1;
						ta |= 1;//for c
						//convert C to T
						full |= 3;
						C[3]++;					

						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}
						range_As_len = 0;
						//END of OLD for 'C'
						
					}//if 'C' in genome
					else if(next_char == 'G' || next_char == 'g'){
						//OLD for 'G'
						cg |= 1;
						full |= 2; 
						C[2]++;						

						consecutive_a = false;
						if(range_As_len >= half_SEED){
							if(max_range_As < range_As_len)
								max_range_As = range_As_len;
							ends_ranges_of_As.push_back(count_bits);
						}
						range_As_len = 0;	
	
						//END of OLD for 'G'
					}//if 'G' in genome
					else if(next_char == 'A' || next_char == 'a'){

						//OLD for 'A'
						C[0]++;
						range_As_len++;//length is in 
						consecutive_a = true;	
						//END OLD for 'A'

					}//if 'A' in genome
					else{
						//OLD for 'N'
						ns |= 1;//N is in genome at this position
						C[0]++;//convert 'N' to 'A' in the original genome
						range_As_len++;//length is in 
						consecutive_a = true;
						//END of OLD for 'N'
					}//if 'N' in genome
		
					full_seed = full & MASK_LMER;
					ht_sizes[full_seed]++;

					count_bits++;
					full_count_bits += 2;

					if((count_bits & mask[two_in_Size]) == 0){//BUG possible: change byte_size = 4
						ta_genome[index] = ta;
						cg_genome[index] = cg;
						ns_genome[index] = ns;
						index++;
					}//if
					if((count_bits & mask[two_in_half_word]) == 0){
						full_genome[full_ind] = full;
						full_ind++;
					}

		}//for ii reversing original chrom

		

		//Delete current chromosome: release space
		delete [] cur_chromosome;
	}//for u all chroms



/**********NOVA: End of ConcatenationOf( ComplementOf(genome)WithCsConvertedToTs, ReverseOf(genome)WithCsConvertedToTs
***************************************************************************************************************/

	if((count_bits & mask[two_in_Size]) > 0){
		ta <<= (SIZE - (count_bits & mask[two_in_Size]));
		cg <<= (SIZE - (count_bits & mask[two_in_Size]));
		ns <<= (SIZE - (count_bits & mask[two_in_Size]));
		ta_genome[index] = ta;
		cg_genome[index] = cg;
		ns_genome[index] = ns;
		index++;
	 }//if
	if((count_bits & mask[two_in_half_word]) > 0){
		full <<= (SIZE - ((count_bits & mask[two_in_half_word]) << 1));
		full_genome[full_ind] = full;
		full_ind++;
	}
	cout << "count_bits before shift = " << count_bits << endl;
	count_bits <<= 1;//corresponds to 2-bit representation
	cout << "count_bits after shift = " << count_bits << endl;
	cout << "full_binary_total = " << full_binary_total << endl;
	_uli i;
	for( i = count_bits - SEED + 2; i < full_binary_total; i += 2){//starting with SEED - 1 character of 2-bit representaion
		//2 abover is because 2-bits per character, and Hash Table is still on 24bits-SEED (i.e. 12 characters)
		_uli res = get_ta_lmer(i, len_mer);//using full-genome, i is in 2-bit repres
		ht_sizes[res]++;
	}//for i
	if(consecutive_a == true){
		range_As_len += 11;//last bases are added to process the last bases of genome to make for SEED starting with last bases of genome
		if(max_range_As < range_As_len)
			max_range_As = range_As_len;
		if(range_As_len >= half_SEED){
			ends_ranges_of_As.push_back(count_bits >> 1 );
		}
	}
	cout << "Max range of As = " << max_range_As << endl;
	shift_ht_lengths = ceil(log(max_range_As + 0.0)/log(2.0));
	if(shift_ht_lengths <= SEED)
		shift_ht_lengths = 0;

	cout << "shift_ht_lengths = " << shift_ht_lengths << endl;
	delete [] buffer;
	cout << "Finished reading genome" << endl;
}//read_genome


void Anchor::collect_ranges(_ui st, _ui en, _ui two_in_bucket, vector< pair<_ui, _ui> > &buckets, vector<_ui> &buckets_sum,
	vector<_uli> &ht_sizes_dif)
{
        _uli cur_st = st, cur_en = en;
		_uli cur_sum = 0;
        _uli i = st;
        _uli bucket_size = binary_total >> two_in_bucket;//don't need to change this: it is count of entries in bucket
        if(binary_total <= (1 << two_in_bucket)){
                for(i = st; i <= en; i++){
                        cur_sum += ht_sizes_dif[i];
                }
                buckets.push_back(make_pair(st, en));
                buckets_sum.push_back(cur_sum);
                return;
        }//

	for(i = st; i <= en; i++){
		if(cur_sum >= bucket_size){
			cur_en = i - 1;
			buckets.push_back(make_pair(cur_st, cur_en));
			cur_st = i; cur_en = i;
			buckets_sum.push_back(cur_sum);
			cur_sum = 0;
		}
		cur_sum += ht_sizes_dif[i];
	}//for
	if(cur_sum > 0){
		cur_en = i - 1;
		buckets.push_back(make_pair(cur_st, cur_en));
		buckets_sum.push_back(cur_sum);
	}
}//collect_ranges

void Anchor::collect_subbuckets(vector<_ui> &bucket_sum)
{
	if(!sub_bucket.empty()){
		sub_bucket.clear();
	}
	if(!sub_bucket_sum.empty())
		sub_bucket_sum.clear();

	_uli cur_st = 0, cur_en = 0;
	_uli cur_sum = 0;
	_uli i = 0;
	_uli asize = bucket_sum.size();
	for(i = 0; i < asize; i++){
		cur_st = cur_sum;
		cur_en = cur_sum + bucket_sum[i] - 1;
		sub_bucket.push_back(make_pair(cur_st, cur_en));
		sub_bucket_sum.push_back(bucket_sum[i]);
		cur_sum += bucket_sum[i];
	}

}//collect_subbuckets

void Anchor::collect_bucket_indices(_ui st, _ui en, _ui asize)
{
	if(ht_indices != 0){
		delete [] ht_indices;
		ht_indices = 0;
	}

	ht_indices = new _uli[asize];
	_uli i = 0;
	_uli cur_total = 0;
	vector<_uli> cur_ptr(en - st + 1);
	for(i = st; i <= en; i++){
		cur_ptr[i - st] = cur_total;
		cur_total += ht_sizes[i];
	}//for
	_uli cur_ind = 0;
	//scan entire genome and if seeds are in the given range, then keep them
	i = 1;
	_uli prev = full_genome[0];
	while(i < full_size){
		_uli next = full_genome[i];
		for(_ui j = 0; j <= SIZE - SEED; j += 2){
			_uli left = (prev >> (SIZE -j - SEED)) & mask[SEED];
			if((left >= st) && (left <= en)){
				ht_indices[cur_ptr[left - st]] = cur_ind;
				cur_ptr[left - st]++;
			}//if
			cur_ind++;
		}//for j
		for(_ui j = SIZE - SEED + 2; j < SIZE; j +=2 ){
			_uli left = prev & mask[SIZE - j];
			_uli right = (next >> (SIZE - (SEED - (SIZE - j)))) & mask[SEED - (SIZE - j)];
			left = ((left << (SEED - (SIZE - j))) | right ) & mask[SEED];
			if((left >= st) && (left <= en)){
				ht_indices[cur_ptr[left - st]] = cur_ind;
				cur_ptr[left - st]++;
			}//if
			cur_ind++;
		}//for j
		prev = next;
		i++;
	}//while
	//process last 
	_uli double_cur = cur_ind;
	for(_uli j = (double_cur << 1); j < full_binary_total; j +=2 ){
		_uli res = get_ta_lmer(j, SEED);
		if((res >= st) && (res <= en)){
			ht_indices[cur_ptr[res - st]] = cur_ind;
			cur_ptr[res - st]++;
		}//if
		cur_ind++;
	}//for j
}//collect_bucket_indices

void Anchor::collect_bucket_indices_first(_ui st, _ui en, _ui asize, vector<_uli> &count_lengths)
{
	if(ht_indices != 0){
		delete [] ht_indices;
		ht_indices = 0;
	}
	if(!lengths.empty()){
		lengths.clear();
	}
	if(!bc_indices.empty())
		bc_indices.clear();

	bc_indices.resize(asize, 0);
	for(_ui j = 0; j < asize; j++)
		bc_indices[j] = j;

	ht_indices = new _uli[asize];
	lengths.resize(asize, 0);

	_uli i = 0;
	_uli cur_total = 0;
	vector<_uli> cur_ptr(en - st + 1);
	for(i = st; i <= en; i++){
		cur_ptr[i - st] = cur_total;
		cur_total += ht_sizes_lengths[i];
	}//for
	_uli cur_ind = 0;
	//scan entire genome and if seeds are in the given range, then keep them
	i = 1;

	_uli prev = full_genome[0];
	while(i < full_size){
		_uli next = full_genome[i];
		for(_ui j = 0; j <= SIZE - SEED; j +=2 ){
			_uli left = (prev >> (SIZE -j - SEED)) & mask[SEED];
			if(left == 0){

				long int range_a_ind = lower_bound(ends_ranges_of_As.begin(), ends_ranges_of_As.end(), cur_ind) - ends_ranges_of_As.begin();
				_uli alength = (ends_ranges_of_As[range_a_ind] - cur_ind) ; 
				_uli len_ind = alength >> shift_ht_lengths;
	
				if((len_ind >= st) && (len_ind <= en) && (ends_ranges_of_As[range_a_ind] != binary_total)){
					ht_indices[cur_ptr[len_ind - st]] = cur_ind;
					lengths[cur_ptr[len_ind - st]] = alength;
					count_lengths[alength]++;
					cur_ptr[len_ind - st]++;
				}//if
			}//if 
			cur_ind++;
		}//for j
		for(_ui j = SIZE - SEED + 2; j < SIZE; j +=2){
			_uli left = prev & mask[SIZE - j];
			_uli right = (next >> (SIZE - (SEED - (SIZE - j)))) & mask[SEED - (SIZE - j)];
			left = ((left << (SEED - (SIZE - j))) | right ) & mask[SEED];
			if(left == 0){

				long int range_a_ind = lower_bound(ends_ranges_of_As.begin(), ends_ranges_of_As.end(), cur_ind) - ends_ranges_of_As.begin();
				_uli alength = (ends_ranges_of_As[range_a_ind] - cur_ind); 
				_uli len_ind = alength >> shift_ht_lengths;
				if((len_ind >= st) && (len_ind <= en) && (ends_ranges_of_As[range_a_ind] != binary_total)){
					ht_indices[cur_ptr[len_ind - st]] = cur_ind;
					lengths[cur_ptr[len_ind - st]] = alength;
					count_lengths[alength]++;
					cur_ptr[len_ind - st]++;
				}//if
			}
			cur_ind++;
		}//for j
		prev = next;
		i++;
	}//while
	//process last 
	_uli double_cur = cur_ind;

	for(_uli j = (double_cur << 1); j < full_binary_total; j += 2){
		_uli res = get_ta_lmer(j, SEED);
		if(res == 0){

			long int range_a_ind = lower_bound(ends_ranges_of_As.begin(), ends_ranges_of_As.end(), cur_ind) - ends_ranges_of_As.begin();
			_uli alength = (ends_ranges_of_As[range_a_ind] - cur_ind); 
			_uli len_ind = alength >> shift_ht_lengths;
			if((len_ind >= st) && (len_ind <= en) && (ends_ranges_of_As[range_a_ind] != binary_total)){
				ht_indices[cur_ptr[len_ind - st]] = cur_ind;
				lengths[cur_ptr[len_ind - st]] = alength;
				count_lengths[alength]++;
				cur_ptr[len_ind - st]++;
			}//if
		}//if
		cur_ind++;
	}//for j
}//collect_bucket_indices_first

void Anchor::collect_indices_first_eos(_uli &count_eos)
{
	//corrects sizes of entries of ht_sizes_lengths excluding those having EOS (end of suffix char)
	if(ht_indices != 0){
		delete [] ht_indices;
		ht_indices = 0;
	}
	if(!lengths.empty()){
		lengths.clear();
	}
	if(!bc_indices.empty())
		bc_indices.clear();
	
	vector<_uli> temp;

	_uli i = 0;
	_uli cur_total = 0;
	_uli cur_ind = 0;
	//scan entire genome and if seeds are in the given range, then keep them
	i = 1;
	count_eos = 0;
	long int ends_ranges_of_As_size =  ends_ranges_of_As.size();
cout << "Starting collect_indices_first_eos mask.size = " << mask.size() << endl;
	_uli prev = full_genome[0];
	while(i < full_size){
		_uli next = full_genome[i];
		for(_uli j = 0; j <= SIZE - SEED; j += 2){
			_uli left = (prev >> (SIZE -j - SEED)) & mask[SEED];
			if(left == 0){

				long int range_a_ind = lower_bound(ends_ranges_of_As.begin(), ends_ranges_of_As.end(), cur_ind) - ends_ranges_of_As.begin();

				_uli alength = (ends_ranges_of_As[range_a_ind] - cur_ind) ; 
				if(ends_ranges_of_As[range_a_ind] == binary_total){
					count_eos++;
					lengths.push_back(alength);
					temp.push_back(cur_ind);
				}//if
				else{
					_uli len_ind = alength >> shift_ht_lengths;
					if(len_ind < small_ht_size)
						ht_sizes_lengths[len_ind]++;//need this to sort by length		
					else{

						cout << "ERROR: alength = " << alength << ", len_ind = " << len_ind << " >= " << small_ht_size << endl;
cout << "Case1: cur_ind = " << cur_ind << ", range_a_ind = " << range_a_ind << ", ends_ranges_of_As.size = " << ends_ranges_of_As.size() << endl;
						space_release();
						exit(0);
					}
				}
			}//if 
			cur_ind++;
		}//for j
		for(_uli j = SIZE - SEED + 2; j < SIZE; j += 2){
			_uli left = prev & mask[SIZE - j];
			_uli right = (next >> (SIZE - (SEED - (SIZE - j)))) & mask[SEED - (SIZE - j)];
			left = ((left << (SEED - (SIZE - j))) | right ) & mask[SEED];
			if(left == 0){

				long int range_a_ind = lower_bound(ends_ranges_of_As.begin(), ends_ranges_of_As.end(), cur_ind) - ends_ranges_of_As.begin();
				_uli alength = (ends_ranges_of_As[range_a_ind] - cur_ind) ; 
				if(ends_ranges_of_As[range_a_ind] == binary_total){
					count_eos++;
					lengths.push_back(alength);
					temp.push_back(cur_ind);
				}//if
				else{
					_uli len_ind = alength >> shift_ht_lengths;
					if(len_ind < small_ht_size)
						ht_sizes_lengths[len_ind]++;//need this to sort by length		
					else{
						cout << "ERROR: alength = " << alength << ", len_ind = " << len_ind << " >= " << small_ht_size << endl;
cout << "Case2: cur_ind = " << cur_ind << ", range_a_ind = " << range_a_ind << ", ends_ranges_of_As.size = " << ends_ranges_of_As.size() << endl;
						space_release();
						exit(0);
					}				
				}
			}
			cur_ind++;
		}//for j
		prev = next;
		i++;
	}//while
	//process last 
	_uli double_cur = cur_ind;

	for(_uli j = (double_cur << 1); j < full_binary_total; j += 2){
		_uli res = get_ta_lmer(j, SEED);
		if(res == 0){

			long int range_a_ind = lower_bound(ends_ranges_of_As.begin(), ends_ranges_of_As.end(), cur_ind) - ends_ranges_of_As.begin();
			_uli alength = (ends_ranges_of_As[range_a_ind] - cur_ind); 
	
			if(ends_ranges_of_As[range_a_ind] == binary_total){
				count_eos++;
				lengths.push_back(alength);
				temp.push_back(cur_ind);
			}//if
			else{
				_uli len_ind = alength >> shift_ht_lengths;
					if(len_ind < small_ht_size)
						ht_sizes_lengths[len_ind]++;//need this to sort by length		
					else{
						cout << "ERROR: alength = " << alength << ", len_ind = " << len_ind << " >= " << small_ht_size << endl;
cout << "Case3: cur_ind = " << cur_ind << ", range_a_ind = " << range_a_ind << ", ends_ranges_of_As.size = " << ends_ranges_of_As.size() << endl;
/*
ERROR: alength = 4293376903, len_ind = 4293376903 >= 16777216
cur_ind = 1590393, range_a_ind = 1238
*/

						space_release();
						exit(0);
					}						
			}
		}//if
		cur_ind++;
	}//for j

	if(count_eos > 0){
		ht_indices = new _uli[count_eos];
		bc_indices.resize(count_eos, 0);		
		for(_uli y = 0 ; y < count_eos; y++){
			ht_indices[y] = temp[y];
			bc_indices[y] = y;
		}//for
	}//if
}//collect_indices_first_eos
////////////////////////// LONG READS
_uli Anchor::get_ta_lmer(_uli cur_ind, _uli asize)//asize is in bits here, for 2-bit genome must be even value
{
	//cur_ind is given by 2-bit representation here
	//asize is given in bits not bases
	_uli zero = 0;
	_uli ta_ind = cur_ind >> two_in_Size;//divide by 64
	if(ta_ind >= full_size){
		return zero;
	}
	_uli prev = full_genome[ta_ind];
	_uli next = (ta_ind + 1 >= full_size ? 0 : full_genome[ta_ind + 1]);
	_uli inside_ind = cur_ind & mask[two_in_Size];//inside 64-bit word
	_uli left_bits = SIZE - inside_ind;
	_uli left = prev & mask[left_bits];
	if(left_bits >= asize){
		left >>= left_bits - asize;
		return (left & mask[asize]);
	}//if
	else{
		_uli right = next >> (SIZE - (asize - left_bits));
		left = (left << (asize - left_bits)) | right;
		return (left & mask[asize]);
	}//else
}

_uli Anchor::get_ns_lmer(_uli cur_ind, _uli asize)
{
	//cur_ind is given in 1-bit representation
	_uli zero = 0;
	_ui ta_ind = cur_ind >> two_in_Size;//divide by 64
	if(ta_ind >= ta_size){
		return zero;
	}
	_uli prev = ns_genome[ta_ind];
	_uli next = (ta_ind + 1 >= ta_size ? 0 : ns_genome[ta_ind + 1]);
	_ui inside_ind = cur_ind & mask[two_in_Size];//inside 64-bit word
	_ui left_bits = SIZE - inside_ind;
	_uli left = prev & mask[left_bits];
	if(left_bits >= asize){
		left >>= left_bits - asize;
		return (left & mask[asize]);
	}//if
	else{
		_uli right = next >> (SIZE - (asize - left_bits));
		left = (left << (asize - left_bits)) | right;
		return (left & mask[asize]);
	}//else
}//get_ns_lmer


void Anchor::sort_first_entry()
{
//DEBUG
cout << "SORT_FIRST_ENTRY: starting, max length of range As = " << max_range_As << endl;
	vector<_uli> count_lengths(max_range_As + 1, 0);
	vector<_uli> ptr_lengths(max_range_As + 1, 0);

	//if there is a range of consecutive As at the end of genome, then sort suffixes of this range first
	_uli count_eos = 0;

	collect_indices_first_eos(count_eos);//ht_indices are filled with count_eos and lengths are filled in, ht_sizes_lengths corrected
	//this collects and sorts indices of suffixes of consecutive As that end with $ (i.e. at the end of genome)
	//In our case, it is the end of Reverse(orig genome), i.e. the start of chromosome 1 that has at the beginning 49 Ns
	
cout << "collect_indices_first_eos is done count_eos = " << count_eos << endl;
	if(count_eos > 0){
		for(_uli i = 0; i < count_eos; i++)
			bc_indices[i] = i;
		_uli_iter st_eos = bc_indices.begin();
		_uli_iter en_eos = bc_indices.begin() + count_eos - 1;
		quick_multikeys(lengths, st_eos, en_eos);
		update_bwt(0);
	}


cout << "sort_first_entry: total EOS suffixes = " << count_eos << endl;

	vector< pair<_ui, _ui> > temp_buckets;
	vector<_ui> temp_sum;
		//NOVA_NEXT ->
	collect_ranges(0, small_ht_size - 1, two_in_block_32, temp_buckets, temp_sum, ht_sizes_lengths);//two_in_block_64 in 14min 268 of 298
	_uli total_ranges = temp_sum.size();
	cout << "total_ranges for first entry = " << total_ranges << endl;

	for(long int i = total_ranges - 1; i >= 0; i--){//i has to be int because we decrement and if i is _ui this loop will run forever
		//initializes ht_indices, lengths, bc_indices

		collect_bucket_indices_first(temp_buckets[i].first, temp_buckets[i].second, temp_sum[i], count_lengths);
cout << "collect_bucket_indices_first completed for i = " << i << ", and total_sum = " << temp_sum[i] << endl;
		//sort in increasing order
/*		_uli_iter st_len = bc_indices.begin();
		_uli_iter en_len = bc_indices.end() - 1;
		quick_multikeys(lengths, st_len, en_len);
*/		
		int size_equal_ranges = 0;
			
		for(_uli y = 1; y <= max_range_As; y++){
			ptr_lengths[y] = ptr_lengths[y -1 ] + count_lengths[y - 1];
			if(count_lengths[y] > 1)
				size_equal_ranges++;
		}//for y
		for(_uli y = 0; y < temp_sum[i]; y++){
			_uli aptr = ptr_lengths[lengths[y]];
			bc_indices[aptr] = y;
			ptr_lengths[lengths[y]] += 1;
		}
		for(_uli y = 0; y <= max_range_As; y++){
			count_lengths[y] = 0;
			ptr_lengths[y] = 0;
		}

		//reverse so that sorting is in decreasing order
		reverse(bc_indices.begin(), bc_indices.end());

		vector<_uli> ranges_lengths(size_equal_ranges);
		//now we need to sort all the equal-length ranges
		vector< pair<_uli_iter, _uli_iter> > equal_ranges(size_equal_ranges);
		_uli_iter prev_ind = bc_indices.begin() ;
		_uli prev_len = lengths[*prev_ind];	

		_uli_iter cur_ind = prev_ind + 1;
		_uli ranges_ind = 0;
		for(; cur_ind < bc_indices.end(); cur_ind++){
			if(lengths[*cur_ind] != prev_len){
				_uli eq_size = cur_ind - prev_ind;
				if(eq_size > 1){
					equal_ranges[ranges_ind] = make_pair(prev_ind, cur_ind - 1);
					ranges_lengths[ranges_ind] = prev_len;
					ranges_ind++;
				}
				prev_ind = cur_ind;
				prev_len = lengths[*cur_ind];
			}//if
		}//for
		if(cur_ind - prev_ind > 1){
			ranges_lengths[ranges_ind] = prev_len;
			equal_ranges[ranges_ind] = make_pair(prev_ind, cur_ind - 1);
		}
		//initialize bucket_components
		bucket_components.resize(temp_sum[i], 0);
		end_of_suffix.resize(temp_sum[i], half_Size);
		prev_same_comp_ids.resize(temp_sum[i], 0);
		_uli sub_bucket_st = 0;

		_uli total_size = equal_ranges.size();
		for(_uli j = 0; j < total_size; j++){
			_ui cur_comp = (ranges_lengths[j] << 1) >> two_in_Size;//ranges_lengths are in bases, each base takes 2 bits
			update_component_first(equal_ranges[j].first, equal_ranges[j].second, cur_comp, sub_bucket_st);
			_ui total_comp = cur_comp + total_components;
			quick_multikeys(cur_comp, equal_ranges[j].first, equal_ranges[j].second, sub_bucket_st, total_comp); 
		}//for j
		update_bwt(0);
//DEBUG
cout << "sort_first_entry: update_bwt i = " << i << " of " << total_ranges << endl;
	}//fori 

}//sort_first_entry


void Anchor::init_comp_indices(_ui st, _ui en, _ui sub_bucket_size)
{
	if(!bucket_components.empty()){
		bucket_components.clear();
	}
	if(!bc_indices.empty())
		bc_indices.clear();
	if(!end_of_suffix.empty())
		end_of_suffix.clear();
	if(!prev_same_comp_ids.empty())
		prev_same_comp_ids.clear();

	bucket_components.resize(sub_bucket_size);
	bc_indices.resize(sub_bucket_size);
	end_of_suffix.resize(sub_bucket_size, half_Size);
	prev_same_comp_ids.resize(sub_bucket_size, 0);//the first component id is 0

	_ui j = 0;
	for(_ui i = st; i <= en ; i++){
		_uli cur_ind = ht_indices[i];//1-bit representation
		_uli bits_before_suffix = half_Size;//max number of bits before end of suffix inside the component
		_uli acomp = get_component(cur_ind, SIZE, 0, bits_before_suffix);//64-bit component
		bucket_components[j] = acomp;
		end_of_suffix[j] = bits_before_suffix;//max number of bits in component preceeding suffix
		bc_indices[j] = j;
		j++;
	}
	
}//init_comp_indices

_uli Anchor::get_component(_uli cur_ind, _uli asize, _uli comp_id, _uli &bits_before_suffix)
{
	//cur_ind is given in 1-bit representation
	cur_ind <<= 1;//*2 for full genome repres
	cur_ind = cur_ind + (comp_id << two_in_Size);
	
	_uli ta_ind = cur_ind >> two_in_Size;//divide by 64
	if(ta_ind >= full_size){
		bits_before_suffix = 0;
		return 0;
	}

	bits_before_suffix = (full_binary_total - cur_ind) >> 1;

	_uli prev = full_genome[ta_ind];
	_uli next = ((ta_ind + 1) >= full_size ? 0 : full_genome[ta_ind + 1]);
	_ui inside_ind = cur_ind & mask[two_in_Size];//inside 64-bit word
	_ui left_bits = SIZE - inside_ind;
	_uli left = prev & mask[left_bits];
	if(left_bits >= asize){
		left >>= left_bits - asize;
		return (left & mask[asize]);
	}//if
	else{
		_uli right = next >> (SIZE - (asize - left_bits));
		left = (left << (asize - left_bits)) | right;
		return (left & mask[asize]);
	}//else

}//get_component

void Anchor::get_pivot_value(_uli_iter st, _uli_iter en, _uli &value, _ui &total_left){

	_uli asize = en - st + 1;
	_uli half_asize = asize >> 1;	

		_uli_iter mid = st + half_asize;
		value = bucket_components[*mid];


}//get_pivot_value

_uli_iter Anchor::partition_multikeys(_uli_iter st, _uli_iter en, _ui &total_left)
{
	//this function is called only if number of values in range [st, en] > 1
	//find pivot value around which we partition the values in range [st, en]
	

	_uli asize = en - st + 1;
	_uli half_asize = asize >> 1;	

	_uli_iter mid = st + half_asize;
	_uli value = bucket_components[*mid];

/* Values in the range [st, en] inclusive are subdivided into two
   halfs of approximately same size, with partitioning values chosen
   to be maximum out of all values in the first half
*/
		_uli_iter orig_st = st;
	
		//first partition into values <= value and values > value
		while(st < en){
			if(bucket_components[*st] > value ){
			  if(bucket_components[*en] <= value){
					swap(*st, *en);
					st++;
					en--;
			  }
			  else
				en--;
			}
			else
			  st++;
        }//while
		if(bucket_components[*st] > value)
			st--;//last pos where elements <= value
		total_left = st - orig_st + 1;
		st = orig_st;


		_uli_iter ahalf = st + (total_left - 1);
		_uli_iter eq = ahalf;

	 //partition the first half into values < value and values = value
		while(st <= eq){
		   if(bucket_components[*st] == value){
			if(bucket_components[*eq] == value){
 			  if(eq > orig_st)
					eq--;
			  else
				st++;
			}//if
			else{
			  swap(*eq, *st);
			  st++;
			  eq--;
			}//else
		   }//if st = val
		   else
			 st++;
		}//while
		if(bucket_components[*eq] == bucket_components[*ahalf]){//CHECK this condition
			return eq;
		}
		else{
			return eq + 1;
		}
}//partition_multikeys

void Anchor::quick_multikeys(_ui component, _uli_iter start, 
	_uli_iter fin, const _ui &sub_bucket_st, _ui total_comp)
{
	if(start >= fin)
	 return;
	if(component >= total_comp){
	  not_sorted.push(make_pair(make_pair(start, fin), make_pair(component, sub_bucket_st)));
	  return;
	}


		//total_left is the number of values in [st, en] whose components <= pivot value
		_ui total_left = 0;
		//left points to the start of range in which the components are equal, the end of this range is start+ total_left -1
		_uli_iter left = partition_multikeys(start, fin, total_left);
		//after partition if flag_sorted
		quick_multikeys(component, start, left - 1, sub_bucket_st, total_components);
		update_end_of_suffix(left, start + (total_left - 1));//left moves to first comp that has no EOS, $, char
		_ui equal_comp = component + 1;
		update_component(left, start + (total_left - 1), equal_comp, sub_bucket_st);
		_ui equal_total_comp = total_components ;//NOVA on h4 server this addition caused core dump Seg fault (+ equal_comp);
		quick_multikeys(equal_comp,  left, start + (total_left - 1), sub_bucket_st, equal_total_comp);
		quick_multikeys(component,  start + total_left, fin, sub_bucket_st, total_components);
	
}//quick_multikeys

void Anchor::quick_multikeys(vector<_uli> &lengths, _uli_iter start, _uli_iter fin){
	if(start >= fin)
	 return;

		//total_left is the number of values in [st, en] whose components <= pivot value
		_ui total_left = 0;
		//left points to the start of range in which the components are equal, the end of this range is start+ total_left -1
		_uli_iter left = partition_multikeys(lengths, start, fin, total_left);
		//after partition if flag_sorted
		quick_multikeys(lengths, start, left - 1);
		quick_multikeys(lengths,  start + total_left, fin);
}//quick_multikeys lengths

_uli_iter Anchor::partition_multikeys(vector<_uli> &lengths, _uli_iter st, _uli_iter en, _ui &total_left)
{
	_uli asize = en - st + 1;
	_uli half_asize = asize >> 1;	

	_uli_iter mid = st + half_asize;
	_uli value = lengths[*mid];

/* Values in the range [st, en] inclusive are subdivided into two
   halfs of approximately same size, with partitioning values chosen
   to be maximum out of all values in the first half
*/
		_uli_iter orig_st = st;
	
		//first partition into values <= value and values > value
		while(st < en){
			if(lengths[*st] > value ){
			  if(lengths[*en] <= value){
					swap(*st, *en);
					st++;
					en--;
			  }
			  else
				en--;
			}
			else
			  st++;
        }//while
		if(lengths[*st] > value)
			st--;//last pos where elements <= value
		total_left = st - orig_st + 1;
		st = orig_st;


		_uli_iter ahalf = st + (total_left - 1);
		_uli_iter eq = ahalf;

	 //partition the first half into values < value and values = value
		while(st <= eq){
		   if(lengths[*st] == value){
			if(lengths[*eq] == value){
 			  if(eq > orig_st)
					eq--;
			  else
				st++;
			}//if
			else{
			  swap(*eq, *st);
			  st++;
			  eq--;
			}//else
		   }//if st = val
		   else
			 st++;
		}//while
		if(lengths[*eq] == lengths[*ahalf]){//CHECK this condition
			return eq;
		}
		else{
			return eq + 1;
		}
}//partition_multikeys lengths

void Anchor::update_end_of_suffix(_uli_iter &st, _uli_iter en){
	/***********************************************************
	components in this range [st, en] are equivalent, since 
	the end of suffix $ is represented as 0, some of components containing
	$ need to be sorted. Components that have $ are done after this function
	***********************************************************/
	if(st >= en)
		return;
	_ui asize = en - st + 1;
	_uli_iter orig_st = st;
	_uli_iter orig_en = en;
	//partition values into those having end of suffix and those that do not
	_ui count_suffix = 0;
	//---------------------------
	while(st <= en){
		if(end_of_suffix[*st] < half_Size)
			count_suffix++;
		st++;
	}//while
	st = orig_st;
	if(count_suffix > 1){
		//sort by number of bits before end of suffix in the bucket	
		vector< pair<_uli, _uli_iter> > dummy(count_suffix);
		int rest_size = en - st + 1 - count_suffix;
		vector< _uli_iter > atemp(rest_size);
		_ui j = 0;
		_ui y = 0;
		for(; st <= en; st++){
			if(end_of_suffix[*st] < half_Size){
				dummy[j].first = end_of_suffix[*st];
				dummy[j].second = st;
				j++;
			}
			else{
				atemp[y] = st;
				y++;
			}
		}//for
		sort(dummy.begin(), dummy.end());
		vector<_uli> dum(count_suffix, 0);
		vector< pair<_uli, _uli_iter> >::iterator st2 = dummy.begin();
		vector< pair<_uli, _uli_iter> >::iterator en2 = dummy.end();
		j = 0;
		for(; st2 < en2; st2++){
			dum[j] = *(st2->second);
			j++;
		}//for

		vector<_uli> dum2(rest_size, 0);
		vector< _uli_iter >::iterator st3 = atemp.begin();
		vector< _uli_iter >::iterator en3 = atemp.end();
		j = 0;
		for(; st3 < en3; st3++){
			dum2[j] = *(*st3);
			j++;
		}//for


		st = orig_st;
		j = 0;
		_uli_iter stop = st + count_suffix;
		for(; st < stop; st++){
			*st = dum[j];
			j++;
		}
		st = orig_st + count_suffix;
		j = 0;
		for(; st <= en; st++){
			*st = dum2[j];
			j++;
		}
		st = orig_st + count_suffix;
	}//if there are end of suffix
	else if(count_suffix == 1){
		st = orig_st;
		for(; st <= en; st++){
			if(end_of_suffix[*st] < half_Size)
				break;
		}
		//st points to EOS
		for(; st > orig_st; st--){
			_uli_iter prev = st - 1;
			swap(*prev, *st);
		}
		st = orig_st + 1;
	}
	else
		st = orig_st;
}//update_end_of_suffix

void Anchor::update_component(_uli_iter st, _uli_iter en, _ui &component, const _ui &sub_bucket_st)
{
	if(st >= en)//already sorted
		return;
	//for each in the range find component
	long int asize = en - st + 1;

	_uli_iter orig_st = st;
	//update prev_same_comp_ids
	for(; st <= en; st++){
		if(prev_same_comp_ids[*st] <= component){
			_uli cur_comp = component;
			_uli bucket_ind = *st;
			_uli ht_ind = sub_bucket_st + bucket_ind;
			_uli genome_ind = ht_indices[ht_ind];
			while(true){
				_uli bits_before_suffix = half_Size;
				//genome_ind is passed as 1-bit repres correctly
				if((genome_ind << 1) + (cur_comp << two_in_Size) >= full_binary_total){
					prev_same_comp_ids[*st] = cur_comp;
					break;
				}
				_uli cur_bucket_comp = get_component(genome_ind, SIZE, cur_comp, bits_before_suffix);
				if((bucket_components[*st] != cur_bucket_comp) || (bits_before_suffix < half_Size)){
					prev_same_comp_ids[*st] = cur_comp;
					break;
				}//if
				cur_comp++;
			
			}//while
		}//if need to change comp id

	}//for
	st = orig_st;
	_ui amin_comp = prev_same_comp_ids[*st];
	st++;
	for(; st <= en ; st++){
		if(amin_comp > prev_same_comp_ids[*st])
			amin_comp = prev_same_comp_ids[*st];
	}
	component = amin_comp;
	st = orig_st;
	for(; st <= en; st++){
		if(prev_same_comp_ids[*st] <= component){
			_uli bucket_ind = *st;
			_uli ht_ind = sub_bucket_st + bucket_ind;
			_uli genome_ind = ht_indices[ht_ind];
			_uli bits_before_suffix = half_Size;
			bucket_components[*st] = get_component(genome_ind, SIZE, component, bits_before_suffix);
			end_of_suffix[bucket_ind] = bits_before_suffix;
		}
	}//for
}//update_component

void Anchor::update_component_first(_uli_iter st, _uli_iter en, _ui component, const _ui &sub_bucket_st)
{
	if(st >= en)//already sorted
		return;

	for(; st <= en; st++){
		_uli bucket_ind = *st;
		_uli ht_ind = sub_bucket_st + bucket_ind;
		_uli genome_ind = ht_indices[ht_ind];
		_uli bits_before_suffix = half_Size;
		bucket_components[*st] = get_component(genome_ind, SIZE, component, bits_before_suffix);
		end_of_suffix[bucket_ind] = bits_before_suffix;
		prev_same_comp_ids[bucket_ind] = component;
	}//for
}//update_component_first

int Anchor::parse_options(int argc, char* argv[])
{ //Genome option -G is depricated
	int res = 0;

	long int i;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-')
			return -1;
		switch(argv[i][1]){
		case 'r': genome_names = argv[++i]; break;
		case 'P': prefix = argv[++i]; break;
		case 'S': STEP = atoi(argv[++i]); break;
		case 'h': usage(); break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl; usage(); space_release(); exit(0); break;
		}
	}//for i
	return res;
}//parse_opt()

void Anchor::flash_genome(ofstream &out, _uli *agenome, _uli asize){

	char abuffer[8];

	for(_uli i = 0; i < asize; i++){
		_uli cur = agenome[i];
		uli_to_char(cur, abuffer);
		out.write(&abuffer[0], 8);
	}

}//flash_genome

void Anchor::flash_bwt()
{
	char abuffer[8];

	for(_uli i = 0; i < cur_ind_bwt; i++){
		_uli cur = bwt[i];
		uli_to_char(cur, abuffer);
		out_bwt.write(&abuffer[0], 8);
	}
}//flash_bwt


void Anchor::flash_marked_ns()
{
	char abuffer[4];

	for(_uli i = 0; i < cur_ind_bwt; i++){
		_ui cur = marked_ns[i];
		ui_to_char(cur, abuffer);
		out_ns_genome.write(&abuffer[0], 4);
	}	
	for(_uli i = 0; i < cur_ind_bwt; i++){
		_ui cur = bwt_marked[i];
		ui_to_char(cur, abuffer);
		out_bwt_marked.write(&abuffer[0], 4);
	}	
	//nova: output pos_strand such that pos_strand[i] = 1 (marked like this only when bwt_marked[i] = 1)
	//then it means that corresponding genomic position is in the second half of concatenated two genome versions
	//During mapping, we need to add orig_genome_size for such genomic positions to know real value
	for(_uli i = 0; i < cur_ind_bwt; i++){
		_ui cur = pos_strand[i];
		ui_to_char(cur, abuffer);
		out_pos_strand.write(&abuffer[0], 4);
	}//for
}//flash_marked

void Anchor::flash_pos_indices(_uli stop)
{
	//count the size of this bucket
	char abuffer2[4];
	//output keys into out_pos_keys and indices into out_pos_indices
	for(_uli i = 0; i < stop; i++){
			_ui cur_ind = pos_indices[i];
			ui_to_char(cur_ind, abuffer2);
			out_pos_indices.write(&abuffer2[0], 4);
	}//fori
	pos_indices_ind = 0;
}//flash_pos_indices()

void Anchor::flash_num_occ(_uli upto)
{
	char abuffer[2];
	char abuffer2[4];

	ui_to_char(num_occ_sqr[0], abuffer2);
	out_num_sqr.write(&abuffer2[0], 4);
	ui_to_char(num_occ_sqr[1], abuffer2);
	out_num_sqr.write(&abuffer2[0], 4);
	ui_to_char(num_occ_sqr[2], abuffer2);
	out_num_sqr.write(&abuffer2[0], 4);
	ui_to_char(num_occ_sqr[3], abuffer2);
	out_num_sqr.write(&abuffer2[0], 4);

	for(_uli i = 0; i < upto; i++){
		_usi cur = num_occ[i][0];
		usi_to_char(cur, abuffer);
		out_num_occ.write(&abuffer[0], 2);

		//Ts
		cur = num_occ[i][1];
		usi_to_char(cur, abuffer);
		out_num_occ.write(&abuffer[0], 2);

		cur = num_occ[i][2];
		usi_to_char(cur, abuffer);
		out_num_occ.write(&abuffer[0], 2);

		cur = num_occ[i][3];
		usi_to_char(cur, abuffer);
		out_num_occ.write(&abuffer[0], 2);
	}//for i

	//genome positions number occ:
	ui_to_char(pos_num_sqr, abuffer2);
	out_pos_num_sqr.write(&abuffer2[0], 4);
	for(_uli i = 0; i < upto; i++){
		_usi cur = pos_num_occ[i];
		usi_to_char(cur, abuffer);
		out_pos_num_occ.write(&abuffer[0], 2);
	}
}//flash_num_occ

void Anchor::usi_to_char(_usi cur, char abuffer[])
{
	_usi x = (cur >> abyte) & MASK_BYTE;
	abuffer[0] = char(x);
	x = cur & MASK_BYTE;
	abuffer[1] = char(x);
	
}//usi_to_char
void Anchor::uli_to_char(_uli cur, char abuffer2[])
{
		unsigned int x = (cur >> seven_bytes) & MASK_BYTE;
		abuffer2[0] = char(x);
		x = (cur >> six_bytes) & MASK_BYTE;
		abuffer2[1] = char(x);
		x = (cur >> five_bytes) & MASK_BYTE;
		abuffer2[2] = char(x);
		x = (cur >> four_bytes) & MASK_BYTE;
		abuffer2[3] = char(x);
		x = (cur >> three_bytes) & MASK_BYTE;
		abuffer2[4] = char(x);
		x = (cur >> two_bytes) & MASK_BYTE;
		abuffer2[5] = char(x);
		x = (cur >> abyte) & MASK_BYTE;
		abuffer2[6] = char(x);
		x = cur & MASK_BYTE;
		abuffer2[7] = char(x);
}//uli_to_char

void Anchor::ui_to_char(_ui cur, char abuffer2[])
{
		unsigned int x = (cur >> three_bytes) & MASK_BYTE;
		abuffer2[0] = char(x);
		x = (cur >> two_bytes) & MASK_BYTE;
		abuffer2[1] = char(x);
		x = (cur >> abyte) & MASK_BYTE;
		abuffer2[2] = char(x);
		x = cur & MASK_BYTE;
		abuffer2[3] = char(x);
}//ui_to_char

void Anchor::output_info(){

	out_chrom << binary_total << endl;
	out_chrom << ta_size << endl;
	out_chrom << full_binary_total << endl;
	out_chrom << full_size << endl;
	//full_binary total calculate from binary total and same full_size

	out_chrom << num_chroms << endl;


	out_chrom << C[0] << " " << C[1] << " " << C[2] << " " << C[3] <<  endl;
	out_chrom << orig_gen_inL << endl;//we need this number for all genomes
	out_chrom << STEP << endl;
	out_chrom << pos_indices_total << endl;

	map<unsigned int, string>::iterator beg = chroms_names.begin();
	map<unsigned int, string>::iterator fin = chroms_names.end();

	//index and name of chroms
	for(; beg != fin; beg++){
		out_chrom << beg->first << "\t" << beg->second << endl;
	}

	//files with references
	_uli orig_num_chroms = num_chroms >> 1;
	for(_uli i = 0; i < orig_num_chroms; i++)
		out_chrom << names[i] << endl;
	//sizes of chroms in bits: each bit per char
	for(_ui i = 0; i < num_chroms; i++){
		out_chrom << size_chrom[i] << endl;
	}
	//starts of chroms in bits
	for(_ui i = 0; i < num_chroms; i++)
		out_chrom << starts_size_chrom[i] << endl;

	out_chrom << orig_genome_size << endl;//nova

}//output_info

void Anchor::precalculate_seeds()
{
	out_bwt.close();
	out_num_sqr.close();
	out_num_occ.close();

	delete [] ns_genome;
	ns_genome = 0;
	delete [] full_genome;
	full_genome = 0;

	if(ht_indices != 0){
		delete [] ht_indices;
		ht_indices = 0;
	}

	vector< vector<_ui> > num_occ_sqr;
	read_num_occ_sqr(num_occ_sqr);
	num_occ.clear();
	read_num_occ(num_occ);
	cout << "ta_size = " << ta_size << ", binary total =  "<< binary_total 
		<< ", C: " << endl;
	for(int y = 0; y < 4; y++)
		cout << C[y] << endl;

	vector<_uli> temp(5, 0);
	temp[4] = C[0] + C[1] + C[2] + C[3] + 1;
	temp[3] = C[0] + C[1] + C[2] + 1;
	temp[2] = C[0] + C[1] + 1;
	temp[1] = C[0] + 1;
	temp[0] = 1;
	C.clear();
	C.resize(5);
	C[0] = temp[0];
	C[1] = temp[1];
	C[2] = temp[2];
	C[3] = temp[3];
	C[4] = temp[4];

	orig_num_sqr = orig_gen_inL >> two_in_block_64;
	orig_num = orig_gen_inL >> two_in_Size;

	string genome_repres = "CT";

	string rem = (prefix[prefix.length() - 1] == '/' ? "" : "/");

	string input_bwt = prefix + rem + genome_repres + ".bwt";

	in_bwt.open(input_bwt.c_str(), ios::in);
	read_bwt_genome();
	cout << "bwt is read" << endl;
	for(_uli i = 0; i < hash_table_size; i++){
		_uli sp = 1;
		_uli ep = 0;
		exact_match(i, sp, ep, num_occ_sqr);
		char buffer[8];
		uli_to_char(sp, buffer);
		out_seeds.write(&buffer[0], 8);
		uli_to_char(ep, buffer);//sp and ep is count of char X in bwt from [0...sp] 1-bit bwt
		out_seeds.write(&buffer[0], 8);
//cout << "seed i = " << i << endl;
	}//for i

	in_bwt.close();
	out_seeds.close();

	
}//precalculate_seeds

void Anchor::exact_match(_uli ta_read, _uli &sp, _uli &ep, vector< vector<_ui> > &num_occ_sqr)//sp=1, ep = 0, remainer = 0: init
{
	_uli ch = ta_read & mask[2];
	ta_read >>= 2;
	int i = SEED - 2;//has to be unsigned
	sp = C[ch];
	ep = C[ch + 1] - 1;
	while( (sp <= ep) && (i > 0)){
		ch = ta_read & mask[2];
		ta_read >>= 2;
		sp = C[ch] + Occ(ch, sp -1, num_occ_sqr);
		ep = C[ch] + Occ(ch, ep, num_occ_sqr) - 1;
		i -= 2;
	}//while

}//exact_match

_uli Anchor::Occ(_uli ch, _uli sp, vector< vector<_ui> > &num_occ_sqr){

	_uli block_ind_sqr = sp >> two_in_block_64;
	_uli asize = num_occ_sqr.size();
	if(block_ind_sqr >= asize)
	{
	cout << "ERROR block_ind_sqr in Occ = " << block_ind_sqr << " > num_occ_sqr.size = " << asize
		<< ", sp = " << sp << endl;
		space_release();
	exit(0);
	}
	_uli num_a_sqr = (block_ind_sqr == 0 ? 0 : num_occ_sqr[block_ind_sqr - 1][ch]);

	_uli bl_ind = sp >> two_in_Size;
    _uli bl_inside_ind =  bl_ind & mask[two_in_Size];
	if(bl_ind >= num_occ.size())
	{
	cout << "ERROR bl_ind in Occ = " << bl_ind << " > num_occ.size "  << endl;
			space_release();
	exit(0);
	}
	_uli num_prev_bl = (bl_inside_ind == 0  ? 0 : num_occ[bl_ind - 1][ch]);
	
	_uli full_ind = (sp << 1) >> two_in_Size;
	_uli full_bits = (sp << 1) & mask[two_in_Size];
	_uli num_cur = get_bwt_lmer(full_ind, full_bits + 2, ch);

	_uli sum_bl_ind = bl_ind << two_in_Size;//not same as sp because of remainder
	_uli full_sum_bl_ind = (sum_bl_ind << 1) + SIZE;
	if((sp << 1) >= full_sum_bl_ind){
		num_cur += get_bwt_lmer(full_ind - 1, SIZE, ch);
	}

	//this is true only for A character, becaus EOS($) char also was represented as A
	//and because this was not counted
	if((ch == 0) && (orig_num == bl_ind) && (orig_gen_inL <= sp))
		num_cur--;

	num_a_sqr += num_prev_bl + num_cur;	

	return num_a_sqr;
}//Occ

_uli Anchor::get_bwt_lmer(_uli cur_byte, _uli aread_len, _uli ch){
	//cur_byte is index of 64-word in full genome and aread_len is in bits
        _uli count_ones = 0;	
	if(cur_byte >= full_size){
	cout << "ERROR: get_bwt_lmer cur_byte is >= full_size " << cur_byte << endl;
			space_release();
	exit(0);
		return count_ones;

	}
        _uli bwt_lmer = bwt_genome[cur_byte] >> (SIZE - aread_len);
        for(int i = 0; i < 4; i++){
                count_ones += num_occ_bytes[bwt_lmer & mask[16]][ch];//64-bit word is subdivided into 16bit sub-words
                bwt_lmer >>= 16;
        }
		if(ch == 0)
			count_ones -= ((SIZE - aread_len) >> 1);//filled in zeros at init bwt_lmer, divide by two = number of A/N/G
        return count_ones;
}//get_bwt_lmer

_uli Anchor::count_mism(unsigned long int bits){

/////////////////Andrew Smith author of the following code:
  bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
  bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
  bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
  bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
  bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
  return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);
/////////////////////END Andrew Smith

}//count_mism()

void Anchor::read_bwt_genome()
{
	long int bwt_size = full_size + 1 ;
	bwt_genome = new _uli[bwt_size];
	assert(bwt_genome != 0);

	const _uli buffer_size = block_256 << two_in_byte;//*8
	char *buffer = new char[buffer_size];
	assert(buffer != 0);

	in_bwt.read(buffer, buffer_size);
	_uli real_size = in_bwt.gcount();

	long int bwt_ind = 0;

	while(true){
	  _uli buf_ind = 0;
	  while(buf_ind < real_size){
		_uli res = 0;
		for(int y = 0; y < 8 ; y++){
			res <<= byte_size;
			res |= int(buffer[buf_ind]) & MASK_BYTE;
			buf_ind++;//buffer size is multiple of 8
		}//for
		bwt_genome[bwt_ind] = res;
		bwt_ind++;
	  }//while buffer size
	  if(in_bwt.eof())
		break;
	  in_bwt.read(buffer, buffer_size);
	  real_size = in_bwt.gcount();
	  if(real_size == 0)
		  break;
	}//while true
	delete [] buffer;
	cout << "finished reading bwt" << endl;
}//read_bwt_genome

void Anchor::read_num_occ_sqr(vector< vector<_ui> > &num_occ_sqr)
{
	ifstream in_num_sqr;
	string rem = (prefix[prefix.length() - 1] == '/' ? "" : "/");
	string genome_repres = "CT";

	string input_num_sqr = prefix + rem + genome_repres + ".nosqr";
	in_num_sqr.open(input_num_sqr.c_str(), ios::in);

	//this is correct because num_occ are calculated in character counts, not bit counts
	_uli full_binary = ta_size  << two_in_Size;//number of bits in 64-word size genome
	_uli total_occ = full_binary >> two_in_block_64;
	if((total_occ & mask[two_in_block_64]) > 0)
		total_occ++;

	const _uli buffer_size = total_occ << 4; //*4 for _ui *4 for four char
	char *buffer = new char[buffer_size];
	in_num_sqr.read(buffer, buffer_size);
	_uli real_size = static_cast<_uli>(in_num_sqr.gcount());

	cout << "total bytes for num_occ_sqr = " << real_size << endl;
	num_occ_sqr.clear();
	vector<_ui> dum(4, 0);
	num_occ_sqr.resize(total_occ, dum);
	_uli buffer_ind = 0;
	for(_uli i = 0; i < total_occ && buffer_ind < real_size; i++){
		_uli res = 0;
		for(_uli y = 0; y < 4; y++){//4 is bytes in unsigned int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		num_occ_sqr[i][0] = res;
		res = 0;
		for(_uli y = 0; y < 4; y++){//4 is bytes in unsigned int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		num_occ_sqr[i][1] = res;
		res = 0;
		for(_uli y = 0; y < 4; y++){//4 is bytes in unsigned int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		num_occ_sqr[i][2] = res;
		res = 0;
		for(_uli y = 0; y < 4; y++){//4 is bytes in unsigned int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		num_occ_sqr[i][3] = res;
	}//for i

	delete [] buffer;
	in_num_sqr.close();
}//read_num_occ_sqr

void Anchor::read_num_occ(vector< vector<_usi> > &num_occ)
{
	string rem = (prefix[prefix.length() - 1] == '/' ? "" : "/");
	string genome_repres ;
	genome_repres = "CT";
	
	string input_num_occ = prefix + rem + genome_repres + ".no";
	ifstream  in_num_occ;
	in_num_occ.open(input_num_occ.c_str(), ios::in);

	_uli total_occ = ta_size ;

	const _uli buffer_size = total_occ << 3; //*2 * 4 size of unsigned short int for four char
	char *buffer = new char[buffer_size];
	in_num_occ.read(buffer, buffer_size);
	_uli real_size = in_num_occ.gcount();
		cout << "total bytes for num_occ = " << real_size << endl;
	num_occ.clear();
	vector<_usi> dum(4, 0);
	num_occ.resize(total_occ, dum);
	_uli buffer_ind = 0;
	for(_uli i = 0; i < total_occ && buffer_ind < real_size; i++){
		_usi res = 0;
		for(_uli y = 0; y < 2; y++){//2 is bytes in unsigned short int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		num_occ[i][0] = res;
		res = 0;
		for(_uli y = 0; y < 2; y++){//2 is bytes in unsigned short int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		num_occ[i][1] = res;
		res = 0;
		for(_uli y = 0; y < 2; y++){//2 is bytes in unsigned short int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		num_occ[i][2] = res;
		res = 0;
		for(_uli y = 0; y < 2; y++){//2 is bytes in unsigned short int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		num_occ[i][3] = res;
	}//for i

	delete [] buffer;
	in_num_occ.close();
}//read_num_occ


int main(int argc, char* argv[])
{
	/**************************************
	This program constructs BWT structures
  *******************************************/
	
	const int SIZE = 64;

	if(argc < 3){
		usage();
		return 0;
	}//if	

	Anchor anc(argc, argv);//reads genome builds ta-, cg- binary representations once
	anc.build_bwt();

	return 0;
}//main


