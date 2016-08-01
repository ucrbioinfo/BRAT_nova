//  BRAT is a whole genome bisulfite sequence mapping program
//	New logic of alignment is used: collect ranges (sp < ep), then collect valid genomic positions, sort them and align full-length only to distinct genomic positions: speeded up 1.25 times
//	Alignment of short and long reads is combined (the same functions are used)
//	Local alignment is supported
//	Supports indels
//	Optimized map_entire_length: MD are calculated only for the winner, if(read-len - mism)> entire_score, then update_singles is called 
//	(mism_threshold is set in init_params to read-len - entire_score)
//	Local alignment score: (mappable length - mismatches inside mappable length) + (mism * mism_penalty)
//	Indel alignment score: (mappable length - mismatches inside mappable length) + (mism * mism_penalty) + gap_open_penalt + gap_penalty*gap_length
//	MD and CIGAR strings are calculated according to SAM format (except for deletion in genome, instead of genome sequence "AAAA..." of length = gap_length is reported
//	update_singles now compares scores and not mismatches
//	NEW OPTIONS: -m option now is for mism_penalty not for the number of mismatches
//	mism_threshold = read_len - entire_score is calculated inside init_params for each read
//	NEW score-system is complete
//	Paired-end alignment is supported (print_pairs is also finished, print_singles is complete), calculation of field FLAG for SAM format is done
//	SAM format is supported, Input is in FASTQ format with st and en placed inside read_name after ? characters
//  Copyright (C) 2015  Elena Yavorska Harris
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

/************************************************
Bugs found:
1. Function map_reads in pairend, istr1 is used instead of istr2 to read in read_name2 (no read name for mate 2 in the output)
2. In print_pairs (or pair_read) read_name is passed as an argument instead of read_name2
3. In update_singles update did not occur for the same position: changed for map_local and map_indels to work correctly:
positions could be equal for local and indel alignment, but indel is done later on, and might have a better score

****************************************************/

#include<iomanip>
#include<map>
#include<bitset>
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

//lookup8.c, by Bob Jenkins, January 4 1997, Public Domain.
#define mix64(a,b) \
{ \
  a -= b;  a ^= (b>>43); \
  b -= a; b ^= (a<<9); \
  a -= b; a ^= (b>>38); \
  b -= a; b ^= (a<<23); \
  a -= b; a ^= (b>>35); \
  b -= a; b ^= (a<<49); \
  a -= b; a ^= (b>>12); \
  b -= a; b ^= (a<<18); \
}

/*
#define mix64(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>43); \
  b -= c; b -= a; b ^= (a<<9); \
  c -= a; c -= b; c ^= (b>>8); \
  a -= b; a -= c; a ^= (c>>38); \
  b -= c; b -= a; b ^= (a<<23); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>35); \
  b -= c; b -= a; b ^= (a<<49); \
  c -= a; c -= b; c ^= (b>>11); \
  a -= b; a -= c; a ^= (c>>12); \
  b -= c; b -= a; b ^= (a<<18); \
  c -= a; c -= b; c ^= (b>>22); \
}
*/
//End lookup8.c, by Bob Jenkins, January 4 1997, Public Domain.


class Stat{
public:


	long int mapped_pairs;
	long int ambiguous;
	long int unmapped;
	long int single_mapped;
	long int single_ambiguous;
	long int broken_dif_chroms;
	long int broken_wrong_insert;
	long int broken_inverse;
	long int broken_mate_unmapped;
	long int broken_mate_amb; //if pair is unmapped and mate is ambiguous, means something else is wrong
	//dif chroms, inverse or wrong insert, otherwise it would be counted as amb pair
	long int broken_both_unmapped;
	long int invalid;
	long int mapped_first_Arich;
	Stat();
};
Stat::Stat(){

	mapped_pairs = 0; 
	ambiguous = 0; 
	unmapped = 0; 
	single_mapped = 0;
	single_ambiguous = 0;
	 broken_dif_chroms = 0;
	 broken_wrong_insert = 0;
	 broken_inverse = 0;
	 broken_mate_unmapped = 0;
	 broken_mate_amb = 0; //if pair is unmapped and mate is ambiguous, means something else is wrong
	//dif chroms, inverse or wrong insert, otherwise it would be counted as amb pair
	 broken_both_unmapped = 0;
	invalid = 0;
	mapped_first_Arich = 0;
}
unsigned int MAX_MISM = 1000000;

class Mates{
public:
	Mates(unsigned int x, unsigned long int y,  short int mark) : cigar1(""), MD1(""), chrom_id(x), pos1(y),  mark_unique(mark), mism1(MAX_MISM),  A_rich(false), clip_left1(0), clip_right1(0), local_score1(0) {}; 
	Mates(): flag(-1), read_len(0), chrom_id(-1), pos1(0), mark_unique(0), mism1(MAX_MISM),   A_rich(false), clip_left1(0), clip_right1(0), MD1(""), cigar1(""), local_score1(0) {};
	void operator =(Mates a){cigar1 = a.cigar1; type_alignment = a.type_alignment; read_len = a.read_len; flag = a.flag; A_rich = a.A_rich; mism1 = a.mism1; chrom_id = a.chrom_id; pos1 = a.pos1 ;  mark_unique = a.mark_unique; clip_left1 = a.clip_left1; clip_right1 = a.clip_right1; MD1 = a.MD1; local_score1 = a.local_score1; };
	unsigned int chrom_id;
	unsigned long int pos1;
	short int mark_unique;
	unsigned int mism1;
	bool A_rich;
	_ui clip_left1;
	_ui clip_right1;
	string MD1;
	string cigar1;
	int local_score1;
	_usi type_alignment; //0 for entire-length, 1 for local, 2 for indels
	_ui read_len;//for indels this is adjusted read length
	int flag;//this is for paired end alignment mostly; for singles values: 0 if mapped to +, 16 if to -
};

class Triple{
public:
	Triple(_uli x, _uli y, _ui z, _ui sh) : sp(x), ep(y), remainder(z), shift(sh) {};
	Triple() : sp(0), ep(0), remainder(0), shift(0) {} ;
	_uli sp;
	_uli ep;
	_ui remainder;
	_ui shift;
};

class AlignInfo{
public:
	AlignInfo(_uli gen_pos, _ui rem, _ui shift, _ui sc, string md) : pos(gen_pos), remainder(rem), shift_bases(shift), score(sc), MD(md), clip_left(0), clip_right(0) {};
	AlignInfo(_uli gen_pos) : pos(gen_pos), remainder(MAX_MISM), shift_bases(MAX_MISM), score(0), MD(""), clip_left(0), clip_right(0) {};
	AlignInfo() : pos(0), remainder(MAX_MISM), shift_bases(MAX_MISM), score(0), MD(""), clip_left(0), clip_right(0) {};
	AlignInfo(_uli gen_pos, _ui rem, _ui shift) : score(0),  pos(gen_pos), remainder(rem), shift_bases(shift) {};
	_uli pos;
	_ui remainder;//from the start of alignment bases not aligned
	_ui shift_bases;//from the end of alignment bases not aligned
	int score;//total number of mapped bases - mism (for local alignment it is based on max_sub_array's score of MD string that identifies best-score subarray, then score = length_sub_array - mism inside it)
	string MD;
	string cigar;
	_ui clip_left;
	_ui clip_right;
};

long int Reverse(unsigned long int i, int bits){//by Nicholas Allen

i = ((i >> 1) & 0x5555555555555555) | ((i & 0x5555555555555555) << 1);
i = ((i >> 2) & 0x3333333333333333) | ((i & 0x3333333333333333) << 2);
i = ((i >> 4) & 0x0F0F0F0F0F0F0F0F) | ((i & 0x0F0F0F0F0F0F0F0F) << 4);
i = ((i >> 8) & 0x00FF00FF00FF00FF) | ((i & 0x00FF00FF00FF00FF) << 8);
i = ((i >> 16) & 0x0000FFFF0000FFFF) | ((i & 0x0000FFFF0000FFFF) << 16);
i = ( i >> 32 ) | ( i << 32);
i = i >> (64 - bits);

return i;
}
void convert_lmer_ta_long(string lmer, vector<unsigned long int> &ta_read, int asize){
	long int j, i;
	j = 0;
	while( j < asize){
		unsigned long int gb = 0; 
		for(i = 0; i < 64 && j < asize; i++){
			gb <<= 1;
			if(lmer[j] == 'T' || lmer[j] == 't' ||  lmer[j] == 'C' || lmer[j] == 'c')
				gb |= 1;
			j++;
		}//for i
		
		ta_read.push_back(gb);

	}//while j
}//convert_ta

void convert_lmer_c_long(string current_chrom, vector<unsigned long int> &cg_read,
						 int num_iter)
{//this function is used only once when there is a beginning of a chromosome
	//aread is an lmer in a genome
	//gb is binary representative of this genome lmer
	long int j =0, i = 0;

	while(i < num_iter){
		unsigned long int gb = 0;

		for(j = 0; j < 64 && i < num_iter; j++){
			gb <<= 1;
			if(current_chrom[i] == 'C' || current_chrom[i] == 'c' ){
				gb |= 1;
			}
			else if(current_chrom[i] == 'G' || current_chrom[i] == 'g'){
				gb |= 1;
			}
			i++;
		}//for j
		cg_read.push_back(gb);
	}//while

}//connvert_lmer_c()

void convert_lmer_c(string current_chrom, unsigned long int &gb, int num_iter)
{//this function is used only once when there is a beginning of a chromosome
	//aread is an lmer in a genome
	//gb is binary representative of this genome lmer
	long int j ;

	gb = 0;
		for(j = 0; j < num_iter; j++){
			gb <<= 1;
			if(current_chrom[j] == 'C' || current_chrom[j] == 'c' ){
				gb |= 1;
			}
			else if(current_chrom[j] == 'G' || current_chrom[j] == 'g'){
				gb |= 1;
			}

		}//for j

}//connvert_lmer_c()



void convert_lmer_ta(string lmer, unsigned long int &gb, unsigned long int &ns, int asize){
	long int j;
	gb = 0; 
	ns = 0;
	for(j = 0; j < asize; j++){
		gb <<= 1;
		ns <<= 1;

		if(lmer[j] == 'T' || lmer[j] == 't' ||  lmer[j] == 'C' || lmer[j] == 'c')
			gb |= 1;
		else if(lmer[j] == 'N' || lmer[j] == 'n')
			ns |= 1;

	}//for j
}//convert_ta

void next_lmer_ta(char ch, unsigned long int &gb, unsigned long int &ns, unsigned long int MASK_READ){

	gb <<= 1;
	ns <<= 1;

	if(ch == 'T' || ch == 't'  || ch == 'C' || ch == 'c')
		gb |= 1;
	else if(ch == 'N' || ch == 'n')
		ns |= 1;

	gb &= MASK_READ;
	ns &= MASK_READ;

}//next_lmer_ta()


long int str_to_int(string s)
{
	istringstream is(s);
	long int result = -1;
	is >> result;
	return result;
}

string int_to_str(long int x)
{
	stringstream outstr;
	outstr << x;
	return outstr.str();

}
string double_to_str(double x)
{
	stringstream outstr;
	outstr << x;
	return outstr.str();

}

long int max(long int x, long int y)
{	
	if(x > y)
		return x;
	else
		return y;

}
unsigned long int min(unsigned long int x, unsigned long int y){ 
	if(x < y)
		return x;
	else
		return y;
}




void next_lmer_cg(char ch, unsigned long int &gb,  unsigned long int MASK_READ){

	gb <<= 1;

	if(ch == 'G' || ch == 'g'  || ch == 'C' || ch == 'c')
		gb |= 1;

	gb &= MASK_READ;

}//next_lmer_ta()
void next_lmer(char ch, unsigned long int &gb, unsigned long int &ns, 
			   unsigned long int &cg, unsigned long int MASK_READ)
{
	gb <<= 1;
	ns <<= 1;
	cg <<=1;

	if(ch == 'T' || ch == 't')  
		gb |= 1;
	else if(ch == 'G' || ch == 'g'){
		cg |= 1;
	}
	else if(ch == 'C' || ch == 'c'){
		gb |= 1;
		cg |= 1;
	}
	else if(ch == 'N' || ch == 'n')
		ns |= 1;

	gb &= MASK_READ;
	cg &= MASK_READ;

	ns &= MASK_READ;

}//next_lmer()

void usage(){

	string mes("\nUSAGE: brat_bw -P <directory-bw-index> -s <input-reads> -o <output> [OPTIONS]\n\nOptions:\n");
	string gen_ref_opt("  -P <directory-bw-index>    directory with BWT index");
	string four_strands(" -W                         sets mapping to four PCR-product strands (default: mapping to two original strands)");
	string singles(    "  -s <input-reads-file>      file with single reads (or queries)");
	string output(     "  -o <output-file>");
	string is_pair(    "  -pe                        set Pair-End mapping (if not specified,\n                          mapping for single reads is done)");
	string is_bs(      "  -bs                        set BS-mapping (if not specified,\n                          normal mapping is done)");
	string first(      "  -1 <paired-ends-file1>     mates 1");
	string second(     "  -2 <paired-ends-file2>     mates 2");
	string min_insert( "  -i <min>                   min insert");
	string max_insert( "  -a <max>                   max insert");
	string arich(      "  -A                         singles that are pair2 file,\n                          that are mapped to an A-rich strand");
	string out_amb(    "  -M                         output ambiguous reads/pairs");
	string allow_C_to_T ("  -C                       allow C to T matches (default: don't allow C in a read to map to T in a genome, count this as a mism)");
	string mism(	   "  -m <non-negative integer>  the number of non-BS-mismatches allowed");
//	string first_bits( "  -f <integer 32..64>        the first <int> bases, within which only one non-BS-mismatch is allowed (default: min(read_len, 48))");
	string mism_first( "  -F <integer 0 or 1>        the number of non-BS-mismatches in the first bits (first bits are defined with -f), default: 0");
//	string multiple_shifts("   -k <int 1 or 2>		 improves accuracy when set to 2, default: 2");
	string shift_delta("  -D <integer>			     interval width, by which multi-seed is shifted (default 16 bases)");
	string total_shifts("  -K <integer>              maximum number of multi-seeds considered (default: 10)");//it is read_len - BIG_seed divided by shift_delta
	cout << mes << endl;
	cout << gen_ref_opt << endl << four_strands << endl << singles << endl
        	<< output << endl
		<< is_pair << endl << is_bs << endl
		<< first << endl << second << endl << min_insert << endl << max_insert
		<< endl << arich << endl <<  out_amb << endl 
		<< mism << endl <<  mism_first << endl  << allow_C_to_T  
		<< shift_delta << endl << total_shifts << endl;

	string mes2(      "Limitations:\n");
	string ref_size(  "  a single reference size <= 2^45");
	string chrom_ids( "  number of references <= 2^18");
	string space(     "  space <= (total references sizes)*4.5*8 + 16*(2^24)");
	
	cout << mes2 << endl;
	cout << ref_size << endl << chrom_ids << endl << space << endl;

	string mes3(      "Output Format is SAM format\n");
	string ri(        "  read id <int>           the order of a read/pair in the input file (starts\n                          with 0)");
	string r1(        "  read1 <string>          mate 1 as in the input file");
	string r2(        "  read2 <string>          (pair-end only) mate 2 as in the input file");
	string chr(       "  reference name <string>");
	string str(       "  strand                  (+) if mate 1 is found in forward strand,\n                          (-) if mate 1 is in reverse strand");
	string pos1(      "  position1 <int>         mate 1 leftmost position in the reference");
	string pos2(      "  position2 <int>         (pair-end only) mate 2 leftmost pos. in the reference");
	string mism1(     "  mism1 <int>             number of mismatches in mate 1, 5' mate");
	string mism2(     "  mism2 <int>             (pair-end only) number of mismatches in mate 2, 3' mate");	
	cout << endl << mes3 << endl;
	cout << ri << endl << r1 << endl << r2 << endl << chr << endl << str << endl
		<< pos1 << endl << pos2 << endl << mism1 << endl << mism2 << endl;
}//usage()


struct gen_info
{
	gen_info() : pos(0), chrom(0) {}
	~gen_info(){}
	gen_info(unsigned int p, unsigned short int c) : pos(p), chrom(c) {}
	void operator =(gen_info a){ chrom = a.chrom; pos = a.pos; };

	unsigned int pos ;
	unsigned short int chrom ;
};

class Anchor{
public:
	Anchor(int argc, char* argv[]);//1 for CT, 2 for GA
	~Anchor();

//NOVA
		void map_entire_length(vector<AlignInfo> &gen_pos, const vector<_uli> &ta_read1, const _uli &read_len1, const vector<_uli> &c_read1, 
		const bool &forward, Mates &mates_pairs,
		const vector<_uli> &ns_read1, bool mism_map);

		void collect_gen_pos(vector<Triple> &sp_ep, vector<AlignInfo> &gen_pos, const _uli & read_len1);

		void map_local(vector<AlignInfo> &gen_pos, const vector<_uli> &ta_read1, const _uli &read_len1, const vector<_uli> &c_read1, 
		const bool &forward, Mates &mates_pairs,
		const vector<_uli> &ns_read1, bool mism_map);

		void map_indels(vector<AlignInfo> &gen_pos, const vector<_uli> &ta_read1, const _uli &read_len1, const vector<_uli> &c_read1, 
			const bool &forward, Mates &mates_pairs,
			const vector<_uli> &ns_read1, bool mism_map);

		void calculate_indel_score(const AlignInfo &leftmost , const AlignInfo &rightmost , const _uli &read_len1, 
			string &new_md, int &indel_score, _ui &mappable_score, _ui &mappable_bases,
			_ui & new_clip_left, _ui & new_clip_right, _ui & new_read_len, string &cigar);
	

		void calculate_MD(const vector<_uli> &ta_gen_mer, const vector<_uli> &cg_gen_mer, 
			const vector<_uli> &ta_check, const _ui &read_len, string &md, _ui &clip_left, _ui &clip_right, int &cur_local_score);
	
		void scores_from_md_left(vector<int> &scores, vector<int> & backtrace, const _uli& read_len1, const string &md);
		void scores_from_md_right(vector<int> &scores, vector<int> & backtrace, const _uli& read_len1, const string &md);

		void concatenate_md(const string & right_md, const string &left_md, string &new_md);

		void clip_md(string &md, const _ui &clip_left, const _ui &clip_right);
		void reverse_md(string &md, const _usi & type_alignment);
		void reverse_simple_md(string &s);
		void reverse_cigar(string &s);

		vector<char> bits_to_char;
		int entire_score;
		int local_score;//threshold on local alignment
		double alignment_score;
		double alignment_length;
		bool is_local_alignment;//user-defined if local alignment is allowed
		bool is_indel_alignment;//user-defined if indels are allowed on the top of local
		bool is_only_entire;//determined from is_local_alignment

		//maximum number of pairs(gen_pos1, gen_pos2) that we want to consider, if greater than this threshold, don't find at all
		_ui indel_threshold;//(similar to total= ep-sp+1 threshold)
		_ui gap_threshold; //max distance between two positions pos2 - pos1 <= gap_threshold to be considered for indels
		int mism_penalty ;//should be Anchor's threshold
		int gap_open_penalty;
		int gap_cont_penalty;
		int match_penalty;//not user-defined: too many functions depend on this to be equal to 1
		_ui sp_ep_range;

//END_NOVA

	vector<string> references_files;
	vector<unsigned long int> chrom_blocks;
	vector<unsigned long int> size_chrom;
	vector<unsigned long int> starts_size_chrom;
	map<long int, string> chroms_names;

	vector< vector<_ui> > num_occ_sqr;
	vector< unsigned short int* > num_occ;//indexed by char
	vector< _usi > num_occ_bytes;//hash table for the number of 1s in a byte

	unsigned long int *seeds;
	unsigned long int *seeds_ep;

	unsigned int *pos_indices;
	unsigned int *pos_num_sqr;//to hold # of pos occ in square blocks 64x64 [0...up to cur]
	unsigned short int *pos_num_occ;
	unsigned long int pos_indices_total;

	vector< _uli > C;//to store total number of occ of As and Ts in ta-genome
	vector<_uli> mask;
	vector< _uli > mask_opp;//
	vector<_uli> mask_opp_zero;

	_uli binary_total;//size of binary representation of the genome
	_uli ta_size;//64-word size of genome
	_uli genome_byte_size;
	int parse_options(int argc, char* argv[]);

	void make_mask();
	void make_mask_opp();

	_uli orig_genome_size;//nova: size of original genome (half of concatenation)
	void input_info();//chrom info, C, pos_sizes (sizes of lists of keys for genome positions identification from given index in BWT)
	void read_num_occ_sqr();
	void read_num_occ();
	void read_pos_sqr();
	void read_pos_occ();
	_uli get_pos_lmer(_uli &cur_byte, _uli aread_len);

	void print(Stat &astat);

	long int convert_lmer_ta_cg(string lmer, _uli &ta, _uli &cg, _uli &ns, vector<_uli>  &full, vector<_uli>  &full_rc, _ui asize, 
		Mates &mates_pairs, _uli &ta_compl, _uli &cg_compl, _uli &ns_compl);
	long int convert_lmer_ta_cg(string lmer, vector<_uli> &ta_read, vector<_uli> &cg_read, vector<_uli> &ns_read,
		vector<_uli> &full, vector<_uli> &full_rc, _ui read_len, Mates &mates_pairs, vector<_uli> &ta_compl_read, vector<_uli> &cg_compl_read, vector<_uli> &ns_compl_read);
	void revcompl_lmer_ta_cg(const vector<_uli> &ta_read, const vector<_uli> &cg_read, const vector<_uli> &ns_read, 
		const _ui &read_len, vector<_uli> &ta_rc_read, vector<_uli> &cg_rc_read, vector<_uli> &ns_rc_read);
	void shift_back(vector<_uli> &ta_read, int asize, const _ui &read_len);

	void convert_gen_mer_ta(vector<_uli> &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int aread_len);
	void convert_gen_mer_cg(vector<_uli> &ta_gen_lmer,
					 unsigned long int gen_pos,
					 unsigned long int aread_len);
	
	int count_mism(unsigned long int bm);
	int count_mism_long(unsigned long int bm);


	void print_singles(ofstream &out, bool second_mates, const _uli &read_id, 
				const string & aread1, const _uli &st1, const _uli &en1,  
				Mates &mates_pairs, Stat &astat, const _uli &read_len, 
				const string &read_name, const string & qual );

	void update_singles(Mates &mates_pairs, const _uli &pos1, const _usi &mism1, 
		bool map_mism, const string &md, const _ui &clip_left, const _ui &clip_right, 
		const int &cur_local_score, const _ui &read_len, const string &cigar, _usi type_map);

	void count_ag_mism(const vector<_uli> &ta_read1, const vector<_uli> &c_read1, const vector<_uli> &cg_gen_mer1,
	vector<_uli> &ag_mism1);
	void count_tc_mism(const vector<_uli> &ta_read1, const vector<_uli> &c_read1, const vector<_uli> &cg_gen_mer1,
	vector<_uli> &ag_mism1);
	int cg_match(const _uli &read_len1, const _uli &pos1, const vector<_uli> &ta_read1, const vector<_uli> &c_read1,
		vector<_uli> &ta_check, const bool &forward, bool &cg_mism1, 
		vector<_uli> &cg_gen_mer1);
	void ta_match(const _uli &read_len1, const _uli &pos1, const vector<_uli> &ta_read1, vector<_uli> &ta_check, 
		const vector<_uli> &ns_read1, int &mism1, vector<_uli> &ta_gen_mer1);

	void retreive_seed(_ui ind, _uli &sp , _uli &ep);

	void align_linear(vector<_uli> ta, _uli &sp, _uli &ep, _uli read_len1, _uli &remainder, _uli y);//overloaded for long reads
	void fill_ta_rest(vector<_uli> &ta, vector<long int> &lengths, long int &bl, _uli &remainder, _uli &ta_rest, int &asize);
	void map_pos(const _uli &sp, const _uli &ep, const vector<_uli> &ta_read1, const _uli &read_len1, const vector<_uli> &c_read1, 
		const bool &forward, const _uli &remainder, Mates &mates_pairs,
		const vector<_uli> &ns_read1, bool mism_map);

	void print_pairs(const _uli &read_id, const string &aread1, const _uli &read_len1, const _uli &st1,
		const _uli & en1, const string &read_name, const string &qual, Mates &mates_pairs1, const string &aread2, const _uli &read_len2, 
		const _uli &st2, const _uli &en2, const string &read_name2, const string &qual2,  Mates &mates_pairs2, Stat &astat);

	_uli pos_Occ(_uli &sp);
	_uli Occ(_uli ch, _uli *sp);	
	void exact_match(_uli ta_read, _uli &sp, _uli &ep, _uli read_len, _uli &remainder, _uli &valid_sp, _uli &valid_ep);//sp=1, ep = 0, remainer = 0: init
	void exact_match_seed(_uli ta_read, _uli &sp, _uli &ep, _uli read_len, _uli &remainder);
	_uli get_bwt_lmer(_uli sp, _uli asize, _uli ch);
	long int find_gen_pos(_uli sp, _uli remainder);
	void lf_map(_uli * &sp);
	void read_pos_index();//index of pos_index in counts not Bytes
	_uli get_marked_char(_uli *sp);
	_uli get_bwt_char(_uli *sp);

	_uli get_pos_strand_lmer(_uli ind);
	void read_genome(ifstream &in, _uli * &agenome);
	void read_seeds();
	void read_bwt_genome();
	void count_char(_uli kmer, vector<_usi> & characters, int asize);
	void init_params(long int read_len, Mates& mates_pairs);
	void space_release();
	void read_rev_complement(string &aread, const _uli &read_len);

	_uli *bwt_genome;
	_uli *ta_genome;
	_uli *cg_genome;

	_uli *bwt_marked;//1 if bwt[i] has gen index
	_uli *pos_strand;//nova 1 for each marked bwt_marked if genomic pos comes from the second half of concatenated genome (positive strand)
	string prefix;
	_ui num_chroms;
	unsigned long int len_mer;
	unsigned long int read_len;
	_uli hash_table_size;
	long int min_insert;
	long int max_insert;
	unsigned long int MASK_BYTE;
	unsigned long int MASK_READ;
	unsigned long int MASK_LMER;
	unsigned long int MASK_HALF;
	unsigned long int MASK_SEED;
	unsigned long int byte_size;
	unsigned long int abyte ;
	unsigned long int SIZE ;//genome char per an unsigned long int
	unsigned long int SEED;
	unsigned long int BIG_SEED;
	unsigned long int MIN_READ;
	unsigned long int STEP;
	unsigned long int two_in_Step;
	unsigned long int two_in_Size; //log_2(SIZE)
	unsigned long int two_in_Min_Read;
	unsigned long int MASK_SIZE;
	unsigned long int block_size ;
	unsigned long int double_block ;//2 end of lines
	unsigned long int two_in_byte ;//log_2(byte_size)
	unsigned long int block_byte_prod_divide ;//=log_2(256*8)
	unsigned long int block_divide ;//log_2(256)
	unsigned long int bytes;//number of bytes per unsigned long int
	unsigned long int half_ui;
	unsigned long int two_in_half_ui;
	_uli block_64;//64*64=4096
	_uli block_32;//32*32 = 1024
	_uli block_256; //256*256
	_uli two_in_block_64;//2^12=4096
	_uli two_in_block_32;//2^10
	_uli two_in_block_256;//2^16 = 256^2
	_uli two_in_32;
	_uli min_read_len;

	int mism_thresh;
	_ui first_bits;//within these first bits only one non-BS-mism is allowed
	_ui all_shifts;//improves accuracy at expence of slower alignment

	unsigned long int MASK_chrom;
	ofstream out_pairs;

	ofstream out_amb1;
	ofstream out_amb2;
	ofstream out_unmapped1;
	ofstream out_unmapped2;

	string output_pairs;
	string reads_file1;
	string reads_file2;

	string is_pe;
	unsigned long int cur_block;



	string is_bs;
	bool cascade_amb;
	bool A_rich;
	bool unmapped;


	long int total_reads;
	unsigned long int dif_len_mer;

	bool out_amb;

	ifstream in_bwt;
	ifstream in_pos_strand;
	ifstream in_num_occ;
	ifstream in_num_sqr;	
	ifstream in_ta_genome;
	ifstream in_cg_genome;

	ifstream in_pos_indices;
	ifstream in_pos_sqr;
	ifstream in_pos_occ;
	ifstream in_chrom;
	ifstream in_marked;
	ifstream in_seeds;

	void error_input(string inp);
	_uli current_read_id;
	_uli orig_gen_inL;
	_uli orig_num_sqr;
	_uli orig_num;

	_uli debug;
	//to support 2-bit BWT
	_uli full_binary_total;
	_uli full_size;


	vector< vector<_usi> > num_occ_char;//hash table for the number of char in 64-bit, 32-bases word

	int allow_C_to_T;//when set allows Cs map to Ts for forward, and G to A to rev-compl
	_uli max_num_shifts;
	_uli shift_delta;
	_uli two_in_shift_delta;
	_uli half_word;
	_uli half_SEED;
	int max_shift_delta;
	int two_in_max_shift_delta;

	bool four_strands;//to support 4-strands PCR-product mapping

};	

void Anchor::init_params(long int read_len, Mates & mates_pairs)
{
	//NOVA
		mates_pairs.A_rich = A_rich;
		entire_score = alignment_score * read_len;
		local_score = alignment_score * alignment_length * read_len;//minimum threshold on local alignment 
		mism_thresh = read_len - entire_score;
		if(mism_thresh == 0)
			max_num_shifts = 1;

        if(read_len <= 55){
                shift_delta = 8;//in bits
                two_in_shift_delta = 2;//log_base2(4) 4 is in bases
				if(first_bits > 0)
					BIG_SEED = read_len;
        }
        else if(read_len <= 64){
                shift_delta = 16;//in bits
                two_in_shift_delta = 3;//log_base2(8) 4 is in bases
  				if(first_bits > 0 && read_len < 60)
					BIG_SEED = read_len;
				else if(first_bits > 0)
					BIG_SEED = 32;
        }
        else if(read_len <= 100){
				if(first_bits > 0)
					BIG_SEED = 32;

                shift_delta = 16;//in bits
                two_in_shift_delta = 3;//log_base2(8) 4 is in bases
        }
        else if(read_len <= 200){
			if(first_bits > 0)
                BIG_SEED = 32;
            shift_delta = 32;//in bits
            two_in_shift_delta = 4;//log_base2(16) 4 is in bases
        }
		else{
			if(first_bits > 0)
                BIG_SEED = 32;
            shift_delta = max_shift_delta;//in bits
            two_in_shift_delta = two_in_max_shift_delta;//log_base2(16) 4 is in bases		
		
		}


}//init_params


void Anchor::space_release()
{
	if(pos_strand != 0)
		delete [] pos_strand;
	pos_strand = 0;
	if(bwt_genome != 0)
		delete [] bwt_genome;
	bwt_genome = 0;


	if(ta_genome != 0)
		delete [] ta_genome;
	ta_genome = 0;
	if(cg_genome != 0)
		delete [] cg_genome;
	cg_genome = 0;
	if(bwt_marked != 0)
		delete [] bwt_marked;
	bwt_marked = 0;
	for(int y = 0; y < 4; y++){
		if(num_occ[y] != 0)
			delete [] num_occ[y];
		num_occ[y] = 0;
	}
	if(pos_indices != 0)
		delete [] pos_indices;
	pos_indices = 0;
	if(pos_num_occ != 0)
		delete [] pos_num_occ;
	pos_num_occ = 0;
	if(pos_num_sqr != 0)
		delete [] pos_num_sqr;
	pos_num_sqr = 0;
	if(seeds != 0)
		delete [] seeds;
	seeds = 0;
	if(seeds_ep != 0)
		delete [] seeds_ep;
	seeds_ep = 0;
	in_pos_strand.close();
	in_bwt.close();
	in_num_occ.close();
	in_num_sqr.close();
	in_pos_sqr.close();
	in_pos_occ.close();
	in_ta_genome.close();
	in_cg_genome.close();

	in_pos_indices.close();
	in_chrom.close();
	in_marked.close();

	out_unmapped1.close();
	out_unmapped2.close();
}//space_release

Anchor::~Anchor(){

	if(bwt_genome != 0)
		delete [] bwt_genome;
	bwt_genome = 0;


	if(ta_genome != 0)
		delete [] ta_genome;
	ta_genome = 0;
	if(cg_genome != 0)
		delete [] cg_genome;
	cg_genome = 0;
	if(bwt_marked != 0)
		delete [] bwt_marked;
	bwt_marked = 0;
	if(pos_strand != 0)
		delete [] pos_strand;
	pos_strand = 0;
	for(int y = 0; y < 4; y++){
		if(num_occ[y] != 0)
			delete [] num_occ[y];
		num_occ[y] = 0;
	}
	if(pos_indices != 0)
		delete [] pos_indices;
	pos_indices = 0;
	if(pos_num_occ != 0)
		delete [] pos_num_occ;
	pos_num_occ = 0;
	if(pos_num_sqr != 0)
		delete [] pos_num_sqr;
	pos_num_sqr = 0;
	if(seeds != 0)
		delete [] seeds;
	seeds = 0;
	if(seeds_ep != 0)
		delete [] seeds_ep;
	seeds_ep = 0;
	in_pos_strand.close();
	in_bwt.close();
	in_num_occ.close();
	in_num_sqr.close();
	in_pos_sqr.close();
	in_pos_occ.close();
	in_ta_genome.close();
	in_cg_genome.close();

	in_pos_indices.close();
	in_chrom.close();
	in_marked.close();

	out_unmapped1.close();
	out_unmapped2.close();
//	cout << "Anchor is finished" << endl;
}
void Anchor::error_input(string inp)
{
	cout << "Could not open " << inp << endl;
	space_release();
	exit(0);
}
int Anchor::count_mism(unsigned long int bits){

	int count_ones = num_occ_bytes[bits & mask[SEED]];
	_uli bits2 = bits >> SEED;
	_uli bits3 = bits2 >> SEED;
	count_ones += num_occ_bytes[bits2 & mask[SEED]];
	count_ones += num_occ_bytes[bits3];

        return count_ones;

}//count_mism()

int Anchor::count_mism_long(unsigned long int bits){

/////////////////Andrew Smith author of the following code:
  bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
  bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
  bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
  bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
  bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
  return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);
/////////////////////END Andrew Smith

}//count_mism_long()

Anchor::Anchor(int argc, char* argv[])
{
//NOVA
	bits_to_char.resize(4, 0);
	bits_to_char[0] = 'A';
	bits_to_char[1] = 'G'; //ta-cg-representation is in binary 01 = 1 in decimal
	bits_to_char[2] = 'T'; //ta-cg-representation is in binary 10 = 2 in decimal
	bits_to_char[3] = 'C';//this is used to convert ta-cg-representation bits to character in genome, C is 1 in ta- and 1 in cg- representation, so put ta, followed by cg, we'll get 11 in binary, which is 3 
	
	is_local_alignment = true;//user-defined if local alignment is allowed
	is_indel_alignment = true;//user-defined if indels are allowed	
	is_only_entire = false;

	//maximum number of pairs(gen_pos1, gen_pos2) that we want to consider, if greater than this threshold, don't find at all
	indel_threshold = 100;// 20;//(similar to total= ep-sp+1 threshold)
	gap_threshold = 50; //max distance between two positions pos2 - pos1 <= gap_threshold to be considered for indels
	mism_penalty = -3;//should be Anchor's threshold
	gap_open_penalty = -5;
	gap_cont_penalty = -3;
	match_penalty = 1;//do not change this: not user-defined: too many functions depend on this to be 1
	sp_ep_range = 1000;
	alignment_score = 90;//this multiplied by mappable length of the read is threshold on the score, i.e. this is quality of alignment
	alignment_length = 30;//in local alignment: min mappable length of the read (without clipped ends)
	//for longer reads with indels, use higher alignment_length: affects mapping accuracy

	//END_NOVA
	four_strands = false;

	max_shift_delta = 64;//in bits
	two_in_max_shift_delta = 5;

	half_word = 32;
	max_num_shifts = 10;
	shift_delta = 32;//in bits
	allow_C_to_T = 0;

	pos_indices_total = 0;
	all_shifts = 2;
	debug = 0;
	SEED = 24;//const
	half_SEED = SEED >> 1;
	BIG_SEED = 32;
	MIN_READ = 32;//const 
	first_bits = 0;//here it is the number of mism in SEED 
	out_amb = false;
	unmapped = false;
	A_rich = false;
	mism_thresh = 0;
	min_insert = 100;
	max_insert = 300;
	is_bs = "1";
	is_pe = "0";
	current_read_id = 0;


	int res = parse_options(argc, argv);


	if(mism_penalty < -10)
		mism_penalty = -10;
	if(gap_threshold > 1000)
		gap_threshold = 1000;
	if(sp_ep_range > 1000)
		sp_ep_range = 1000;
	if(max_num_shifts > 20)
		max_num_shifts = 20;

	if(is_local_alignment == false)
		is_only_entire = true;

	if(res == -1){
		usage();
		space_release();
		exit(0);
	}
	is_bs = "1";
	if(shift_delta != 32){
		//was changed
		if(shift_delta == 0)
			shift_delta = 32;
		else
			shift_delta = shift_delta << 1;//multiply by 2, convert from bases to bits
	}
	two_in_shift_delta = floor(log(static_cast<double>(shift_delta >> 1))/log(2.0));
	two_in_max_shift_delta = floor(log(static_cast<double>(max_shift_delta >> 1))/log(2.0));

	if(all_shifts != 2)
		all_shifts = 2;

	if(reads_file1.empty()){
		cout << "ERROR: file with reads is missing (option -s or -1 or -2)" << endl;
		space_release();
		exit(0);
	}
	if(is_pe == "1" && reads_file2.empty())
	{
		cout << "ERROR: file with mates 2 is missing for paired-end alignment (option -2)." << endl;
		space_release();
		exit(0);
	}
	if(output_pairs.empty()){
		cout << "ERROR: output file is not provided (option -o)" << endl;
		space_release();
		exit(0);
	}
	if(prefix.empty()){
		cout << "ERROR: path to BW-index is not provided (option -P)" << endl;
		space_release();
		exit(0);
	}
	if(BIG_SEED > 64)
		BIG_SEED = 64;

	if(first_bits > 1)
		first_bits = 1;

	if(res == -1 || prefix.empty()){
		usage();
		space_release();
		exit(0);
	}



	out_pairs.open(output_pairs.c_str(), ios::out);

	string unm1 =  reads_file1 + ".unm";
	string unm2 = reads_file2 + ".unm";
	if(unmapped == true ){
		out_unmapped1.open(unm1.c_str(), ios::out);
		if(is_pe == "1"){
			out_unmapped2.open(unm2.c_str(), ios::out);
		}
	}//if unm

	string amb1 = reads_file1 + ".amb";
	string amb2 = reads_file2 + ".amb";
	if(out_amb == true ){
		out_amb1.open(amb1.c_str(), ios::out);
		if(is_pe == "1")
			out_amb2.open(amb2.c_str(), ios::out);
	}
		ofstream out_pairs;

	string rem = (prefix[prefix.length() - 1] == '/' ? "" : "/");

	string genome_repres = "CT"; 

	string input_ta_genome = prefix + rem + "ta.ta";
	string input_cg_genome = prefix + rem +  "cg.cg";

	string input_chrom = prefix + rem + genome_repres + ".chr";
	string input_pos_indices = prefix + rem + genome_repres + ".pos_ind";
	string input_pos_sqr = prefix + rem + genome_repres + ".pos_num_sqr";
	string input_pos_occ = prefix + rem + genome_repres + ".pos_num_occ";
	string input_num_occ = prefix + rem + genome_repres + ".no";
	string input_num_sqr = prefix + rem + genome_repres + ".nosqr";
	string input_bwt = prefix + rem + genome_repres + ".bwt";
	string input_marked = prefix + rem + genome_repres + ".bwt_marked";
	string input_seeds = prefix + rem + genome_repres + ".seeds";
	string input_pos_strand = prefix + rem + "pos_strand.txt";

	in_bwt.open(input_bwt.c_str(), ios::in);
	if(!in_bwt){
		error_input(input_bwt);
	}
	in_pos_strand.open(input_pos_strand.c_str(), ios::in);
	if(!in_pos_strand)
		error_input(input_pos_strand);
	in_num_occ.open(input_num_occ.c_str(), ios::in);
	if(!in_num_occ)
		error_input(input_num_occ);
	in_num_sqr.open(input_num_sqr.c_str(), ios::in);
	if(!in_num_sqr)
		error_input(input_num_sqr);
	
		in_ta_genome.open(input_ta_genome.c_str(), ios::in);
		if(!in_ta_genome)
			error_input(input_ta_genome);
		in_cg_genome.open(input_cg_genome.c_str(), ios::in);
		if(!in_cg_genome)
			error_input(input_cg_genome);


	in_pos_indices.open(input_pos_indices.c_str(), ios::in);
	if(!in_pos_indices)
		error_input(input_pos_indices);
	in_pos_sqr.open(input_pos_sqr.c_str(), ios::in);
	if(!in_pos_sqr)
		error_input(input_pos_sqr);
	in_pos_occ.open(input_pos_occ.c_str(), ios::in);
	if(!in_pos_occ)
		error_input(input_pos_occ);
	in_chrom.open(input_chrom.c_str(), ios::in);	
	if(!in_chrom)
		error_input(input_chrom);
	in_marked.open(input_marked.c_str(), ios::in);
	if(!in_marked)
		error_input(input_marked);
	in_seeds.open(input_seeds.c_str(), ios::in);
	if(!in_seeds)
		error_input(input_seeds);

	_uli one = 1;
	min_read_len = 24;
	half_ui = 16;
	two_in_half_ui = 4;
	abyte = 8;
	SIZE = 64;//const: genome char per an unsigned long int
	two_in_Size = 6;
	MASK_SIZE = (one << two_in_Size) - 1;
	MASK_HALF = (one << 32) -1 ;
	block_size = 256;
	double_block = 256 + 256 + 2;//2 end of lines
	byte_size = 8;//byte size is characters(genome) per a byte
	two_in_byte = 3;//log_2(byte_size)

	block_byte_prod_divide = 11;//=log_2(256*8)
	block_divide = 8;//log_2(256)
	bytes = SIZE >> two_in_byte;//number of bytes per unsigned long int

	MASK_SEED = (one << SEED) -1;
	len_mer = SEED;
	block_32 = 1024;
	block_64 = SIZE << two_in_Size;//64*64
	two_in_Min_Read = 5;

	two_in_block_32 = 10;
	two_in_block_64 = 12;
	two_in_32 = 5;

	MASK_BYTE = (one << (byte_size)) -1;
	if(len_mer == SIZE)
		MASK_LMER = 0xFFFFFFFFFFFFFFFF;
	else if(len_mer < 32)
		MASK_LMER = (one << len_mer) - 1;//seed mask
	else{
		int dist = len_mer - 32;
		MASK_LMER = (((one << dist) - 1) << 32) | 0x00000000FFFFFFFF;	
	}

	dif_len_mer = SIZE - len_mer;
	MASK_READ = MASK_LMER;//for initialization only
	make_mask();
	make_mask_opp();

	input_info();//reads chroms info, C, pos_sizes and values for binary_total and ta-size

	two_in_Step = log(static_cast<float>(STEP))/log(static_cast<float>(2));
	block_256 = pow(2.0, 24.0);
	two_in_block_256 = 24;

	read_pos_sqr();

	read_pos_occ();

	read_pos_index();


	read_num_occ();

	read_num_occ_sqr();


	hash_table_size = pow(2.0, 24.0);
	num_occ_bytes.resize(hash_table_size);//2^16	
	for(long int i = 0; i < hash_table_size; i++){
        num_occ_bytes[i] = count_mism_long(i);
	}//for

	read_seeds();


	genome_byte_size = ((binary_total + 1) >> two_in_byte) + 1;


	read_genome(in_marked, bwt_marked);
	read_genome(in_pos_strand, pos_strand);


	read_genome(in_ta_genome, ta_genome);
	read_genome(in_cg_genome, cg_genome);


	read_bwt_genome();


	int apow_seed = 16;
	int alphabet_size = 4;
	int hash_small_size = pow(2.0, apow_seed + 0.0);
	num_occ_char.resize(hash_small_size);//2^SEED	
	for(_uli i = 0; i < hash_small_size; i++){
		vector< _usi > characters(alphabet_size, 0);
		count_char(i, characters, apow_seed);
        num_occ_char[i] = characters;

	}//for

}//Anchor()

void Anchor::count_char(_uli kmer, vector<_usi> & characters, int asize){

	for(int i = 0; i < asize; i += 2){
		_uli cur_char = kmer & mask[2];
		characters[cur_char]++;
		kmer >>= 2;
	}
}//count_char

void Anchor::read_genome(ifstream &in_genome, _uli * &agenome)
{

	agenome = new _uli[ta_size];
	assert(agenome != 0);

	long int buffer_size = block_256 << two_in_byte;//*8
	char *buffer = new char[buffer_size];
	assert(buffer != 0);

	in_genome.read(buffer, buffer_size);
	long int real_size = in_genome.gcount();

	long int gen_ind = 0;

	while(true){
		long int buf_ind = 0;
		while(buf_ind < real_size){
			_uli res = 0;
			for(int y = 0; y < 8 ; y++){
				res <<= byte_size;
				res |= int(buffer[buf_ind]) & MASK_BYTE;
				buf_ind++;
			}//for
			agenome[gen_ind] = res;
			gen_ind++;
		}//while buffer size
		if(in_genome.eof())
			break;
		in_genome.read(buffer, buffer_size);
		real_size = in_genome.gcount();
		if(real_size == 0)
			break;

	}//while true
	delete [] buffer;
//	cout << "finished reading bwt" << endl;
}//read_genome

void Anchor::read_bwt_genome()
{
	long int bwt_size = full_size + 1 ;
	bwt_genome = new _uli[bwt_size];
	assert(bwt_genome != 0);

	const _ui buffer_size = block_256 << two_in_byte;//*8
	char *buffer = new char[buffer_size];
	assert(buffer != 0);

	in_bwt.read(buffer, buffer_size);
	_uli real_size = in_bwt.gcount();

	long int bwt_ind = 0;

	while(true){
	  long int buf_ind = 0;
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

}//read_bwt_genome

void Anchor::make_mask(){
	_uli one = 1;
	mask.push_back(0);
	for(int i = 1; i <= SIZE; i++){
		_uli res = (one << i) - 1;
		mask.push_back(res);
	}//for
	mask[SIZE] = 0xFFFFFFFFFFFFFFFF;
}//make mask

void Anchor::make_mask_opp(){
	_uli one = 1;
	mask_opp.resize(SIZE);
	mask_opp[0] = 1;
	mask_opp_zero.resize(SIZE);
	mask_opp_zero[0] = ~mask_opp[0];
	for(int i = 1; i < SIZE; i++){
		mask_opp[i] = one << i;
		mask_opp_zero[i] = ~mask_opp[i];
	}//for

}//make mask

void Anchor::input_info(){

	in_chrom >> binary_total ;
	in_chrom >> ta_size;
	in_chrom >> full_binary_total;
	in_chrom >> full_size;
	in_chrom >> num_chroms ;	

	_uli zero = 0;
	C.resize(6, zero);
	in_chrom >> C[0] >>  C[1] >> C[2] >> C[3] ;

	C[5] = C[0] + C[1] + C[2] + C[3] + C[4] + 1;
	C[4] = C[0] + C[1] + C[2] + C[3] + 1;
	C[3] = C[0] + C[1] + C[2] + 1;
	C[2] = C[0] + C[1] + 1;
	C[1] = C[0] + 1;
	C[0] = 1;

	in_chrom >> orig_gen_inL;
	in_chrom >> STEP ;
	in_chrom >> pos_indices_total;

	orig_num_sqr = orig_gen_inL >> two_in_block_64;
	orig_num = orig_gen_inL >> two_in_Size;
/*
//NOVA_DEBUG
cout << binary_total << endl
	<< ta_size << endl
	<< full_binary_total << endl
	<< full_size << endl
	<< num_chroms << endl
	<< C[0] << " " << C[1] << " " << C[2] << " " << C[3] << " " << C[4] << " " << C[5] << endl	
	<< orig_gen_inL << endl
	<< STEP << endl
	<< pos_indices_total << endl;
*/	

	map<long int, string>::iterator beg = chroms_names.begin();
	map<long int, string>::iterator fin = chroms_names.end();

	//index and name of chroms

	for(_ui i = 0; i < num_chroms; i++){
		_ui chr_id;
		string achr_name;
		in_chrom >> chr_id >> achr_name ;
		chroms_names.insert(make_pair(chr_id, achr_name));
	}

	//files with references
	_uli orig_num_chroms = num_chroms >> 1;
	for(_uli i = 0; i < orig_num_chroms; i++){
		string afile;
		in_chrom >> afile;
		references_files.push_back(afile);
	}//for i
	//sizes of chroms in bits: each bit per char //over entire concatenation
	for(_ui i = 0; i < num_chroms; i++){
		_uli achr_size;
		in_chrom >> achr_size;
		size_chrom.push_back(achr_size);
	}

	//starts of chroms in bits
	for(_ui i = 0; i < num_chroms; i++){
		_uli achr_start ;
		in_chrom >> achr_start;
//cout << achr_start << endl;
		starts_size_chrom.push_back(achr_start);
	}
	in_chrom >> orig_genome_size;//nova: reads original genome's size (half of the concatination that is used here)

//cout << orig_genome_size << endl;

}//input_info
//lookup8.c, by Bob Jenkins, January 4 1997, Public Domain.

void Anchor::retreive_seed(_ui ind, _uli &sp , _uli &ep){

	sp = seeds[ind];
	ep = seeds_ep[ind];

}//retreive_seed

void Anchor::read_seeds()
{
	seeds = new _uli[hash_table_size];//stores sp
	seeds_ep = new _uli[hash_table_size]; // stores ep

	_uli buffer_size = hash_table_size << two_in_half_ui;// *16 8byts for sp and 8bytes for ep
	char *buffer = new char[buffer_size];
	in_seeds.read(buffer, buffer_size);
	_uli read_char = static_cast<_uli>(in_seeds.gcount());
	if(read_char < buffer_size){
		cout << "\nERROR: Unexpected size of seeds file" << endl;
		space_release();
		exit(0);
	}
	_uli buf_ind = 0;

	for(_uli i = 0; i < hash_table_size; i++){
		_uli res = 0;
		for(_ui j = 0; j < 8; j++){
			//sp 8 bytes and ep 8 bytes
			res <<= 8;
			res |= (int(buffer[buf_ind]) & MASK_BYTE);
			buf_ind++;
		}//for j
		seeds[i] = res;
		res = 0;
		for(_ui j = 0; j < 8; j++){
			//sp 8 bytes and ep 8 bytes
			res <<= 8;
			res |= (int(buffer[buf_ind]) & MASK_BYTE);
			buf_ind++;
		}//for j
		seeds_ep[i] = res;
	}//for i
	delete [] buffer;
}//read_seeds

void Anchor::read_pos_sqr(){
	_uli full_binary = ta_size << two_in_Size;
	_uli total_occ = full_binary >> two_in_block_64;
	if((full_binary & mask[two_in_block_64]) > 0)
		total_occ++;

	const _uli buffer_size = total_occ << 2; //*4 for unsigned int
	pos_num_sqr = new _ui[total_occ];
	assert(pos_num_sqr != 0);

	char *buffer = new char[buffer_size];
	assert(buffer != 0);

	in_pos_sqr.read(buffer, buffer_size);
	_uli read_char = static_cast<_uli>(in_pos_sqr.gcount());
	if(read_char < buffer_size){
		cout << "\nERROR: Unexpected size of positions' sqr occ file" << endl;
		space_release();
		exit(0);
	}
	_uli buffer_ind = 0;
	for(_uli i = 0; i < total_occ; i++){
		_uli res = 0;
		for(_uli y = 0; y < 4; y++){//4 is bytes in unsigned int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		pos_num_sqr[i] = res;
	}//for i

	delete [] buffer;
	in_pos_sqr.close();
}//read_pos_sqr

void Anchor::read_pos_occ(){
	_uli total_occ = ta_size;
	const _uli buffer_size = total_occ << 1; //*2 size of unsigned short int
	char *buffer = new char[buffer_size];
	in_pos_occ.read(buffer, buffer_size);
	_uli read_char = in_pos_occ.gcount();
	if(read_char < buffer_size){
		cout << "\nERROR: Unexpected size of positions' keys file" << endl;
		space_release();
		exit(0);
	}

	pos_num_occ = new _usi[total_occ];
	assert(pos_num_occ != 0);

	_uli buffer_ind = 0;
	for(_uli i = 0; i < total_occ; i++){
		_uli res = 0;
		for(_uli y = 0; y < 2; y++){//2 is bytes in unsigned short int
			res <<= byte_size;
			char x = buffer[buffer_ind];
			buffer_ind++;
			res |= (int(x) & MASK_BYTE);
		}//for y	
		pos_num_occ[i] = res;
	}//for i

	delete [] buffer;
	in_pos_occ.close();
}

void Anchor::read_num_occ_sqr()
{
	_uli full_binary = ta_size << two_in_Size;
	_uli total_occ = full_binary >> two_in_block_64;
	if((full_binary & mask[two_in_block_64]) > 0)
		total_occ++;

	const _uli buffer_size = total_occ << 4; //*4 *4 for 4 char each 4 bytes
	char *buffer = new char[buffer_size];
	assert(buffer != 0);

	in_num_sqr.read(buffer, buffer_size);
	_uli read_char = static_cast<_uli>(in_num_sqr.gcount());
	if(read_char < buffer_size){
		cout << "\nERROR: Unexpected size of positions' keys file" << endl;
		space_release();
		exit(0);
	}
	vector<_ui> dum(4, 0);
	num_occ_sqr.resize(total_occ, dum);
	_uli buffer_ind = 0;
	for(_uli i = 0; i < total_occ; i++){
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
}//read_num_occ_sqr

void Anchor::read_num_occ()
{
	//ta_size because we count number of occ of a char in bwt of N bases length
	//the count is in bases not in bits: there are ta_size of 64-bases words
	_uli total_occ = ta_size ;//*4 for four characters
	num_occ.resize(4);
	for(int y = 0; y < 4; y++){
		num_occ[y] = new unsigned short int[total_occ];
		assert(num_occ[y] != 0);
	}//for y

	const _ui buffer_size = block_256 << two_in_byte;//*8: 4char each 2 bytes
	char *buffer = new char[buffer_size];
	assert(buffer != 0);

	in_num_occ.read(buffer, buffer_size);
	_uli real_size = in_num_occ.gcount();
 
	_uli num_occ_ind = 0;

	while(true){//while reading entire file

		_uli buffer_ind = 0;
		while(buffer_ind < real_size){
			//As
			_ui res = 0;
			for(_uli y = 0; y < 2; y++){//2 is bytes in unsigned short int
				res <<= byte_size;
				char x = buffer[buffer_ind];
				buffer_ind++;
				res |= (int(x) & MASK_BYTE);
			}//for y	
			num_occ[0][num_occ_ind] = res;

			//Cs
			res = 0;
			for(_uli y = 0; y < 2; y++){//2 is bytes in unsigned short int
				res <<= byte_size;
				char x = buffer[buffer_ind];
				buffer_ind++;
				res |= (int(x) & MASK_BYTE);
			}//for y	
			num_occ[1][num_occ_ind] = res;

			//Gs
			res = 0;
			for(_uli y = 0; y < 2; y++){//2 is bytes in unsigned short int
				res <<= byte_size;
				char x = buffer[buffer_ind];
				buffer_ind++;
				res |= (int(x) & MASK_BYTE);
			}//for y	
			num_occ[2][num_occ_ind] = res;

			//Ts
			res = 0;
			for(_uli y = 0; y < 2; y++){//2 is bytes in unsigned short int
				res <<= byte_size;
				char x = buffer[buffer_ind];
				buffer_ind++;
				res |= (int(x) & MASK_BYTE);
			}//for y	
			num_occ[3][num_occ_ind] = res;
			num_occ_ind++;

		}//while buffer size
		if(in_num_occ.eof())
			break;
		in_num_occ.read(buffer, buffer_size);
		real_size = in_num_occ.gcount();
		if(real_size == 0)
			break;
	}//while entire file

	delete [] buffer;
}//read_num_occ

void Anchor::read_pos_index()//index of pos_index in counts not Bytes
{
	_uli total_keys = pos_indices_total;

	//total bytes to read = total_keys * 4 for _ui
	const _ui buffer_size = block_256 << 2;//*4
	char *buffer = new char[buffer_size];
	assert(buffer != 0);

	in_pos_indices.read(buffer, buffer_size);
	_uli real_size = in_pos_indices.gcount();

	pos_indices = new _ui[total_keys];
	assert(pos_indices != 0);

	long int last_lmer_char = 0;
	
	long int pos_indices_ind = 0;
	while(true){
		long int buf_ind = 0;
		while( buf_ind < real_size){
			_ui res = 0;
			for(_ui y = 0; y < 4; y++){//4 is bytes in unsigned int
				res <<= byte_size;
				char x = buffer[buf_ind];
				buf_ind++;
				res |= (int(x) & MASK_BYTE);
			}//for y
			pos_indices[pos_indices_ind] = res;
			pos_indices_ind++;
		}//while buffer size

		if(in_pos_indices.eof())
			break;
		in_pos_indices.read(buffer, buffer_size);
		real_size = in_pos_indices.gcount();

		if(real_size == 0)
				break;	
	}//while
	delete [] buffer;

}//read_pos_index


long int Anchor::find_gen_pos(_uli sp, _uli remainder)
{
	//remainder is in bases

	_uli *ptr_sp = &sp;
	if(STEP == 1){
		_uli gen_ind = pos_indices[(*ptr_sp) - 1];
		_uli is_second_half = get_pos_strand_lmer(sp);
		if(is_second_half > 0){
			gen_ind = gen_ind + orig_genome_size;
		}
		return (gen_ind - remainder);
	}
	else{
		_uli gen_ind = 0;

		_uli count_lfmap = 0;
		while(true){
			_uli is_marked = get_marked_char(ptr_sp);
			if(is_marked > 0){
				gen_ind = pos_indices[pos_Occ(*ptr_sp) - 1];//pos_Occ returns count, index will be 1 less
				break;
			}//marked
			else{
				count_lfmap++;
				lf_map(ptr_sp);
			}//else
		}//while
		_uli is_second_half = get_pos_strand_lmer(sp);//must be updated with lf_map sp of marked bw_marked (only for marked bw_marked we set is_second_half bit
		if(is_second_half > 0){
			gen_ind = gen_ind + orig_genome_size;
		}
		//remainder is in bases
		return (gen_ind + count_lfmap - remainder);
	}
}//find_gen_pos

void Anchor::exact_match(_uli ta_read, _uli &sp, _uli &ep, _uli read_len, _uli &remainder, 
	_uli &valid_sp, _uli &valid_ep)//sp=1, ep = 0, remainer = 0: init
{
	//this function continues reading a read: either a long read
	//or after matching the seed
	//remainder is the rest of the read after seed

	//read_len is in bases 
	//remainder is in bases

	int i = read_len ;
	while( (sp < ep) && (i > 0)){
		_uli ch = ta_read & mask[2];
		ta_read >>= 2;
		_uli prev_sp = sp -1;
		_uli *ptr_sp = &prev_sp;
		_uli prev_ep = ep;
		_uli *ptr_ep =  &prev_ep;
		sp = C[ch] + Occ(ch, ptr_sp);
		ep = C[ch] + Occ(ch, ptr_ep) - 1;
		if(sp <= ep){
			valid_sp = sp;
			valid_ep = ep;
		}
		i--;
	}//while
	remainder = remainder - (read_len - i);
}//exact_match

_uli Anchor::Occ(_uli ch, _uli *sp){
	
	_uli double_sp = (*sp) << 1;
	_uli full_ind = (double_sp) >> two_in_Size;
	_uli full_bits = (double_sp) & mask[two_in_Size];
	_uli num_cur = get_bwt_lmer(full_ind, full_bits + 2, ch);
	//is full_ind in the first or second half of the original 64-bit sp range?
	//number of sp (1-bit count) before cur block
	_uli bl_ind = (*sp) >> two_in_Size;
	_uli sum_bl_ind = bl_ind << two_in_Size;//not same as sp because of remainder
	_uli full_sum_bl_ind = (sum_bl_ind << 1) + SIZE;
	if(double_sp >= full_sum_bl_ind){
		num_cur += get_bwt_lmer(full_ind - 1, SIZE, ch);
	}

	_uli block_ind_sqr = (*sp) >> two_in_block_64;
	_uli num_a_sqr = (block_ind_sqr == 0 ? 0 : num_occ_sqr[block_ind_sqr - 1][ch]);


    _uli bl_inside_ind =  bl_ind & mask[two_in_Size];//index inside of a conceptual block
	_uli num_prev_bl = (bl_inside_ind == 0 ? 0 : num_occ[ch][bl_ind - 1]);


	//this is true only for A character, becaus EOS($) char also was represented as A
	//and because this was not counted
	if((ch == 0) && (orig_num == bl_ind) && (orig_gen_inL <= (*sp)))
		num_cur--;

	num_a_sqr += num_prev_bl + num_cur;	
	

	return num_a_sqr;
}//Occ

_uli Anchor::get_bwt_lmer(_uli cur_byte, _uli aread_len, _uli ch){
	//aread_len here is in bits not bases

        _uli bwt_lmer = bwt_genome[cur_byte];
        _uli count_ones = 0;
		_uli ashift = (SIZE - aread_len);
		bwt_lmer >>= ashift; 
        for(int i = 0; i < 4; i++){
                count_ones += num_occ_char[bwt_lmer & mask[16]][ch];//64-bit word is subdivided into 16bit sub-words
                bwt_lmer >>= 16;
        }
		if(ch == 0)
			count_ones -= ((SIZE - aread_len) >> 1);//filled in zeros at init bwt_lmer, divide by two = number of A/N/G
        return count_ones;
}//get_bwt_lmer

_uli Anchor::pos_Occ(_uli &sp){

	_uli bl_ind = (sp) >> two_in_Size;
	_uli ind_cur = (sp) & mask[two_in_Size];
	_uli num_cur = get_pos_lmer(bl_ind, ind_cur + 1);//for Ts returns number of ones

	_uli block_ind_sqr = (sp) >> two_in_block_64;
	_uli num_a_sqr = (block_ind_sqr == 0 ? 0 : pos_num_sqr[block_ind_sqr - 1]);
    _uli bl_inside_ind =  bl_ind & mask[two_in_Size];//index inside of a conceptual block
	_uli num_prev_bl = (bl_inside_ind == 0 ? 0 : pos_num_occ[bl_ind - 1]);

	num_a_sqr += num_prev_bl + num_cur;	
	
	return num_a_sqr;
}//pos_Occ

_uli Anchor::get_pos_lmer(_uli &cur_byte, _uli aread_len)
{
    _uli bwt_lmer = bwt_marked[cur_byte];
	_uli count_ones = 0;
	_uli ashift = SIZE - aread_len;
	bwt_lmer >>= ashift;
    count_ones += num_occ_bytes[bwt_lmer & mask[SEED]];
    bwt_lmer >>= SEED;
	count_ones += num_occ_bytes[bwt_lmer & mask[SEED]];
	count_ones += num_occ_bytes[bwt_lmer >> SEED];
        
        return count_ones;
}//get_pos_lmer


_uli Anchor::get_marked_char(_uli *ind){

	if((*ind) > binary_total)//bwt size = binary_total + 1
		return 0;
	_uli ta_ind = (*ind) >> two_in_Size;//ind inside bwt_marked
	_uli bwt_lmer = bwt_marked[ta_ind];
	_uli ind_inside_Size = (*ind) & mask[two_in_Size];
	_uli ashift = SIZE - ind_inside_Size - 1;
	_uli res = (bwt_lmer >> ashift) & mask[1];
	return res;
}//get_marked_char

void Anchor::lf_map(_uli * &sp)
{
	_uli ch = get_bwt_char(sp);
	(*sp) = C[ch] + Occ(ch, sp) - 1;
}//lf_map

_uli Anchor::get_bwt_char(_uli *sp){

	_uli full_ind = ((*sp) << 1) >> two_in_Size;
	_uli bwt_lmer = bwt_genome[full_ind];
	_uli full_bits = (((*sp) << 1) & mask[two_in_Size]) + 2;

	if((full_ind) >= full_size)
		return 0;
	_uli ashift = SIZE - full_bits;
	_uli res = (bwt_lmer >> ashift) & mask[2];
	return res;
}//get_bwt_char

//only marked positions = binary_total/STEP have information in pos_strand
_uli Anchor::get_pos_strand_lmer(_uli ind){

	if(ind > binary_total)//bwt size = binary_total + 1
		return 0;
	_uli ta_ind = ind >> two_in_Size;//ind inside bwt
	_uli ns_lmer = pos_strand[ta_ind];
	_uli ind_inside_Size = ind & mask[two_in_Size];
	_uli ashift = (SIZE - ind_inside_Size - 1);
	_uli res = (ns_lmer >> ashift ) & mask[1];
	return res;
}//get_pos_strand_lmer
////////////////////////// LONG READS


void Anchor::convert_gen_mer_ta(vector<_uli> &ta_gen_read,
					 unsigned long int pos,
					 unsigned long int aread_len)
{
	if(aread_len <= SIZE){
	
		_uli ta_gen_lmer = 0;
		if(pos >= binary_total){
			cout << "ERROR reading from ta-genome: pos index is out of boundary" << endl;
			space_release();
			exit(0);
			return;
		}
		//current binary char inside current block: there 8 characters fit in one binary characters
		//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
		_uli gen_ind = (pos >> two_in_Size) ;//index of 64-bit word

		_uli prev = ta_genome[gen_ind];
		_uli next = (gen_ind + 1 >= ta_size ? 0 : ta_genome[gen_ind + 1]);

		_uli ind_inside_byte = pos & mask[two_in_Size];//pos of genome char inside 64-bit word	
		_uli bits_in_left = SIZE - ind_inside_byte;

		_uli left = prev & mask[bits_in_left];
		if(bits_in_left >= aread_len){
			ta_gen_lmer = left >> (bits_in_left - aread_len);
	
		}
		else{
		_uli right = next >> bits_in_left;
		ta_gen_lmer = ((left << (SIZE - bits_in_left)) | right) >> (SIZE - aread_len);
		}
		ta_gen_read.resize(1, 0);
		ta_gen_read[0] = ta_gen_lmer;

	}//short
	else{ //long
		_uli ta_gen_lmer = 0;

		//current binary char inside current block: there 8 characters fit in one binary characters
		//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
		_uli cur_byte = (pos >> two_in_Size) ;//index of the first byte
		_uli ind_inside_byte = pos & mask[two_in_Size];//pos of genome char inside byte-block	
		_uli bits_in_first_byte = SIZE - ind_inside_byte;
		_uli new_len = aread_len - bits_in_first_byte;
		_uli bits_in_last_byte = new_len & mask[two_in_Size];
		_uli bytes_rest = (new_len >> two_in_Size) + ( bits_in_last_byte > 0 ? 1 : 0);

		_uli num_bytes = bytes_rest + 1;
		if((cur_byte + num_bytes - 1) >= ta_size){
			ta_gen_read.resize((aread_len >> two_in_Size) + ((aread_len & mask[two_in_Size]) > 0 ? 1 : 0), 0);
			return;
		}
		ta_gen_read.resize(num_bytes, 0);
		ta_gen_read[0] = ta_genome[cur_byte] & mask[bits_in_first_byte];

		int asize = 1, j = cur_byte + 1;
		while(asize < num_bytes){
			ta_gen_read[asize] = ta_genome[j];
			asize++;
			j++;
		}//while
		//shift the last byte;
		if(bits_in_last_byte > 0){
			ta_gen_read[asize - 1] >>= (SIZE - bits_in_last_byte);
			shift_back(ta_gen_read, asize, bits_in_last_byte);
		}
		reverse(ta_gen_read.begin(), ta_gen_read.end());
		int correct_size = (aread_len >> two_in_Size) + ((aread_len & mask[two_in_Size]) > 0 ? 1 : 0);
		if(asize > correct_size)
			ta_gen_read.resize(correct_size);
	}//long reads
}//convert_gen_mer_ta

void Anchor::convert_gen_mer_cg(vector<_uli> &ta_gen_read,
					 unsigned long int pos,
					 unsigned long int aread_len)
{

	if(aread_len <= SIZE){
		_uli ta_gen_lmer = 0;
		if(pos >= binary_total){
			cout << "ERROR reading from ta-genome: pos index is out of boundary" << endl;
			space_release();
			exit(0);
			return;
		}
		//current binary char inside current block: there 8 characters fit in one binary characters
		//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
		_uli gen_ind = (pos >> two_in_Size) ;//index of 64-bit word

		_uli prev = cg_genome[gen_ind];
		_uli next = (gen_ind + 1 >= ta_size ? 0 : cg_genome[gen_ind + 1]);

		_uli ind_inside_byte = pos & mask[two_in_Size];//pos of genome char inside 64-bit word	
		_uli bits_in_left = SIZE - ind_inside_byte;

		_uli left = prev & mask[bits_in_left];
		if(bits_in_left >= aread_len){
			ta_gen_lmer = left >> (bits_in_left - aread_len);

		}
		else{
		_uli right = next >> bits_in_left;
		ta_gen_lmer = ((left << (SIZE - bits_in_left)) | right) >> (SIZE - aread_len);	
		}
		ta_gen_read.resize(1, 0);
		ta_gen_read[0] = ta_gen_lmer;

	}//short
	else{
_uli ta_gen_lmer = 0;

	//current binary char inside current block: there 8 characters fit in one binary characters
	//(pos >> two_in_byte) = pos / byte_size = number of binary characters needed to cover this pos
	_uli cur_byte = (pos >> two_in_Size) ;//index of the first byte
	_uli ind_inside_byte = pos & mask[two_in_Size];//pos of genome char inside byte-block	
	_uli bits_in_first_byte = SIZE - ind_inside_byte;
	_uli new_len = aread_len - bits_in_first_byte;
	_uli bits_in_last_byte = new_len & mask[two_in_Size];
	_uli bytes_rest = (new_len >> two_in_Size) + ( bits_in_last_byte > 0 ? 1 : 0);

	_uli num_bytes = bytes_rest + 1;
	if((cur_byte + num_bytes - 1) >= ta_size){
		ta_gen_read.resize((aread_len >> two_in_Size) + ((aread_len & mask[two_in_Size]) > 0 ? 1 : 0), 0);
		return;
	}

	ta_gen_read.resize(num_bytes, 0);
	ta_gen_read[0] = cg_genome[cur_byte] & mask[bits_in_first_byte];

	int asize = 1, j = cur_byte + 1;
	while(asize < num_bytes){
		ta_gen_read[asize] = cg_genome[j];
		asize++;
		j++;
	}//while

	//shift the last byte;
    if(bits_in_last_byte > 0){
	    ta_gen_read[asize - 1] >>= (SIZE - bits_in_last_byte);
		shift_back(ta_gen_read, asize, bits_in_last_byte);
	}
	reverse(ta_gen_read.begin(), ta_gen_read.end());
	int correct_size = (aread_len >> two_in_Size) + ((aread_len & mask[two_in_Size]) > 0 ? 1 : 0);
	if(asize > correct_size)
		ta_gen_read.resize(correct_size);
	}//else long reads
}//convert_gen_mer_cg

void Anchor::align_linear(vector<_uli> ta, _uli &sp, _uli &ep, _uli read_len1, _uli &remainder, _uli total_bits)
{
	//ta must be not reference but a copy because we need it later in program and are going to change it here
	//first entries of ta contain full 64-bit words
	_uli prev_sp = 0, prev_ep = 0;
	int size_ta = ta.size();
	vector<long int> lengths(size_ta, half_word);
	int read_rest = read_len1 - ((size_ta - 1) << two_in_32);

	lengths[size_ta - 1] = read_rest;

	long int bl = total_bits >> two_in_Size;//current 64-bit block
	_uli ashift = total_bits & mask[two_in_Size];
	sp = 0 ;
	ep = 0;
	int full_remainder = read_len1 - (total_bits >> 1);//read_len1 in bases, total_bits in bits
	
	_uli ta_rest = (ta[bl] >> (ashift));
	int rest_len = lengths[bl] - (ashift >> 1);
	remainder = rest_len;
	bl++;//points to next not fully used
	if(remainder < half_SEED){
		fill_ta_rest(ta, lengths, bl, remainder, ta_rest, size_ta);
		rest_len = remainder;
	}

	while(true){
		if(sp >= ep){

			_uli ta_seed = ta_rest & MASK_SEED;
			ta_rest >>= SEED;
			rest_len -= half_SEED;
			full_remainder -= half_SEED;
			remainder -= half_SEED;
			retreive_seed(ta_seed, sp, ep);


		}//if
		exact_match(ta_rest, sp, ep, rest_len, remainder, prev_sp, prev_ep);
		int mapped = rest_len - remainder;
		if(mapped == half_word)
			ta_rest = 0;
		else
			ta_rest >>= (mapped << 1);//mapped in bases
		full_remainder -= mapped;
		rest_len -= mapped;
		if(sp == ep){//mapped uniquelly
			break;
		}
		if(sp > ep){//found mism, for k=2, it means, do next shift
			break;
		}
		if(full_remainder == 0)
			break;
		if(remainder < half_SEED){
			fill_ta_rest(ta, lengths, bl, remainder, ta_rest, size_ta);
			rest_len = remainder;
		}

	}//while
	remainder = full_remainder;
	if((sp > ep) && (read_len1 - remainder - total_bits + 1 >= BIG_SEED)){
		sp = prev_sp;
		ep = prev_ep;
		remainder++;
	}
}//align_linear for long reads
void Anchor::fill_ta_rest(vector<_uli> &ta, vector<long int> &lengths, long int &bl, _uli &remainder, _uli &ta_rest, int &asize)
{

	long int unfilled = SIZE - (remainder << 1);//remainder is in bases
	while((bl < asize) && (unfilled > 0)){
		int next_len = lengths[bl];
		_uli ta_next = ta[bl];
		ta_next <<= (remainder << 1);
		ta_rest |= ta_next;
		if((next_len << 1) <= unfilled ){
			unfilled -= (next_len << 1);
			remainder += next_len;
			lengths[bl] = 0;
			ta[bl] >>= (next_len << 1);
			bl++;
		}
		else{
			remainder += (unfilled >> 1);
			ta[bl] >>= unfilled;
			lengths[bl] -= (unfilled >> 1);
			unfilled = 0;
		}
		
	}//while
}//fill_ta_rest called from align_linear for long reads

//is called only on case when sp = ep
void Anchor::map_pos(const _uli &sp, const _uli &ep, const vector<_uli> &ta_read1, const _uli &read_len1, const vector<_uli> &c_read1, 
		const bool &forward, const _uli &remainder, Mates &mates_pairs,
		const vector<_uli> &ns_read1, bool map_mism)
{
	_uli start1 = sp;
	_uli end1 = ep + 1;//must be so, in case sp = ep , don;t do twice
	int asize = ta_read1.size();

	long int total = ep - sp + 1;
	if(total > sp_ep_range){
		return;
	}

	while(start1 < end1){


			long int pos1 = find_gen_pos(start1, remainder);
			if((pos1 > 0) && ((pos1 + read_len1) <= binary_total)){ 
				bool cg_mism1 = false;
				vector<_uli> ta_check(asize, 0);
				int mism1 = 0;
				
				vector<_uli> cg_gen_mer;
				mism1 = cg_match(read_len1, pos1, ta_read1, c_read1, ta_check, forward, cg_mism1, cg_gen_mer);

					vector<_uli> ta_gen_mer;				
				if(mism1 <= mism_thresh ){
					//ta-check: need only for partially mapped reads

					ta_match(read_len1, pos1, ta_read1, ta_check, ns_read1, mism1, ta_gen_mer);
				}

				if(mism1 <= mism_thresh){
					int ascore =  read_len1 - mism1 ;//this function is called only from entire-length sp=ep

					_ui clip_left = 0, clip_right = 0; //from original read, not reverse read
					int cur_local_score = 0;
					string md("");
					calculate_MD(ta_gen_mer, cg_gen_mer, ta_check, read_len1, md , clip_left, clip_right, cur_local_score);
					//reset score: number of matched bases - (mism * mism_penalty)
					ascore = (read_len1 - mism1) + (mism1 * mism_penalty);
					stringstream ss_cigar;
					ss_cigar << read_len1 << 'M';
					string cigar1 = ss_cigar.str();
					update_singles(mates_pairs, pos1, mism1, map_mism, md, 0, 0, ascore, read_len1, cigar1, 0);		
				}//both reads passed cg-check
		
			}//if pos is within the genome size

		start1++;
	}//for forward read's list
}//map_pos

void pair_reads_singles_long(_uli read_id, string aread1,
				Stat &astat, _uli read_len1, _uli st1, _uli en1, Mates &mates_pairs , Anchor &anc_ct)
{
	//initialization
	anc_ct.init_params(read_len1, mates_pairs);


	unsigned long int one = 1, mism1, ta_lmer1, cg_lmer1;
	unsigned long int asize,   pos1, dist, i;
	char strand;

	vector<unsigned long int> c_compl_read, c_read1;//cg-representation
	vector<unsigned long int>  ta_read1, ta_compl_read;//ta-representation
	vector<_uli> ns_read1, ns_compl_read;

	unsigned long int start1, end1;

	//get ta-read, cg-read, ns-read: vectors must be empty, we use push_back here
	vector<_uli> acgt_read;
	vector<_uli> acgt_rc_read;
	long int num_ns = anc_ct.convert_lmer_ta_cg(aread1, ta_read1, c_read1, ns_read1, acgt_read, acgt_rc_read, read_len1, mates_pairs, ta_compl_read, c_compl_read, ns_compl_read);
	if(num_ns > anc_ct.mism_thresh)
		return;

	_uli min_shift = (anc_ct.SEED << 1);
	if(anc_ct.MIN_READ < read_len1)
		min_shift = (read_len1 - anc_ct.MIN_READ) << 1;
	_uli cur_shift_delta = min(anc_ct.shift_delta, min_shift);//num_shifts calculated from this

	_uli anew_read_len = read_len1 - anc_ct.SEED;
	if(read_len1 >= anc_ct.MIN_READ )
		anew_read_len = read_len1 - anc_ct.MIN_READ;//was SEED, but we don't want short reads < 48
	_uli num_shifts = min((anew_read_len >> anc_ct.two_in_shift_delta) + 1, anc_ct.max_num_shifts);
	
	//try to match read exactly (fully or partially)
	vector<Triple> sp_ep;

	_uli total_bits = 0;
	
	for(_uli y = 0; y < num_shifts; y++){
		_uli sp = 0, ep = 0; //retreive_seed shifts these values: have to be 0	
	
		_uli remainder = read_len1;
		if(anc_ct.A_rich == false)
			anc_ct.align_linear(acgt_read, sp, ep, read_len1, remainder, total_bits);//only ta-alignment
		else
			anc_ct.align_linear(acgt_rc_read, sp, ep, read_len1, remainder, total_bits);//only ta-alignment

		if(sp == ep){//mapped exactly, even if sp is the same, remainders might be different, hence positions different

			//process forward read1 
			//even though we align it here for sp = ep, we need the alignment information 
			//for local and indel alignments later one, so we are going to re-align to these positions sp=ep
			_ui shift_bases = total_bits >> 1;
			Triple atriple(sp, ep, remainder, shift_bases);	
			sp_ep.push_back(atriple);			
			
			bool forward = true;
			if(anc_ct.A_rich == false){
				
				anc_ct.map_pos(sp, ep, ta_read1, read_len1, c_read1, 
					forward,  remainder, mates_pairs, ns_read1, false);
			}
			else{
				anc_ct.map_pos(sp, ep, ta_compl_read, read_len1, c_compl_read, 
					forward,  remainder, mates_pairs, ns_compl_read, false);			
			}
	
		}//if mapped exactly to forward strand
		else if(sp < ep){
			_ui shift_bases = total_bits >> 1;
			Triple atriple(sp, ep, remainder, shift_bases);	
			sp_ep.push_back(atriple);
		}//else

		//if a read is mapped with 0, then at other location it will map with mism > 0, and we'll select best anyway
		//if a read is not mapped with 0 so far, then min mism is 1, need to make sure it is not ambiguous, if it is amb and mism 1, then this is best posible alignment
		//third condition is for ambiguous that is aligned perfectly
	
		//can read be mapped with 0 mism at one location and then with 0 mismatches at another location? (
		//because in this case, it will be considered ambiguous (and we'll miss ambiguous)
		//if read maps to two locations with 0 mismatches (C or not C option), then sp < ep, and mism > 0
		//If read can be mapped with 0, then condition mism = 0 will occur before mism = 1
			if((mates_pairs.mism1 == 0) || ((mates_pairs.mark_unique > 1) && (mates_pairs.mism1 == 1))){
				//set FLAG for SAM format
				if(mates_pairs.A_rich && mates_pairs.pos1 >= anc_ct.orig_genome_size)
					mates_pairs.flag = 16;
				else if( (!mates_pairs.A_rich) && mates_pairs.pos1 < anc_ct.orig_genome_size )
					mates_pairs.flag = 16;
				else
					mates_pairs.flag = 0;
							return;
			}
			if((total_bits == 0) && (remainder == 0) && (sp < ep) && (anc_ct.allow_C_to_T == 1)){
				mates_pairs.mark_unique = 2;//ambiguous 
				mates_pairs.mism1 = 0;
				return;
			}

		total_bits += cur_shift_delta ;
		if(total_bits > ((read_len1 - anc_ct.MIN_READ ) << 1))
			break;
	}//for y


	//NOVA:
	//set mates_pairs informatin back, so we could re-align sp=ep reads and build MD string for local alignment
	mates_pairs.mism1 = MAX_MISM;
	mates_pairs.mark_unique = 0;//treat it as unmapped yet


	//otherwise, align full-length to positions stored at sp_ep
	vector<AlignInfo> gen_pos;
	anc_ct.collect_gen_pos(sp_ep, gen_pos, read_len1);//gen pos are calculated fro sp,ep, then sorted, then only distinct positions are selected in gen_pos

	//Entire-length alignment to all collected distinct sorted positions
	bool forward = true;
	if(anc_ct.A_rich == false){
		anc_ct.map_entire_length(gen_pos, ta_read1, read_len1, c_read1, 
			forward,   mates_pairs, ns_read1, false);
	}
	else{
		anc_ct.map_entire_length(gen_pos, ta_compl_read, read_len1, c_compl_read, 
			forward,  mates_pairs, ns_compl_read, false);	
	}

	//if not mapped with entire-length alignment, map with local
	if(mates_pairs.mark_unique == 0 && anc_ct.is_local_alignment){
		if(anc_ct.A_rich == false){
			anc_ct.map_local(gen_pos, ta_read1, read_len1, c_read1, 
				forward, mates_pairs, ns_read1, false);
		}
		else{
			anc_ct.map_local(gen_pos, ta_compl_read, read_len1, c_compl_read, 
				forward, mates_pairs, ns_compl_read, false);	
		}	
	}//if not mapped with entire-length alignment, map with local

	//if not mapped with entire-length alignment, and did not map with local alignment
	//Do alignments with indels
	if( (read_len1 >= 75) && anc_ct.is_indel_alignment){
		if(anc_ct.A_rich == false){
			anc_ct.map_indels(gen_pos, ta_read1, read_len1, c_read1, 
				forward,  mates_pairs, ns_read1, false);
		}
		else{
			anc_ct.map_indels(gen_pos, ta_compl_read, read_len1, c_compl_read, 
				forward, mates_pairs, ns_compl_read, false);	
		}	
	}//if not mapped with entire-length alignment, and did not map with local alignment

	//set FLAG for SAM format
	if(mates_pairs.A_rich && mates_pairs.pos1 >= anc_ct.orig_genome_size)
		mates_pairs.flag = 16;
	else if( (!mates_pairs.A_rich) && mates_pairs.pos1 < anc_ct.orig_genome_size )
		mates_pairs.flag = 16;
	else
		mates_pairs.flag = 0;
}//pair_reads_singles_long new



//PAIRED_END reads
void pair_reads(_uli read_id, string aread1,
				Stat &astat, _uli read_len1, _uli st1, _uli en1, string read_name, string qual,
				string aread2, _uli read_len2, _uli st2, _uli en2 , string read_name2, string qual2, Anchor &anc_ct)
{
	//initialization
	Mates mates_pairs1, mates_pairs2;
	
	anc_ct.A_rich = false; 
	if(read_len1 >= anc_ct.MIN_READ)
		pair_reads_singles_long(read_id, aread1, astat, read_len1, st1, en1, mates_pairs1, anc_ct);
	
	anc_ct.A_rich = true;
	if(read_len2 >= anc_ct.MIN_READ){
		pair_reads_singles_long(read_id, aread2, astat, read_len2, st2, en2, mates_pairs2, anc_ct);

	}

	anc_ct.A_rich = false; 
	anc_ct.print_pairs(read_id, aread1, read_len1, st1, en1, read_name, qual, mates_pairs1, aread2, read_len2, st2, en2, read_name2, qual2, mates_pairs2, astat);

}//pair_reads

//END PAIRED END READS

void Anchor::ta_match(const _uli &read_len1, const _uli &pos1, const vector<_uli> &ta_read1, vector<_uli> &ta_check, 
	const vector<_uli> &ns_read1, int &mism1,  vector<_uli> &ta_gen_mer1)
{
	//ta_check must be declared and initialized with zeros outside this function before it's call

	convert_gen_mer_ta(ta_gen_mer1, pos1, read_len1);
	int asize = ta_read1.size();
	mism1 = 0;
	for(int j =0; j < asize; j++){
		ta_check[j] |= (ta_read1[j] ^ ta_gen_mer1[j]) | ns_read1[j];
		mism1 += count_mism(ta_check[j]);
		if(is_only_entire && (mism1 > mism_thresh))
			  break;
	}

}//ta_match long

int Anchor::cg_match(const _uli &read_len1, const _uli &pos1, const vector<_uli> &ta_read1, const vector<_uli> &c_read1,
	vector<_uli> &ta_check, const bool &forward, bool &cg_mism1,
	vector<_uli> &cg_gen_mer1)
{
	int asize = ta_read1.size();


	convert_gen_mer_cg(cg_gen_mer1, pos1, read_len1);

	//since in BRAT-nova we map all reads to positive and negative strand simaltaneously,
	//forward strand is always true
	//Moreover, rev-complements of 3'-mates map with T->C mismatches (to positive and negative strands matching strands of their 5'mates)
	count_ag_mism(ta_read1, c_read1, cg_gen_mer1, ta_check);

	_uli mism1 = 0;	
	for(int j = 0; j < asize; j++){
		if(allow_C_to_T == 0)
			ta_check[j] |= (c_read1[j] ^ (c_read1[j] & cg_gen_mer1[j])) ;
		else{
			_uli c_only_read = ta_read1[j] & c_read1[j];
			_uli g_only_read = (~(ta_read1[j])) & c_read1[j];
			ta_check[j] |= g_only_read ^ (g_only_read & cg_gen_mer1[j]);//checked up to here NOVA

		}//allow_C_to_T is true
	  mism1 += count_mism(ta_check[j]);
	  if(mism1 > mism_thresh){
		  cg_mism1 = true;
		  if(is_only_entire)
			  break;
	  }
	}//for j
	return mism1;
}//cg_match long 

void Anchor::count_ag_mism(const vector<_uli> &ta_read1, const vector<_uli> &c_read1, const vector<_uli> &cg_gen_mer1,
	vector<_uli> &ta_check)
{
	int asize = ta_read1.size();
	for(int j = 0; j < asize; j++)
		ta_check[j] |= (~(ta_read1[j] | c_read1[j])) & cg_gen_mer1[j];

}//count_ag_mism



void Anchor::count_tc_mism(const vector<_uli> &ta_read1, const vector<_uli> &c_read1, const vector<_uli> &cg_gen_mer1,
	vector<_uli> &ta_check)
{
	int asize = ta_read1.size();
	for(int j = 0; j < asize; j++)
		ta_check[j] |= (ta_read1[j] &(~c_read1[j])) & cg_gen_mer1[j];
}//count_ag_mism

//end BISULFITE reads FIRST mates

void Anchor::update_singles(Mates &mates_pairs, 
	const _uli &pos1, const _usi &mism1, bool map_mism, const string &md, 
	const _ui &clip_left, const _ui &clip_right, const int & cur_local_score, const _ui & read_len, const string &cigar, _usi type_map)
{
//NOVA: if type_map = 0 then entire-length, if type-map = 1, then local alignment, and 2 for indels
	if(mates_pairs.mark_unique == 0){
		mates_pairs.mark_unique++;
		mates_pairs.pos1 = pos1;
		mates_pairs.read_len = read_len;
		mates_pairs.mism1 = mism1;
		mates_pairs.A_rich = A_rich;
		mates_pairs.type_alignment = type_map;
		mates_pairs.local_score1 = cur_local_score;		
			
		if(!md.empty()){
			mates_pairs.MD1 = md;
			mates_pairs.clip_left1 = clip_left;
			mates_pairs.clip_right1 = clip_right;
			mates_pairs.cigar1 = cigar;
		}
	}//if
	else if(mates_pairs.pos1 != pos1 ){//|| mates_pairs.A_rich != A_rich){
		if(cur_local_score > mates_pairs.local_score1){
			//unique
			mates_pairs.mark_unique = 1;
			mates_pairs.pos1 = pos1;
			mates_pairs.read_len = read_len;		
			mates_pairs.mism1 = mism1;			
			mates_pairs.A_rich = A_rich;
			mates_pairs.type_alignment = type_map;
			mates_pairs.local_score1 = cur_local_score;
			if(!md.empty()){
				mates_pairs.MD1 = md;
				mates_pairs.clip_left1 = clip_left;
				mates_pairs.clip_right1 = clip_right;
				mates_pairs.cigar1 = cigar;
			}
		}
		else if(cur_local_score == mates_pairs.local_score1){
			mates_pairs.mark_unique++;
			}
	}//else	
	else if(cur_local_score > mates_pairs.local_score1){
                        //unique
                        mates_pairs.mark_unique = 1;
                        mates_pairs.pos1 = pos1;
                        mates_pairs.read_len = read_len;
                        mates_pairs.mism1 = mism1;
                        mates_pairs.A_rich = A_rich;
                        mates_pairs.type_alignment = type_map;
                        mates_pairs.local_score1 = cur_local_score;
                        if(!md.empty()){
                                mates_pairs.MD1 = md;
                                mates_pairs.clip_left1 = clip_left;
                                mates_pairs.clip_right1 = clip_right;
                                mates_pairs.cigar1 = cigar;
                        }
                }
	
	
}//update_singles


long int Anchor::convert_lmer_ta_cg(string lmer, vector<_uli> &ta_read, vector<_uli> &cg_read, vector<_uli> &ns_read, 
	vector<_uli> &full, vector<_uli> &full_rc, _ui read_len, Mates &mates_pairs, vector<_uli> &ta_compl_read, vector<_uli> &cg_compl_read, vector<_uli> &ns_compl_read){

	long int num_ns = 0;

	int asize = (read_len >> two_in_Size) + ((read_len & mask[two_in_Size]) > 0 ? 1 : 0) ;
	int full_size = (read_len >> two_in_32) + ((read_len & mask[two_in_32]) > 0 ? 1 : 0);

	full.resize(full_size, 0);
	full_rc.resize(full_size, 0);

	ta_read.resize(asize, 0);
	cg_read.resize(asize, 0);
	ns_read.resize(asize, 0);
		
	ta_compl_read.resize(asize, 0);
	cg_compl_read.resize(asize, 0);
	ns_compl_read.resize(asize, 0);

	int full_ind = 0;
	int ind = 0;
	int i = read_len - 1;
	int count_char = 0;
	while(i >= 0 ){ //Take reverse of read and convert C->T
		_uli ta = 0, cg = 0, ns = 0;

		for(int j = 0; (j < SIZE) && (i >= 0 ); j++){
			ta <<= 1;
			ns <<= 1;
			cg <<= 1;
			full[full_ind] <<= 2;

			if(lmer[i] == 'T' || lmer[i] == 't') {
				ta |= 1;
				full[full_ind] |= 3;
			}
			else if(lmer[i] == 'C' || lmer[i] == 'c'){
				ta |= 1;
				cg |= 1;
				full[full_ind] |= 3;//Convert C to T

			}
			else if(lmer[i] == 'G' || lmer[i] == 'g'){
				cg |= 1;

				full[full_ind] |= 2;
			}
			else if(lmer[i] == 'N' || lmer[i] == 'n'){
				ns |= 1;
				num_ns++;
			}	
			i--;
			count_char++;
			if((count_char & mask[two_in_32]) == 0)
				full_ind++;
		}//for j
		ta_read[ind] = ta;
		cg_read[ind] = cg;
		ns_read[ind] = ns;
		ind++;

	}//while

	int rest_bits = read_len & mask[two_in_Size];//= (read_len1 mod SIZE) 
	if(rest_bits > 0){
		shift_back(ta_read, asize, rest_bits);
		shift_back(cg_read, asize, rest_bits);
		shift_back(ns_read, asize, rest_bits);
	}
	reverse(ta_read.begin(), ta_read.end());
	reverse(cg_read.begin(), cg_read.end());
	reverse(ns_read.begin(), ns_read.end());

	
	//take complement of original string = reverse(rev-complement of ogiginal string), and change Cs to Ts for full_rc
	i = 0;
	full_ind = 0;
	ind = 0;
	count_char = 0;
	while(i < read_len ){ //
		_uli ta = 0, cg = 0, ns = 0;

		for(int j = 0; (j < SIZE) && (i < read_len ); j++){
			ta <<= 1;
			ns <<= 1;
			cg <<= 1;
			full_rc[full_ind] <<= 2;

			if(lmer[i] == 'A' || lmer[i] == 'a') {//convert to T
				ta |= 1;
				full_rc[full_ind] |= 3;
			}
			else if(lmer[i] == 'C' || lmer[i] == 'c'){//convert to G
				
				cg |= 1;
				full_rc[full_ind] |= 2;//

			}
			else if(lmer[i] == 'G' || lmer[i] == 'g'){//convert to C, then to T for full_rc genome only
				cg |= 1;
				ta |= 1;//for C
				full_rc[full_ind] |= 3;
			}
			else if(lmer[i] == 'N' || lmer[i] == 'n'){//complement of N is T
				ns |= 1;
				ta |= 1;
				full_rc[full_ind] |= 3;
			}	
			i++;
			count_char++;
			if((count_char & mask[two_in_32]) == 0)
				full_ind++;
		}//for j
		ta_compl_read[ind] = ta;
		cg_compl_read[ind] = cg;
		ns_compl_read[ind] = ns;
		ind++;

	}//while

	rest_bits = read_len & mask[two_in_Size];//= (read_len1 mod SIZE) 
	if(rest_bits > 0){
		shift_back(ta_compl_read, asize, rest_bits);
		shift_back(cg_compl_read, asize, rest_bits);
		shift_back(ns_compl_read, asize, rest_bits);
	}
	reverse(ta_compl_read.begin(), ta_compl_read.end());
	reverse(cg_compl_read.begin(), cg_compl_read.end());
	reverse(ns_compl_read.begin(), ns_compl_read.end());
//END of COMPLEMENT
	rest_bits = (read_len & mask[two_in_32]) << 1;
	if(rest_bits > 0){
		shift_back(full, full_size, rest_bits);
		shift_back(full_rc, full_size, rest_bits);
	}
	reverse(full.begin(), full.end());
	reverse(full_rc.begin(), full_rc.end());

	return num_ns;
}//convert_ta_cg
void Anchor::shift_back(vector<_uli> &ta_read, int asize, const _ui &rest_bits)
{
	asize--;//last index, asize is size of ta
	//rest_bits is rest of bits of a read, when taking multiple of 64
	for(int i = asize; i > 0; i--){
		_uli prev_ta = ta_read[i - 1] << rest_bits;
		ta_read[i] |= prev_ta;
		ta_read[i-1] >>= SIZE - rest_bits;

	}//for i

	//now full first entries correspond to the end of read
}//shift_back

void Anchor::revcompl_lmer_ta_cg(const vector<_uli> &ta_read, const vector<_uli> &cg_read, const vector<_uli> &ns_read, const _ui &read_len,
	vector<_uli> &ta_rc_read, vector<_uli> &cg_rc_read, vector<_uli> &ns_rc_read)
{	
/*
1111111011111110000101101101101111101101100101111111100000100010
000100010101010000010101010101010101
*/
	int len1 = read_len & mask[two_in_Size];//= (read_len1 mod SIZE) 
	if(len1 == 0)
		len1 = SIZE;
	_uli one = 1;
	unsigned long int MASK_READ;
	if(len1 < SIZE){
		MASK_READ = ( one << len1) - 1;
	}
	else
		MASK_READ = 0xFFFFFFFFFFFFFFFF;

	int asize = ta_read.size();
	ta_rc_read.resize(asize);
	cg_rc_read.resize(asize);
	ns_rc_read.resize(asize);

	int stop = asize - 1;
	for(int i = 0; i < stop; i++){
		ta_rc_read[i] = ~Reverse(ta_read[i], SIZE);
		cg_rc_read[i] =  Reverse(cg_read[i], SIZE);//shows positions of Gs in reverse strand
		ns_rc_read[i] = Reverse(ns_read[i], SIZE);
	}
	ta_rc_read[stop] = MASK_READ & (~Reverse(ta_read[stop], len1));
	cg_rc_read[stop] =  MASK_READ & Reverse(cg_read[stop], len1);//shows positions of Gs in reverse strand
	ns_rc_read[stop] = Reverse(ns_read[stop], len1) & MASK_READ;

	//shift to make rev-compl in same order, i.e. first entry contain last bits of rev-compl
	shift_back(ta_rc_read, asize, len1);
	shift_back(cg_rc_read, asize, len1);
	shift_back(ns_rc_read, asize, len1);

	reverse(ta_rc_read.begin(), ta_rc_read.end());
	reverse(cg_rc_read.begin(), cg_rc_read.end());
	reverse(ns_rc_read.begin(), ns_rc_read.end());

}//revcompl_lmer_ta_cg

int Anchor::parse_options(int argc, char* argv[])
{
	int res = 0;

	long int i, j;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-')
			return -1;
		switch(argv[i][1]){
		case 'P': prefix = argv[++i]; break;
		case 's': reads_file1 = argv[++i]; break;
		case '1': reads_file1 = argv[++i]; break;
		case '2': reads_file2 = argv[++i]; break;
		case 'o': output_pairs = argv[++i]; break;
		case 'p': is_pe = "1"; break;
		case 'i': min_insert = atoi(argv[++i]); break;
		case 'a': max_insert = atoi(argv[++i]); break;
		case 'A': A_rich = true; break;
		case 'C': allow_C_to_T = 1; break;
		case 'K': max_num_shifts = atoi(argv[++i]); break;
		case 'W': four_strands = true; break;
		case 'm': mism_penalty = atoi(argv[++i]); break;//can keep this as mism-penalty
		case 'L': is_local_alignment = true; break;
		case 'E' : is_only_entire = true; break;
		case 'G': is_indel_alignment = true; break;
		case 'g': ((argv[i][2] == 'o') ? (gap_open_penalty = atoi(argv[++i])) : ((argv[i][2] == 'e') ? (gap_cont_penalty = atoi(argv[++i])) : (gap_threshold = atoi(argv[++i])))); break;
		case 'l' : alignment_length = atoi(argv[++i]); break;//need to divide by 100 to get proportion of read length instead of percentage
		case 'q' : alignment_score = atoi(argv[++i]); break;//need to divide by 100
		case 't' : sp_ep_range = atoi(argv[++i]); break;//max allowable is 1000 set by the program
		case 'h': usage(); break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl; usage(); space_release(); exit(0); break;
//		case 'F': first_bits = atoi(argv[++i]); break; //depricated
//		case 'D': max_shift_delta = atoi(argv[++i]); break;
		}
	}//for i
	mism_penalty = 0 - (abs(mism_penalty));
	gap_open_penalty = 0 - (abs(gap_open_penalty));
	gap_cont_penalty = 0 - (abs(gap_cont_penalty));
	alignment_length = alignment_length / (100.0);
	alignment_score = alignment_score / (100.0);
	if(is_only_entire){
		is_local_alignment = false;
		is_indel_alignment = false;
	}
//Depricated:
//		case 'u': unmapped = true; break;
//		case 'k': all_shifts = atoi(argv[++i]); break;
//		case 'b': is_bs = "1"; break;
//		case 'M': out_amb = true; break;
//		case 'f': BIG_SEED = atoi(argv[++i]); break;

	return res;
}//parse_opt()
void Anchor::read_rev_complement(string &aread, const _uli  &read_len){
	reverse(aread.begin(), aread.end());

	for(_uli i = 0; i < read_len; i++){
		if(aread[i] == 'T' || aread[i] == 't')
			aread[i] = 'A';
		else if(aread[i] == 'A' || aread[i] == 'a')
			aread[i] = 'T';
		else if(aread[i] == 'G' || aread[i] == 'g')
			aread[i] = 'C';
		else if(aread[i] == 'C' || aread[i] == 'c')
			aread[i] = 'G';
		else
			aread[i] = 'N';
	}//for i
}//read_rev_complement

void Anchor::print_pairs(const _uli &read_id, const string &aread1, const _uli &read_len1, const _uli &st1,
		const _uli & en1, const string & read_name, const string &qual, Mates &mates_pairs1, const string &aread2, const _uli &read_len2, 
		const _uli &st2, const _uli &en2, const string & read_name2, const string &qual2, Mates &mates_pairs2, Stat &astat)
{


	long int orig_start1, orig_start2;

	if((mates_pairs1.mark_unique == 1) && (mates_pairs2.mark_unique == 1)){

		//determine chrom id for this read
		long int chrom_id1 = 0, chrom_id2 = 0;

		//find the first start of chrom that is greater than pos
		long int pos1 = mates_pairs1.pos1; 
		long int asize = size_chrom.size();
		for(long int i = 0; i < asize; i++){
			if(pos1 >= starts_size_chrom[i])
				chrom_id1++;
			else
				break;
		}//for i
		chrom_id1--;
//		pos1 = pos1 - starts_size_chrom[chrom_id1];

		long int pos2 = mates_pairs2.pos1;//must by signed
		for(long int i = 0; i < asize; i++){
			if(pos2 >= starts_size_chrom[i])
				chrom_id2++;
			else
				break;
		}//for i
		chrom_id2--;
//		pos2 = pos2 - starts_size_chrom[chrom_id2];

		if(chrom_id1 == chrom_id2){ //same chromosomes

			/* In pair_reads_singles_long, we map a mate according to its A-rich
			If A-rich = False, then original read is mapped
			If A-rich = True, then reverse-complement of the read is mapped
			Whatever mapping position, it is stored in mates_pairs.pos1
			Strands are depricated in Mates
			Both mates must map to the same strand (one mate = original read, another is rev-complement of orig)
			*/

			if(((pos1 < orig_genome_size) && (pos2 < orig_genome_size)) || ((pos1 >= orig_genome_size) && (pos2 >= orig_genome_size))){
				//check if insert is correct
				long int dist = 0;
				if(pos2 >= pos1)
					dist = pos2 - pos1;
				else
					dist = pos1 - pos2;

				//no need to adjust for mates1 and A-rich: it is already taken care of by strand = '-'
				if(dist >= min_insert && dist <= max_insert ){
					astat.mapped_pairs++;
					//note flag is already set to 16 if sequence is rev-complemented inside pair_reads_singles_long
					//also note that mates_pairs1 is for 5'-mate
					mates_pairs1.flag += 3 + 64 + ((mates_pairs1.pos1 >= orig_genome_size) ? 32 : 0);
					//note that mates_pairs2 are for 3'-mate
					mates_pairs2.flag += 3 + 128 + ((mates_pairs2.pos1 < orig_genome_size) ? 32 : 0);//16 is already set
					print_singles(out_pairs, false, read_id, aread1, st1, en1, mates_pairs1, astat, read_len1, read_name, qual);
					print_singles(out_pairs, true, read_id, aread2, st2, en2, mates_pairs2, astat, read_len2, read_name2, qual2);

				}//correct insert
				else{
					astat.broken_wrong_insert++;
				}
			}//if correct strands
			else{
				astat.broken_inverse++;
			}
		}//if correct_chrom
		else{
			astat.broken_dif_chroms++;
		}
	}//if both mates are mapped uniquely
	else{
		if((mates_pairs1.mark_unique == 1 && mates_pairs2.mark_unique == 0) ||
			(mates_pairs1.mark_unique == 0 && mates_pairs2.mark_unique == 1)){
				astat.broken_mate_unmapped++;
		}
		else if((mates_pairs1.mark_unique == 1 && mates_pairs2.mark_unique > 1) ||
			(mates_pairs1.mark_unique > 1 && mates_pairs2.mark_unique == 1)){
			astat.broken_mate_amb++;
		}
		else if(mates_pairs1.mark_unique == 0 && mates_pairs2.mark_unique == 0){
			astat.broken_both_unmapped++;
		}
		else
			astat.ambiguous++;//both reads are amb

		//set flag that mate is unmapped
		if(mates_pairs1.mark_unique == 1){
			//sequence rev-complemented is set in pair_reads_singles_long to flag = 16
			mates_pairs1.flag += 9 +  64;//mates_pairs1 refers to 5'-mate
			print_singles(out_pairs, false, read_id, aread1, st1, en1, mates_pairs1, astat, read_len1, read_name, qual);

		}
		else if(mates_pairs2.mark_unique == 1){
			//sequence rev-complemented is set in pair_reads_singles_long to flag = 16
			mates_pairs2.flag += 9 + 128;//mates_pairs2 refers to 3'-mate (hence, + 128)
			print_singles(out_pairs, true, read_id, aread2, st2, en2, mates_pairs2, astat, read_len2, read_name2, qual);
		}
	}
	
}//print_pairs

void Anchor::print_singles(ofstream &out, bool second_mates, const _uli &read_id, 
	const string & aread1, const _uli &st1, const _uli &en1,  Mates &mates_pairs, Stat &astat
	, const _uli &read_len, const string & read_name, const string &qual)
{

	long int orig_start1;
	if(mates_pairs.mark_unique == 1){
		long int chrom_id = 0;
		long int asize = size_chrom.size();
		//find the first start of chrom that is greater than pos
		long int pos1 = mates_pairs.pos1; 
		char strand = '-';
		if(pos1 >= orig_genome_size)//in the second half
			strand = '+';

		for(long int i = 0; i < asize; i++){
			if(pos1 >= starts_size_chrom[i])
				chrom_id++;
			else
				break;
		}//for i
		chrom_id--;
		pos1 = pos1 - starts_size_chrom[chrom_id];
		if(strand == '+'){
			pos1 = size_chrom[chrom_id] - pos1 - mates_pairs.read_len;//in the second half original chromosome direction is from right-to-left
			//since we need the leftmost point on the original chromosome
			if(mates_pairs.type_alignment == 2){ //type-alignment 0(entire), 1(local), 2(indels)
				//for deletion, we need adjust position by subtracting deleted bases
				string cigar = mates_pairs.cigar1;
				int len = cigar.length();
				bool is_deletion = false;
				int j = 0;
				for(j = 0; j < len; j++){
					if(cigar[j] == 'D')
					{
						is_deletion = true;
						break;
					}//if
				}//for
				if(is_deletion){
					//find the leftmost not digit: it must be true because we have indels in between two seeds
					j--;
					int sub_len = 0;
					while(j >= 0 && isdigit(cigar[j])){
						j--;
						sub_len++;
					}//while
					string sub_del = cigar.substr(j + 1, sub_len);
					istringstream istr(sub_del);
					int del_len;
					istr >> del_len;
					pos1 = pos1 - del_len;
					
				}//if
			}//if indel check if deletion
			
		}//if indels

		map<long int, string>::iterator afind = chroms_names.find(chrom_id);
		string chr_name = afind->second;
		if(strand == '+')
			orig_start1 = pos1 - st1;
		else
			orig_start1 = pos1 - en1;
		if(orig_start1 < 0)
			orig_start1 = pos1;
		

		string out_read = aread1;
		string out_qual = qual;
		int new_read_len = read_len;
		if(mates_pairs.type_alignment == 1 || mates_pairs.type_alignment == 2){ //type-alignment 0(entire), 1(local), 2(indels)
			//clip the original read: important: here keep read_len not mates_pairs.read_len, we need original read length
			new_read_len = read_len - mates_pairs.clip_left1 - mates_pairs.clip_right1;
			out_read = aread1.substr(mates_pairs.clip_left1, new_read_len);
			out_qual = qual.substr(mates_pairs.clip_left1, new_read_len);

			//DEBUG_HERE:

			//adjust position
			int adjust_pos =  mates_pairs.clip_left1;
			if(strand == '-'){
				adjust_pos = mates_pairs.clip_right1;
			}
			pos1 += adjust_pos;

		}//local

		
		if((mates_pairs.A_rich == true) && (strand == '+')){//this is not reverse-complement for negative strand
			//when counting methylation, if strand is +, then we count T->C toward positive strand's cytocines
			//if strand is -, then we count A->G toward negative strand's cytocines
			//case: rev-complement of 3'-mate is mapped to positive strand
			read_rev_complement(out_read, new_read_len);
			reverse(out_qual.begin(), out_qual.end());
			reverse_md(mates_pairs.MD1, mates_pairs.type_alignment);//for negative strand 
			reverse_cigar(mates_pairs.cigar1);
			/* In BRAT-bw, we swapped strands, but now, in BRAT-nova, strand is determined by pos >= genome_orig_size (if so, then strand = +)
			so, we do need to take reverse complement of the read if A-rich is true, but we do not swap strands
			strand = (strand == '+' ? strand = '-' : strand = '+');	
			*/
		}
		if((mates_pairs.A_rich == false) && (strand == '-')){//5'-mate from negative strand
			//in SAM format, reads must be ready to positive strand alignment

			read_rev_complement(out_read, new_read_len);
			reverse(out_qual.begin(), out_qual.end());
			reverse_md(mates_pairs.MD1, mates_pairs.type_alignment);//for negative strand 
			reverse_cigar(mates_pairs.cigar1);

		}

		/*//BRAT-bw output
		out << read_id << "\t" << out_read << "\t" 
			<< chr_name << "\t" << strand <<"\t"
			<< pos1 << "\t" << mates_pairs.mism1 << "\t" << orig_start1 << endl;
		*/
		
		//DEBUG only: Remove after testing
		if(mates_pairs.cigar1.empty() || mates_pairs.MD1.empty()){
			cout << "ERROR: cigar or MD strings are empty, read_id: " << read_id << "\t" << aread1 
				<< "\t" << chr_name << "\t" << strand << "\t" << (pos1 + 1) << endl;
		}

		//If XB:A:+, then count methylation toward positive (top) strand (C->C, and T->C)
		//If XB:A:-, then count methylation toward negative (bottom) strand (G->G, and A->G)
		out << read_name << "\t" << mates_pairs.flag  << "\t" << chr_name
			<< "\t" << (pos1 + 1) << "\t" << 255 << "\t" << mates_pairs.cigar1 
			<< "\t" << "*" << "\t" << 0  << "\t" << 0 << "\t" 
			<< out_read << "\t" << out_qual << "\t" << "AS:i:" << mates_pairs.local_score1 
			<< "\t" << "MD:Z:" << mates_pairs.MD1 
			<< "\t" << "XB:A:" << strand << "\t" << "XC:i:" << (orig_start1 + 1) << "\t" << endl;
	}//if

}//print_singles

void find_st_en(string &name, _uli &st, _uli &en){

	//find first character ?
	int i = 0, first_ind = -1, second_ind = -1;
	int asize = name.length();
	while(i < asize){
		if(name[i] == '?'){
			first_ind = i;
			break;
		}
		i++;
	}//while
	if(first_ind == -1)
	{
		st = 0;
		en = 0;
		return;
	}
	i = first_ind + 1;
	second_ind = asize;
	while(i < asize){
		if(name[i] == '?'){
			second_ind = i;
			break;
		}
		i++;
	}//while
	string sub_st = name.substr(first_ind + 1, second_ind - first_ind - 1 );
	istringstream istr1(sub_st);
	istr1 >> st;
	string sub_en = name.substr(second_ind + 1, asize - second_ind - 1);
	istringstream istr2(sub_en);
	istr2 >> en;
	name = name.substr(0, first_ind);//as it was originally

}//find_st_en

void map_reads(Anchor &anc_ct)
{
	//pair_reads_singles(unsigned long int read_id, string aread1, long int &contain_Ns, Stat &astat, short int len1, int t, int rs);
	_uli read_id = 0;
   Stat astat;
   string read, read2;
   _uli st1, en1, st2, en2;

	_uli read_len1, read_len2;
	if(anc_ct.is_pe == "1"){
		read_id = 0;

		ifstream in, in2;
		in.open(anc_ct.reads_file1.c_str(), ios::in);
		if(!in){
			cout << "can't open " << anc_ct.reads_file1 << " reads_file " << endl;
			anc_ct.space_release();
			exit(0);
		}
		in2.open(anc_ct.reads_file2.c_str(), ios::in);
		if(!in2){
			cout << "can't open " << anc_ct.reads_file2 << " reads_file2" << endl;
			anc_ct.space_release();
			exit(0);
		}
		string qual, qual2, id, id2, dum, dum2;
		getline(in, id);
		getline(in2, id2);
		
//		in >> read;
//		in2 >> read2;

		while(!in.eof()){
			//mate1
			istringstream istr1(id);
			string read_name;
			istr1 >> read_name;
			if(read_name[0] == '@')
				read_name = read_name.substr(1, read_name.length() - 1);
			//find st and en (attached to read's id after ? character)
			find_st_en(read_name, st1, en1);

			//mate2
			istringstream istr2(id2);
			string read_name2;
			istr2 >> read_name2;
			if(read_name2[0] == '@')
				read_name2 = read_name2.substr(1, read_name2.length() - 1);
			find_st_en(read_name2, st2, en2);

			//read the second line of fastq
			getline(in, read);
			getline(in2, read2);
	
			//read the third line of fastq
			getline(in, dum);
			getline(in2, dum2);

			//read the fourth line with quals of fastq
			getline(in, qual);
			getline(in2, qual2);

			read_len1 = read.length();
			read_len2 = read2.length();

			pair_reads(read_id, read, astat, read_len1, st1, en1, read_name, qual,
				read2, read_len2, st2, en2, read_name2, qual2 , anc_ct);

			read_id++;
			getline(in, id);
			getline(in2, id2);
		}//while

		cout << "Mapped uniquely pairs (pair count): " << astat.mapped_pairs << endl;
		cout << "Of them, with first mate being A-rich (from not original strands) (pair count): " << astat.mapped_first_Arich << endl;
		cout << "Broken pair: one mate is uniquely mapped, the other is unmapped (single-read count): " << astat.broken_mate_unmapped << endl;
		cout << "Broken pair: one mate is uniquely mapped, the other is ambiguos (single-read count): " << astat.broken_mate_amb << endl;
		cout << "Broken pair: mates mapped to different chromosomes (pair count): " << astat.broken_dif_chroms << endl;
		cout << "Broken pair: mates have wrong orientation (pair count): " << astat.broken_inverse << endl;
		cout << "Broken pair: mates have wrong inster size (pair count): " << astat.broken_wrong_insert << endl;
		cout << "Unmapped pair: both mates are unmapped (pair count): " << astat.unmapped << endl;
		cout << "Ambiguous pair: both mates are ambiguous (pair count): " << astat.ambiguous << endl;
		
	}//paired-end mapping
	else{

		ifstream in;
		in.open(anc_ct.reads_file1.c_str(), ios::in);
		if(!in){
			cout << "can't open " << anc_ct.reads_file1 << " reads_file (option -s) " << endl;
			anc_ct.space_release();
			exit(0);
		}
			string qual, dum, id;
			getline(in, id);
			read_id = 0;

			while(!in.eof()){

				istringstream istr1(id);
				string read_name;
				istr1 >> read_name;
				if(read_name[0] == '@')
					read_name = read_name.substr(1, read_name.length() - 1);
				//find st and en (attached to read's id after ? character)
				find_st_en(read_name, st1, en1);

				getline(in, read);
				getline(in, dum);
				getline(in, qual);

				read_len1 = read.length();
				Mates mates_pairs;	
				if(read_len1 >= anc_ct.MIN_READ ){
					pair_reads_singles_long(read_id, read, astat, read_len1, st1, en1, mates_pairs, anc_ct);
				}

				anc_ct.print_singles(anc_ct.out_pairs, false, read_id, read, st1, en1, mates_pairs, astat, read_len1, read_name, qual);	

if(mates_pairs.mark_unique > 1)
	astat.ambiguous++;
else if(mates_pairs.mark_unique == 0)
	astat.unmapped++;
else
	astat.mapped_pairs++;

				read_id++;

				getline(in, id); //in >> read;
			}//while
			in.close(); in.clear();
cout << "Mapped uniquely: " << astat.mapped_pairs << endl;
cout << "Unmapped: " << astat.unmapped << endl;
cout << "Ambiguous: " << astat.ambiguous << endl;

	}//else singles

}//map_reads


//NOVA
bool is_less_gen_pos(const AlignInfo &x, const AlignInfo &y){
	if(x.pos < y.pos)
		return true;
	else
		return false;

}



//from ta-check vector of bits
void Anchor::calculate_MD(const vector<_uli> &ta_gen_mer, const vector<_uli> &cg_gen_mer, const vector<_uli> &ta_check, 
	const _ui &read_len, string &md, _ui &clip_left, _ui &clip_right, int &cur_local_score){
	
	int prev_score = mism_penalty - 1;
	int prev_ind = 0;
	int prev_len = 0;
	int max_score = 0;
	int max_ind = 0;
	int max_len = 0;
	int prev_mism = 0;
	int max_mism = 0;


	int asize = ta_check.size();
	stringstream ss;
	_uli one = 1;
	int cur_score = 0;//number of matched exactly bases between two consecutive mismatches
	_uli total_bits = read_len - ((asize - 1) << two_in_Size);
	if(read_len <= SIZE)
		total_bits = read_len;

	for(int i = 0; i <= asize - 1; i++){
		_uli cur_word = ta_check[i];
		int stop = (i < asize - 1) ? SIZE : total_bits;
		for(int j = 0; j < stop; j++){
			_uli abit = cur_word & mask[1];
			cur_word >>= 1;

			if(abit == one){
				//mismatch is found
				//calculate total mapped bases from previous mismatch to now
				if(!(i == 0 && j == 0)){//not beginning of the MD string (no 0-scores are at the end of MD string)
					ss << cur_score;

					//*****************calculate score for current entry of array with value = cur_score
					if(prev_score + cur_score >= cur_score){
						prev_score += cur_score;
						//prev_ind does not change
						prev_len += cur_score;
					}
					else{
						prev_score = cur_score;
						prev_ind = ((i << two_in_Size) + j ) - cur_score;//index in original string
						prev_len = cur_score;
						prev_mism = 0;
					}//else

					//set correct max
					if(prev_score > max_score){
						max_score = prev_score;
						max_ind = prev_ind;
						max_len = prev_len;
						max_mism = prev_mism;
					}//if max
				}//if mismatch is not at the beginning of the string

				//find out the character in the genome
				_uli ta_cg_bits = (((ta_gen_mer[i] >> j) & mask[1]) << 1) | ((cg_gen_mer[i] >> j) & mask[1]);
				ss << (bits_to_char[ta_cg_bits]);

				//update mismatch change to the score
				if(prev_score + mism_penalty >= mism_penalty){
					prev_score += mism_penalty;
					//index does not change
					prev_len++;//1 mism
					prev_mism++;
				}
				else{
					prev_score = mism_penalty;
					prev_ind = (i << two_in_Size) + j;
					prev_len = 1;
					prev_mism = 1;
				}//eles
					//set correct max
				if(prev_score > max_score){
						max_score = prev_score;
						max_ind = prev_ind;
						max_len = prev_len;
						max_mism = prev_mism;
				}//if max

				cur_score = 0;
			}//if 
			else
				cur_score++;
		}//for j
		
	}//for i

	//if mismatch is not at the end of the read
	if(cur_score > 0){
		ss << cur_score;
					//*****************calculate score for current entry of array with value = cur_score
					if(prev_score + cur_score >= cur_score){
						prev_score += cur_score;
						//prev_ind does not change
						prev_len += cur_score;
					}
					else{
						prev_score = cur_score;
						prev_ind = read_len - cur_score;//index in original string
						prev_len = cur_score;
						prev_mism = 0;
					}//else

					//set correct max
					if(prev_score > max_score){
						max_score = prev_score;
						max_ind = prev_ind;
						max_len = prev_len;
						max_mism = prev_mism;
					}//if max
	}//if mismatch is not at the end of string
	
	md = ss.str();//direction is as in the original string (since we started from the end of the reverse of original string)

	clip_left = max_ind;
	clip_right = read_len - (clip_left + max_len);
	cur_local_score = (read_len - clip_left - clip_right) - max_mism;

	//dynamic programming is used to calculate overall max-score and corresponding max-sub-array,
	//hence, clip_left and clip_right.
	//cur_local_score is not dynamic-programming score
}//calculate_MD

void Anchor::reverse_simple_md(string &s){
	int size = s.length();
	stringstream ss;
       for(int i = size - 1; i >= 0; i--){
                if(isdigit(s[i])){
                        int j = i;
                        while(j >= 0 && isdigit(s[j])){
                                j--;
                        }
                        ss << s.substr(j+1, i - j);
                        i = j + 1;
                }//if
                else{
			ss << s[i];
                        
                }//char
        }
	s = ss.str();
}
void Anchor::reverse_md(string &s, const _usi & type_alignment){
	int size = s.length();
	_usi two = 2;
	bool is_deletion = false;
	int del_ind = -1;
	if(type_alignment == two){
		//check if there is a deletiong
		for(int i = 0; i < size; i++){
			if(s[i] == '^')
			{
				is_deletion = true;
				del_ind = i;
				break;
			}
		}//for
	}
	if(!is_deletion){
		reverse_simple_md(s);
	}//if entire-length or local alignments
	else{//deletion
		//find next digit
		int next_digit = del_ind + 1;
		for(int i = del_ind; i < size; i++){
			if(isdigit(s[i])){
				next_digit = i;
				break;
			}//if
		}//for
		string first = s.substr(0, del_ind);
		reverse_simple_md(first);
		if(!isdigit(first[0]))
			first = "0" + first;
		string second = s.substr(next_digit, size - next_digit);
		reverse_simple_md(second);
		string del_substr = s.substr(del_ind + 1, next_digit - del_ind - 1);
		read_rev_complement(del_substr, del_substr.length());
		stringstream ss;
		ss << second << '^' << del_substr << first;
		s = ss.str();
	}//indels
}//reverse_md

void Anchor::reverse_cigar(string &s){
	int size = s.length();
	stringstream ss;
	vector<string> chuncks;
	int i = 0, prev = 0;
	while(i < size){
		
		while(i < size && isdigit(s[i])){
			i++;
		}//while
		i++;
		chuncks.push_back(s.substr(prev, i - prev));
		prev = i;

	}//while
	size = chuncks.size();
      for(i = size - 1; i >= 0; i--){
            ss << chuncks[i] ;   
      }//for
	s = ss.str();
}
void Anchor::clip_md(string &s, const _ui &left, const _ui & right){
	if(left == 0 && right == 0)
		return;

	int size = s.length();
	_ui count = 0;
	int i = 0;
	if(left > 0){
	for(i = 0; i < size; i++){
		if(isdigit(s[i])){
			int j = i;
			stringstream ss;
			while(j < size && isdigit(s[j])){
				ss << s[j];
				j++;
			}
			istringstream istr(ss.str());
			int x;
			istr >> x;
			count += x;
			i = j - 1;
			//check if clip_left is reached
			if(count > left){
				//split matched bases into two parts: clipped and kept
				int dif = count - left;
				stringstream read_new_score;
				read_new_score << dif;
				string rem = (j < size) ? s.substr(j, size - j) : "";
				s = read_new_score.str() + rem;
				break;
			}//if count > left
		}//if
		else{
			count++;
		}//char
		if(count == left){
			int next = i + 1;
			s = s.substr(next, size - next);
			break;
		}//if count=left

	}//for
	
	}//if left > 0

	if(right > 0){
		count = 0;
		size = s.length();
		for(i = size - 1; i >= 0; i--){
			if(isdigit(s[i])){
				int j = i;
				stringstream ss;
				while(j >= 0 && isdigit(s[j])){                    
					j--;
				}//while
				ss << s.substr(j+1, i - j);
				istringstream istr(ss.str());
				int x;
				istr >> x;
				count += x;
				i = j + 1;
				//check if clip_right is reached
				if(count > right){
					//split matched bases into two parts: clipped and kept
					int dif = count - right;
					stringstream read_new_score;
					read_new_score << dif;
					string rem = (j >= 0) ? s.substr(0, i) : "";
					s = rem + read_new_score.str() ;
					break;
				}//if count > left
            }//if
            else{
				count++;
            }//char
            if(count == right){
				s = s.substr(0, i);       
				break;
			}//if count = right
		}//for
        
	}//if right > 0
	size = s.length();
	//check the ends 
	if(s[0] == '0')
		s = s.substr(1, size - 1);
	if(size >= 2){
		if((s[size - 1] == '0') &&(!(isdigit(s[size-2]))))
			s = s.substr(0, size - 1);
	}//if length allows to look at len-2 pos
}//clip_md

void Anchor::collect_gen_pos(vector<Triple> &sp_ep, vector<AlignInfo> &gen_pos, const _uli & read_len1){
	
	int sp_ep_size = sp_ep.size();
	if(sp_ep_size == 0)
		return;
	//collect all valid positions
	gen_pos.reserve(sp_ep_size << 10); //multiply by 1024 closest to total threshold 1000

	for(_ui i = 0; i < sp_ep_size; i++){
		_uli start1 = sp_ep[i].sp;
		_uli end1 = sp_ep[i].ep + 1;
		_ui remainder = sp_ep[i].remainder;
		_ui shift_bases = sp_ep[i].shift;
		_uli total = end1 - start1;
		if(total <= sp_ep_range){

			while(start1 < end1){

					long int pos1 = find_gen_pos(start1, remainder);
					if((pos1 >= 0) && ((pos1 + read_len1) <= binary_total)){ 
						AlignInfo al(pos1, remainder, shift_bases);
						gen_pos.push_back(al);

					}

				start1++;
			}//while
		}//if sp-ep ranges is small
	}//for all (sp, ep) ranges

	int gen_pos_size = gen_pos.size();
	if(gen_pos_size > 0){
	//sort positions to do full-length alignment only to distinct positions
	sort(gen_pos.begin(), gen_pos.end(), is_less_gen_pos);//also should help in retreival of ta-cg- lmers

	////////////////////////////////////////Collect DISTINCT only////////////////////////

		vector<AlignInfo> distinct_gen_pos;
		distinct_gen_pos.reserve(gen_pos_size);

		_uli prev = gen_pos[0].pos;
		distinct_gen_pos.push_back(gen_pos[0]);

		int prev_ind = 0;
	    for(_ui i = 1; i < gen_pos_size; i++){
			_uli next = gen_pos[i].pos;
			if(next != prev){
//NOVA				//map previous position including gen_pos[0]
				distinct_gen_pos.push_back(gen_pos[i]);
				prev = next;

			}//if distinct
	
		}//for i

		gen_pos.clear();
		gen_pos = distinct_gen_pos;
		//distinct_gen_pos.clear(); //local variable will disappear anyway
		//distinct_gen_pos.resize(0);
	}//if any pos
}//collect_gen_pos

void Anchor::map_entire_length(vector<AlignInfo> &gen_pos, const vector<_uli> &ta_read1, const _uli &read_len1, const vector<_uli> &c_read1, 
		const bool &forward, Mates &mates_pairs,
		const vector<_uli> &ns_read1, bool mism_map){		

	int gen_pos_size = gen_pos.size();

	int asize = ta_read1.size();

	if(gen_pos_size > 0){
		
		for(_ui i = 0; i < gen_pos_size; i++){
	
			_uli pos1 = gen_pos[i].pos;
			bool cg_mism1 = false;
			vector<_uli> ta_check(asize, 0);
			int mism1 = 0;
				
			vector<_uli> cg_gen_mer;
				
			mism1 = cg_match(read_len1, pos1, ta_read1, c_read1, ta_check, forward, cg_mism1, cg_gen_mer);

			vector<_uli> ta_gen_mer;
			//ta-check: need only for partially mapped reads
			ta_match(read_len1, pos1, ta_read1, ta_check, ns_read1, mism1, ta_gen_mer);
			_ui clip_left = 0, clip_right = 0;
			int cur_local_score = 0;//from original read, not reverse read
			int ascore = read_len1 - mism1;
				
			if(ascore >= entire_score){
					//if only entire-len alignment is allowed, then we don't need to calculate MD for all candidates, 
					//we can calculate MD once if any entire-length valid match is found
					//if not only entire-len, i.e. local is allowed, then we need to calculate MD here to try best local 
					
					update_singles(mates_pairs, pos1, mism1, mism_map, "", 0, 0, ascore, read_len1, "", 0);					
					if(mates_pairs.mark_unique > 1 && mates_pairs.mism1 == 1)//no need for further alignment			
						return;
			}//if ascore
			else if(ascore >= local_score && is_local_alignment){
					calculate_MD(ta_gen_mer, cg_gen_mer, ta_check, read_len1, gen_pos[i].MD, clip_left, clip_right, cur_local_score);
					gen_pos[i].score = cur_local_score;
					gen_pos[i].clip_left = clip_left;
					gen_pos[i].clip_right = clip_right;			
			}//if local is allowed
			else if(is_indel_alignment){
					calculate_MD(ta_gen_mer, cg_gen_mer, ta_check, read_len1, gen_pos[i].MD, clip_left, clip_right, cur_local_score);
					gen_pos[i].score = cur_local_score;
					gen_pos[i].clip_left = clip_left;
					gen_pos[i].clip_right = clip_right;				
			}//if indels are allowed
			
//END_NOVA
	
		}//for i

		// If any match for entire length, then calculate MD for the winner
		//This will save time of calculating MD for all entire-length candidates and only re-calculates mism once
		if(mates_pairs.mark_unique == 1){
			//calculate MD for the winner
			_uli pos1 = mates_pairs.pos1;
			bool cg_mism1 = false;
			vector<_uli> ta_check(asize, 0);
			int mism1 = 0;
			vector<_uli> cg_gen_mer;
			mism1 = cg_match(read_len1, pos1, ta_read1, c_read1, ta_check, forward, cg_mism1, cg_gen_mer);
			vector<_uli> ta_gen_mer;
			//ta-check: need only for partially mapped reads
			ta_match(read_len1, pos1, ta_read1, ta_check, ns_read1, mism1, ta_gen_mer);
			_ui clip_left = 0, clip_right = 0; //from original read, not reverse read
			int cur_local_score = 0;
			calculate_MD(ta_gen_mer, cg_gen_mer, ta_check, read_len1, mates_pairs.MD1 , clip_left, clip_right, cur_local_score);
			//reset score: number of matched bases - (mism * mism_penalty)
			mates_pairs.local_score1 = (read_len1 - mates_pairs.mism1) + (mates_pairs.mism1 * mism_penalty);
			stringstream ss_cigar;
			ss_cigar << read_len1 << 'M';
			mates_pairs.cigar1 = ss_cigar.str();	
		}
	
	}//if any positions

	//if read is aligned entirely with the number of mismatches <= threshold on mismatches, then we are done

}//map_entire_length

void Anchor::map_local(vector<AlignInfo> &gen_pos, const vector<_uli> &ta_read1, const _uli &read_len1, const vector<_uli> &c_read1, 
		const bool &forward, Mates &mates_pairs,
		const vector<_uli> &ns_read1, bool mism_map){

	
	int gen_pos_size = gen_pos.size();	
	if(gen_pos_size == 0)
		return;

	//if read is aligned entirely with the number of mismatches <= threshold on mismatches, then we are done
	if(mates_pairs.mark_unique == 0){
		//scan gen_pos to find MD strings that not empty => we can compare local scores
		_ui max_score = 0;
		int ind_max_score = 0;

		bool amb = false;
		for(int i = 0; i < gen_pos_size; i++){
			//gen_pos.score = read_len - clip_left - clip_right - mism (mappable score)
			if(gen_pos[i].score >= local_score){//need this since if indels, then some scores < local_score 
				//if score is of passing quality
				_ui current_score = alignment_score *(read_len1 - gen_pos[i].clip_left - gen_pos[i].clip_right);
				if(gen_pos[i].score >= current_score){ //choose maximum only among passing scores
					if(max_score < gen_pos[i].score){
						max_score = gen_pos[i].score;
						ind_max_score = i;
						amb = false;
					}
					else if(max_score == gen_pos[i].score)
						amb = true;//only distinct positions have score != 0
				}//if passes quality score
			}//if only for good score 
		}//for i

		//best local alignment is found
		//we want to have read-quality same or better of alignment_score*read_lenght(new read len)
		if(max_score > 0){//found some

				if(!amb){
					//mismatches inside clipped read
					int mappable_bases = read_len1 - gen_pos[ind_max_score].clip_left - gen_pos[ind_max_score].clip_right;
					int mism1 = mappable_bases - gen_pos[ind_max_score].score ;
					clip_md(gen_pos[ind_max_score].MD, gen_pos[ind_max_score].clip_left , gen_pos[ind_max_score].clip_right);
					//md now contains only the mappable bases
					string cigar;
					stringstream ss;
					ss << gen_pos[ind_max_score].clip_left << 'H' << mappable_bases << 'M' << gen_pos[ind_max_score].clip_right << 'H' ;
					cigar = ss.str();
					//update_singles will be called only once for local_alignment, since max-score already being chosen
					//and ambiguous alread being determined over all distinct genomic positions
					//so comparison of mismatches inside update_singles happens only for entire-length alignment, which is ok
					//since in entire-length alignment mappable score directly depends on mism, and read length is same for entire-length alignment 
					int new_local_score = (mappable_bases - mism1) + (mism1 * mism_penalty);

					update_singles(mates_pairs, gen_pos[ind_max_score].pos, mism1, mism_map, gen_pos[ind_max_score].MD, gen_pos[ind_max_score].clip_left , gen_pos[ind_max_score].clip_right, new_local_score, read_len1, cigar, 1); //local alignment	
					//prepare MD string for output: clip bad ends
				
					//update local score:

				}
				else{
					//ambiguous read is found with local alignment
					mates_pairs.mark_unique = 2;
				}

		}//if found

	}//read is not mapped using entire-length-alignment

}//map_local


void Anchor::map_indels(vector<AlignInfo> &gen_pos, const vector<_uli> &ta_read1, const _uli &read_len1, const vector<_uli> &c_read1, 
		const bool &forward, Mates &mates_pairs,
		const vector<_uli> &ns_read1, bool mism_map){

	
	int gen_pos_size = gen_pos.size();	
	if(gen_pos_size == 0)
		return;

	int MIN_MAX_SCORE = 0 - MAX_MISM;

	//if read is aligned entirely with the number of mismatches <= threshold on mismatches, then we are done
//	if(mates_pairs.mark_unique == 0){
		_ui zero = 0;
		//scan gen_pos to find MD strings that not empty => we can compare local scores
		vector< pair< AlignInfo, AlignInfo > > seed_pairs;
		seed_pairs.reserve(gen_pos_size);

		for(int i = 0; i < gen_pos_size - 1; i++){
			//st1 is within the read
			_ui seed_st1 = gen_pos[i].remainder ;
			_ui seed_en1 = read_len1 - 1 - gen_pos[i].shift_bases;//inclusive
			for(int j = i+1; j < gen_pos_size; j++){
				_ui delta_pos = gen_pos[j].pos - gen_pos[i].pos;
				if(delta_pos > gap_threshold)
					break;//too far from the position i
				_ui seed_st2 = gen_pos[j].remainder ;
				_ui seed_en2 = read_len1 - 1 - gen_pos[j].shift_bases;//inclusive	

				if(seed_st1 > seed_en2){//insertion pos1 < pos2, but st1 > en2
					int new_gap_threshold = seed_st1 - seed_en2 - 1;
					if(delta_pos <= new_gap_threshold){
						seed_pairs.push_back( make_pair(gen_pos[i], gen_pos[j]));
						//break;//one pair per gen_pos: break is not good for mapping accuracy
					}
				}//if insertion
				else if(seed_st2 > seed_en1){
					if(delta_pos <= gap_threshold)
						seed_pairs.push_back( make_pair(gen_pos[i], gen_pos[j]));
						//break;
				}//else if deletion

			}//for j
		}//for i
			
		int seed_pairs_size = seed_pairs.size();		
		if(seed_pairs_size > 0 && seed_pairs_size <= indel_threshold){

			int max_score = MIN_MAX_SCORE;
			_uli max_ind = 0;
			_ui max_mism = 0;
			_ui max_clip_left = 0;
			_ui max_clip_right = 0;
			_uli max_read_len = read_len1;
			string max_md("");
			string max_cigar("");
			bool amb = false;
			_uli max_pos = 0;
			
			for(int i = 0; i < seed_pairs_size; i++){
				_ui seed_st1 = seed_pairs[i].first.remainder ;
				_ui seed_st2 = seed_pairs[i].second.remainder ;
				if(seed_st1 < seed_st2){
				//calculate scores array by bases from left to right of the aligned string for leftmost seed from MD string
					string new_md("");
					int new_score = 0;//can be negative
					_ui mappable_score = 0;
					_ui mappable_bases = 0;
					_ui new_clip_left = 0;
					_ui new_clip_right = 0;
					_ui new_read_len = read_len1;
					string cigar("");
					calculate_indel_score(seed_pairs[i].first, seed_pairs[i].second, read_len1, new_md, new_score, mappable_score, mappable_bases, new_clip_left, new_clip_right, new_read_len, cigar);
					//if mappable_score did not pass threshold entire-score or local-score, then it is set to 0
					
					if(mappable_score > zero){
						if(max_score < new_score){
							max_score = new_score;
							max_ind = i;
							max_md = new_md;
							max_mism = mappable_bases - mappable_score;
							max_clip_left = new_clip_left;
							max_clip_right = new_clip_right;
							max_pos = seed_pairs[i].first.pos;//of the leftmost seed
							max_cigar = cigar;
							amb = false;
						}
						else if(max_score > MIN_MAX_SCORE && max_score == new_score)
							amb = true;				
					}//if alignment is of very high qualit outside of indel
				}//if
				else{
					string new_md("");
					int new_score = 0;
					_ui mappable_score = 0;//read length minus mismatches in mappable bases (not including indels)
					_ui mappable_bases = 0;
					_ui new_clip_left = 0, new_clip_right = 0;
					_ui new_read_len = read_len1;
					string cigar("");
					//new_score is alignment score that includes gap-opening and gap-continue penalties
					//current_score is the number of mapped bases in the portion of aligned bases to genome (i.e. inserts are not included here)
					//this allows to check for quality of alignment that includes aligned bases to genome only
					calculate_indel_score(seed_pairs[i].second, seed_pairs[i].first, read_len1, new_md, new_score, mappable_score, mappable_bases, new_clip_left, new_clip_right, new_read_len, cigar);//for insertions read length will be shorter than real length
					//mappable_score already passed all necessary thresholds in calculate_indel_score
					//if mappable_score did not pass thresholds it was set to 0
					if(mappable_score > zero){
						if(max_score < new_score){
							max_score = new_score;
							max_ind = i;
							max_md = new_md;
							max_mism = mappable_bases - mappable_score;
							max_clip_left = new_clip_left;
							max_clip_right = new_clip_right;
							max_read_len = new_read_len;
							max_pos = seed_pairs[i].second.pos;//of the leftmost seed
							max_cigar = cigar;
							amb = false;
						}
						else if(max_score > MIN_MAX_SCORE && max_score == new_score)
							amb = true;				
					}
				}//else
			}//for

			if(max_score > MIN_MAX_SCORE){//found some, but it could be negative if indel is too long

					if(!amb){
						//mismatches inside clipped read
						int mism1 = max_mism;
						//updates_singles is done only once for all indels-alignment candidates, so comparisons
						//inside of update_singles will work
						update_singles(mates_pairs, max_pos, mism1, mism_map, max_md, max_clip_left , max_clip_right, max_score, max_read_len, max_cigar, 2); //indel alignment	
					}
					else{
					    if(mates_pairs.mark_unique == 0)
						mates_pairs.mark_unique = 2;
					    else{
						if(max_score > mates_pairs.local_score1)
							mates_pairs.mark_unique = 2;
					    }//else
					}//else
			}//if found
		}//if not empty
//	}//read is not mapped using entire-length-alignment

}//map_indels

void Anchor::calculate_indel_score(const AlignInfo &leftmost , const AlignInfo &rightmost, 
	const _uli &read_len1,  string &new_md, int &indel_score, _ui &mappable_score, _ui &mappable_bases,
	_ui & new_clip_left, _ui & new_clip_right, _ui & new_read_len, string &cigar)
{
	new_clip_left = 0;
	new_clip_right = 0;
	//leftmost contains the lefmost seed's info, rightmost of the rightmost seed 
	//new_md will be calculated from old MDs of leftmost and rightmost
	//indel_score includes gap_opening and gap_continue penalty and will be in the output
	//mappable_score is the number of exactly matched bases in the alignment to the genome (insertion is not included) (this is also equal to mappable_bases - mismatches inside them, here mappable is all reads' bases except for an insert portion; in case of deletion, all reads' bases)
	//mappable bases = read_len for case of deletion, and = read_len - insertion length
	//for insertion case, indel_threshold must be reset everytime here (current_indel_thresh = min(indel_thresh, st2 - en1 -1);
	bool insertion = false;
	if(leftmost.pos > rightmost.pos)
		insertion = true;

	int MIN_MAX_SCORE = 0 - MAX_MISM;

	//calculate base-by-base score for leftmost from left to right of the string, which is from end of the MD (since MD is for original string, not reversed)
	vector<int> scores_left;
	vector<int> backtrace_left;
	scores_from_md_left(scores_left, backtrace_left, read_len1, leftmost.MD);
	
	//calculate base-by-base score for rightmost from right to left
	vector<int> scores_right;
	vector<int> backtrace_right;
	scores_from_md_right(scores_right, backtrace_right, read_len1, rightmost.MD);

	//calculate the sum and choose index after which there is an insertion/deletion from the max_sum
	//if(leftmost.pos > rightmost.pos, then insertion) and sum is over i and j = i + (leftmost.pos - rightmost.pos) + 1
	//for deletions sum is over i and j = i + 1;
	int max_score = MIN_MAX_SCORE;// - 1000000
	//first possible base after which there is an indel is last base of leftmost seed
	int max_ind = read_len1 - leftmost.shift_bases - 1;
	//last base after which there might be an indel, this is a previous position of the start of the rightmost seed
	int stop = rightmost.remainder;//not inclusive
	_ui delta = 0;
	if(insertion){
		delta = leftmost.pos - rightmost.pos ;//positive, this condition is checked above
		stop = stop - delta ;
		_ui concat_pos = max_ind + delta + 1;//correct
		for(int i = max_ind; i < stop; i++){
			int sum = scores_left[i] + scores_right[concat_pos];//NOVA_NEXT
			concat_pos++;

			if(max_score < sum){
				max_score = sum;
				max_ind = i;
			}//set max_score
		}//for
	
	}else{
		delta = rightmost.pos - leftmost.pos;
		_ui concat_pos = max_ind + 1;
		for(int i = max_ind; i < stop; i++){
			int sum = scores_left[i] + scores_right[concat_pos];//NOVA_NEXT
			concat_pos++;

			if(max_score < sum){
				max_score = sum;
				max_ind = i;
			}//set max_score
		}//for	
	}//deletion


			//calculates max-score subarray, given backtraces of scores 

					
	//now we know that indel occurs right after max_ind and is of length delta (pos-max - pos-min)
	string left_md = leftmost.MD; //MD is for original read: from start to end of the read
	_ui clip_left = read_len1 - (max_ind + 1), clip_right = 0;
	//since max_score was found for these specific substrings, we need to clip right too
	clip_right = backtrace_left[max_ind];//clip right of left_md, i.e. the end of original string
	new_clip_right = clip_right;
	clip_md(left_md, clip_left, clip_right );
	int mappable_bases_left_md = read_len1 - clip_left - clip_right;

	string right_md = rightmost.MD;
	_ui clip_right_delta = delta;
	if(!insertion)
		clip_right_delta = 0;	
	clip_left = 0;
	//according to max_score trim the end of the right_md (not just mergeable part, but also the other end)
	clip_left = read_len1 - backtrace_right[max_ind + clip_right_delta + 1] - 1;//clip left of right_md, i.e. the start of the original read
	new_clip_left = clip_left;
	clip_right = max_ind + clip_right_delta + 1;
	clip_md(right_md, clip_left, clip_right);
	int mappable_bases_right_md = read_len1 - clip_left - clip_right;
		
	_ui mism_inside_mappable = 0;
	stop = left_md.length();
	for(int i = 0; i < stop; i++){
		if(!isdigit(left_md[i]))
			mism_inside_mappable++;
	}
	stop = right_md.length();
	for(int i = 0; i < stop; i++){
		if(!isdigit(right_md[i]))
			mism_inside_mappable++;
	}

	//need to add insert's length, i.e. mappable score is read-len - mism inside mappable bases
	mappable_bases = mappable_bases_left_md + mappable_bases_right_md;
	mappable_score = mappable_bases - mism_inside_mappable;//does not include insert bases, only alignable to genome
	_ui adjusted_entire_score = (read_len1 - clip_right_delta) * alignment_score;
	_ui adjusted_local_score = adjusted_entire_score * alignment_length;//since length of read is adjusted
	//is mappabel_score good enough to make entire-length alignment of mappable reads?
	if(mappable_score >= adjusted_entire_score){
		//do nothing

	}//if mappable bases can align with entire-length alignment quality
	else if(mappable_score >= adjusted_local_score){//local alignment: trim not-mappable-ends of the read
		//calculate max-subarray-score from concatenation of right_md and left_md strings

		_ui current_thresh = mappable_bases * alignment_score;
		if(mappable_score < current_thresh){
			mappable_score = 0;
			indel_score = 0 - MAX_MISM;
			return;
		}


	}//else if mappable bases can align with local alignment quality
	else{
		mappable_score = 0;//should not pick this alignment: too much uncertainty
		indel_score = 0 - MAX_MISM;
		return;
	}//mappable score did not pass

	//calculate new MD: i12i means insertion of 12bp
	//mappable_score = mappable_bases - mismatches inside mappable bases, so this is number of matches in mappable bases
	indel_score = mappable_score + (mism_inside_mappable * mism_penalty) + gap_open_penalty + (gap_cont_penalty * delta);
	if(insertion){
	
		//if right_md ends with an integer and left_md starts with an integer, we'll get longer integer: bug
		concatenate_md(right_md, left_md, new_md);
		stringstream cigar_ss;

		if(new_clip_left > 0)
			cigar_ss << new_clip_left << 'H';
		cigar_ss << mappable_bases_right_md << 'M' << delta << 'I' << mappable_bases_left_md << 'M';
		if(new_clip_right > 0)
			cigar_ss << new_clip_right << 'H';
		cigar = cigar_ss.str();
	}
	else{
		stringstream ss;
		string del("");
		//this needs to be changed: retreived from genome
		for(int i = 0 ; i < delta; i++){
			del = del + 'A';
		}
		if(!isdigit(left_md[0]))
			left_md = "0" + left_md;
		ss << right_md << "^" << del <<  left_md;//left_md must start with digit
		new_md = ss.str();

		stringstream cigar_ss;
		if(new_clip_left > 0)
			cigar_ss << new_clip_left << 'H';
		cigar_ss << mappable_bases_right_md << 'M' << delta << 'D' << mappable_bases_left_md << 'M';
		if(new_clip_right > 0)
			cigar_ss << new_clip_right << 'H';
		cigar = cigar_ss.str();	
	}
	if(insertion)
		new_read_len = read_len1 - delta;
	else
		new_read_len = read_len1;
}//calculate_indel_score

void Anchor::concatenate_md(const string & right_md, const string &left_md, string &new_md){
	//if right_md ends with an integer and left_md starts with an integer, then we can just concatenate the strings
	int lenR = right_md.length();
	int lenL = left_md.length();

	if(!(isdigit(right_md[lenR - 1])))
		new_md = right_md + left_md;
	else if(!(isdigit(left_md[0])))
		new_md = right_md + left_md;
	else{
		//find index of the first mismatch/char from the end of right_md
		int endR = lenR - 1;
		while(endR >= 0 && isdigit(right_md[endR])){
			endR--;
		}
		endR++;
		string last_intR = right_md.substr(endR, lenR - endR);
		istringstream istrR(last_intR);
		int sum;
		istrR >> sum;

		//find index of the firstm mismatch from the start of left_md
		int startL = 0;
		while(startL < lenL && isdigit(left_md[startL])){
			startL++;
		}
		string first_intL = left_md.substr(0, startL);
		istringstream istrL(first_intL);
		int x ;
		istrL >> x;
		sum += x;

		stringstream ss;
		string mdR("");
		if(endR > 0)
			mdR = right_md.substr(0, endR);
		string mdL("");
		if(startL < lenL)
			mdL = left_md.substr(startL, lenL - startL);
		ss << mdR << sum << mdL;
		new_md = ss.str();
	}//both end of right_md and beginning of left_md are digits
}//concatenate_md

/******************************************************************************************
MD is calculated for original string: from start of original read to end
Here we need to calculated score for the right seed in the aligned read (that is a reverse of original, not rev-compl)
Calculations should proceed from the right end of the aligned read toward the leftmost end of alignment
Hence, it should follow MD from start to the end of MD
********************************************************************************************/
void Anchor::scores_from_md_right(vector<int> &scores, vector< int > &backtr, const _uli & read_len1, const string &md)
{

	scores.resize(read_len1, 0);
	backtr.resize(read_len1, 0);
	int size = md.length();
	_ui count = 0;
	//ind is for filling in the vector scores of size read_len1 from end of vector to the beginning
	int ind = read_len1 - 1, i = 0;
	int prev_score = 0;
	backtr[read_len1 - 1] = read_len1 - 1;

	for(i = 0; i < size; i++){
		if(ind < 0){
			cout << "ERROR: scores_from_md_right: ind < 0, but MD is not processed." << endl;
			space_release();
			exit(0);
		}
		if(isdigit(md[i])){
			int j = i;
			stringstream ss;
			while(j < size && isdigit(md[j]) ){
				ss << md[j];
				j++;
			}
			istringstream istr(ss.str());
			int x;
			istr >> x;
			int stop = ind - x;
			for(int e = ind; e > stop; e--){
				//update for these bases score
				int from_left = prev_score + match_penalty;
				if(from_left >= match_penalty){
					scores[e] = from_left;
					if(e != read_len1 - 1)
						backtr[e] = backtr[e + 1];
				}//
				else{
					scores[e] = match_penalty;
					backtr[e] = e;
				}
				prev_score = scores[e];
			}//for
			ind = stop;
			i = j - 1;
		}//if
		else{
			int from_left = prev_score + mism_penalty;
			if(from_left >= mism_penalty){
				scores[ind] = from_left;
				if(ind != read_len1 - 1)
					backtr[ind] = backtr[ind + 1];
			}//
			else{
				scores[ind] = mism_penalty;
				backtr[ind] = ind;
			}
			prev_score = scores[ind];
			ind--;
		}//char
		
	}//for i


}//scores_from_md_right to left

void Anchor::scores_from_md_left(vector<int> &scores, vector< int > &backtr, const _uli& read_len1, const string &md)
{

	scores.resize(read_len1, 0);
	backtr.resize(read_len1, 0);

	int size = md.length();
	_ui count = 0;
	int ind = 0, i = size - 1;
	int prev_score = 0;
	
	for(i = size - 1; i >= 0; i--){
		if(ind < 0){
			cout << "ERROR: scores_from_md_right: ind = read_len1, but MD is not processed." << endl;
			space_release();
			exit(0);
		}
		if(isdigit(md[i])){
			int j = i;
			stringstream ss;
			while(j >= 0 && isdigit(md[j])){
				j--;
			}
			ss << md.substr(j+1, i-j);
			istringstream istr(ss.str());
			int x;
			istr >> x;
			int stop = ind + x ;
			for(int e = ind; e < stop; e++){
				//update for these bases score
				int from_left = prev_score + match_penalty;
				if(from_left >= match_penalty){
					scores[e] = from_left;
					if(e != 0)
						backtr[e] = backtr[e - 1];
				}//
				else{
					scores[e] = match_penalty;
					backtr[e] = e;
				}
				prev_score = scores[e];
			}//for
			ind = stop;
			i = j + 1;
		}//if
		else{
			int from_left = prev_score + mism_penalty;
			if(from_left >= mism_penalty){
				scores[ind] = from_left;
				if(ind != 0)
					backtr[ind] = backtr[ind - 1];
			}//
			else{
				scores[ind] = mism_penalty;
				backtr[ind] = ind;
			}
			prev_score = scores[ind];
			ind++;
		}//char
		
	}//for i


}//scores_from_md_left to right of the aligned reads (i.e. reverse of the read)

int main(int argc, char* argv[])
{
	/**************************************
	This program calculates distribution of different fingerprints in the genome given
	*******************************************/
		/**************************************
		bitset [11000] is translated into a string as 00011 and into long int as 3
		so, least significant bits in the bitset are the leftmost
  *******************************************/
	
	if(argc < 6){//prefix, reads, output
		usage();
		return 0;
	}//if
	Anchor anc_ct(argc, argv);//reads genome builds ta-, cg- binary representations once	
	map_reads(anc_ct);
	cout << "map_reads is done" << endl;
	return 0;
}//main

