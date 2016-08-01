//  BRAT is a whole genome bisulfite sequence mapping program
//  Copyright (C) 2009  Elena Yavorska Harris
//
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

#include<iomanip>
#include<map>
#include<bitset>
#include<ctime>
#include<list>
#include<cmath>
#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
using namespace std;

class Point
{

public:
	//constructors of this class
	Point(long int xx, long int yy, long int zz) : x(xx), y(yy), z(zz) {};
	Point() : x(0), y(0), z(0) {};
	
	
	void operator =(Point a) { x = a.x; y = a.y; z = a.z; };
	//x is index, y is tails
	bool operator <(const Point& a) const
	{ 
		if(y < a.y) return true; 
		else return false;
	}


	long int x;
	long int y;
	long int z;

};
class Cov{
public:
	Cov(unsigned short int x, unsigned short int y, unsigned short int z, unsigned short int w) : a(x), c(y), g(z), t(w) {};
	Cov() : a(0), c(0), g(0), t(0) {};

	void operator =(Cov s){ a = s.a; c = s.c; g = s.g; t = s.t; };
	unsigned short int a;
	unsigned short int c;
	unsigned short int g;
	unsigned short int t;
};
long int str_to_int(string s)
{
	istringstream is(s);
	long int result = -1;
	is >> result;
	return result;
}
double str_to_double(string s)
{
	istringstream is(s);
	double result = -1;
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
long int min(long int x, long int y){ 
	if(x < y)
		return x;
	else
		return y;
}

char rev(char x)
{
	if(x == 'a')
		return 'c';
	else if(x == 'c')
		return 'g';
	else
		return 'a';

}//rev



string rev_complement(string s)
{
	string res("");
	int alen = s.length() - 1; 
	for(int j = alen; j >= 0; j--)
{
		if(s[j] == 'A')
			res = res + 'T';
		else if(s[j] == 'T')
			res = res + 'A';
		else if(s[j] == 't')
			res = res + 'a';
		else if(s[j] == 'a')
			res = res + 't';
		else if(s[j] == 'C')
			res = res + 'G';
		else if(s[j] == 'G')
			res = res + 'C';
		else if(s[j] == 'c')
			res = res + 'g';
		else if(s[j] == 'g')
			res = res + 'c';
		else
			res = res + 'N';
		
	}
	return res;

}//rev_comple

void fill_names_sizes(string references, map<string, long int> &chrom_names, vector<string> &chroms_names, vector<long int> &size_chrom)
{	

	ifstream in;
	in.open(references.c_str(), ios::in);
	if(!in){
		cerr << "\nERROR: cannot open " << references << endl;
		exit(0);
	}
	vector<string> gen_names;
	string afile;
	in >> afile;
	while(!in.eof()){
		gen_names.push_back(afile);
		in >> afile;
	}//while
	in.close(); in.clear();

	
	
	long int asize = gen_names.size();
	string line;

		string chrom_name ;
	long int i;
	for(i = 0; i < asize; i++){
	
		in.open(gen_names[i].c_str(), ios::in);

		if(!in){
			cerr << "\nERROR: cannot open " << gen_names[i] << endl;
			exit(0);
		}//if
cout << "opening " << gen_names[i] << endl;
		getline(in, line);//first line of fasta file
		istringstream istr(line);
		char ch;
		istr >> ch;
		istr >> chrom_name;
		chrom_names.insert(make_pair(chrom_name, i));
		chroms_names.push_back(chrom_name);

		getline(in, line);
		long int cur_size = 0;
		while(!in.eof()){
			cur_size += line.length();
			getline(in, line);
		}//while
		in.close(); in.clear();

		size_chrom.push_back(cur_size);
cout << "size of " << chrom_name << " is " << cur_size << endl;
	}//for

}//fill_names_sizes()

void fill_ref_seq(vector< map<string, Point> > &ref_seq, vector<string> &chroms_names, vector<long int> &size_chrom,
	const long int thresh_size, long int num_chroms )
{
	map<string, Point> first_map;
	ref_seq.push_back(first_map);
	long int last_map_ind = 0;
	long int cur_index = 0;
	long int cur_size = 0;

cout << "current map: " ;
	long int i = 0;
	for(i = 0; i < num_chroms; i++){
		cur_size += size_chrom[i];
		ref_seq[last_map_ind].insert(make_pair(chroms_names[i], Point(cur_index, size_chrom[i], i)));
		cur_index++;
		cout << chroms_names[i] << ", " ;
		if(cur_size >= thresh_size && i < num_chroms-1){
			map<string, Point> amap;
			ref_seq.push_back(amap);
			last_map_ind++;
			cur_size = 0;
			cur_index = 0;
cout << endl;
cout << "current map: " ;
		}//if
	
	}//for i
cout << endl;
}//fill_ref_seq()

void read_file_names(string from_file, vector<string> &file_names)
{		ifstream in;
		string aname;
		in.open(from_file.c_str(), ios::in);
		if(!in){
			cerr << "ERROR: can't open " << from_file << endl;
			exit(0);
		}
		in >> aname;
		while(!in.eof()){
			file_names.push_back(aname);
			in >> aname;
		}//while
		in.close(); in.clear();
}//read_file_names

void count_copies(vector< vector< int> > &genome, vector< vector<int> > &rev_genome, 
	vector< vector<int> > &pairs_forw, vector< vector<int> > &pairs_rev,
	vector<string> names, 
	vector< map<string, Point> > &ref_seq, long int t, long int chrom_num)
{
	string read1, read2;
	char strand;
	string id, chrom_name ;
	long int   asize,  y, chrom_id,  pos2 ;

	ifstream in;
		map<string, Point>::iterator end = ref_seq[t].end();
			long int asize_names = names.size();
			for( y = 0; y < asize_names; y++){
				in.open(names[y].c_str(), ios::in);
				if(!in){
					cerr << "ERROR: can't open " << names[y] << endl;
					exit(0);
				}
				else
					cout << "opening " << names[y] << endl;

				long int count_lines = 0;
				in >> id;
				while(!in.eof()){
					count_lines++;
					long int flag, pos;
					string dum, cigar1, qual1, md1, xb1, xc1, chr;

					in >> flag >> chr >> pos >> dum >> cigar1 ;
					in >> dum >> dum >> dum ;
					in >> read1 >> qual1 >> dum >> md1 >> xb1 >> xc1;

					istringstream istr(xc1.substr(5, xc1.length() - 5));
					long int orig_start1;
					istr >> orig_start1;
					orig_start1--;
					strand = xb1[5];				

					map<string, Point>::iterator afind = ref_seq[t].find(chrom_name);
					if(afind != end){
						chrom_id = (afind->second).x ;
						asize = (afind->second).y;
						if(chrom_id < 0 || chrom_id >= chrom_num)
						{
							cout << "chrom_id is out of boundary at line " << count_lines << endl;
							exit(0);
						}
						if((flag == 99) || (flag == 83)){
							string read_name2, chr2, cigar2, qual2, md2, xb2, xc2;
							long int flag2;
							in >> read_name2 >> flag2 >> chr2 >> pos2 >> dum >> cigar2;
							in >> dum >> dum >> dum ;
							in >> read2 >> qual2 >> dum >> md2 >> xb2 >> xc2;
							istringstream istr2(xc2.substr(5, xc2.length() - 5));//XC:i:30300
							long int orig_start2 ;
							istr2 >> orig_start2;
							orig_start2--;//since original start is 1-base count in SAM output
							long int forw_pos = orig_start1;
							if(orig_start2 < orig_start1)//by leftmost read of a pair
								forw_pos = orig_start2;
							
							if(strand == '+'){
								pairs_forw[chrom_id][forw_pos]++;
							}
							else
								pairs_rev[chrom_id][forw_pos]++;

						}else{//single or an unpaired mate
							if(strand == '+')
								genome[chrom_id][orig_start1]++;
							else
								rev_genome[chrom_id][orig_start1]++;
						
						}//unpaired mate or single read

					}//if there is chromosome in a current chunk
					in >> id;
				}//while
				in.close(); in.clear();
			}//for y all files with results
  
}//count_copies()

void choose_random(vector< vector< int> > &genome, vector< vector< int> > &rev_genome,
	vector< vector<int> > &pairs_forw, vector< vector<int> > &pairs_rev,
	vector<long int> &sizes_chroms_cur_map, long int amap_size)
{
	long int i, j;
		for(i = 0; i < amap_size; i++){
			long int cur_size = sizes_chroms_cur_map[i];
			for(j = 0; j < cur_size; j++){
				if(pairs_forw[i][j] > 1){
						long int res = rand() % pairs_forw[i][j] + 1;//random pair
						pairs_forw[i][j] = res;
				}//if pairs
				if(pairs_rev[i][j] > 1){
						long int res = rand() % pairs_rev[i][j] + 1;//random pair
						pairs_rev[i][j] = res;
				}//if pairs
				if(genome[i][j] > 1){
						long int res = rand() % genome[i][j] + 1;//random single
						genome[i][j] = res;
				}//only singles
				if(rev_genome[i][j] > 1){
						long int res = rand() % rev_genome[i][j] + 1;//random single
						rev_genome[i][j] = res;
				}//only singles
			}//for j
		}//for i
}//choose_random

void print(vector< vector< int> > &genome, vector< vector< int> > &rev_genome, 
	vector< vector<int> > &pairs_forw, vector< vector<int> > &pairs_rev,
	vector<string> names, vector< map<string, Point> > &ref_seq, long int t, long int chrom_num)
{
	string id, read1, read2;
	char strand;
	string chrom_name ;
	long int   asize,  y, chrom_id;

	ifstream in;
	ofstream out;
		map<string, Point>::iterator end = ref_seq[t].end();
			long int asize_names = names.size();
			for( y = 0; y < asize_names; y++){
				in.open(names[y].c_str(), ios::in);
				if(!in){
					cerr << "ERROR: can't open " << names[y] << endl;
					exit(0);
				}
				else
					cout << "opening " << names[y] << endl;

				string output = names[y] + ".nodupl";
				out.open(output.c_str(), ios::app);

				long int count_lines = 0;
				in >> id;


				while(!in.eof()){
					count_lines++;
					long int flag, pos;
					string dum1, dum2, dum3, dum4, dum5, qual1, md1, xb1, xc1, cigar1, chr;

					in >> flag >> chr >> pos >> dum1 >> cigar1 ;
					in >> dum2 >> dum3 >> dum4 ;
					in >> read1 >> qual1 >> dum5 >> md1 >> xb1 >> xc1;

					istringstream istr(xc1.substr(5, xc1.length() - 5));
					long int orig_start1;
					istr >> orig_start1;
					orig_start1--;
					strand = xb1[5];			

					map<string, Point>::iterator afind = ref_seq[t].find(chrom_name);

					if(afind != end){

						chrom_id = (afind->second).x ;
						asize = (afind->second).y;
						if(chrom_id < 0 || chrom_id >= chrom_num)
						{
							cout << "chrom_id is out of boundary at line " << count_lines << endl;
							exit(0);
						}
						if((flag == 99) || (flag == 83)){
							string dum6, dum7, dum8, dum9, dum10, dum11;
							string read_name2, chr2, cigar2, qual2, md2, xb2, xc2;
							long int flag2, pos2;
							in >> read_name2 >> flag2 >> chr2 >> pos2 >> dum6 >> cigar2;
							in >> dum7 >> dum8 >> dum9 ;
							in >> read2 >> qual2 >> dum10 >> md2 >> xb2 >> xc2;
							istringstream istr2(xc2.substr(5, xc2.length() - 5));//XC:i:30300
							long int orig_start2 ;
							istr2 >> orig_start2;
							orig_start2--;//since original start is 1-base count in SAM output
							long int forw_pos = orig_start1;
							if(orig_start2 < orig_start1)//by leftmost read of a pair
								forw_pos = orig_start2;
							
							if(strand == '+'){
								
								pairs_forw[chrom_id][forw_pos]--;
								if(pairs_forw[chrom_id][forw_pos] == 0){
									out << id << '\t' <<  flag << '\t' << chr << '\t' 
									<< pos << '\t' << dum1 << '\t' << cigar1 << '\t'  ;
									out << dum2 << '\t' << dum3 << '\t' << dum4 << '\t'  ;
									out << read1 << '\t' << qual1 << '\t'  
									<< dum5 << '\t' << md1 << '\t' << xb1 << '\t' << xc1 << '\t' << endl;

									out << read_name2 << '\t' <<  flag2 << '\t' << chr2 << '\t' 
									<< pos2 << '\t' << dum6 << '\t' << cigar2 << '\t'  ;
									out << dum7 << '\t' << dum8 << '\t' << dum9 << '\t'  ;
									out << read2 << '\t' << qual2 << '\t'  
									<< dum10 << '\t' << md2 << '\t' << xb2 << '\t' << xc2 << '\t' << endl;
								}//if random choice
							}
							else{
								pairs_rev[chrom_id][forw_pos]--;
								if(pairs_rev[chrom_id][forw_pos] == 0){
									out << id << '\t' <<  flag << '\t' << chr << '\t' 
									<< pos << '\t' << dum1 << '\t' << cigar1 << '\t'  ;
									out << dum2 << '\t' << dum3 << '\t' << dum4 << '\t'  ;
									out << read1 << '\t' << qual1 << '\t'  
									<< dum5 << '\t' << md1 << '\t' << xb1 << '\t' << xc1 << '\t' << endl;

									out << read_name2 << '\t' <<  flag2 << '\t' << chr2 << '\t' 
									<< pos2 << '\t' << dum6 << '\t' << cigar2 << '\t'  ;
									out << dum7 << '\t' << dum8 << '\t' << dum9 << '\t'  ;
									out << read2 << '\t' << qual2 << '\t'  
									<< dum10 << '\t' << md2 << '\t' << xb2 << '\t' << xc2 << '\t' << endl;
								}//if random choice
							}
						}else{//single or an unpaired mate
							if(strand == '+'){
								genome[chrom_id][orig_start1]--;
								if(genome[chrom_id][orig_start1] == 0){
									out << id << '\t' <<  flag << '\t' << chr << '\t' 
									<< pos << '\t' << dum1 << '\t' << cigar1 << '\t'  ;
									out << dum2 << '\t' << dum3 << '\t' << dum4 << '\t'  ;
									out << read1 << '\t' << qual1 << '\t'  
									<< dum5 << '\t' << md1 << '\t' << xb1 << '\t' << xc1 << '\t' << endl;
								}
							}
							else{
								rev_genome[chrom_id][orig_start1]--;
								if(rev_genome[chrom_id][orig_start1] == 0){
									out << id << '\t' <<  flag << '\t' << chr << '\t' 
									<< pos << '\t' << dum1 << '\t' << cigar1 << '\t'  ;
									out << dum2 << '\t' << dum3 << '\t' << dum4 << '\t'  ;
									out << read1 << '\t' << qual1 << '\t'  
									<< dum5 << '\t' << md1 << '\t' << xb1 << '\t' << xc1 << '\t' << endl;
								}//if
							}//else
						
						}//unpaired mate or single read

					}//if there is chromosome in a current chunk
					in >> id;
				}//while
				in.close(); in.clear();
				out.close(); out.clear();
			}//for y all files with results
		
}//print()

void parse_options(int argc, char* argv[], string &input_files, string &genome)
{
	long int i;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-'){
			cerr << "\nERROR: Invalid char=" << argv[i][0] << "=this" << endl ;		
			return;
		}
		switch(argv[i][1]){
		case 's': input_files = argv[++i]; break;
		case 'r': genome = argv[++i]; break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl ;  exit(0); break;
		}
	}//for i

}//parse_opt()


int main(int argc, char* argv[]){

	srand( time(0) );

	string read1, read2, strand;
	long int  i, j;

	long int count = 0;
	bool no_dupl = false;
	bool is_pe = false;
	string input_files(""), references("");//out_type default is matrix of ACGT-counts
	parse_options(argc, argv, input_files, references);
	if(input_files.empty()){
		cout << "ERROR: please provide file with the names of files with the mapped reads" << endl;
		exit(0);
	}
	if(references.empty())
	{
		cout << "ERROR: please provide file with the names of fasta files with references" << endl;
		exit(0);
	}
	ifstream in;
	vector<long int> size_chrom;
	map<string, long int> chrom_names;
	vector<string> chroms_names;	
	fill_names_sizes(references, chrom_names, chroms_names, size_chrom);
	
	long int num_chroms = size_chrom.size();

	//to save space, we process references sequentially
	long int thresh_size = 15000000;
	vector< map<string, Point> > ref_seq;
	fill_ref_seq(ref_seq, chroms_names, size_chrom,	thresh_size, num_chroms);

	vector<string> names;

		read_file_names(input_files, names);
		cout << "the number of files with results is " << names.size() << endl;


	long int chrom_num = size_chrom.size();
	string aname;

	long int t;
	long int cur_size = ref_seq.size();
	long int amap_size ;
	for(t = 0; t < cur_size; t++){

		amap_size = ref_seq[t].size();

		vector<long int> sizes_chroms_cur_map(amap_size);
		map<string, Point>::iterator start = ref_seq[t].begin();
		map<string, Point>::iterator end = ref_seq[t].end();
		for(; start != end; start++){
			long int ind = (start->second).x;
			long int asize_cur = (start->second).y;
			sizes_chroms_cur_map[ind] = asize_cur;
		}//for

		vector< vector<int> > pairs_forw;
		vector< vector<int> > pairs_rev;
		vector< vector<int> > genome;
		vector< vector<int> > rev_genome;//for single-end reads to count the number of reads mapped to each genome pos
		for(i = 0; i < amap_size; i++){
			vector<int> dum(sizes_chroms_cur_map[i]);
			for(j = 0; j < sizes_chroms_cur_map[i]; j++)
				dum[j] = 0;
			genome.push_back(dum);
			rev_genome.push_back(dum);
			pairs_forw.push_back(dum);
			pairs_rev.push_back(dum);
		}//for i

		//count copy-duplicates for each location
		count_copies(genome, rev_genome, pairs_forw, pairs_rev, names, ref_seq, t, chrom_num);

		//if count > 1 at some location, choose random representative
		choose_random(genome, rev_genome, pairs_forw, pairs_rev, sizes_chroms_cur_map, amap_size);

		//print out the mapped results without copy-duplicates
		print(genome, rev_genome, pairs_forw, pairs_rev, names, ref_seq, t, chrom_num);

	}//for ref_seq size

	return 0;

}

