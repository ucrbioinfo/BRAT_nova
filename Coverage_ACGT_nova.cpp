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
long int rand_res(long int total_len, double i)
{
	srand( time(0)*i );
	long int res = 0 + rand() % (total_len - 1);
	return res;
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

double choose(double top, double bottom){
	//given top and bottom it returns Log_10(choose(top, bottom)), where choose is binomial coef
	if(bottom == 0 || bottom == top)
		return 0; //choose = 1, log(choose)=0
	else if(bottom == 1)
		return log10(top);
	double sum = 0.0, current;
	long int i, j;

	long int stop = top - bottom;
	for(i = top; i > stop; i--){
		current = log10(top) - log10(bottom);
		sum += current;
		top -= 1.0;
		bottom -= 1.0;
		
	}//for
	
	 return sum; 
	
}//choose
double p_value(double total_genes, double cluster, double colors, double total_all){
	/*
	Given a set of total_genes' sequences, 
	given a set of chosen genes (cluster) sequences,
	given # occurrences of motif in cluster,
	given # occurrences of motif in all_genes
	
	Find hypergeometric log_10(p_value)
	
	*/
	double stop = total_all;
	if(cluster < total_all)
		stop = cluster;
	double N = total_genes;
	double N_cluster = N - cluster;
	
	long int i;
	double j;//total_all - i
	double sum = 0.0, cluster_choose_i, N_cluster_choose_j , N_choose_total_all, current;
	for(i = colors; i <= stop; i++){
		j = total_all - i;
		cluster_choose_i = choose(cluster, i);
		N_cluster_choose_j = choose(N_cluster, j);
		N_choose_total_all = choose(N, total_all);
		current = (cluster_choose_i + N_cluster_choose_j - N_choose_total_all);// 
		sum += pow(10.0, current);
	}
	
	
	return (-log10(sum));
//return sum;	
}

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

void fill_line(vector<unsigned int> &mark, string &prev, string &cur, string &next, long int alen_prev, 
	long int alen, long int alen_next, long int &cur_base)
{
	long int j = 0;
			for(j = 0; j < alen; j++){					
				if(cur[j] == 'C' || cur[j] == 'c'){
					//determine CHH, CHG or CG
					if(j < alen - 2){//two following bases exist
						if(cur[j+1] == 'G' || cur[j+1] == 'g')
							mark[cur_base]=1;
						else{
							if(cur[j+2] == 'G' || cur[j + 2] == 'g')
								mark[cur_base] = 2;
							else
								mark[cur_base] = 3;
						}//else first consecutive is not G
					}//if fits into a line
					else if(j == alen - 2){
						if(cur[j+1] == 'G' || cur[j+1] == 'g')
							mark[cur_base]=1;
						else{
							if(!next.empty()){
								if(next[0] == 'G' || next[0] == 'g')
									mark[cur_base] = 2;
								else
									mark[cur_base] = 3;
							}
							else
								mark[cur_base] = 2;//CH? preference is to CHG
						}//else
					}//else if 
					else if(j == alen - 1){
						if(!next.empty()){
							if(next[0] == 'G' || next[0] == 'g')
								mark[cur_base] = 1;
							else{
								if(alen_next > 1){
									if(next[1] == 'G' || next[1] == 'g')
										mark[cur_base] = 2;
									else
										mark[cur_base] = 3;
								}
								else
									mark[cur_base] = 2;//CH? preference is to CHG
							}//else
						}//next
						else
							mark[cur_base] = 1;//we don't know, but preference is to CG
					}//alen - 1
				}//if C							
				else if(cur[j] == 'G' || cur[j] == 'g'){
					if(j > 1){
						if(cur[j-1] == 'C' || cur[j-1] == 'c')
							mark[cur_base] = 4;
						else{
							if(cur[j-2] == 'C' || cur[j-2] == 'c')
								mark[cur_base] = 5;
							else
								mark[cur_base] = 6;
						}//not C

					}//if fits into line
					else if(j == 1){
						if(cur[0] == 'C' || cur[0] == 'c')
							mark[cur_base] = 4;
						else{
							if(alen_prev > 0){
								if(prev[alen_prev - 1] == 'C' || prev[alen_prev - 1] == 'c')
									mark[cur_base] = 5;
								else
									mark[cur_base] = 6;
							}//prev not empty
							else
								mark[cur_base] = 5;
						}//else CH
					}//j = 1
					else if(j == 0){
						if(alen_prev > 0){
							if(prev[alen_prev - 1] == 'C' || prev[alen_prev - 1] == 'c')
								mark[cur_base] = 4;
							else{//CH
								if(alen_prev > 1){
									if(prev[alen_prev - 2] == 'C' || prev[alen_prev - 2] == 'c')
										mark[cur_base] = 5;
									else
										mark[cur_base] = 6;
								}
								else
									mark[cur_base] = 5;
							}//else CH
						}//prev
						else
							mark[cur_base] = 4;
					}//if j = 0
				}//if base is G						
												
				cur_base++;						
			}//for j				
}
void fill_content(vector<string> &gen_names, vector<unsigned int> &mark, long int gen_ind)
{
	ifstream in;
	string line, prev, next, cur;
	long int j;
	in.open(gen_names[gen_ind].c_str(), ios::in);					
	getline(in, line);//first FASTA line					
	getline(in, prev);//at least one line in the reference will exist
	getline(in, cur);//current line
	//NEED to process previous	
	long int alen_prev = 0, alen = 0;
	alen_prev = prev.length();
	alen = cur.length();
	long int cur_base = 0;		
	line="";
	if(!prev.empty()){
		//mark previous line's characters
		if(in.eof()){
			cur ="";
			fill_line(mark, line, prev, cur, 0, alen_prev, 0, cur_base);
			return;
		}//if second line is eof
		else
			fill_line(mark, line, prev, cur, 0, alen_prev, alen, cur_base);//current exists
	}//process prev
	else
		return;//no lines in reference

	getline(in, line);//next
	if(in.eof()){//next is eof, but current is not
		next="";
		fill_line(mark, prev, cur, next, alen_prev, alen, 0, cur_base);//process current
		return;
	}//if

		long int alen_next;
		while(!in.eof()){//if reference is more than three lines		
			next = line;
			alen_next = next.length();		
			fill_line(mark, prev, cur, next, alen_prev, alen, alen_next, cur_base);

			prev = cur;
			cur = next;
			alen_prev = alen;
			alen = alen_next;

			getline(in, line);	
		}//while

	
	in.close(); in.clear();	
	next = "";
	//process the last line
	fill_line(mark, prev, cur, next, alen_prev, alen, 0, cur_base);	
}//fill_content()

// prints in BED format
void print(vector< map<string, Point> > &ref_seq, vector< vector<Cov> > &genome,
	vector< vector<Cov> > &rev_genome, vector<string> &gen_names, string pref, 
	ofstream &outf,  long int t){

		string line;
		long int i, j;

		ofstream out;
		ifstream in;

		map<string, Point>::iterator start = ref_seq[t].begin();
		map<string, Point>::iterator end = ref_seq[t].end();


//NEW format: chr, start, end, coverage=C+T (forw) or G+A(rev), meth level C/(C+T) or G/(G+A)

				start = ref_seq[t].begin();

				for(; start != end; start++){
					string achrom =  start->first;

					long int ind = (start->second).x;
					long int asize_of_chrom = (start->second).y;
					long int gen_ind = (start->second).z;

					vector<unsigned int> mark(asize_of_chrom);

					for(j = 0; j < asize_of_chrom; j++){
						mark[j] = 0;
					}
					fill_content(gen_names, mark, gen_ind);

					for(j = 0; j < asize_of_chrom; j++){
						
						if(mark[j] >= 1 && mark[j] < 4){
							double total = 0.0 + genome[ind][j].c + genome[ind][j].t;
							double meth_level = 0;
							if(total > 0)
								meth_level = (genome[ind][j].c + 0.0)/total;
							string acontent("CpG:");
							if(mark[j] == 2)
								acontent = "CHG:";
							else if(mark[j] == 3)
								acontent = "CHH:";

							outf << achrom << "\t" << j << "\t" << j << "\t" << acontent << total << "\t" << meth_level << "\t" 
								<< "+" <<  endl;

						}//if C
						else if(mark[j] >= 4){
							double total = 0.0 + rev_genome[ind][j].g + rev_genome[ind][j].a;
							double meth_level = 0;
							string acontent("CpG:");
							if(mark[j] == 5)
								acontent = "CHG:";
							else if(mark[j] == 6)
								acontent = "CHH:";

							if(total > 0)
								meth_level = (rev_genome[ind][j].g + 0.0)/total;
							outf << achrom << "\t" << j << "\t" << j << "\t" << acontent << total << "\t" << meth_level << "\t" 
								<< "-" <<  endl;
						}//else if G
					}//for j

				}//for i

}//print()

bool process_cigar(string cigar, vector<int> &parts, char &type){
	//if there is no indel, then parts has only one integer, else three integers
//int count_mism(const string & cigar, bool &is_indel, const string &read){
	bool is_indel = false;
	//break cigar into chunks ending with character
	vector<int> indices;
	int asize = cigar.length();
	for(int i = 0; i < asize; i++){
		if(!isdigit(cigar[i])){
			indices.push_back(i);
			if(cigar[i] == 'D' || cigar[i] == 'I'){
				is_indel = true;
				type = cigar[i];
			}
		}//if not digit
	}//for	

	if(!is_indel)
		return is_indel;

	//else find match, indel, match
	int ind_size = indices.size();
	//check if there is D, deletion 

	for(int i = 0; i < ind_size; i++){
		int ind = indices[i];
		if(cigar[ind] == 'M' || cigar[ind] == 'D' || cigar[ind] == 'I'){
			int ind2 = 0;//start of number
			if(i > 0)
				ind2 = indices[i-1] + 1;
			int len2 = ind - ind2;
			string sub2 = cigar.substr(ind2, len2);
			istringstream istr1(sub2);
			int match ;
			istr1 >> match;
			parts.push_back(match);
		}

	}//for

	return is_indel;

}//process_cigar

//increment_forw(genome, chrom_id, asize, pos2, stop2, j, read2);
void increment_forw(vector< vector<Cov> > &genome, long int chrom_id, long int asize, 
	long int forw_pos,  string forw_read, string cigar)
{
	//long int stop, long int read_pos,
	vector<int> parts;
	char indel_type = 'M';
	bool is_indel = process_cigar(cigar, parts, indel_type);

	long int i, j = 0;
	if(!is_indel){
		long int stop = forw_pos + forw_read.length();
		for(i = forw_pos; i < stop; i++){

			if(i >= 0 && i < asize){
					if(forw_read[j] == 'A' || forw_read[j] == 'a')
						genome[chrom_id][i].a++;
					else if(forw_read[j] == 'C' || forw_read[j] == 'c')
						genome[chrom_id][i].c++;
					else if(forw_read[j] == 'G' || forw_read[j] == 'g')
						genome[chrom_id][i].g++;
					else if(forw_read[j] == 'T' || forw_read[j] == 't')
						genome[chrom_id][i].t++;
			}//if valid position
			j++;
		}//for i
	}//no indels
	else if(indel_type == 'D'){//deletion
		//align the first match
		i = forw_pos;
		j = 0;
		long int stop = forw_pos + parts[0];//match1
		for(i = forw_pos; i < stop; i++){

			if(i >= 0 && i < asize){
					if(forw_read[j] == 'A' || forw_read[j] == 'a')
						genome[chrom_id][i].a++;
					else if(forw_read[j] == 'C' || forw_read[j] == 'c')
						genome[chrom_id][i].c++;
					else if(forw_read[j] == 'G' || forw_read[j] == 'g')
						genome[chrom_id][i].g++;
					else if(forw_read[j] == 'T' || forw_read[j] == 't')
						genome[chrom_id][i].t++;
			}//if valid position
			j++;
		}//for i
		//align the second match: j has the correct position within the read now
		i = stop + parts[1];//indel
		stop = i + parts[2];//match2
		for(; i < stop; i++){

			if(i >= 0 && i < asize){
					if(forw_read[j] == 'A' || forw_read[j] == 'a')
						genome[chrom_id][i].a++;
					else if(forw_read[j] == 'C' || forw_read[j] == 'c')
						genome[chrom_id][i].c++;
					else if(forw_read[j] == 'G' || forw_read[j] == 'g')
						genome[chrom_id][i].g++;
					else if(forw_read[j] == 'T' || forw_read[j] == 't')
						genome[chrom_id][i].t++;
			}//if valid position
			j++;
		}//for i

	}//deletion
	else if(indel_type == 'I'){
		//align the first match
		i = forw_pos;
		j = 0;
		long int stop = forw_pos + parts[0];//match1
		for(i = forw_pos; i < stop; i++){

			if(i >= 0 && i < asize){
					if(forw_read[j] == 'A' || forw_read[j] == 'a')
						genome[chrom_id][i].a++;
					else if(forw_read[j] == 'C' || forw_read[j] == 'c')
						genome[chrom_id][i].c++;
					else if(forw_read[j] == 'G' || forw_read[j] == 'g')
						genome[chrom_id][i].g++;
					else if(forw_read[j] == 'T' || forw_read[j] == 't')
						genome[chrom_id][i].t++;
			}//if valid position
			j++;
		}//for i
		//align the second match: i has the correct position within the read now	
		j += parts[1];//indel
		stop += parts[2];
		for(; i < stop; i++){

			if(i >= 0 && i < asize){
					if(forw_read[j] == 'A' || forw_read[j] == 'a')
						genome[chrom_id][i].a++;
					else if(forw_read[j] == 'C' || forw_read[j] == 'c')
						genome[chrom_id][i].c++;
					else if(forw_read[j] == 'G' || forw_read[j] == 'g')
						genome[chrom_id][i].g++;
					else if(forw_read[j] == 'T' || forw_read[j] == 't')
						genome[chrom_id][i].t++;
			}//if valid position
			j++;
		}//for i
	}//insertion
	
}//increment_forw()


void parse_options(int argc, char* argv[], string &singles, string &pref, string &genome)
{
	/*
	pairs or singles are in SAM format
	*/
	long int i, j;
	for(i = 1; i < argc; i++){
		if(argv[i][0] != '-'){
			cerr << "\nERROR: Invalid char=" << argv[i][0] << "=this" << endl ;		
			return;
		}
		switch(argv[i][1]){
		case 's': singles = argv[++i]; break;//SAM format any reads: singles or pairs
		case 'P': pref = argv[++i]; break;
		case 'r': genome = argv[++i]; break;
		default: cerr << "ERROR: Invalid char=" << argv[i][1] << "=this" << endl ;  exit(0); break;
		}
	}//for i

}//parse_opt()



void find_indel_type(string cigar, char &indel_type, bool & is_indel, int &indel){
	indel = 0;
	int asize = cigar.length();
	
	int indel_ind = 0;
	for(int i = 0; i < asize; i++){
		if(!isdigit(cigar[i])){
			if(cigar[i] == 'D' || cigar[i] == 'I'){
				is_indel = true;
				if(cigar[i] == 'D')
					indel_type = 'D';
				else 
					indel_type = 'I';
				indel_ind = i;
				break;
			}//if indel
		}//if
	}//for

	if(is_indel){
		int base_ind = indel_ind - 1;
		while(base_ind >= 0){
			if(isdigit(cigar[base_ind]))
				base_ind--;
			else
				break;
		}//while
		if(base_ind < 0)
			base_ind = 0;
		else
			base_ind++;
		string sub = cigar.substr(base_ind, indel_ind - base_ind);
		istringstream istr(sub);

		istr >> indel;

	}//if indel
}//find_indel_type

long int find_mate1_stop(long int pos1, int len1, string cigar){
	//finds the first position in the genome that is on the right of mate1 and not covered by mate1
	//BRAT-nova supports only single insertion or deletion
	int asize = cigar.length();
	char indel_type = 'M';
	bool is_indel = false;
	int indel_ind = 0;
	int indel = 0;
	find_indel_type(cigar, indel_type,is_indel, indel);

	if(!is_indel)
		return (pos1 + len1);
	if(is_indel){
		if(indel_type == 'D')
			return (pos1 + len1 + indel);
		else
			return (pos1 + len1 - indel);
	}//if indel


}//find_mate1_stop

void reduce_match1(string &cigar, int delta){

	//find index of the first occurrence of M
	int asize = cigar.length();
	int i = 0;
	while(i < asize){
		if(cigar[i] == 'M')
			break;
		i++;
	}//while
	if(i == asize)
	{
		cout << "ERROR: no character M is found in CIGAR string: " << cigar << endl;
		exit(0);
	}
	//find index of the number of bases referring to the first M
	int ind_M = i;
	i--;
	while(i >= 0){
		if(isdigit(cigar[i]))
			i--;
		else{
			i++;
			break;
		}
	}//while
	if(i < 0)
		i = 0;
	string sub = cigar.substr(i, ind_M - i);
	istringstream istr(sub);
	int match ;
	istr >> match;//length of match
	match = match - delta;
	stringstream ss;
	ss << match << (cigar.substr(ind_M, asize - ind_M));
	cigar = ss.str();

}//reduce_match1

void clip_reduce_match2(string &cigar, int delta){
	//clips everything before the second match, for example 10H15M2D20M12H with delta 2, will be
	// 18M12H (10H15M have been clipped and 20M reduced by 2 resulting in 18M)

	//find index of the first occurrence of M
	int asize = cigar.length();
	int i = 0;
	while(i < asize){
		if(cigar[i] == 'M')
			break;
		i++;
	}//while
	if(i == asize)
	{
		cout << "ERROR: no character M is found in CIGAR string: " << cigar << endl;
		exit(0);
	}
	i++;
	//find second M
	while(i < asize){
		if(cigar[i] == 'M')
			break;
		i++;
	}//while
	if(i == asize)
	{
		cout << "ERROR: no second character M is found in CIGAR string with indels: " << cigar << endl;
		exit(0);
	}
	//find index of the number of bases referring to the first M
	int ind_M = i;
	i--;
	while(i >= 0){
		if(isdigit(cigar[i]))
			i--;
		else{
			i++;
			break;
		}
	}//while
	if(i < 0)
		i = 0;
	string sub = cigar.substr(i, ind_M - i);
	istringstream istr(sub);
	int match ;
	istr >> match;//length of match
	match = match - delta;
	stringstream ss;
	ss << match << 'M';
	cigar = ss.str();
}//clip_reduce_match2

void reduce_match2(string &cigar, int delta){
	// keep everything before match2, and reduce match2
	//for example 10H15M2D20M12H with delta 2, will be 10H15M2D18M 

	//find index of the first occurrence of M
	int asize = cigar.length();
	int i = 0;
	while(i < asize){
		if(cigar[i] == 'M')
			break;
		i++;
	}//while
	if(i == asize)
	{
		cout << "ERROR: no character M is found in CIGAR string: " << cigar << endl;
		exit(0);
	}
	i++;
	//find second M
	while(i < asize){
		if(cigar[i] == 'M')
			break;
		i++;
	}//while
	if(i == asize)
	{
		cout << "ERROR: no second character M is found in CIGAR string with indels: " << cigar << endl;
		exit(0);
	}
	//find index of the number of bases referring to the first M
	int ind_M = i;
	i--;
	while(i >= 0){
		if(isdigit(cigar[i]))
			i--;
		else{
			i++;
			break;
		}
	}//while
	if(i < 0)
		i = 0;
	string sub = cigar.substr(i, ind_M - i);
	string leftmost = cigar.substr(0, i);
	istringstream istr(sub);
	int match ;
	istr >> match;//length of match
	match = match - delta;
	stringstream ss;
	ss << leftmost << match << 'M';
	cigar = ss.str();
}//reduce_match2

void reduce_match1_clip_right(string &cigar, int delta){

	//clips everything after the first match, for example 10H15M2D20M12H with delta 2, will be
	// 10H13M (2D20M12H have been clipped and 15M reduced by 2 resulting in 13M)

	//find index of the first occurrence of M
	int asize = cigar.length();
	int i = 0;
	while(i < asize){
		if(cigar[i] == 'M')
			break;
		i++;
	}//while
	if(i == asize)
	{
		cout << "ERROR: no character M is found in CIGAR string: " << cigar << endl;
		exit(0);
	}

	//find index of the number of bases referring to the first M
	int ind_M = i;
	i--;
	while(i >= 0){
		if(isdigit(cigar[i]))
			i--;
		else{
			i++;
			break;
		}
	}//while
	if(i < 0)
		i = 0;
	string sub = cigar.substr(i, ind_M - i);
	
	istringstream istr(sub);
	int match ;
	istr >> match;//length of match
	match = match - delta;
	stringstream ss;
	ss << match << 'M';
	cigar = ss.str();

}//reduce_match1_clip_right

void find_match1_match2(const string & cigar, int & match1, int & match2){
	//finds the length of match1 and match2, for example 10H15M2D20M12H 
	// match1 is 15, match2 is 20 (integers preceding M)

	//find index of the first occurrence of M
	int asize = cigar.length();
	int i = 0;
	while(i < asize){
		if(cigar[i] == 'M')
			break;
		i++;
	}//while
	if(i == asize)
	{
		cout << "ERROR: no character M is found in CIGAR string: " << cigar << endl;
		exit(0);
	}
	int ind_M = i;
	i--;
	while(i >= 0){
		if(isdigit(cigar[i]))
			i--;
		else{
			i++;
			break;
		}
	}//while
	if(i < 0)
		i = 0;
	string sub = cigar.substr(i, ind_M - i);
	istringstream istr(sub);
	
	istr >> match1;//length of match


	i = ind_M + 1;
	//find second M
	while(i < asize){
		if(cigar[i] == 'M')
			break;
		i++;
	}//while
	if(i == asize)
	{
		cout << "ERROR: no second character M is found in CIGAR string with indels: " << cigar << endl;
		exit(0);
	}
	//find index of the number of bases referring to the first M
	ind_M = i;
	i--;
	while(i >= 0){
		if(isdigit(cigar[i]))
			i--;
		else{
			i++;
			break;
		}
	}//while
	if(i < 0)
		i = 0;
	sub = cigar.substr(i, ind_M - i);
	istringstream istr2(sub);
	istr2 >> match2;//length of match

}//find_match1_match2

void trim_first_mappable(long int pos, long int stop, long int & pos2, long int stop2, 
	string & read2, string & cigar2){
	//finds mappable bases of mate2 and trims overlapping bases of mate2

	//when mate2 is entirely inside mate1
	if((pos2 >= pos) && (stop2 <= stop)){
		cigar2 = "0M";
		read2 = "";
		return;
	}
	//find out if there is an indel, if so then type of indel
	char indel_type = 'M';//no indels
	int indel = 0; //bases of indel
	bool is_indel = false;
	find_indel_type(cigar2, indel_type, is_indel, indel);

	int fm = 0;//first mappable overlapping with the end of mate1
	int len2 = read2.length();
	//Case 1: mate2 overlaps with stop, then take the rightmost portion of mate2 starting with 
	//the base overlapping with stop
	if(stop2 > stop){
		//Subcase1: no indels
		if(!is_indel){
			fm = stop - pos2;
			read2 = read2.substr(fm, len2 - fm);
			reduce_match1(cigar2, fm);
			return;
		}//no indels
		else{
			if(indel_type == 'D'){
				int match1 = 0, match2 = 0;
				find_match1_match2(cigar2, match1, match2);
				if((pos2 + match1 - 1) >= stop){
					fm = stop - pos2;
					pos2 = stop;
					read2 = read2.substr(fm, len2 - fm);
					reduce_match1(cigar2, fm);
				}//if stop overlaps with mate2 within match1 
				else if((pos2 + match1 + indel - 1) >= stop){
					fm = match1;
					pos2 = pos2 + match1 + indel;
					read2 = read2.substr(match1, len2 - match1);
					clip_reduce_match2(cigar2, 0);

				}//else if stop overlaps with mate2 within deletion
				else{
					fm = stop - pos2 - indel;
					pos2 = stop;
					read2 = read2.substr(fm, len2 - fm);
					clip_reduce_match2(cigar2, (fm - match1));
				}//else stop overlaps with mate2 within match2
			}//if deletion
			else{
				int match1 = 0, match2 = 0;
				find_match1_match2(cigar2, match1, match2);
				if((pos2 + match1 - 1) >= stop){
					fm = stop - pos2;
					pos2 = stop;
					read2 = read2.substr(fm, len2 - fm);
					reduce_match1(cigar2, fm);
				}//stop overlaps with mate2 within match1 with an insertion
				else{
					fm = stop - pos2 + indel;
					pos2 = stop;
					read2 = read2.substr(fm, len2 - fm);
					clip_reduce_match2(cigar2, (fm - match1 - indel));
				}//else case insertion and stop overlaps mate2 within match2
			}//insertion
		}//indels
	}//if
	else{
		//Case2: mate2 does not overlap with stop, but overlaps with pos, then
		//take the leftmost portion of mate2 just before pos1	
		if(!is_indel){
				long int delta = len2 - (pos - pos2);
				read2 = read2.substr(0, pos - pos2);
				reduce_match1(cigar2, delta);		
		}//no indels
		else{
			int match1 = 0, match2 = 0;
			find_match1_match2(cigar2, match1, match2);
			if(indel_type == 'D'){
				if((pos2 + match1 - 1) >= (pos - 1)){
					read2 = read2.substr(0, pos - pos2);
					int delta = match1 - (pos - pos2);
					reduce_match1_clip_right(cigar2, delta);
				}//if mate2 overlaps with pos1 within match1
				else if((pos2 + match1 + indel - 1) >= (pos - 1)){
					read2 = read2.substr(0, match1);
					reduce_match1_clip_right(cigar2, 0);
				}//if mate2 overlaps with pos1 within deletion
				else{
					read2 = read2.substr(0, pos - pos2 - indel);
					int delta = match2 - ( pos - pos2 - indel - match1);
					reduce_match2(cigar2, delta);
				}//if mate2 overlaps with pos1 within the second match

			}//deletion
			else{
				if((pos2 + match1 - 1) >= (pos - 1)){
					read2 = read2.substr(0, pos - pos2);
					int delta = match1 - (pos - pos2);
					reduce_match1_clip_right(cigar2, delta);
				}//overlaps within first match
				else{
					read2 = read2.substr(0, pos - pos2 + indel);
					int delta = match2 - (pos - pos2 - match1);
					reduce_match2(cigar2, delta);
				}//overlaps within second match
			}//insertion
		}//indels
	}//Case2


}//trim_first_mappable

int main(int argc, char* argv[]){

	srand( time(0) );

	string read;
	char strand;
	long int id, i, j, asize, u, y, chrom_id, pos, mism, mism1, mism2, st1, st2, en1, en2;
	long int orig_start1, orig_start2;
	long int count = 0;
	bool no_dupl = false;
	bool is_pe = false;
	string singles(""), pref, references;//out_type default is meth level only for Cs
	parse_options(argc, argv, singles, pref, references);

	ofstream outf;
	outf.open(pref.c_str(), ios::out);

	ifstream in;
	in.open(references.c_str(), ios::in);
	if(!in){
		cerr << "\nERROR: cannot open " << references << endl;
		exit(0);
	}
	vector<string> gen_names;
	in >> read;
	while(!in.eof()){
		gen_names.push_back(read);
		in >> read;
	}//while
	in.close(); in.clear();

	vector<long int> size_chrom;
	
	asize = gen_names.size();
	string line;
	map<string, long int> chrom_names;
	vector<string> chroms_names;
		string chrom_name ;

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

	long int num_chroms = size_chrom.size();

	//to save space, we process references sequentially
	long int thresh_size = 150000000;
	vector< map<string, Point> > ref_seq;
	map<string, Point> first_map;
	ref_seq.push_back(first_map);
	long int last_map_ind = 0;
	long int cur_index = 0;
	long int cur_size = 0;

cout << "current map: " ;
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


	vector<string> names;

		string aname;
		in.open(singles.c_str(), ios::in);
		if(!in){
			cerr << "ERROR: can't open " << singles << endl;
			exit(0);
		}
		getline(in, aname);
		while(!in.eof()){
			names.push_back(aname);
			getline(in, aname);
		}//while
		in.close(); in.clear();


	long int chrom_num = size_chrom.size();
	
 	long int t;
	cur_size = ref_seq.size();
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

		vector< vector<Cov> > genome;
		vector< vector<Cov> > rev_genome;
		for(i = 0; i < amap_size; i++){
			vector<Cov> dum(sizes_chroms_cur_map[i]);
			genome.push_back(dum);
			rev_genome.push_back(dum);
		}//for i


			for( y = 0; y < names.size(); y++){


				in.open(names[y].c_str(), ios::in);
				if(!in){
					cerr << "ERROR: can't open " << names[y] << endl;
					exit(0);
				}
				//else
				//	cout << "opening " << names[y] << endl;

				long int count_lines = 0, flag;
				string read_name, qual, dum, rest, cigar, md, xb, xc, xc2, md2, cigar2, read_name2;
				string qual2, read2;
				//xb if + then align to positive, if - then align to negative
				//xc shows original starting point

				in >> read_name;
				string prev_name("");

				while(!in.eof()){
					if(read_name[0] != '@'){

					count_lines++;

					in >> flag >> chrom_name >> pos >> dum >> cigar;
					in >> dum >> dum >> dum;
					in >> read >> qual >> dum >> md >> xb >> xc;
					pos--;		

					xc = xc.substr(5, xc.length() - 5);//XC:i:NUMBER, NUMBER starts at 5-th pos
					istringstream istr_xc(xc);
					istr_xc >> orig_start1;
					
					strand = '+';
					if(xb[5] == '-')
						strand = '-';


					map<string, Point>::iterator afind = ref_seq[t].find(chrom_name);

					if(afind != end){

						chrom_id = (afind->second).x ;
						asize = (afind->second).y;
						if(chrom_id < 0 || chrom_id >= chrom_num)
						{
							cout << "chrom_id is out of boundary at line " << count_lines << endl;
							exit(0);
						}



						long int len_mer1 = read.length();
						
						j = 0;
						long int stop = pos + len_mer1;
						if(strand == '+'){//ADD md and cigar here to support indels
							increment_forw(genome, chrom_id, asize, pos, read, cigar);
						}//if first mate is mapped to positive strand 
						else{
							increment_forw(rev_genome, chrom_id, asize, pos, read, cigar);
						}//else rev-complement of first map is mapped to positive strand
						//process reverse read
						//check if there is an overlap between two mates

						read2 = "";
						//check if there is a mate
						long int flag2, pos2;
						string chr2;						
						if((flag == 99) || (flag == 83)){
							in >> read_name2 >> flag2 >> chr2 >> pos2 >> dum >> cigar2;
							in >> dum >> dum >> dum ;
							in >> read2 >> qual2 >> dum >> md2 >> dum >> xc2;
							pos2--;
							//identify is there overlap
							long int len_mer2 = read2.length();
							long int stop2 = pos2 + len_mer2;
							j = 0;
							//find stop for mate1
							stop = find_mate1_stop(pos, len_mer1, cigar);
							//find stop2 for mate2 (the first rightmost not covered genome base)
							stop2 = find_mate1_stop(pos2, len_mer2, cigar2);
							bool overlap = false;
							long int start2 = pos2;
							if((pos <= pos2) && (pos2 < stop)){
								
									overlap = true;

							}//if overlapping with pisitive mate1
							else if((pos2 <= pos) && (pos < stop2) ){
								
									overlap = true;

							}//if overlapping with negative mate1
							if(overlap)
								trim_first_mappable(pos, stop, pos2, stop2, read2, cigar2);

							if(strand == '+'){//ADD md and cigar here to support indels
								increment_forw(genome, chrom_id, asize, pos2,  read2, cigar2);
							}//if first mate is mapped to positive strand 
							else{
								increment_forw(rev_genome, chrom_id, asize, pos2, read2, cigar2);
							}//else rev-complement of first map is mapped to positive strand
						}//mate exists
						
					}//if there is chromosome in a current chunk
					}//if does not start with @
					else
						getline(in, rest);//read the rest of line

					in >> read_name;
				}//while
				in.close(); in.clear();

			}//for y all files with results
			cout << "Finished calculating ACGT-count from paired-end reads." << endl;

			print(ref_seq, genome, rev_genome, gen_names, pref, outf, t);				

	}//for ref_seq size

	outf.close(); 


	return 0;

}
