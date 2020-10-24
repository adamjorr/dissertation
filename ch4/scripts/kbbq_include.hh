
// FILE:/home/adam/code/kbbq/include/kbbq/kseq.hh 


#ifndef __KBBQ_KSEQ_H
#define __KBBQ_KSEQ_H

#include <htslib/bgzf.h>
#include <htslib/kseq.h>

namespace kseq{
	KSEQ_INIT(BGZF*, bgzf_read)
}

#endif
// FILE:/home/adam/code/kbbq/include/kbbq/readutils.hh 


#ifndef READUTILS_H
#define READUTILS_H

#include <vector>
#include <unordered_map>
#include <string>
#include <errno.h>
//
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <htslib/bgzf.h>
#include "bloom.hh"
#include "covariateutils.hh"
#include "kseq.hh"

#define INFER_ERROR_BAD_QUAL 2

//fwd declare
namespace covariateutils{
	struct dq_t;
}

namespace readutils{

	static int8_t complement[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

	// int correction_len(const bloom::bloomary_t& t, int k);

	//get read sequence as a string in the forward orientation
	inline std::string bam_seq_str(bam1_t* bamrecord){
		std::string seq;
		char* s = (char*)bam_get_seq(bamrecord);
		for(size_t i = 0; i < bamrecord->core.l_qseq; ++i){
			seq.push_back(bam_is_rev(bamrecord) ?
				seq_nt16_str[seq_nt16_table['0' + 3-seq_nt16_int[bam_seqi(s, i)]]] :
				seq_nt16_str[bam_seqi(s, i)]);
		}
		if(bam_is_rev(bamrecord)){
			std::reverse(seq.begin(), seq.end());
		}
		return seq;
	}

	class CReadData{
		public:
			static std::unordered_map<std::string, std::string> rg_to_pu;
			static std::unordered_map<std::string, int> rg_to_int;
			CReadData(){}
			CReadData(bam1_t* bamrecord, bool use_oq = false);
			// hello?
			CReadData(kseq::kseq_t* fastqrecord, std::string rg = "", int second = 2, std::string namedelimiter = "_");
			std::string seq;
			std::vector<uint8_t> qual;
			std::vector<bool> skips;
			std::string name;
			std::string rg;
			bool second;
			std::vector<bool> errors;
			std::string str_qual();
			std::string canonical_name();
			inline int get_rg_int() const{return this->rg_to_int[this->rg];}
			inline std::string get_pu() const{return this->rg_to_pu[this->rg];}
			std::vector<bool> not_skipped_errors() const;
			//fill errors attribute given sampled kmers and thresholds.
			void infer_read_errors(const bloom::Bloom& b, const std::vector<int>& thresholds, int k);
			//fix one error and return the index of the fixed base; std::string::npos if no fixes are found
			size_t correct_one(const bloom::Bloom& t, int k);
			static void load_rgs_from_bamfile(bam_hdr_t* header);
			//fill errors attribute given trusted kmers
			std::vector<bool> get_errors(const bloom::Bloom& trusted, int k, int minqual = 6, bool first_call = true);
			std::vector<uint8_t> recalibrate(const covariateutils::dq_t& dqs, int minqual = 6) const;
			CReadData substr(size_t pos = 0, size_t count = std::string::npos) const;

	};
}

#endif

// FILE:/home/adam/code/kbbq/include/kbbq/gatkreport.hh 


#ifndef KBBQ_GATKREPORT_HH
#define KBBQ_GATKREPORT_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>

namespace gatkreport{

enum column_type{STRING = 0, FLOAT, INT};

class TableValue{
private:
	column_type type;
	union{
		long double f;
		unsigned long long i;
	}
	std::string s;
public:
	TableValue(): type(), f(), s() {}
	TableValue(unsigned long long i): TableValue() {this->set(i);}
	TableValue(long double f): TableValue() {this->set(f);}
	TableValue(std::string s): TableValue() {this->set(s);}

	inline column_type current_type(){return type;}
	template<typename T> T get();
	template<>
	unsigned long long get<unsigned long long>(){
		return type == INT ? i : 0;
	}
	template<>
	long double get<long double>(){
		return type == FLOAT ? f : 0;
	}
	template<>
	std::string get<std::string>(){
		return type == STRING ? s : "";
	}
	inline TableValue& set(unsigned long long i){this->i = i; this->type = INT; return this;}
	inline TableValue& set(long double f){this->f = f; this->type = FLOAT; return this;}
	inline TableValue& set(std::string s){this->s = s; this->type = STRING; return this;}
};

class TableRow{
private:
	std::vector<TableValue> columns;
	std::vector<column_type> types;
	std::vector<std::string> headers;
public:
	TableRow(): columns(), types() {}
	TableRow(const std::string& line, const std::vector<column_type>& t, const std::vector<std::string>& h):
		types(t), columns(), headers(h)
	{
		std::istringstream is(line);
		std::string token;
		for(column_type t : types){
			is >> token;
			switch(t){
				case STRING: columns.emplace_back(token);
					break;
				case INT: columns.emplace_back(std::stoull(token));
					break;
				case FLOAT: columns.emplace_back(std::stold(token));
					break;
			}
		}
	}
	//space-delimited values in the row.
	inline explicit operator std::string() const {
		std::string ret;
		for(size_t i = 0; i < columns.size(); ++i){
			std::string column_val = "";
			switch(types[i]){
				case FLOAT: long double val = columns[i].get<long double>();
					std::string format = headers[i] == "Errors" ? "%.2f" : "%.4f"
					int outsize = std::snprintf(nullptr, 0, format.c_str(), val);
					column_val.resize(outsize + 1);
					std::snprintf(&column_val[0], column_val.size(), format.c_str(), val);
					break;
				case INT: column_val = std::to_string(columns[i].get<unsigned long long>());
					break
				case STRING:
				default:
					column_val = columns[i].get<std::string>;
			}
			ret += column_val + " ";
		}
		ret.pop_back(); //remove trailing space
		return ret;
	}
	template <typename T>
	T get(size_t i){return columns[i].get<T>();}
}

class GATKTable{
private:
	std::vector<TableRow> rows;
	std::string title;
	std::string description;
	std::vector<std::string> headers;
	std::vector<column_type> types;
	std::vector<size_t> widths;
	/*
	const std::map<std::string, std::string> precision = {
		{"EmpiricalQuality", ".4"},
		{"EstimatedQReported",".4"},
		{"Errors",".2"}
	};
	*/
public:
	GATKTable(): rows(), title(), description(), headers(), types() {}
	GATKTable(std::string tablestr):
		rows(), title(), description(), headers(), types()
	{
		std::istringstream in(tablestr);
		std::string headerline;
		std::string headerprefix = "#:GATKTable:"
		std::getline(in, headerline); //process first line
		if(headerline.substr(0,headerprefix.length()) != headerprefix){
			throw std::invalid_argument("Error: Unable to parse first line of table. " +
				"Ensure input is a valid GATKTable.");
		}
		headerline.erase(0,headerprefix.length());
		std::istringstream hdrin(headerline);
		std::string token;
		std::getline(hdrin, token, ':');
		size_t ncols = std::stoull(token);
		std::getline(hdrin, token, ':');
		size_t nrows = std::stoull(token);

		while(std::getline(hdrin, token, ':')){
			switch(token.pop_back()){
				case 'd': types.push_back(INT);
					break;
				case 'f': types.push_back(FLOAT);
					break;
				case ';': //do nothing
					break;
				case 's':
				default: types.push_back(STRING);
			}
		}
		if(types.size() != ncols){
			throw std::invalid_argument("Error: Number of types doesn't match header!" +
				"Ensure the type string includes all types.\n");
		}

		std::getline(in, headerline); //process second line; it's the title and description
		if(headerline.substr(0,headerprefix.length()) != headerprefix){
			throw std::invalid_argument("Error: Unable to parse second line of table. " +
				"Ensure input is a valid GATKTable.");
		}
		headerline.erase(0,headerprefix.length());
		hdrin.str(headerline);
		std::getline(hdrin, title, ':');
		std::getline(hdrin, description, ':');
		std::getline(in, headerline); //process the 3rd line; it's the column headers
		hdrin.str(headerline);
		std::copy(std::istream_iterator<std::string>(hdrin),
				std::istream_iterator<std::string>(),
				std::back_inserter(headers));
		if(types.size() != headers.size()){
			throw std::invalid_argument(
				"Error: Number of headers doesn't match the number of columns.\n");
		}
		while(std::getline(in, token)){
			rows.emplace_back(token, types, headers);
		}
		if(rows.size() != nrows){
			throw std::invalid_argument(
				"Error: Number of stated rows doesn't match number of actual rows.\n");
		}
	}
	//TODO: ctor, default and std::string
	//return the header string without a trailing newline
	inline std::string headerstring() const{
		std::string headerstr = "#:GATKTable:" + std::to_string(types.size()) +
			":" + std::to_string(rows.size()) + ":";
		for(size_t i = 0; i < types.size(); ++i){
			column_type t = types[i];
			std::string h = headers[i];
			std::string colcode = "";
			switch(t){
				case STRING: colcode = "%s";
					break;
				case INT: colcode = "%d";
					break;
				case FLOAT: colcode = h == "Errors" ? "%.2f" : "%.4f";
					break;
			}
			headerstr += colcode + ":";
		}
		return headerstr + ";";
	}

	//the title string without a trailing newline
	inline std::string titlestring() const{
		return "#:GATKTable:" + title + ":" + description;
	}

	inline std::vector<size_t> get_col_widths() const{
		std::vector<size_t> widths{headers.size(),0};
		std::transform(headers.begin(), headers.end(), widths.begin(),
			[](std::string str) -> size_t {return str.length();});
		for(size_t i = 0; i < rows.size(); ++i){
			std::istringstream rowin = rows[i];
			for(size_t j = 0; j < widths.size(); ++j){
				std::string colstr;
				rowin >> colstr;
				if(colstr.length() > widths[j]){
					widths[j] = colstr.length();
				}
			}
		}
		return widths;
	}

	inline explicit operator std::string() const {
		std::string ret = this->headerstring() + "\n";
		ret += this->titlestring() + "\n";
		std::ostringstream os();
		std::vector<size_t> widths = this->get_col_widths();
		for(size_t i = 0; i < widths.size(); ++i){
			os << std::setw(widths[i]) << headers[i] << "  "; 
		}
		ret += os.str();
		ret.erase(ret.end()-2);
		os.str(""); //reset output stream
		std::istringstream is("");
		std::string token;
		for(const std::string& s : rows){
			is.str(s);
			is >> token;
			os << std::setw(widths[0]) << token;
			for(size_t i = 1; i < widths.size(); ++i){
				is >> token;
				os << "  " << std::setw(widths[i]) << token;
			}
			os << "\n";
		}
		ret += os.str();
		return ret;
	}
	TableRow& get(size_t i){return &rows[i];}
};

class GATKReport{
private:
	std::vector<GATKTable> tables;
	std::string version;
public:
	GATKReport(){}
	GATKReport(const std::vector<GATKTable>& tables): tables(tables) {}
	GATKReport(const std::string& filename): tables(), version(){
		std::ifstream inf(filename);
		std::string tablestr;
		std::string headerline;
		std::string versionprefix = "#:GATKReport.v";
		std::getline(inf, headerline);
		if(headerline.substr(0,versionprefix.length()) != versionprefix){
			throw std::invalid_argument("Error: Unable to parse first line of input file " +
				filename + " ; Ensure it is a valid GATKReport.");
		}
		headerline.erase(0,versionprefix.length());
		size_t colon_pos = headerline.find(':');
		version = headerline.substr(0, colon_pos);
		size_t ntables = std::stoull(headerline.substr(colon_pos+1, std::string::npos));
		for(std::string line; std::getline(inf, line);){
			if(line.empty()){ //blank line delimits tables
				tablestr.pop_back(); //get rid of ending newline
				tables.emplace_back(tablestr);
				tablestr.clear();
			} else {
				tablestr += line + '\n';
			}
		}
		if(ntables != tables.size()){
			throw std::invalid_argument("Found " + std::to_string(ntables) + " tables in " +
				filename + ", but " + std::to_string(tables.size()) " were declared." +
				"Ensure the file is not truncated and adjust the declared number of tables.");
		}
	}
	//return the header string without a trailing newline
	inline std::string headerstring() const{
		return "#:GATKReport.v" + version + ":" + std::to_string(tables.size());
	}
	inline explicit operator std::string() const {
		std::string ret = this->headerstring() + "\n";
		for(const std::string& s : tables){ret += s + "\n";}
		return ret;
	}
};

}

#endif
// FILE:/home/adam/code/kbbq/include/kbbq/bloom.hh 


#ifndef KBBQ_BLOOM_HH
#define KBBQ_BLOOM_HH
#include <cstdint>
#include <utility>
#include <htslib/hts.h>
#include <minion.hpp>
#include <memory>
#include <functional>
#include <iostream>
#include "bloom_filter.hpp"
#include <stdexcept>
#include <immintrin.h>

#define PREFIXBITS 10
#define KBBQ_MAX_KMER 32

//the number of prefix hashes is 1<<PREFIXBITS - 1

namespace bloom{

class blocked_bloom_filter: public bloom_filter
{
protected:
	typedef unsigned char v32uqi __attribute__ ((__vector_size__ (32))); //activate gcc vectorization
	typedef long long base_type;
	typedef base_type v4di __attribute__ ((__vector_size__ (32)));
	typedef v4di cell_type;
	static constexpr cell_type cell_type_zero = {0ULL,0ULL,0ULL,0ULL}; //used to initialize
	typedef std::unique_ptr<cell_type, std::function<void(cell_type*)>> table_type;
	// typedef std::unique_ptr<unsigned char, std::function<void(unsigned char*)>> table_type;
public:
	static const size_t block_size = 512; //512 bits = 64 bytes
	table_type bit_table_;
	//TODO: ensure table size is a multiple of block_size
	blocked_bloom_filter(): bloom_filter(){}
	blocked_bloom_filter(const bloom_parameters& p){
		projected_element_count_ = p.projected_element_count;
     	inserted_element_count_ = 0;
     	random_seed_ = (p.random_seed * 0xA5A5A5A5) + 1 ;
     	desired_false_positive_probability_ = p.false_positive_probability;
    	salt_count_ = std::max(p.optimal_parameters.number_of_hashes, 2u);
		table_size_ = p.optimal_parameters.table_size;
		//ensure table fits a full block
		table_size_ += (table_size_ % block_size) != 0 ? block_size - (table_size_ % block_size) : 0;
		generate_unique_salt();
		void* ptr = 0;
		int ret = posix_memalign(&ptr, block_size / bits_per_char,
			table_size_ / bits_per_char);
		if(ret != 0){
			throw std::bad_alloc();
		}
		bit_table_ = table_type(static_cast<cell_type*>(ptr),
			[](cell_type* x){free(x);});
		std::uninitialized_fill_n(bit_table_.get(), table_size_ / bits_per_char / sizeof(cell_type),
			cell_type_zero);
	}
	//delete copy ctor
	blocked_bloom_filter(const blocked_bloom_filter&) = delete;
	//delete copy assignment
	blocked_bloom_filter& operator=(const blocked_bloom_filter&) = delete;
	//move ctor
	blocked_bloom_filter(blocked_bloom_filter&& o){
		salt_count_ = std::move(o.salt_count_);
		table_size_ = std::move(o.table_size_);
		bit_table_ = std::move(o.bit_table_);
		salt_ = std::move(o.salt_);
		projected_element_count_ = std::move(o.projected_element_count_);
		inserted_element_count_ = std::move(o.inserted_element_count_);
		random_seed_ = std::move(o.random_seed_);
		desired_false_positive_probability_ = std::move(o.desired_false_positive_probability_);
	}

	//move function
	inline blocked_bloom_filter& operator=(blocked_bloom_filter&& o){
		if(this != &o){
			salt_count_ = std::move(o.salt_count_);
			table_size_ = std::move(o.table_size_);
			bit_table_ = std::move(o.bit_table_);
			salt_ = std::move(o.salt_);
			projected_element_count_ = std::move(o.projected_element_count_);
			inserted_element_count_ = std::move(o.inserted_element_count_);
			random_seed_ = std::move(o.random_seed_);
			desired_false_positive_probability_ = std::move(o.desired_false_positive_probability_);
		}
		return *this;
	}


	inline virtual size_t num_blocks() const{
		assert((table_size_ % block_size) == 0);
		return table_size_ / block_size;
	}

	// inline virtual void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const {
	// 	bit_index = hash % block_size; //which bit in the block?
 //      	bit       = bit_index % bits_per_char; // 
	// }

	inline virtual size_t get_block(const bloom_type& hash) const{
		// how lighter does it
		// return (hash % (table_size_ - block_size + 1)) / bits_per_char; 
		// we intead pick an aligned block.
		// hash index * sizeof(block) / sizeof(cell)
		return (hash % num_blocks()) * (block_size / bits_per_char) / (sizeof(cell_type));
		
	}

	//pick the correct vector inside the block and then the correct
	//base type inside the vector.
	inline virtual std::pair<size_t,size_t> get_vector_unit(const size_t& bit_index) const{
		return std::make_pair((bit_index / bits_per_char) / sizeof(cell_type),
			(bit_index / bits_per_char) % (sizeof(cell_type)/sizeof(base_type)));
	}

	inline virtual void insert(const unsigned char* key_begin, const size_t& length){
		size_t bit_index = 0;
		size_t bit = 0;
		//index in table with first byte of block
		size_t block_idx = get_block(hash_ap(key_begin, length, salt_[0]));
		cell_type* block = bit_table_.get() + block_idx;
		size_t vec, unit;
		for(size_t i = 1; i < salt_.size(); ++i){
			compute_indices(hash_ap(key_begin, length, salt_[i]), bit_index, bit);
			//the bit index is out of 512, the bit is out of 8.
			//we advance to the proper vector, then subscript to the proper byte
			std::tie(vec, unit) = get_vector_unit(bit_index);
			(*(block + vec))[unit] |= static_cast<base_type>(1) << bit;
		}
		++inserted_element_count_;
	}

	template <typename T>
	inline void insert(const T& t){
		insert(reinterpret_cast<const unsigned char*>(&t), sizeof(T));
	}

	inline virtual bool contains(const unsigned char* key_begin, const std::size_t length) const {
		size_t bit_index = 0;
		size_t bit = 0;
		//index in table with first byte of block
		size_t block_idx = get_block(hash_ap(key_begin, length, salt_[0])); 
		cell_type* block = bit_table_.get() + block_idx;
		size_t vec, unit;
		for(size_t i = 1; i < salt_.size(); ++i){
			compute_indices(hash_ap(key_begin, length, salt_[i]), bit_index, bit);
			std::tie(vec, unit) = get_vector_unit(bit_index);
			if( ((*(block + vec))[unit] & static_cast<base_type>(1) << bit) !=
				static_cast<base_type>(1) << bit)
			{
				return false;
			}
		}
		return true;
	}

	template <typename T>
	inline bool contains(const T& t) const {
		return contains(reinterpret_cast<const unsigned char*>(&t),static_cast<std::size_t>(sizeof(T)));
	}

	inline double effective_fpp() const {
		long double c = size() / element_count() ;
		long double lambda = block_size / c;
		long double fpp = 0;
		for(int i = 0; i < 3 * lambda; ++i){ //var = lambda, so 3*lambda should include most anything
			long double k = i;
			long double p_block = std::pow(lambda, k) * std::exp(-lambda) / std::tgammal(k+1);
			long double fpr_inner = std::pow(1.0 - std::exp(-1.0 * salt_.size() * i / block_size), 1.0 * salt_.size());
			fpp += p_block * fpr_inner;
		}
		return fpp;
	}

protected:
	inline virtual void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const {
		bit_index = hash % block_size; //which bit in the block?
      	bit       = bit_index % (bits_per_char * sizeof(base_type)); //which bit in the base type?
	}
};

class pattern_blocked_bf: public blocked_bloom_filter
{
protected:
	typedef std::unique_ptr<cell_type, std::function<void(cell_type*)>> pattern_type;
	static const size_t num_patterns = 65536; //4MiB
	pattern_type patterns;
public:
	pattern_blocked_bf(): blocked_bloom_filter(){}
	pattern_blocked_bf(const bloom_parameters& p): blocked_bloom_filter(p){
		void* ptr = 0;
		//we have num_patterns patterns, each with size block_size (in bits)
		int ret = posix_memalign(&ptr, block_size / bits_per_char,
			num_patterns * block_size / bits_per_char);
		if(ret != 0){
			throw std::bad_alloc();
		}
		patterns = pattern_type(static_cast<cell_type*>(ptr),
			[](cell_type* x){free(x);});
		std::uninitialized_fill_n(patterns.get(),
			num_patterns * block_size / bits_per_char / sizeof(cell_type),
			cell_type_zero);
		minion::Random rng;
		rng.Seed(random_seed_); //todo: check that seeding is proper for multiple rng instances
		//ie. the rng in kmersubsampler.
		std::uniform_int_distribution<> d(0, block_size-1);

		std::vector<size_t> possible_bits(block_size);
		std::iota(possible_bits.begin(), possible_bits.end(), 0);
		//begin with a fully shuffled array
		std::shuffle(possible_bits.begin(), possible_bits.end(), rng);
		for(size_t i = 0; i < num_patterns; ++i){
			//first index of the block containing the pattern
			size_t block_start_idx = i * (block_size / bits_per_char) / sizeof(cell_type);
			//sample salt_.size() bits
			for(int j = 0; j < salt_.size(); ++j){
				size_t sampled_bit = d(rng,
					std::uniform_int_distribution<>::param_type{j, block_size - 1});
				std::swap(possible_bits[j], possible_bits[sampled_bit]);
			}
			//set the bits in the appropriate pattern
			size_t vec, unit;
			for(size_t j = 0; j < salt_.size(); ++j){
				size_t sampled_bit_number = possible_bits[j]; //the bit out of [0, 511]
				//the index of the correct vector and byte within the block
				std::tie(vec, unit) = get_vector_unit(sampled_bit_number);
				size_t sampled_bit_idx = block_start_idx + vec;
				//set the index of the bit within the char
				(*(patterns.get() + sampled_bit_idx))[unit] |= (static_cast<base_type>(1) << (sampled_bit_number % (sizeof(base_type) * bits_per_char)));
			}
		}
	}
	//delete copy ctor
	pattern_blocked_bf(const pattern_blocked_bf&) = delete;
	//delete copy assignment
	pattern_blocked_bf& operator=(const pattern_blocked_bf&) = delete;
	//move ctor
	pattern_blocked_bf(pattern_blocked_bf&& o):
		blocked_bloom_filter(std::move(o)), patterns(std::move(o.patterns))
		{}

	//move function
	inline pattern_blocked_bf& operator=(pattern_blocked_bf&& o){
		if(this != &o){
			blocked_bloom_filter::operator=(std::move(o));
			patterns = std::move(o.patterns);
		}
		return *this;
	}

	inline virtual size_t get_pattern(const bloom_type& hash) const{
		size_t pattern_number = hash & (num_patterns - 1); //which block
		return pattern_number * (block_size / bits_per_char) / sizeof(cell_type); // how lighter does it
	}

	inline virtual void insert(const unsigned char* key_begin, const size_t& length){
		//index in table with first byte of block
		size_t block = get_block(hash_ap(key_begin, length, salt_[0]));
		size_t pattern = get_pattern(hash_ap(key_begin, length, salt_[1]));
		static_assert(block_size / bits_per_char / sizeof(cell_type) > 0,
			"Block size must be greater than or equal to size of cell type.");
		cell_type* bit_block = reinterpret_cast<cell_type*>(__builtin_assume_aligned(bit_table_.get() + block, block_size / bits_per_char));
		cell_type* pattern_block = reinterpret_cast<cell_type*>(__builtin_assume_aligned(patterns.get() + pattern, block_size / bits_per_char));
		for(size_t i = 0; i < block_size / bits_per_char / sizeof(cell_type); ++i){
			*(bit_block + i) |= *(pattern_block + i);
		}
		++inserted_element_count_;
	}

	template <typename T>
	inline void insert(const T& t)
	{
		// Note: T must be a C++ POD type.
		insert(reinterpret_cast<const unsigned char*>(&t),sizeof(T));
	}

	inline virtual bool contains(const unsigned char* key_begin, const size_t length) const{
		//index in table with first byte of block
		size_t block = get_block(hash_ap(key_begin, length, salt_[0]));
		size_t pattern = get_pattern(hash_ap(key_begin, length, salt_[1]));
		static_assert(block_size / bits_per_char / sizeof(cell_type) > 0,
			"Block size must be greater than or equal to size of cell type.");
		cell_type* bit_block = reinterpret_cast<cell_type*>(__builtin_assume_aligned(bit_table_.get() + block, block_size / bits_per_char));
		cell_type* pattern_block = reinterpret_cast<cell_type*>(__builtin_assume_aligned(patterns.get() + pattern, block_size / bits_per_char));
		for(size_t i = 0; i < block_size / bits_per_char / sizeof(cell_type); ++i){
			if(!_mm256_testc_si256((*(bit_block + i) & *(pattern_block + i)), //==
				*(pattern_block + i)))
			{
				return false;
			}
		}
		return true;
	}

	template <typename T>
	inline bool contains(const T& t) const
	{
		return contains(reinterpret_cast<const unsigned char*>(&t),static_cast<std::size_t>(sizeof(T)));
	}

	inline virtual bloom_type block_hash(const unsigned char* key_begin, const size_t& length) const{
		return hash_ap(key_begin, length, salt_[0]);
	}

	template <typename T>
	inline bloom_type block_hash(const T& t) const{
		return block_hash(reinterpret_cast<const unsigned char*>(&t),static_cast<std::size_t>(sizeof(T)));
	}

	inline virtual bloom_type pattern_hash(const unsigned char* key_begin, const size_t& length) const{
		return hash_ap(key_begin, length, salt_[1]);
	}

	template <typename T>
	inline bloom_type pattern_hash(const T& t) const{
		return pattern_hash(reinterpret_cast<const unsigned char*>(&t),static_cast<std::size_t>(sizeof(T)));
	}

	inline double effective_fpp() const {
		long double c = size() / element_count() ;
		long double lambda = block_size / c;
		long double fpp = 0;
		for(int i = 0; i < 3 * lambda; ++i){ //var = lambda, so 3*lambda should include most anything
			long double p_block = std::pow(lambda, i) * std::exp(-lambda) / std::tgammal(i+1);
			long double p_collision = 1.0l - std::pow(1.0l - 1.0l / (num_patterns), i);
			long double fpr_inner = std::pow(1.0l - std::exp(-1.0l * salt_.size() * i / block_size), 1.0l * salt_.size());
			fpr_inner = p_collision + (1.0l - p_collision) * fpr_inner;
			fpp += p_block * fpr_inner;
		}
		return fpp;
	}

};





//a class to hold an encoded kmer
class Kmer{
protected:
	size_t s; //num times kmer added to since last reset
	int k;
	uint64_t x[2]; //fwd and reverse
	uint64_t mask;
	uint64_t shift;
public:
	Kmer(int k): k(k), s(0), mask(k < 32 ? (1ULL<<k*2) - 1 : -1), shift((k-1)*2) {x[0] = x[1] = 0;}
	Kmer(const Kmer& o): k(o.k), s(o.s), mask(o.mask), shift(o.shift), x{o.x[0], o.x[1]} {}
	//add a character and return the number of times the kmer has been added to since last reset
	inline size_t push_back(char ch){
		int c = seq_nt16_int[seq_nt16_table[ch]];
		if (c < 4){
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			++s;
		} else this->reset(); // if there is an "N", restart
		return s;
	}
	//get encoded kmer
	inline uint64_t get() const{return x[0] < x[1] ? x[0] : x[1];} //min of x[0] and x[1]
	//get encoded prefix
	inline uint64_t prefix() const{return this->get()&((1<<PREFIXBITS)-1);}
	//empty the kmer and set s to 0
	inline void reset(){s = 0; x[0] = x[1] = 0;}
	//the number of times the kmer has been added to since the last reset
	inline size_t size() const{return s;}
	//whether the kmer has enough bases to be of length k
	inline bool valid() const{return (s >= k);}
	//the length of the kmer
	inline int ksize() const{return k;}
	inline operator std::string() const {
		std::string ret{};
		for(int i = 1; i <= k; ++i){
			ret.push_back( seq_nt16_str[seq_nt16_table['0' + ((x[0] & (3ULL << (2*(k-i)))) >> (2*(k-i)))]] );
		} 
		return ret;
	}
	inline explicit operator bool() const{return this->valid();}
};

//a bloom filter. TODO: make blocked; hold an array of bloom filters and
//delegate insert/query fn's to the appropriate one.
//calculating the right fpr may be a bit involved.
//if we do threads we can also add a mutex for insertion
class Bloom
{
public:
	typedef pattern_blocked_bf bloom_type;
	Bloom(unsigned long long int projected_element_count, double fpr, unsigned long long int seed = 0xA5A5A5A55A5A5A5AULL);
	Bloom(Bloom&& b) noexcept: bloom(std::move(b.bloom)){} //move ctor
	Bloom& operator=(Bloom&& o){bloom = std::move(o.bloom); return *this;} //move assign
	~Bloom();
	bloom_parameters params;
	bloom_type bloom;
	inline void insert(const Kmer& kmer){if(kmer.valid()){bloom.insert(kmer.get());}}
	template <typename T>
	inline void insert(const T& t){bloom.insert(t);}
	inline bool query(const Kmer& kmer) const {return (kmer.valid() && bloom.contains(kmer.get()));}
	inline double fprate() const {return bloom.effective_fpp();}
	inline unsigned long long inserted_elements() const {return bloom.element_count();}
	// inline double fprate() const {return bloom.GetActualFP();}
};

// typedef std::array<Bloom,(1<<PREFIXBITS)> bloomary_t;

std::array<std::vector<size_t>,2> overlapping_kmers_in_bf(std::string seq, const Bloom& b, int k = 31);

//return the total number of kmers in b
int nkmers_in_bf(std::string seq, const Bloom& b, int k);

//given a kmer, get the next character (in ACGT order) that would create a trusted
//kmer when appended and return it. Return 0 if none would be trusted.
//Set the flag to test in TGCA order instead.
char get_next_trusted_char(const bloom::Kmer& kmer, const Bloom& trusted, bool reverse_test_order = false);

//return the INCLUSIVE indices bounding the largest stretch of trusted sequence
//if the first value is -1, there are no trusted kmers.
//if the 2nd value is -1 (== std::string::npos), until the end of the string is trusted.
//thus the whole string being trusted looks like {0, std::string::npos}
//while no part of the string being trusted looks like {std::string::npos, std::string::npos}
std::array<size_t,2> find_longest_trusted_seq(std::string seq, const Bloom& b, int k);

//find the longest possible fix for the kmer at position (k-1) until the end
//return the best character (multiple in case of a tie), the index of the next untrusted base,
// and whether multiple corrections were considered for the fix.
//if the length of the fix character vector is 0, no fix was found and correction should end.
//If the sequence is reverse-complemented, set the flag to test bases in reverse order
// (TGCA) instead of (ACGT)
std::tuple<std::vector<char>, size_t, bool> find_longest_fix(std::string seq, const Bloom& t, int k, bool reverse_test_order = false);

//given the sampling rate, calculate the probability any kmer is in the array.
long double calculate_phit(const Bloom& bf, long double alpha);

//given the number of inserts and the desired fpr, calculate the total size of the hash needed
uint64_t numbits(uint64_t numinserts, long double fpr);

//given the desired fpr, calculate the number of hash fn's needed
//assuming we use the optimal number of bits
int numhashes(long double fpr);

//given a sequence and an anchor, see if the anchor could be improved by moving it to the left.
//this is used in the correction step.
//ensure anchor >= k before this function is called.
//return {anchor, multiple}, the location of the new anchor and whether multiple corrections
//were possible during anchor adjustment.
std::pair<size_t, bool> adjust_right_anchor(size_t anchor, std::string seq, const Bloom& trusted, int k);

//get the biggest consecutive trusted block, for creating a trusted anchor when
//one doesn't exist. If there are too many consecutive misses, the procedure
//will end early. 94518
int biggest_consecutive_trusted_block(std::string seq, const Bloom& trusted, int k, int current_len);

}

inline std::ostream& operator<< (std::ostream& stream, const bloom::Kmer& kmer){
	stream << std::string(kmer);
	return stream;
}

#endif
// FILE:/home/adam/code/kbbq/include/kbbq/covariateutils.hh 


#ifndef COVARIATEUTILS_H
#define COVARIATEUTILS_H
#define KBBQ_MAXQ 93

#include <vector>
#include <utility>
#include <cmath>
#include <random>
#include <array>
#include <limits>
#include "readutils.hh"
#include "recalibrateutils.hh"

//fwd declare
namespace readutils{
	class CReadData;
}

namespace covariateutils{

typedef std::vector<int> meanq_t;
typedef std::vector<int> rgdq_t;
typedef std::vector<std::vector<int>> qscoredq_t;
typedef std::vector<std::vector<std::array<std::vector<int>,2>>> cycledq_t;
typedef std::vector<std::vector<std::vector<int>>> dinucdq_t;

struct dq_t
{
	meanq_t meanq;
	rgdq_t rgdq;
	qscoredq_t qscoredq;
	cycledq_t cycledq;
	dinucdq_t dinucdq;
};

typedef std::vector<int> prior1_t; //1 prior for each rg
typedef std::vector<std::vector<int>> prior2_t; //1 prior for each rg -> q pair

class NormalPrior{
public:
	static std::vector<long double> normal_prior;
	NormalPrior(){}
	~NormalPrior(){}
	static long double get_normal_prior(size_t j);
};

/*
def _logpmf(self, x, n, p):
    k = floor(x)
    combiln = (gamln(n+1) - (gamln(k+1) + gamln(n-k+1)))
    return combiln + special.xlogy(k, p) + special.xlog1py(n-k, -p)
*/

inline long double log_binom_pmf(unsigned long long k, unsigned long long n, long double p){
	// k+=1; //+1/+2 correction
	// n+=2;
	long double coefficient = (std::lgamma(n+1) - (std::lgamma(k+1) + std::lgamma(n-k+1)));
	return coefficient + (long double)k * std::log(p) + (long double)(n-k) * std::log1p(-p);
}

inline std::vector<long double> log_binom_cdf(unsigned long long k, long double p){
	std::vector<long double> ret(k+1);
	ret[0] = log_binom_pmf(0,k,p);
	for(int i = 1; i <= k; ++i){
		ret[i] = std::log(std::exp(ret[i-1]) + std::exp(log_binom_pmf(i,k,p)));
	}
	return ret;
}

//return vector of each k (k < n) 
inline std::vector<int> calculate_thresholds(unsigned long long k, long double p, long double quartile = .995l){
	std::vector<int> threshold(k+1,0);
	threshold[0] = 0;
	for(int i = 1; i <= k; ++i){
		std::vector<long double> cdf = log_binom_cdf(i,p);
		for(int j = 0; j < cdf.size(); ++j){
			long double logp = cdf[j];
			if(logp >= std::log(quartile)){
				threshold[i] = j;
				break;
			}
		}
	}
	return threshold;
}

inline bool nt_is_not_n(char c){
	return seq_nt16_int[seq_nt16_table[c]] < 4;
}

//ensure there's no N!!!!
inline int8_t dinuc_to_int(char first, char second){
	int8_t f = seq_nt16_int[seq_nt16_table[first]];
	int8_t s = seq_nt16_int[seq_nt16_table[second]];
	return 15 & ((f << 2) | s); // 1111 & (xx00 | 00xx)
}

inline std::array<char, 2> int_to_dinuc(int8_t dinuc){
	std::array<char,2> x = {std::to_string(dinuc >> 2)[0], std::to_string(dinuc | 3)[0]};
	return x;
}

typedef std::array<unsigned long long, 2> covariate_t;

class CCovariate: public std::vector<covariate_t>
{
public:
	CCovariate(): std::vector<covariate_t>(){}
	CCovariate(size_t len): std::vector<covariate_t>(len) {}
	void increment(size_t idx, covariate_t value);
	void increment(size_t idx, unsigned long long err, unsigned long long total);
};

class CRGCovariate: public CCovariate
{
public:
	CRGCovariate(){}
	CRGCovariate(size_t len): CCovariate(len){}
	void consume_read(const readutils::CReadData& read);
	rgdq_t delta_q(prior1_t prior);
};

class CQCovariate: public std::vector<CCovariate>
{
public:
	CQCovariate(){}
	CQCovariate(size_t rgs, size_t qlen): std::vector<CCovariate>(rgs, CCovariate(qlen)){}
	void consume_read(const readutils::CReadData& read);
	qscoredq_t delta_q(prior1_t prior);
};

typedef std::array<CCovariate,2> cycle_t;

//The first Covariate is for fwd reads, the 2nd Covariate is for reverse reads.
class CCycleCovariate: public std::vector<std::vector<cycle_t>>
{
public:
	CCycleCovariate(){}
	CCycleCovariate(size_t rgs, size_t qlen, size_t cylen):
		std::vector<std::vector<cycle_t>>(rgs, std::vector<cycle_t>(qlen, cycle_t({CCovariate(cylen),CCovariate(cylen)})))
		{}
	void consume_read(const readutils::CReadData& read);
	cycledq_t delta_q(prior2_t prior);
};

class CDinucCovariate: public std::vector<std::vector<CCovariate>>
{
public:
	CDinucCovariate(){}
	CDinucCovariate(size_t rgs, size_t qlen, size_t dilen):
		std::vector<std::vector<CCovariate>>(rgs, std::vector<CCovariate>(qlen, CCovariate(dilen)))
		{}
	void consume_read(const readutils::CReadData& read, int minscore = 6);
	dinucdq_t delta_q(prior2_t prior);
};

class CCovariateData
{
// protected:
// 	CRGCovariate rgcov;
// 	CQCovariate qcov;
// 	CCycleCovariate cycov;
// 	CDinucCovariate dicov;
public:
	CRGCovariate rgcov;
	CQCovariate qcov;
	CCycleCovariate cycov;
	CDinucCovariate dicov;
	CCovariateData(){};
	void consume_read(readutils::CReadData& read, int minscore = 6);
	dq_t get_dqs();
};

}
#endif
// FILE:/home/adam/code/kbbq/include/kbbq/recalibrateutils.hh 


#ifndef RECALIBRATEUTILS_H
#define RECALIBRATEUTILS_H

#include <vector>
#include <string>
#include <cmath>
#include "bloom.hh"
#include "htsiter.hh"
#include "covariateutils.hh"

//debug
#include <fstream>


//fwd declare covariateutils stuff
namespace covariateutils{
	class CCovariateData;
	struct dq_t;
}

namespace htsiter{
	class KmerSubsampler;
	class HTSFile;
}

namespace recalibrateutils{

//subsample kmers, hash them, and add them to the bloom filter
void subsample_kmers(htsiter::KmerSubsampler& s, bloom::Bloom& sampled);

//get some reads from a file, whether a kmer is trusted and put it in a cache.
//this can probably be parallelized if needed because the reads are independent
void find_trusted_kmers(htsiter::HTSFile* file, bloom::Bloom& trusted,
	const bloom::Bloom& sampled, std::vector<int> thresholds, int k);

inline long double q_to_p(int q){return std::pow(10.0l, -((long double)q / 10.0l));}
inline int p_to_q(long double p, int maxscore = 42){return p > 0 ? (int)(-10 * std::log10(p)) : maxscore;}

//get covariate data using the trusted kmers
//this can be parallelized easily since each read is independent
covariateutils::CCovariateData get_covariatedata(htsiter::HTSFile* file, const bloom::Bloom& trusted, int k);

//recalibrate all reads given the CovariateData
void recalibrate_and_write(htsiter::HTSFile* in, const covariateutils::dq_t& dqs, std::string outfn);

}

#endif
// FILE:/home/adam/code/kbbq/include/kbbq/htsiter.hh 


#ifndef KBBQ_HTSITER_H
#define KBBQ_HTSITER_H

#include <iterator>
#include <string>
#include <random>
#include <algorithm>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/thread_pool.h>
#include <minion.hpp>
#include "readutils.hh"
#include "kseq.hh"
#include "bloom.hh"
#include <cstdlib>

//This is defined starting in version 1.10
#ifndef HTS_VERSION
#define HTS_VERSION 0 
#endif

static_assert(HTS_VERSION >= 101000, "Your version of htslib is out of date. KBBQ requires version >= 1.10.");

#ifndef NDEBUG
#define KBBQ_USE_RAND_SAMPLER
#endif

//fwd declare
namespace readutils{
	class CReadData;
	std::string bam_seq_str(bam1_t* bamrecord);
}

//unsigned char* s = bam_get_seq(bamrecord);
namespace htsiter{

class HTSFile{
public:
	virtual ~HTSFile(){}
	virtual int next()=0;
	virtual std::string next_str()=0;
	virtual readutils::CReadData get()=0;
	virtual void recalibrate(const std::vector<uint8_t>& qual)=0;
	virtual int open_out(std::string filename)=0; //open an output file so it can be written to later.
	virtual int write()=0; //write the current read to the opened file.
};

class BamFile: public HTSFile{
public:
	samFile *sf;
	// hts_idx_t *idx;
	sam_hdr_t *h;
	// hts_itr_t *itr;
	bam1_t *r;
	samFile *of;
	bool use_oq;
	bool set_oq;
	htsThreadPool* tp;
	BamFile(std::string filename, htsThreadPool* tp, bool use_oq = false, bool set_oq = false):
		use_oq(use_oq), set_oq(set_oq), tp(tp){
		r = bam_init1();
		sf = sam_open(filename.c_str(), "r");
		if(tp->pool && hts_set_thread_pool(sf, tp) != 0){
			std::cerr << "Couldn't attach thread pool to file " << filename << std::endl;
		};
	    h = sam_hdr_read(sf);
	    //TODO: support iteration with index?
	    // idx = sam_index_load(sf, filename.c_str());
	   	//TODO: throw when index can't be found
	    // itr = sam_itr_queryi(idx, HTS_IDX_START, 0, 0); //iterate over whole file;
	    of = NULL; //this can be opened later with open_out
	}
	~BamFile(){
		if(r != NULL){bam_destroy1(r);}
		// if(itr != NULL){sam_itr_destroy(itr);}
		if(sf != NULL){sam_close(sf);}
		if(of != NULL){sam_close(of);}
		if(h != NULL){sam_hdr_destroy(h);}
		// if(idx != NULL){hts_idx_destroy(idx);}
	}

	// to use: while (ret = BamFile.next() >= 0){//do something with this->r}
	// >= 0 on success, -1 on EOF, <-1 on error
	int next();
	// return next read sequence as a string. if there are no more, return the empty string.
	std::string next_str();
	//
	readutils::CReadData get();
	//
	void recalibrate(const std::vector<uint8_t>& qual);
	// TODO:: add a PG tag to the header
	int open_out(std::string filename);
	//
	int write();
	//
}; //end of BamFile class

class FastqFile: public HTSFile
{
public:
	BGZF* fh;
	kseq::kseq_t* r;
	BGZF* ofh;
	htsThreadPool* tp;
	FastqFile(std::string filename, htsThreadPool* tp): ofh(NULL), tp(tp){
		fh = bgzf_open(filename.c_str(),"r");
		if(tp->pool && bgzf_thread_pool(fh, tp->pool, tp->qsize) < 0){
			std::cerr << "Couldn't attach thread pool to file " << filename << std::endl;
		}
		r = kseq::kseq_init(fh);
	};
	~FastqFile(){
		if(r != NULL){kseq::kseq_destroy(r);}
		if(fh != NULL){bgzf_close(fh);}
		if(ofh != NULL){bgzf_close(ofh);}
	}
	int next();
	std::string next_str();
	readutils::CReadData get();
	void recalibrate(const std::vector<uint8_t>& qual);
	int open_out(std::string filename);
	int write();
};

class KmerSubsampler{
public:
	HTSFile* file;
	minion::Random rng;
	std::bernoulli_distribution d;
	std::string readseq = "";
	std::vector<bloom::Kmer> kmers;
	bloom::Kmer kmer;
	size_t cur_kmer = 0;
	size_t total_kmers = 0;
	bool not_eof = true;
	int k;
	KmerSubsampler(HTSFile* file): KmerSubsampler(file, KBBQ_MAX_KMER){}
	KmerSubsampler(HTSFile* file, int k): KmerSubsampler(file, k, .15){}
	KmerSubsampler(HTSFile* file, int k, double alpha): KmerSubsampler(file, k, alpha,  minion::create_seed_seq().GenerateOne()){}
	KmerSubsampler(HTSFile* file, int k, double alpha, uint64_t seed): file(file), k(k), kmer(k), d(alpha) {rng.Seed(seed); std::cerr << "p: " << d.p() << std::endl;} //todo remove srand
	
	//return the next kmer
	//once the file is finished iterating and there are no remaining kmers,
	//return an empty kmer. You should check the not_eof flag.
	//the kmer returned is not guaranteed to be valid.
	bloom::Kmer next_kmer();

	//return the next kmer that survives sampling.
	//it is not gauranteed to be valid.
	bloom::Kmer next();
	inline bloom::Kmer operator()() {return this->next();}
	inline explicit operator bool() const{return not_eof;}
};

}

#endif
