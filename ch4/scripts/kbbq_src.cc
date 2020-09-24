
// FILE:/home/adam/code/cbbq/src/kbbq/test.cc 


#include <bitset>

int test(){
	bloom_parameters p{};
	p.projected_element_count = 100;
	p.false_positive_probability = .05;
	p.compute_optimal_parameters();
	bloom::pattern_blocked_bf bf(p);
	std::cerr << "Table size: " << bf.size() << std::endl;
	std::cerr << "Block size: " << bf.block_size << std::endl;
	std::cerr << "Block bytes: " << bf.block_size / bits_per_char << std::endl;
	std::cerr << "Pattern bytes: " << bf.num_patterns * bf.block_size / bits_per_char << std::endl;
	std::cerr << "Num hashes: " << bf.hash_count() << std::endl;
	std::cerr << "Optimal K: " << p.optimal_parameters.number_of_hashes << std::endl;
	for(size_t i = 0; i < 10; ++i){
		// std::cerr << (void*)bf.patterns.get() << std::endl;
		std::cerr << i << "(" << (void*)&bf.patterns.get()[2 * i] << ")" << ": ";
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.patterns.get()[2*i][j]) << ", ";
		}
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.patterns.get()[2*i+1][j]) << ", ";
		}
		std::cerr << std::endl;
	}
	uint64_t kmer = 21345423534512;
	uint64_t gibberish = 123543451243ull;
	bf.insert(kmer);
	std::cerr << "kmer inserted into block: " << bf.get_block(bf.block_hash(kmer)) << std::endl;
	size_t kpattern = bf.get_pattern(bf.pattern_hash(kmer));
	std::cerr << "kmer pattern number: " << kpattern << std::endl;
	std::cerr << "gibberish block: " << bf.get_block(bf.block_hash(gibberish)) << std::endl;
	std::cerr << "gibberish pattern number: " << bf.get_pattern(bf.pattern_hash(gibberish)) << std::endl;
	std::cerr << "kmerp: ";
	for(size_t i = kpattern; i < kpattern + 2;++i){
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.patterns.get()[i][j]) << ", ";	
		}
	}
	std::cerr << std::endl;

	std::cerr << "Table: \n";
	for(size_t i = 0; i < bf.size() / bits_per_char / 64; ++i){
		// if(i % (bf.block_size/bits_per_char) == 0){
		// 	std::cerr << "\nBLOCK" << std::endl;
		// }
		std::cerr << "block: ";
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.bit_table_.get()[2*i][j]) << ", ";
		}
		for(size_t j = 0; j < 4; ++j){
			std::cerr << std::bitset<64>(bf.bit_table_.get()[2*i+1][j]) << ", ";
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
	std::cerr << "bf contains kmer ? " << bf.contains(kmer) << std::endl;
	std::cerr << "bf contains gibberish ? " << bf.contains(gibberish) << std::endl;
	assert(bf.contains(kmer));


	return 0;
}
// FILE:/home/adam/code/cbbq/src/kbbq/htsiter.cc 


#include "htsiter.hh"

namespace htsiter{

int BamFile::next(){return sam_read1(sf, h, r);}//return sam_itr_next(sf, itr, r);
	// return next read as a string. if there are no more, return the empty string.
std::string BamFile::next_str(){return this->next() >= 0 ? readutils::bam_seq_str(r) : "";}
//
readutils::CReadData BamFile::get(){return readutils::CReadData(this->r, use_oq);}
//
void BamFile::recalibrate(const std::vector<uint8_t>& qual){
	uint8_t* q = bam_get_qual(this->r);
	if(set_oq){
		std::string qstr;
		std::transform(q, q + this->r->core.l_qseq, std::back_inserter(qstr),
			[](uint8_t c) -> char {return c + 33;}); //qual value to actual str
		//returns 0 on success, -1 on fail. We should consider throwing if it fails.
		if(bam_aux_update_str(this->r, "OQ", qstr.length()+1, qstr.c_str()) != 0){
			if(errno == ENOMEM){
				std::cerr << "Insufficient memory to expand bam record." << std::endl;
			} else if(errno == EINVAL){
				std::cerr << "Tag data is corrupt. Repair the tags and try again." << std::endl;
			}
			throw std::invalid_argument("Unable to update OQ tag.");
		}
	}
	if(bam_is_rev(this->r)){
		std::reverse_copy(qual.begin(), qual.end(), q);
	} else {
		std::copy(qual.begin(), qual.end(), q);
	}
	// for(int i = 0; i < this->r->core.l_qseq; ++i){
	// 	q[i] = (char)qual[i];
	// }
}
// TODO:: add a PG tag to the header
int BamFile::open_out(std::string filename){this->of = sam_open(filename.c_str(), "wb"); return sam_hdr_write(this->of, this->h);}
//
int BamFile::write(){return sam_write1(this->of, this->h, this->r);}

// FastqFile class

int FastqFile::next(){
	return kseq_read(r);
}

std::string FastqFile::next_str(){
	return this->next() >= 0? std::string(this->r->seq.s): "";
}

readutils::CReadData FastqFile::get(){
	return readutils::CReadData(this->r);
}

void FastqFile::recalibrate(const std::vector<uint8_t>& qual){
	for(int i = 0; i < this->r->qual.l; ++i){
		this->r->qual.s[i] = (char)(qual[i]+33);
	}
}

int FastqFile::open_out(std::string filename){
	ofh = bgzf_open(filename.c_str(),"w"); //mode should be "wu" if uncompressed output is desired.
	return ofh == 0 ? -1 : 0;
}

int FastqFile::write(){
	std::string name = ks_c_str(&this->r->name);
	std::string seq = ks_c_str(&this->r->seq);
	std::string comment = ks_c_str(&this->r->comment);
	std::string qual = ks_c_str(&this->r->qual);
	// std::string name = this->r->name.s ? std::string(this->r->name.s) : "";
	// std::string seq = this->r->seq.s ? std::string(this->r->seq.s) : "";
	// std::string comment = this->r->comment.s ? std::string(this->r->comment.s) : "";
	// std::string qual = this->r->qual.s ? std::string(this->r->qual.s) : "";
	std::string s("@" + name + "\n" + seq + "\n+" + comment + "\n" + qual + "\n");
	return bgzf_write(ofh, s.c_str(), s.length());
}

// KmerSubsampler
bloom::Kmer KmerSubsampler::next_kmer(){
	if(cur_kmer < kmers.size()){
		return kmers[cur_kmer++]; //return the current kmer and advance
	} else {
		readseq = file->next_str();
		kmer.reset();
		if(readseq.empty()){
			this->not_eof = false;
			return kmer; //no more sequences
		} else {
			kmers.clear();
			for(size_t i = 0; i < readseq.length(); ++i){
				kmer.push_back(readseq[i]);
				if(i >= k-1){
					kmers.push_back(kmer);
				}
			}
			cur_kmer = 0; //reset current kmer
			total_kmers += kmers.size();
			return this->next_kmer(); //try again
		}
	}
}

bloom::Kmer KmerSubsampler::next(){
	bloom::Kmer kmer = this->next_kmer();
	if(this->not_eof){
#ifdef KBBQ_USE_RAND_SAMPLER
		double p = std::rand() / (double)RAND_MAX;
		if(p < this->d.p()){
#else
		if(d(rng)){ //sampled
#endif
			return kmer;
		}
		else{ //try again
			return this->next();
		}
	}
	return kmer; // empty
}








}

// FILE:/home/adam/code/cbbq/src/kbbq/cbbq.cc 


#include "htsiter.hh"
#include "bloom.hh"
#include "covariateutils.hh"
#include "recalibrateutils.hh"
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <htslib/hfile.h>
#include <getopt.h>
#include <cassert>
#include <functional>
#include <iomanip>
#include <ctime>

#ifndef NDEBUG
#define KBBQ_USE_RAND_SAMPLER
#endif

//opens file filename and returns a unique_ptr to the result.
std::unique_ptr<htsiter::HTSFile> open_file(std::string filename, bool is_bam = true, bool use_oq = false, bool set_oq = false){
	std::unique_ptr<htsiter::HTSFile> f(nullptr);
	if(is_bam){
		f = std::move(std::unique_ptr<htsiter::BamFile>(new htsiter::BamFile(filename, use_oq, set_oq)));
		// f.reset(new htsiter::BamFile(filename));
	} else {
		f = std::move(std::unique_ptr<htsiter::FastqFile>(new htsiter::FastqFile(filename)));
		// f.reset(new htsiter::FastqFile(filename));
	}
	return f;
}



template<typename T>
std::ostream& operator<< (std::ostream& stream, const std::vector<T>& v){
	stream << "[";
	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(stream, ", "));
	stream << "]";
	return stream;
}

template<typename T>
void print_vec(const std::vector<T>& v){
	std::cerr << v;
}

std::ostream& put_now(std::ostream& os){
	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);
	return os << std::put_time(&tm, "[%F %T %Z]");
}

int check_args(int argc, char* argv[]){
	if(argc < 2){
		std::cerr << put_now << " Usage: " << argv[0] << " input.[bam,fq]" << std::endl;
		return 1;
	} else {
		std::cerr << put_now << "  Selected file: " << std::string(argv[1]) << std::endl;
		return 0;
	}
}

//long option, required arg?, flag, value
struct option long_options[] = {
	{"ksize",required_argument,0,'k'}, //default: 31
	{"use-oq",no_argument,0,'u'}, //default: off
	{"set-oq",no_argument,0,'s'}, //default: off
	{"genomelen",required_argument,0,'g'}, //estimated for bam input, required for fastq input
	{"coverage",required_argument,0,'c'}, //default: estimated
	{"fixed",required_argument,0,'f'}, //default: none
	{"alpha",required_argument,0,'a'}, //default: 7 / coverage
#ifndef NDEBUG
	{"debug",required_argument,0,'d'},
#endif
	{0, 0, 0, 0}
};

int main(int argc, char* argv[]){
	int k = 32;
	long double alpha = 0;
	uint64_t genomelen = 0; //est w/ index with bam, w/ fq estimate w/ coverage
	uint coverage = 0; //if not given, will be estimated.
	uint32_t seed = 0; 
	bool set_oq = false;
	bool use_oq = false;
	std::string fixedinput = "";

	int opt = 0;
	int opt_idx = 0;
#ifndef NDEBUG
	std::string kmerlist("");
	std::string trustedlist("");
#endif
	while((opt = getopt_long(argc,argv,"k:usg:c:f:a:d:",long_options, &opt_idx)) != -1){
		switch(opt){
			case 'k':
				k = std::stoi(std::string(optarg));
				if(k <= 0 || k > KBBQ_MAX_KMER){
					std::cerr << put_now << "  Error: k must be <= " << KBBQ_MAX_KMER << " and > 0." << std::endl;
				}
				break;
			case 'u':
				use_oq = true;
				break;
			case 's':
				set_oq = true;
				break;
			case 'g':
				genomelen = std::stoull(std::string(optarg));
				break;
			case 'c':
				coverage = std::stoul(std::string(optarg));
				break;
			case 'f':
				fixedinput = std::string(optarg);
				break;
			case 'a':
				alpha = std::stold(std::string(optarg));
				break;
#ifndef NDEBUG
			case 'd': {
				std::string optstr(optarg);
				std::istringstream stream(optstr);
				std::getline(stream, kmerlist, ',');
				std::getline(stream, trustedlist, ',');
				break;
			}
#endif
			case '?':
			default:
				std::cerr << put_now << "  Unknown argument " << (char)opt << std::endl;
				return 1;
				break;
		}
	}

	std::string filename("-");
	if(optind < argc){
		filename = std::string(argv[optind]);
		while(++optind < argc){
			std::cerr << put_now << " Warning: Extra argument " << argv[optind] << " ignored." << std::endl;
		}
	}

	long double sampler_desiredfpr = 0.01; //Lighter uses .01
	long double trusted_desiredfpr = 0.0005; // and .0005

	//see if we have a bam
	htsFormat fmt;
	hFILE* fp = hopen(filename.c_str(), "r");
	if (hts_detect_format(fp, &fmt) < 0) {
		//error
		std::cerr << put_now << " Error opening file " << filename << std::endl;
		hclose_abruptly(fp);
		return 1;
	}
	bool is_bam = true;
	if(fmt.format == bam || fmt.format == cram){
		is_bam = true;
	} else if (fmt.format == fastq_format){
		is_bam = false;
	} else {
		//error
		std::cerr << put_now << " Error: File format must be bam, cram, or fastq." << std::endl;
		hclose_abruptly(fp);
		return 1;
	}
	std::unique_ptr<htsiter::HTSFile> file;
	covariateutils::CCovariateData data;

if(fixedinput == ""){ //no fixed input provided

	if(genomelen == 0){
		if(is_bam){
			std::cerr << put_now << " Estimating genome length" << std::endl;
			samFile* sf = hts_hopen(fp, filename.c_str(), "r");
			sam_hdr_t* h = sam_hdr_read(sf);
			for(int i = 0; i < sam_hdr_nref(h); ++i){
				genomelen += sam_hdr_tid2len(h, i);
			}
			sam_hdr_destroy(h);
			hts_close(sf);
			if(genomelen == 0){
				std::cerr << put_now << " Header does not contain genome information." <<
					" Unable to estimate genome length; please provide it on the command line" <<
					" using the --genomelen option." << std::endl;
				return 1;
			} else {
				std::cerr << put_now << " Genome length is " << genomelen <<" bp." << std::endl;
			}
		} else {
			std::cerr << put_now << " Error: --genomelen must be specified if input is not a bam." << std::endl;
		}
	} else {
		if(hclose(fp) != 0){
			std::cerr << put_now << " Error closing file!" << std::endl;
		}
	}
	
	//alpha not provided, coverage not provided
	if(alpha == 0){
		std::cerr << put_now << " Estimating alpha." << std::endl;
		if(coverage == 0){
			std::cerr << put_now << " Estimating coverage." << std::endl;
			uint64_t seqlen = 0;
			file = std::move(open_file(filename, is_bam, use_oq, set_oq));
			std::string seq("");
			while((seq = file->next_str()) != ""){
				seqlen += seq.length();
			}
			if (seqlen == 0){
				std::cerr << put_now << " Error: total sequence length in file " << filename <<
					" is 0. Check that the file isn't empty." << std::endl;
				return 1;
			}
			std::cerr << put_now << " Total Sequence length: " << seqlen << std::endl;
			std::cerr << put_now << " Genome length: " << genomelen << std::endl;
			coverage = seqlen/genomelen;
			std::cerr << put_now << " Estimated coverage: " << coverage << std::endl;
			if(coverage == 0){
				std::cerr << put_now << " Error: estimated coverage is 0." << std::endl;
				return 1;
			}
		}
		alpha = 7.0l / (long double)coverage; // recommended by Lighter authors		
	}

	if(coverage == 0){ //coverage hasn't been estimated but alpha is given
		coverage = 7.0l/alpha;
	}

	file = std::move(open_file(filename, is_bam, use_oq, set_oq));

	std::cerr << put_now << " Sampling kmers at rate " << alpha << std::endl;
	recalibrateutils::kmer_cache_t subsampled_hashes;

	//in the worst case, every kmer is unique, so we have genomelen * coverage kmers
	//then we will sample proportion alpha of those.
	unsigned long long int approx_kmers = genomelen*coverage*alpha;
	bloom::Bloom subsampled(approx_kmers, sampler_desiredfpr); //lighter uses 1.5 * genomelen
	bloom::Bloom trusted(approx_kmers, trusted_desiredfpr);

	if(seed == 0){
		seed = minion::create_seed_seq().GenerateOne();
	}
	std::cerr << put_now << " Seed: " << seed << std::endl ;

	//sample kmers here.
#ifdef KBBQ_USE_RAND_SAMPLER
	std::srand(seed); //lighter uses 17
#endif
	htsiter::KmerSubsampler subsampler(file.get(), k, alpha, seed);
	//load subsampled bf.
	//these are hashed kmers.
	subsampled_hashes = recalibrateutils::subsample_kmers(subsampler);

	//report number of sampled kmers
	uint64_t nsampled = 0;
	for(std::vector<uint64_t>& v : subsampled_hashes){
		nsampled += v.size();
	}
	std::cerr << put_now << " Sampled " << nsampled << " valid kmers." << std::endl;
	recalibrateutils::add_kmers_to_bloom(subsampled_hashes, subsampled);

#ifndef NDEBUG
	//ensure kmers are properly sampled
	if(kmerlist != ""){
		std::ifstream kmersin(kmerlist);
		bloom::Kmer kin(k);
		for(std::string line; std::getline(kmersin, line); ){
			kin.reset();
			for(char c: line){
				kin.push_back(c);
			}
			if(kin.valid()){
				assert(subsampled.query(kin));
			}
		}
	}
#endif


	//calculate thresholds
	long double fpr = subsampled.fprate();
	std::cerr << put_now << " Approximate false positive rate: " << fpr << std::endl;
	if(fpr > .15){
		std::cerr << put_now << " Error: false positive rate is too high. " <<
			"Increase genomelen parameter and try again." << std::endl;
		return 1;
	}

	long double p = bloom::calculate_phit(subsampled, alpha);
	std::vector<int> thresholds = covariateutils::calculate_thresholds(k, p);
#ifndef NDEBUG
	std::vector<int> lighter_thresholds = {0, 1, 2, 3, 4, 4, 5, 5, 6, 6,
		7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 12, 13, 13, 14, 14, 15, 15, 15, 16, 16, 17};

	std::cerr << put_now << " Thresholds: [ " ;
	std::copy(thresholds.begin(), thresholds.end(), std::ostream_iterator<int>(std::cerr, " "));
	std::cerr << "]" << std::endl;
	std::cerr << put_now << " Lighter Th: [ " ;
	std::copy(lighter_thresholds.begin(), lighter_thresholds.end(), std::ostream_iterator<int>(std::cerr, " "));
	std::cerr << "]" << std::endl;
	assert(lighter_thresholds == thresholds);
#endif

	std::vector<long double> cdf = covariateutils::log_binom_cdf(k,p);

	std::cerr << put_now << " log CDF: [ " ;
	for(auto c : cdf){std::cerr << c << " ";}
	std::cerr << "]" << std::endl;

	//get trusted kmers bf using subsampled bf
	std::cerr << put_now << " Finding trusted kmers" << std::endl;

	file = std::move(open_file(filename, is_bam, use_oq, set_oq));
	recalibrateutils::kmer_cache_t trusted_hashes = 
		recalibrateutils::find_trusted_kmers(file.get(), subsampled, thresholds, k);
	recalibrateutils::add_kmers_to_bloom(trusted_hashes, trusted);

#ifndef NDEBUG
// check that all kmers in trusted list are actually trusted in our list.
// it seems that lighter has quite a few hash collisions that end up making
// it trust slightly more kmers than it should

if(trustedlist != ""){
	std::ifstream kmersin(trustedlist);
	bloom::Kmer kin(k);
	for(std::string line; std::getline(kmersin, line); ){
		// std::cerr << "Trusted kmer: " << line << std::endl;
		kin.reset();
		for(char c: line){
			kin.push_back(c);
		}
		if(!trusted.query(kin)){
			std::cerr << "Trusted kmer not found!" << std::endl;
			std::cerr << "Line: " << line << std::endl;
			std::cerr << "Kmer: " << kin << std::endl;
		}
		assert(trusted.query(kin));
	}
}
#endif

	//use trusted kmers to find errors
	std::cerr << put_now << " Finding errors" << std::endl;
	file = std::move(open_file(filename, is_bam, use_oq, set_oq));
	data = recalibrateutils::get_covariatedata(file.get(), trusted, k);
} else { //use fixedfile to find errors
	std::cerr << put_now << " Using fixed file to find errors." << std::endl;
	file = std::move(open_file(filename, is_bam, use_oq, set_oq));
	std::unique_ptr<htsiter::HTSFile> fixedfile = std::move(open_file(fixedinput, is_bam, use_oq, set_oq));
	while(file->next() >= 0 && fixedfile->next() >= 0){
		readutils::CReadData read = file->get();
		readutils::CReadData fixedread = fixedfile->get();
		std::transform(read.seq.begin(), read.seq.end(), fixedread.seq.begin(),
			read.errors.begin(), std::not_equal_to<char>{});
		data.consume_read(read);
	}
}



	std::vector<std::string> rgvals(readutils::CReadData::rg_to_int.size(), "");
	for(auto i : readutils::CReadData::rg_to_int){
		rgvals[i.second] = i.first;
	}

	std::cerr << put_now << " Covariate data:" << std::endl;
	std::cerr << "rgcov:";
	for(int i = 0; i < data.rgcov.size(); ++i){ //rgcov[rg][0] = errors
		std::cerr << i << ": " << rgvals[i] << " {" << data.rgcov[i][0] << ", " << data.rgcov[i][1] << "}" << std::endl;
	}
	std::cerr << "qcov:" << "(" << data.qcov.size() << ")" << std::endl;
	for(int i = 0; i < data.qcov.size(); ++i){
		std::cerr << i << " " << rgvals[i] << "(" << data.qcov[i].size() << ")" << ": [";
		for(int j = 0; j < data.qcov[i].size(); ++j){
			if(data.qcov[i][j][1] != 0){
				std::cerr << j << ":{" << data.qcov[i][j][0] << ", " << data.qcov[i][j][1] << "} ";
			}
		}
		std::cerr << "]" << std::endl;
	}


	//recalibrate reads and write to file
	std::cerr << put_now << " Training model" << std::endl;
	covariateutils::dq_t dqs = data.get_dqs();

	std::cerr << put_now << " dqs:\n" << "meanq: ";
	print_vec<int>(dqs.meanq);
	std::cerr << "\nrgdq:" << std::endl;
	for(int i = 0; i < dqs.rgdq.size(); ++i){
		std::cerr << rgvals[i] << ": " << dqs.rgdq[i] << " (" << dqs.meanq[i] + dqs.rgdq[i] << ")" << std::endl;
	}
	std::cerr << "qscoredq:" << std::endl;
	for(int i = 0; i < dqs.qscoredq.size(); ++i){
		for(int j = 0; j < dqs.qscoredq[i].size(); ++j){
			if(data.qcov[i][j][1] != 0){
				std::cerr << rgvals[i] << ", " << "q = " << j << ": " << dqs.qscoredq[i][j] << " (" <<
					dqs.meanq[i] + dqs.rgdq[i] + dqs.qscoredq[i][j] << ") " << 
					data.qcov[i][j][1] << " " << data.qcov[i][j][0] << std::endl;
			}
		}
	}
	std::cerr << "cycledq:" << std::endl;
	for(int i = 0; i < dqs.cycledq.size(); ++i){
		for(int j = 0; j < dqs.cycledq[i].size(); ++j){
			if(data.qcov[i][j][1] != 0){
				for(size_t k = 0; k < dqs.cycledq[i][j].size(); ++k){
					for(int l = 0; l < dqs.cycledq[i][j][k].size(); ++l){
						std::cerr << rgvals[i] << ", " << "q = " << j << ", cycle = " << (k ? -(l+1) : l+1) << ": " << dqs.cycledq[i][j][k][l] << " (" <<
							dqs.meanq[i] + dqs.rgdq[i] + dqs.qscoredq[i][j] + dqs.cycledq[i][j][k][l] << ") " << 
							data.cycov[i][j][k][l][1] << " " << data.cycov[i][j][k][l][0] << std::endl;
					}
				}
			}
		}
	}
	std::cerr << "dinucdq:" << std::endl;
	for(int i = 0; i < dqs.dinucdq.size(); ++i){
		for(int j = 0; j < dqs.dinucdq[i].size(); ++j){
			if(data.qcov[i][j][1] != 0){
				for(size_t k = 0; k < dqs.dinucdq[i][j].size(); ++k){
					std::cerr << rgvals[i] << ", " << "q = " << j << ", dinuc = " << seq_nt16_str[seq_nt16_table['0' + (k >> 2)]] << seq_nt16_str[seq_nt16_table['0' + (k & 3)]] << ": " << dqs.dinucdq[i][j][k] << " (" <<
						dqs.meanq[i] + dqs.rgdq[i] + dqs.qscoredq[i][j] + dqs.dinucdq[i][j][k] << ") " << 
						data.dicov[i][j][k][1] << " " << data.dicov[i][j][k][0] << std::endl;
				}
			}
		}
	}

	std::cerr << put_now << " Recalibrating file" << std::endl;
	file = std::move(open_file(filename, is_bam, use_oq, set_oq));
	recalibrateutils::recalibrate_and_write(file.get(), dqs, "-");
	return 0;
}

// FILE:/home/adam/code/cbbq/src/kbbq/recalibrateutils.cc 


#include "recalibrateutils.hh"

using namespace htsiter;

namespace recalibrateutils{

//if the number of reads doesn't fit in a uint64_t call this twice :)
kmer_cache_t subsample_kmers(KmerSubsampler& s, uint64_t chunksize){
	uint64_t counted = 0;
	kmer_cache_t ret;
	ret.fill(std::vector<uint64_t>());
	for(bloom::Kmer kmer = s.next(); s.not_eof && ++counted < chunksize; kmer = s.next()){
		if(kmer.valid()){
			ret[kmer.prefix()].push_back(kmer.get());
		}
	}
	std::cerr << "Sampled kmers: " << counted << std::endl;
	std::cerr << "Total kmers in dataset: " << s.total_kmers << std::endl;
	return ret;
}

//this is a good target for multithreading ;)
void add_kmers_to_bloom(const kmer_cache_t& kmers, bloom::Bloom& filters){
	for(const std::vector<uint64_t>& v : kmers){
		for(uint64_t kmer : v){
			filters.insert(kmer);
		}
	}
}

kmer_cache_t find_trusted_kmers(HTSFile* file, const bloom::Bloom& sampled, std::vector<int> thresholds, int k, uint64_t chunksize){
	uint64_t counted = 0;
	kmer_cache_t ret;
	ret.fill(std::vector<uint64_t>());
	int n_trusted;
	bloom::Kmer kmer(k);
	//the order here matters since we don't want to advance the iterator if we're chunked out
	while(counted++ < chunksize && file->next() >= 0){
		readutils::CReadData read = file->get();
		read.infer_read_errors(sampled, thresholds, k);
		n_trusted = 0;
		kmer.reset();
		for(int i = 0; i < read.seq.length(); ++i){
			kmer.push_back(read.seq[i]);
			if(!read.errors[i]){
				++n_trusted;
			}
			if(i >= k && !read.errors[i-k]){
				--n_trusted;
			}
			if(kmer.valid() && n_trusted == k){
				ret[kmer.prefix()].push_back(kmer.get());
				//trusted kmer here
			}
		}
	}
	return ret;
}

covariateutils::CCovariateData get_covariatedata(HTSFile* file, const bloom::Bloom& trusted, int k){
	covariateutils::CCovariateData data;
#ifndef NDEBUG
	std::ifstream errorsin("../../adamjorr-Lighter/corrected.txt");
	int linenum = 0;
	std::string line = "";
#endif
	while(file->next() >= 0){
		readutils::CReadData read = file->get();
		read.get_errors(trusted, k, 6);
#ifndef NDEBUG
		//check that errors are same
		std::getline(errorsin, line);
		linenum++;
		std::vector<bool> lighter_errors(line.length(), false);
		std::transform(line.begin(), line.end(), lighter_errors.begin(),
			[](char c) -> bool {return (c == '1');});
		if( lighter_errors != read.errors){
			std::string message("Line num: " + std::to_string(linenum));
			// std::array<size_t,2> anchors = bloom::find_longest_trusted_seq(read.seq, trusted, k);
			// if(anchors[1] - anchors[0] - k + 1 >= k){ //number of trusted kmers >= k
			// 	anchors[1] = bloom::adjust_right_anchor(anchors[1], read.seq, trusted, k);
			// }
			// std::cerr << "Anchors: [" << anchors[0] << ", " << anchors[1] << "]";
			// std::cerr << " (npos is " << std::string::npos << ")\n";
			std::cerr << message << std::endl << "Errors : " ;
			for(const bool& v : read.errors){
				std::cerr << v;
			}
			std::cerr << std::endl;
			std::cerr << "Lighter: " ;
			for(const bool& v : lighter_errors){
				std::cerr << v;
			}
			std::cerr << std::endl;
			std::cerr << "Seq: " << read.seq << std::endl;
			// bloom::Kmer kmer(k);
			// for(const char& c : std::string("CAGAATAGAAAGATTTATAAATTAAATACTC")){
			// 	std::cerr << c << ":" << seq_nt4_table[c] << ":" << kmer.push_back(c) << "," ;
			// }
			// std::cerr << "Last kmer trusted?" << trusted[kmer.hashed_prefix()].query(kmer.get_query()) << std::endl; 
		}
		assert(lighter_errors == read.errors);
#endif
		data.consume_read(read);
	}
	return data;
}

void recalibrate_and_write(HTSFile* in, const covariateutils::dq_t& dqs, std::string outfn){
	if(in->open_out(outfn) < 0){
		//error!! TODO
		return;
	}
	while(in->next() >= 0){
		readutils::CReadData read = in->get();
		std::vector<uint8_t> newquals = read.recalibrate(dqs);
		in->recalibrate(newquals);
		if(in->write() < 0){
			//error! TODO
			return;
		}
	}
}

}

// FILE:/home/adam/code/cbbq/src/kbbq/bloom.cc 


#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <random>
#include "bloom.hh"

namespace bloom
{

	constexpr blocked_bloom_filter::cell_type blocked_bloom_filter::cell_type_zero;

	Bloom::Bloom(unsigned long long int projected_element_count, double fpr,
	unsigned long long int seed): params(){
		// params(projected_element_count, fpr, seed), bloom(params){
		params.projected_element_count = projected_element_count;
		params.false_positive_probability = fpr;
		params.random_seed = seed;
		if(!params){
			throw std::invalid_argument("Error: Invalid bloom filter parameters. \
				Adjust parameters and try again.");
		}
		params.compute_optimal_parameters();
		bloom = bloom_type(params);
	}

	Bloom::~Bloom(){}

	std::array<std::vector<size_t>,2> overlapping_kmers_in_bf(std::string seq, const Bloom& b, int k){
		bloom::Kmer kmer(k);
		std::vector<bool> kmer_present(seq.length()-k+1, false);
		std::vector<size_t> kmers_in(seq.length(), 0);
		std::vector<size_t> kmers_possible(seq.length(), 0);
		size_t incount = 0;
		size_t outcount = 0;
		for(size_t i = 0; i < seq.length(); ++i){
			kmer.push_back(seq[i]);
			if(i >= k-1){
				kmer_present[i-k+1] = kmer.valid() ? b.query(kmer) : false;
			}
#ifndef NDEBUG
			// if(i >= k-1){std::cerr << kmer << " " << kmer_present[i-k+1] << std::endl;}
			if(i >= k-1 && seq == "AAGTGGGTTTCTCAGTATTTTATTCTTTTGATATTATCATACATGATACTATCGTCTTGATTTCTTCTTCAGAGAGTTTATTGTTGTTGTAGAAATACAATTGATTTTTGTGTATTGATTTTGTATCCTGCAGCTTTGCTGAATTTTATTT"){
				std::string kmerstr(seq, i-k+1, k);
				std::cerr << kmerstr << " " << kmer_present[i-k+1] << std::endl;
			}
#endif
		}
		for(size_t i = 0; i < seq.length(); ++i){
			if(i < seq.length() - k + 1){ //add kmers now in our window
				if(kmer_present[i]){
					++incount;
				} else {
					++outcount;
				}
			}
			if(i >= k){ //remove kmers outside our window
				if(kmer_present[i - k]){
					--incount;
				} else {
					--outcount;
				}
			}
			kmers_in[i] = incount;
			kmers_possible[i] = incount + outcount;
		}
		return {kmers_in, kmers_possible};
	}

	int nkmers_in_bf(std::string seq, const Bloom& b, int k){
		Kmer kmer(k);
		int count = 0;
		for (size_t i = 0; i < seq.length(); ++i) {
			kmer.push_back(seq[i]);
			if (kmer.valid()) { //we have a full k-mer
				if(b.query(kmer)){
					count++;
				}
			}
		}
		return count;
	}

	char get_next_trusted_char(const Kmer& kmer, const Bloom& trusted, bool reverse_test_order){
		const std::array<char, 4> test_bases = !reverse_test_order ?
			std::array<char, 4>{'A','C','G','T'}: std::array<char, 4>{'T','G','C','A'};
		for(const char& c: test_bases){
			Kmer extra = kmer;
			extra.push_back(c);
			if(trusted.query(extra)){
				return c;
			}
		}
		return 0;
	}

	std::array<size_t, 2> find_longest_trusted_seq(std::string seq, const Bloom& b, int k){
		Kmer kmer(k);
		size_t anchor_start, anchor_end, anchor_best, anchor_current;
		anchor_start = anchor_end = std::string::npos;
		anchor_best = anchor_current = 0;
		for(size_t i = 0; i < seq.length(); ++i){
			if(kmer.push_back(seq[i]) >= k){
				if(b.query(kmer)){
					anchor_current++; //length of current stretch
				} else { //we had a streak but the kmer is not trusted
					if(anchor_current > anchor_best){
						anchor_best = anchor_current;
						anchor_end = i - 1; // always > 0 because i >= kmer.size() >= k
						anchor_start = i + 1 - k - anchor_current;
					}
					anchor_current = 0;
				}
			} else if(anchor_current != 0){ //we had a streak but ran into a non-ATGC base
				if(anchor_current > anchor_best){
					anchor_best = anchor_current;
					anchor_end = i - 1; // always > 0 because i >= kmer.size() >= k
					anchor_start = i + 1 - k - anchor_current;
				}
				anchor_current = 0;
			}
		} //we got to the end
		if(anchor_current > anchor_best){
			anchor_best = anchor_current;
			anchor_end = std::string::npos;
			anchor_start = seq.length() + 1 - k - anchor_current;
		}
		return std::array<size_t,2>{{anchor_start, anchor_end}};
	}

	std::tuple<std::vector<char>,size_t,bool> find_longest_fix(std::string seq, const Bloom& trusted, int k, bool reverse_test_order){
#ifndef NDEBUG
		std::cerr << seq << std::endl;
#endif
		Kmer kmer(k);
		std::vector<char> best_c{};
		size_t best_i = 0;
		bool single = false; //this pair of flags will be used to determine whether
		bool multiple = false; //multiple corrections were considered
		char unfixed_char = seq[k-1];
		const std::array<char,4> test_bases = !reverse_test_order ?
			std::array<char, 4>{'A','C','G','T'}: std::array<char, 4>{'T','G','C','A'};
		for(const char& c: test_bases){
			if(c == unfixed_char){continue;}
			seq[k-1] = c;
			kmer.reset();
			size_t i;
			size_t i_stop = std::max((size_t)2*k-1, seq.length()); //2k-1 -> 2k?
			for(i = 0; i < i_stop; ++i){ //i goes to max(2*k-1, seq.length())
				char n = i < seq.length() ? seq[i] : get_next_trusted_char(kmer, trusted, reverse_test_order);
				if(n == 0){ // no next trusted kmer
					break;
				}
				kmer.push_back(n);
				if(i >= k-1){
					// std::cerr << kmer.size() << std::endl;
					if(kmer.valid()){
#ifndef NDEBUG
						std::cerr << std::string(kmer) << " " << i << " " << trusted.query(kmer) << std::endl;
#endif
						if(!trusted.query(kmer)){
							break;
						} else { //we have a trusted kmer with this fix
							if(i == k-1){ //the first possible trusted kmer
								if(single){multiple = true;} //if we had one already, set multiple
								single = true;
							}
						}
					} else { //a non-ATCG base must've been added
						break;
					}
				}
			}
			if(i > best_i){
#ifndef NDEBUG
				std::cerr << std::string(kmer) << " L" << std::endl;
#endif
				best_c.clear();
				best_c.push_back(c);
				best_i = i;
			} else if (i == best_i){
#ifndef NDEBUG
				std::cerr << std::string(kmer) << " T" << std::endl;
#endif
				best_c.push_back(c);
			}
		}
		return std::make_tuple(best_c, best_i, multiple);
	}

	long double calculate_phit(const Bloom& bf, long double alpha){
		long double fpr = bf.fprate();
		double exponent = alpha < 0.1 ? 0.2 / alpha : 2;
		long double pa = 1 - pow(1-alpha,exponent);
		return pa + fpr - fpr * pa;
	}

	uint64_t numbits(uint64_t numinserts, long double fpr){
		//m = - n * log2(fpr) / ln2
		return numinserts * (-log2(fpr) / log(2));
	}

	int numhashes(long double fpr){
		//k = -log2(fpr)
		return ceil(-log2(fpr));
	}

	//ensure anchor >= k - 1 before this.
	std::pair<size_t,bool> adjust_right_anchor(size_t anchor, std::string seq, const Bloom& trusted, int k){
		Kmer kmer(k);
		bool multiple = false; //multiple corrections were considered
		size_t modified_idx = anchor + 1;
		assert((anchor >= k - 1));
		for(size_t i = modified_idx - k + 1; i < modified_idx; ++i){
			kmer.push_back(seq[i]);
		}
		for(const char& c : {'A','C','G','T'}){
			if(seq[modified_idx] == c){continue;}
			bloom::Kmer new_kmer = kmer;
			new_kmer.push_back(c);
			//if this fix works, we don't need to adjust the anchor if it fixes all remaining kmers.
			//modified_idx + 1 to modified_idx + 1 + k - 1
			for(size_t i = 0; i <= k; ++i){ //668188
#ifndef NDEBUG
				std::cerr << "i: " << i << " anchor + 2 + i: " << anchor + 2 + i << " len: " << seq.length() << std::endl;
#endif
				if(!trusted.query(new_kmer)){
					break;
				}
				//if we get to the end of the altered kmers and they're all fixed, the anchor is fine.
				if(modified_idx + i == seq.length()-1 || i == k){
#ifndef NDEBUG
					std::cerr << "First Try!" << std::endl;
#endif
					return std::make_pair(anchor, multiple);
				} else { //we haven't gotten to the end yet, so modified_idx+i+1 is valid.
					new_kmer.push_back(seq[modified_idx+i+1]);
				}
			}
		} // if we make it through this loop, we need to adjust the anchor.
		//we will test fixes starting with halfway through the last kmer to the end.
		//how much we're winding back the anchor; anchor-i-k must be > 0.
		for(int i = k/2-1; i >= 0 && anchor > i+k-1; --i){ 
			kmer.reset();
			modified_idx = anchor-i;
			// size_t modified_idx = anchor+1-i+k-2;
			// for(size_t j = anchor + 1 - i - k; j < anchor + 1 - i; ++j){ //start can be -1 from this
			for(size_t j = modified_idx - k + 1; j < modified_idx; ++j){
				kmer.push_back(seq[j]);
			}
			for(const char& c: {'A','C','G','T'}){
				if(seq[modified_idx] == c){continue;}
				bloom::Kmer new_kmer = kmer;
				new_kmer.push_back(c);
#ifndef NDEBUG
				std::cerr << "i: " << i << " modified_idx: " << modified_idx << " kmer:" << seq.substr(modified_idx-k+1,k-1) << c << std::endl;
#endif
				if(new_kmer.valid() && trusted.query(new_kmer)){
#ifndef NDEBUG
					std::cerr << "Trusted!" << std::endl;
#endif
					multiple = true;
					for(size_t j = 0; new_kmer.valid() &&
					trusted.query(new_kmer) &&
					modified_idx+1+j < seq.length() && j <= k/2; ++j){
						new_kmer.push_back(seq[modified_idx+1+j]);
						if(j == k/2 && new_kmer.valid() && trusted.query(new_kmer)){
							//we went the full length
							return std::make_pair(modified_idx-1, multiple); //the idx before the new base that needs fixing
							//only return here if j == k/2, NOT if you run out of sequence!!!
						}
					}
				}
			}
		}
		//couldn't find a better adjustment
		return std::make_pair(anchor, multiple);
	}

	int biggest_consecutive_trusted_block(std::string seq, const Bloom& trusted, int k, int current_len){
		Kmer kmer(k);
		int in = 0;
		int out = 0;
		int len = 0;
		for (size_t i = 0; i < seq.length(); ++i) {
			kmer.push_back(seq[i]);
			if (i >= k-1){ //we have a full k-mer
				if(trusted.query(kmer)){ //if it's in, increment
					++in;
				} else { //otherwise, reset but record the miss
					if(in > len){
						len = in;
					}
					in = 0;
					++out;
					if(k - out < current_len){ //end if we have too many misses
						break;
					}
				}
			}
		}
		if(in > len){
			len = in;
		}
		return len;
	}

//end namespace
}

// FILE:/home/adam/code/cbbq/src/kbbq/covariateutils.cc 


#include "covariateutils.hh"

namespace covariateutils{

	std::vector<long double> NormalPrior::normal_prior{};

	long double NormalPrior::get_normal_prior(size_t j){
		if(j >= normal_prior.size()){
			for(int i = normal_prior.size(); i < j+1; ++i){
				errno = 0;
				// long double prior_linspace = .9l * std::exp(-std::pow((long double)i,2.0l) * 2.0l)
				normal_prior.push_back(std::log(.9l * std::exp(-(std::pow(((long double)i/.5l),2.0l))/2.0l)));
				if(errno != 0){ //if an underflow happens just set the prior to smallest possible #
					normal_prior[i] = std::numeric_limits<long double>::lowest();
				}
			}
		}
		return normal_prior[j];
	}

	void CCovariate::increment(size_t idx, covariate_t value){
		this->increment(idx,value[0],value[1]);
	}

	void CCovariate::increment(size_t idx, unsigned long long err, unsigned long long total){
		this->at(idx)[0] += err;
		this->at(idx)[1] += total;
	}

	void CRGCovariate::consume_read(const readutils::CReadData& read){
		int rg = read.get_rg_int();
		if(this->size() <= rg){this->resize(rg+1);}
		std::vector<bool> nse = read.not_skipped_errors();
		this->increment(rg,
			std::accumulate(nse.begin(), nse.end(), 0),
			read.skips.size() - std::accumulate(read.skips.begin(), read.skips.end(), 0));
		// for(size_t i = 0; i < read.skips.size(); ++i){
		// 	if(!read.skips[i]){
		// 		this->increment(rg, !!read.errors[i], 1);
		// 	}
		// }
	}

	rgdq_t CRGCovariate::delta_q(meanq_t prior){
		rgdq_t dq(this->size());
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			int map_q = 0; //maximum a posteriori q
			long double best_posterior = std::numeric_limits<long double>::lowest();
			for(int possible = 0; possible < KBBQ_MAXQ+1; ++possible){
				int diff = std::abs(prior[i] - possible);
				long double prior_prob = NormalPrior::get_normal_prior(diff);
				long double p = recalibrateutils::q_to_p(possible);
				long double loglike = log_binom_pmf((*this)[i][0] + 1, (*this)[i][1] + 2, p);
				long double posterior = prior_prob + loglike;
				if(posterior > best_posterior){
					map_q = possible;
					best_posterior = posterior;
				}
			}
			dq[i] = map_q - prior[i];
		}
		return dq;
	}

	void CQCovariate::consume_read(const readutils::CReadData& read){
		int rg = read.get_rg_int();
		int q;
		if(this->size() <= rg){this->resize(rg+1);}
		for(size_t i = 0; i < read.skips.size(); ++i){
			if(!read.skips[i]){
				q = read.qual[i];
				if(this->at(rg).size() <= q){this->at(rg).resize(q+1);}
				(*this)[rg].increment(q, std::array<unsigned long long, 2>({read.errors[i], 1}));
			}
		}
	}

	qscoredq_t CQCovariate::delta_q(prior1_t prior){
		qscoredq_t dq(this->size());
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			dq[i].resize((*this)[i].size());
			for(int j = 0; j < (*this)[i].size(); ++j){ //j is q's here
				int map_q = 0; //maximum a posteriori q
				long double best_posterior = std::numeric_limits<long double>::lowest();
				for(int possible = 0; possible < KBBQ_MAXQ+1; possible++){
					int diff = std::abs(prior[i] - possible);
					long double prior_prob = NormalPrior::get_normal_prior(diff);
					long double p = recalibrateutils::q_to_p(possible);
					long double loglike = log_binom_pmf((*this)[i][j][0] + 1, (*this)[i][j][1] + 2, p);
					long double posterior = prior_prob + loglike;
					if(posterior > best_posterior){
						map_q = possible;
						best_posterior = posterior;
					}
				}
				dq[i][j] = map_q - prior[i];
			}
		}
		return dq;
	}

	void CCycleCovariate::consume_read(const readutils::CReadData& read){
		int rg = read.get_rg_int();
		int q;
		int cycle;
		if(this->size() <= rg){this->resize(rg+1);}
		size_t readlen = read.skips.size();
		for(size_t i = 0; i < readlen; ++i){
			if(!read.skips[i]){
				q = read.qual[i];
				if((*this)[rg].size() <= q){(*this)[rg].resize(q+1);}
				if((*this)[rg][q][read.second].size() <= i){(*this)[rg][q][read.second].resize(i+1);}
				(*this)[rg][q][read.second].increment(i, read.errors[i], 1);
			}
		}
	}

	cycledq_t CCycleCovariate::delta_q(prior2_t prior){
		cycledq_t dq(this->size()); //rg -> q -> fwd(0)/rev(1) -> cycle -> values
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			dq[i].resize((*this)[i].size());
			for(int j = 0; j < (*this)[i].size(); ++j){ //j is q's here
				for(int k = 0; k < 2; ++k){ //fwd/rev
					dq[i][j][k].resize((*this)[i][j][k].size());
					for(int l = 0; l < (*this)[i][j][k].size(); ++l){ //cycle value
						int map_q = 0; //maximum a posteriori q
						long double best_posterior = std::numeric_limits<long double>::lowest();
						for(int possible = 0; possible < KBBQ_MAXQ+1; possible++){
							int diff = std::abs(prior[i][j] - possible);
							long double prior_prob = NormalPrior::get_normal_prior(diff);
							long double p = recalibrateutils::q_to_p(possible);
							long double loglike = log_binom_pmf((*this)[i][j][k][l][0] + 1, (*this)[i][j][k][l][1] + 2, p);
							long double posterior = prior_prob + loglike;
							if(posterior > best_posterior){
								map_q = possible;
								best_posterior = posterior;
							}
						}
						dq[i][j][k][l] = map_q - prior[i][j];
					}
				}
			}
		}
		return dq;
	}

	void CDinucCovariate::consume_read(const readutils::CReadData& read, int minscore){
		int rg = read.get_rg_int();
		if(this->size() <= rg){this->resize(rg+1);}
		int q;
		for(size_t i = 1; i < read.seq.length(); ++i){
			q = read.qual[i];
			if(!read.skips[i] && nt_is_not_n(read.seq[i]) &&
				nt_is_not_n(read.seq[i-1]) && read.qual[i] >= minscore)
			{
				if((*this)[rg].size() <= q){(*this)[rg].resize(q+1);}
				if((*this)[rg][q].size() < 16){(*this)[rg][q].resize(16);}
				(*this)[rg][q].increment(dinuc_to_int(read.seq[i-1], read.seq[i]), read.errors[i], 1);
			}
		}
		// seq_nt16_table[256]: char -> 4 bit encoded (1/2/4/8)
		// seq_nt16_str[]: 4 bit -> char
		// seq_nt16_int[]: 4 bit -> 2 bits (0/1/2/3)
	}

	dinucdq_t CDinucCovariate::delta_q(prior2_t prior){
		dinucdq_t dq(this->size()); //rg -> q -> fwd(0)/rev(1) -> cycle -> values
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			dq[i].resize((*this)[i].size());
			for(int j = 0; j < (*this)[i].size(); ++j){ //j is q's here
				dq[i][j].resize((*this)[i][j].size());
				for(int k = 0; k < (*this)[i][j].size(); ++k){ //k is dinuc
					int map_q = 0; //maximum a posteriori q
					long double best_posterior = std::numeric_limits<long double>::lowest();
					for(int possible = 0; possible < KBBQ_MAXQ+1; possible++){
						int diff = std::abs(prior[i][j] - possible);
						long double prior_prob = NormalPrior::get_normal_prior(diff);
						long double p = recalibrateutils::q_to_p(possible);
						long double loglike = log_binom_pmf((*this)[i][j][k][0] + 1, (*this)[i][j][k][1] + 2, p);
						long double posterior = prior_prob + loglike;
						if(posterior > best_posterior){
							map_q = possible;
							best_posterior = posterior;
						}
					}
					dq[i][j][k] = map_q - prior[i][j];
				}
			}
		}
		return dq;
	}

	void CCovariateData::consume_read(readutils::CReadData& read, int minscore){
		//TODO: readd skips once the error correction code works properly
		// for(int i = 0; i < read.seq.length(); ++i){
		// 	read.skips[i] = (read.skips[i] || seq_nt16_int[seq_nt16_table[read.seq[i]]] >= 4 || read.qual[i] < minscore);
		// }
		rgcov.consume_read(read);
		qcov.consume_read(read);
		cycov.consume_read(read);
		dicov.consume_read(read, minscore);
	}

	dq_t CCovariateData::get_dqs(){
		dq_t dq;
		std::vector<long double> expected_errors(this->qcov.size(),0);
		meanq_t meanq(this->qcov.size(),0);
		for(int rg = 0; rg < this->qcov.size(); ++rg){
			for(int q = 0; q < this->qcov[rg].size(); ++q){
				expected_errors[rg] += (recalibrateutils::q_to_p(q) * this->qcov[rg][q][1]);
			}
			meanq[rg] = recalibrateutils::p_to_q(expected_errors[rg] / this->rgcov[rg][1]);
		}
		dq.meanq = meanq;
		dq.rgdq = this->rgcov.delta_q(meanq);
		prior1_t rgprior(this->qcov.size());
		for(int rg = 0; rg < this->qcov.size(); ++rg){
			rgprior[rg] = meanq[rg] + dq.rgdq[rg];
		}
		dq.qscoredq = this->qcov.delta_q(rgprior);
		prior2_t qprior(this->qcov.size());
		for(int rg = 0; rg < this->qcov.size(); ++rg){
			for(int q = 0; q < this->qcov[rg].size(); ++q){
				qprior[rg].push_back(rgprior[rg] + dq.qscoredq[rg][q]);
			}
		}
		dq.cycledq = this->cycov.delta_q(qprior);
		dq.dinucdq = this->dicov.delta_q(qprior);
		return dq;
	}



}
// FILE:/home/adam/code/cbbq/src/kbbq/readutils.cc 


#include "readutils.hh"
#include <algorithm>
#include <iterator>
#include <iostream>

KSEQ_DECLARE(BGZF*)

namespace readutils{

	std::unordered_map<std::string, std::string> CReadData::rg_to_pu{};
	std::unordered_map<std::string, int> CReadData::rg_to_int{};

	CReadData::CReadData(bam1_t* bamrecord, bool use_oq){ //TODO: flag to use OQ
		this->name = bam_get_qname(bamrecord);
		this->seq = bam_seq_str(bamrecord);
		if(use_oq){
			const uint8_t* oqdata = bam_aux_get(bamrecord, "OQ"); // this will be null on error
			//we should throw in that case
			if(oqdata == NULL){
				std::cerr << "Error: --use-oq was specified but unable to read OQ tag " << 
				"on read " << this->name << std::endl;
				if(errno == ENOENT){
					std::cerr << "OQ not found. Try again without the --use-oq option." << std::endl;
				} else if(errno == EINVAL){
					std::cerr << "Tag data is corrupt. Repair the tags and try again." << std::endl;
				}
				throw std::invalid_argument("Unable to read OQ tag.");
			}
			const std::string oq = bam_aux2Z(oqdata);
			std::transform(oq.cbegin(), oq.cend(), std::back_inserter(this->qual),
				[](const char& c) -> uint8_t {return c - 33;});
		} else {
			std::copy(bam_get_qual(bamrecord), bam_get_qual(bamrecord) + bamrecord->core.l_qseq,
				std::back_inserter(this->qual));
		}
		if(bam_is_rev(bamrecord)){
			//seq is already in fwd orientation
			std::reverse(this->qual.begin(), this->qual.end());
		}
		this->skips.resize(bamrecord->core.l_qseq, 0);
		uint8_t* rgdata = bam_aux_get(bamrecord, "RG");
		if(rgdata == NULL){
			std::cerr << "Error: Unable to read RG tag on read " << this->name << std::endl;
			if(errno == ENOENT){
				std::cerr << "RG not found. " <<
				"Every read in the BAM must have an RG tag; add tags with " <<
				"samtools addreplacerg and try again." << std::endl;
			} else if(errno == EINVAL){
				std::cerr << "Tag data is corrupt. Repair the tags and try again." << std::endl;
			}
			throw std::invalid_argument("Unable to read RG tag.");;
		}
		this->rg = bam_aux2Z(bam_aux_get(bamrecord, "RG"));
		if(rg_to_pu.count(this->rg) == 0){
			rg_to_int[this->rg] = rg_to_int.size();
			// we just need something unique here.
			rg_to_pu[this->rg] = rg; //when loaded from the header this is actually a PU 
		}
		this->second = bamrecord->core.flag & BAM_FREAD2; // 0x80
		this->errors.resize(bamrecord->core.l_qseq,0);
	}
	// if second is >1, that means infer.

	CReadData::CReadData(kseq::kseq_t* fastqrecord, std::string rg, int second, std::string namedelimiter):
	seq(fastqrecord->seq.s), skips(seq.length(), false), errors(seq.length(), false)
	{
		// this->seq = std::string(fastqrecord->seq.s);
		// this->qual.assign(fastqrecord->qual.s, fastqrecord->qual.l);
		std::string quals(fastqrecord->qual.s);
		std::transform(quals.begin(), quals.end(), std::back_inserter(this->qual),
			[](char c) -> int {return c - 33;});
		// this->skips.resize(this->seq.length(), false);

		std::string fullname(fastqrecord->name.s);
		size_t current_pos = fullname.find(namedelimiter);
		std::string first_name = fullname.substr(0, current_pos);

		while(rg == "" && current_pos != std::string::npos){
			// if we need to find rg
			fullname = fullname.substr(current_pos+1); //get the right part, excluding the delimiter
			current_pos = fullname.find(namedelimiter); //reset the delimiter; this is npos if no more fields
			if(fullname.substr(0,3) == "RG:"){
				size_t last_colon = fullname.find_last_of(":", current_pos); // current_pos is last char to search
				rg = fullname.substr(last_colon+1, current_pos);
			}
		}
		this->rg = rg;

		if(second > 1){
			std::string tail = first_name.substr(first_name.length() - 2);
			second = (tail == "/2");
			if(second || tail == "/1"){
				first_name = first_name.substr(0, first_name.length() - 2);
			}
		}
		this->name = first_name;
		this->second = second;
		// this->errors.resize(this->seq.length(), false);

		if(rg_to_pu.count(this->rg) == 0){
			rg_to_int[this->rg] = rg_to_int.size();
			rg_to_pu[this->rg] = rg; //when loaded from the header this is actually a PU.
		}
	}

	void CReadData::load_rgs_from_bamfile(bam_hdr_t* header){
		std::string hdrtxt(header->text);
		size_t linedelim = hdrtxt.find('\n');
		// kstring_t val;
		// while (sam_hdr_find_tag_pos(header, "RG", i, "ID", val) == 0){
		// 	std::string rgid(val.s);
		// 	sam_hdr_find_tag_pos(header, "RG", i, "PU", val);
		// 	std::string pu(val.s);
		// 	if(rg_to_pu.count(rgid) == 0){
		// 		rg_to_int[rgid] = rg_to_int.size();
		// 		rg_to_pu[rgid] = pu;
		// 	}
		// }
		while(linedelim != std::string::npos){
			if(hdrtxt.substr(0,3) == "@RG"){
				std::string id("");
				std::string pu("");
				std::string line = hdrtxt.substr(0, linedelim);
				size_t tokendelim = line.find_first_of("\t ");
				while(tokendelim != std::string::npos){
					if(line.substr(0,3) == "ID:"){
						//id
						id = line.substr(4, tokendelim);
					} else if (line.substr(0, 3) == "PU:"){
						//pu
						pu = line.substr(4, tokendelim);
					}
					line = line.substr(tokendelim);
					tokendelim = line.find_first_of("\t ");
				}
				if(id != ""){
					rg_to_int[id] = rg_to_int.size();
					rg_to_pu[id] = pu;
				}
			}
			hdrtxt = hdrtxt.substr(linedelim);
			linedelim = hdrtxt.find('\n');
		}
	}

	std::string CReadData::str_qual(){
		std::string str_qual;
		for(size_t i = 0; i < this->qual.size(); i++){
			str_qual.push_back(this->qual[i] + 33);
		}
		return str_qual;
	}

	std::string CReadData::canonical_name(){
		std::string suffix;
		if (this->second){
			suffix = "/2";
		}
		else{
			suffix = "/1";
		}
		return this->name + suffix;
	}

	std::vector<bool> CReadData::not_skipped_errors() const{
		std::vector<bool> unskipped_errs;
		std::transform(this->skips.cbegin(), this->skips.cend(), this->errors.cbegin(),
			std::back_inserter(unskipped_errs),
			[](const bool &s, const bool &e) -> bool {return (!s) && e;});
		return unskipped_errs;
	}

	void CReadData::infer_read_errors(const bloom::Bloom& b, const std::vector<int>& thresholds, int k){
		std::array<std::vector<size_t>,2> overlapping = bloom::overlapping_kmers_in_bf(this->seq, b, k);
		std::vector<size_t> in = overlapping[0];
		std::vector<size_t> possible = overlapping[1];
		for(size_t i = 0; i < errors.size(); ++i){
			this->errors[i] = (in[i] <= thresholds[possible[i]] || this->qual[i] <= INFER_ERROR_BAD_QUAL);
#ifndef NDEBUG
			if(possible[i] > k){std::cerr << "seq:" << this->seq << " " << possible[i] << "WARNING: Invalid i: " << i << std::endl;}
			// if(i>=k-1 && std::string(this->seq, i-k+1, k) == "CCCCCCCCCTCGCCCCCCCCCCCCCCCCCCC"){
			// 	std::cerr << "seq: " << this->seq << "\n";
			// 	std::cerr << "in: ";
			// 	std::copy(in.begin()+i-k+1, in.begin()+i+1, std::ostream_iterator<size_t>(std::cerr, ", "));
			// 	std::cerr << "\npossible: ";
			// 	std::copy(possible.begin()+i-k+1, possible.begin()+i+1, std::ostream_iterator<size_t>(std::cerr, ", "));
			// 	std::cerr << "\nerrors: ";
			// 	std::copy(this->errors.begin()+i-k+1, this->errors.begin()+i+1, std::ostream_iterator<bool>(std::cerr, ", "));
			// 	std::cerr << std::endl;
			// }			
#endif
		}
	}

	size_t CReadData::correct_one(const bloom::Bloom& t, int k){
		int best_fix_len = 0;
		char best_fix_base;
		size_t best_fix_pos = std::string::npos;
		for(size_t i = 0; i < this->seq.length(); ++i){
			std::string original_seq(this->seq); //copy original
			for(const char& c : {'A','C','G','T'}){
				if(this->seq[i] == c){continue;}
				original_seq[i] = c;
				size_t start = i > k - 1 ? i - k + 1 : 0;
				//test a kmer to see whether its worth counting them all
				//i'm not sure any performance gain is worth it, but this is how Lighter does it
				size_t magic_start = i > k/2 - 1 ? std::min(i - k/2 + 1, original_seq.length()-k) : 0;
				bloom::Kmer magic_kmer(k);
				for(size_t j = magic_start; j <= magic_start + k - 1; ++j){
					magic_kmer.push_back(original_seq[j]);
				}
				//
				if(t.query(magic_kmer)){
					//94518
					int n_in = bloom::biggest_consecutive_trusted_block(
						original_seq.substr(start, 2*k - 1),t,k,best_fix_len);
#ifndef NDEBUG					
					std::cerr << "Found a kmer: " << magic_kmer << " i: " << i << " Fix len: " << n_in << std::endl;
#endif
					if(n_in > best_fix_len){ //94518
						best_fix_base = c;
						best_fix_pos = i;
						best_fix_len = n_in;
					} else if(n_in == best_fix_len && this->qual[i] < this->qual[best_fix_pos]){
						best_fix_base = c;
						best_fix_pos = i;
					}
				}
			}
		}
		if(best_fix_len > 0){
			this->seq[best_fix_pos] = best_fix_base;
		}
		return best_fix_pos;
	}

	//this is a chonky boi
	std::vector<bool> CReadData::get_errors(const bloom::Bloom& trusted, int k, int minqual, bool first_call){
		std::string original_seq(this->seq);
		size_t bad_prefix = 0;
		size_t bad_suffix = std::string::npos;
		bool multiple = false; //whether there were any ties
#ifndef NDEBUG
		std::cerr << "Correcting seq: " << original_seq << std::endl;
#endif
		std::array<size_t,2> anchor = bloom::find_longest_trusted_seq(this->seq, trusted, k);
#ifndef NDEBUG
		std::cerr << "Initial anchors: [" << anchor[0] << ", " << anchor[1] << "]" << std::endl;
#endif
		if(anchor[0] == std::string::npos){ //no trusted kmers in this read.
			multiple = true;
			size_t corrected_idx = this->correct_one(trusted, k);
			if(corrected_idx == std::string::npos){
				return this->errors;
			} else {
				anchor = bloom::find_longest_trusted_seq(this->seq, trusted, k);
#ifndef	NDEBUG
				std::cerr << "Created anchor: [" << anchor[0] << ", " << anchor[1] << "]" << std::endl;
#endif				
				this->errors[corrected_idx] = true;
			}
		}
		if(anchor[0] == 0 && anchor[1] == std::string::npos){ //all kmers are trusted
			return this->errors;
		}
		//we're guaranteed to have a valid anchor now.
		//min is in case anchor[1] is npos.
		size_t anchor_len = std::min(anchor[1], this->seq.length()-1) + 1 - anchor[0];
		bool corrected = false; //whether there were any corrections
		//right side
		if(anchor[1] != std::string::npos){
			//number of trusted kmers >= k
			if(anchor_len - k + 1 >= k){ //number of trusted kmers 
				bool current_multiple;
				std::tie(anchor[1], current_multiple) = bloom::adjust_right_anchor(anchor[1], this->seq, trusted, k);
#ifndef NDEBUG
				std::cerr << "Adjust R Multiple: " << current_multiple << std::endl;
#endif
				multiple = multiple || current_multiple;
			}
			for(size_t i = anchor[1] + 1; i < this->seq.length();){
				size_t start = i - k + 1; //seq containing all kmers that are affected
				std::vector<char> fix;
				size_t fixlen;
				bool current_multiple;
				std::tie(fix, fixlen, current_multiple) = bloom::find_longest_fix(this->seq.substr(start, std::string::npos), trusted, k);
#ifndef NDEBUG
				std::cerr << "R fix Multiple: " << current_multiple << std::endl;
#endif
				multiple = multiple || current_multiple;
				size_t next_untrusted_idx = start + fixlen;
				if(next_untrusted_idx > i){
					if(fix.size() > 1){
						// multiple = true; //we take care of this in function
						// i = index of erroneous base
						//to = ( i + kmerLength - 1 < readLength ) ? i + kmerLength - 1 : readLength - 1 ;
						//to = inclusive biggest index the fix can possibly go to. (= largest_possible_idx)
						//maxto = largest index the fix actually goes to exclusive; =next_untrusted_idx
						//if( maxTo <= to || to - i + 1 < kmerLength ) ...
						//i + k
						multiple = true;
						size_t largest_possible_idx = std::min(i + k - 1, this->seq.length()-1);
#ifndef NDEBUG
						std::cerr << "next_untrusted_idx: " << next_untrusted_idx << std::endl;
						std::cerr << "largest_possible_idx (to): " << largest_possible_idx << std::endl;
						std::cerr << "largest_possible_idx-i+1: " << largest_possible_idx-i+1 << std::endl;
#endif					//20629746
						if(next_untrusted_idx <= largest_possible_idx || largest_possible_idx - i + 1 < k){
							// size_t trimstart = next_untrusted_idx;
							bad_suffix = i; //readlength - trimstart
#ifndef NDEBUG
							//if there's a tie and we haven't gone the max number of kmers, end correction
							std::cerr << "i: " << i << " next_untrusted_idx " << next_untrusted_idx << std::endl;
							std::cerr << "Fixlen is " << fixlen << " Fix: ";
							for(char c : fix){
								std::cerr << c << " ";
							}
							std::cerr << std::endl;
							std::cerr << "Tie and fix not long enough; ending correction early!" << std::endl;
#endif
							break;
						}
					} else {
						this->seq[i] = fix[0];
						this->errors[i] = true;
					}
					corrected = true;
#ifndef NDEBUG
					std::cerr << "Error detected at position " << i << ". Advancing " << fixlen - k + 1 << "." << std::endl;
#endif
					i += fixlen - k + 1; // i = next_untrusted_idx
				} else {
					//couldn't find a fix; skip ahead and try again if long enough
#ifndef NDEBUG
					std::cerr << "Couldn't fix position " << i << " Cutting read and trying again." << std::endl;//". Skipping ahead " << k - 1 << "." << std::endl;
#endif
					bad_suffix = i;
					break;
					// i += k-1; //move ahead and make i = k-1 as the first base in the new kmer
					// if(this->seq.length() - i + k <= (seq.length()/2) || this->seq.length() - i + k <= 2*k ){
					// 	//sequence not long enough. end this side.
					// 	break;
					// }
				}
			}
		}
		//left side
		if(anchor[0] != 0){
			//the bad base is at anchor[0]-1, then include the full kmer for that base.
			// std::string sub = this->seq.substr(0, anchor[0] - 1 + k);
			std::string revcomped(this->seq.length(), 'N');
			std::transform(this->seq.rbegin(), this->seq.rend(), revcomped.begin(),
				[](char c) -> char {return seq_nt16_str[seq_nt16_table[('0' + 3-seq_nt16_int[seq_nt16_table[c]])]];});
			//if num of trusted kmers >= k, see if anchor needs adjusting.
			if(anchor_len - k + 1 >= k){ 
				size_t left_adjust;
				bool current_multiple;
				std::tie(left_adjust, current_multiple) = bloom::adjust_right_anchor(revcomped.length()-anchor[0]-1, revcomped, trusted, k);
				anchor[0] = revcomped.length()-left_adjust-1; //change back to original coordinates
#ifndef NDEBUG
				std::cerr << "Adjust L Multiple: " << current_multiple << std::endl;
#endif
				multiple = multiple || current_multiple;
			}
			//
			for(int i = anchor[0] - 1; i >= 0;){ //index of erroneous base in original seq
				int j = revcomped.length()-i-1; //index of erroneous base in reversed seq
				size_t start = j - k + 1; //seq containing all kmers that are affected
				//but [j -k + 1, npos) in reverse space.
				std::string sub = revcomped.substr(start, std::string::npos); //get the right subsequence
				std::vector<char> fix;
				size_t fixlen;
				bool current_multiple;
				std::tie(fix, fixlen, current_multiple) = bloom::find_longest_fix(sub, trusted, k, true); //155392 TODO: add reverse_test to adjust_right_anchor
#ifndef NDEBUG
				std::cerr << "L Fix Multiple: " << current_multiple << std::endl;
#endif
				multiple = multiple || current_multiple;
				//new i is return value + start // i += r -k + 1
				//new j should be return value + start // j += r - k + 1
				// new j should be fixlen + start; j += fixlen - k + 1
				// j = fixlen + start; j = j - k + 1 + fixlen; j += fixlen -k + 1
				size_t next_untrusted_idx = start + fixlen; // next untrusted idx in j space
				if( next_untrusted_idx > j){
					if(fix.size() > 1){
						multiple = true; // we don't have to do this here because we do it in function
						//this value needs to be fixed i think
						//to = ( i - kmerLength + 1 < 0 ) ? 0 : ( i - kmerLength + 1 ) ; 
						//( minTo >= to || i - to + 1 < kmerLength )
						//Line num: 21537
						// if(next_untrusted_idx <= std::max(start + (size_t)2*k - 1, revcomped.length())){
						size_t largest_possible_idx = std::min(j + (size_t)k - 1, revcomped.length()-1);
						if(next_untrusted_idx <= largest_possible_idx || largest_possible_idx - j + 1 < k){ //595573
							//if there's a tie and we haven't gone the max number of kmers, end correction
							bad_prefix = i;
							break;
						}
					} else {
						revcomped[j] = fix[0];
						this->errors[i] = true;
					}
					corrected = true;
#ifndef NDEBUG
					std::cerr << "next_untrusted_idx: " << next_untrusted_idx << " j: " << j << std::endl;
					std::cerr << "Error detected at position " << i << ". Advancing " << next_untrusted_idx-j << " " << (fixlen - k + 1) << "." << std::endl;
#endif					
					i -= next_untrusted_idx - j; //fixlen - k + 1;
				} else {
					//couldn't find a fix; skip ahead and try again if long enough
#ifndef NDEBUG
					std::cerr << "Couldn't fix position " << i << "Cutting read and trying again." << std::endl;//". Skipping ahead " << k - 1 << "." << std::endl;
#endif			
					bad_prefix = i; //the last position in the new read
					break;
					// i -= k-1;
					// if(i + k <= (seq.length()/2) || i + k <= 2*k ){
					// 	//sequence not long enough. end this side.
					// 	break;
					// }
				}
			}
		}
#ifndef NDEBUG
		std::cerr << "Anchors: [" << anchor[0] << ", " << anchor[1] << "]" << std::endl; 
#endif


		// check for overcorrection and fix it
		if(corrected){
			bool adjust = true;
			//check that no trusted kmers were "fixed"
			bloom::Kmer kmer(k);
			size_t trusted_start = std::string::npos;
			size_t trusted_end = std::string::npos;
			for(size_t i = 0; i < original_seq.length() && adjust == true; ++i){
				kmer.push_back(original_seq[i]);
				if(kmer.valid() && trusted.query(kmer)){
					trusted_start = std::min(trusted_start,i-k+1);
					trusted_end = i;
					// std::cerr << "Trusted: " << trusted_start << " " << trusted_end << std::endl;
				} else {
					if(i > trusted_end){
						//we clear everything from beginning to end
						for(size_t j = trusted_start; j <= trusted_end; ++j){
							if(this->errors[j]){
								adjust = false;
								break;
							}
						}
						trusted_start = std::string::npos;
						trusted_end = std::string::npos;
					}
				}
			}
			//if we get to the end but have a trusted block
			//we don't actually trust that block :|
			// if(trusted_end != std::string::npos){
			// 	for(size_t j = trusted_start; j <= trusted_end; ++j){
			// 		if(this->errors[j]){
			// 			std::cerr << "(2) Problem: j: " << j << " " << trusted_start << " " << trusted_end << std::endl;
			// 			adjust = false;
			// 			break;
			// 		}
			// 	}
			// }
#ifndef NDEBUG
			std::cerr << "Adjust: " << adjust << " Multiple: " << multiple << std::endl;
#endif
			adjust = adjust && !multiple;//made it through the loop and no ties during correction
			// std::cerr << "Read corrected. Adjust threshold? " << adjust << std::endl;
#ifndef NDEBUG
			std::cerr << "Errors before adjustment: ";
			for(const bool& b: this->errors){
				std::cerr << b;
			}
			std::cerr << std::endl;
			std::cerr << "Seq After Correction: " << seq << std::endl;
#endif
			int ocwindow = 20;
			int base_threshold = 4;
			int threshold = base_threshold;
			double occount = 0;
			//check for overcorrection
			std::vector<int> overcorrected_idx; //push_back overcorrected indices in order
			//then from overcorrected.begin() - k to overcorrected.end() + k should all be reset.
			for(int i = 0; i < this->seq.length(); ++i){
				if(this->errors[i] && seq_nt16_int[seq_nt16_table[original_seq[i]]] < 4){ //increment correction count
					if(this->qual[i] <= minqual){
						occount += 0.5;
					} else {
						++occount;
					}
				}
				if(i >= ocwindow && this->errors[i-ocwindow] && seq_nt16_int[seq_nt16_table[original_seq[i-ocwindow]]] < 4){ //decrement count for not in window
					if(this->qual[i-ocwindow] <= minqual){
						occount -= 0.5;
					} else {
						--occount;
					}
				}
				//set threshold
				threshold = adjust && i >= ocwindow && i + ocwindow - 1 < this->seq.length() ? 
					base_threshold + 1 : base_threshold;
				//determine if overcorrected
#ifndef NDEBUG
				std::cerr << "Occount: " << occount << " Threshold: " << threshold;
				std::cerr << " Seq: " << original_seq[i] << " (" << seq_nt16_int[seq_nt16_table[original_seq[i]]];
				std::cerr << ")" << " Q: " << +this->qual[i] << std::endl;
#endif
				if(occount > threshold && this->errors[i]){
					overcorrected_idx.push_back(i);
				}

			}
#ifndef NDEBUG
			std::cerr << "Overcorrected indices (" << overcorrected_idx.size() << "): ";
			std::copy(overcorrected_idx.begin(), overcorrected_idx.end(), std::ostream_iterator<int>(std::cerr, ", "));
			std::cerr << std::endl;
#endif
			//Line num: 23026
			for(int oc_idx : overcorrected_idx){
				if(this->errors[oc_idx]){ //overcorrected idx hasn't been addressed yet
					int start = oc_idx-k+1; //the beginningmost position to check
					start = start >= 0? start : 0;
					int end = oc_idx+k; //the endmost position to check
					end = end < this->seq.length() ? end : this->seq.length();
					//we start iteration but we need to unfix anything within k of an overcorrected position
					//OR within k of one of those fixed positions.
					for(int i = start; i < end; ++i){
						if(this->errors[i]){
							this->errors[i] = false;
							//changint the end must come before the start change because the start
							//change changes i!
							if(i+k > end){ //change the end if we need to
								end = i+k < this->seq.length() ? i+k : this->seq.length();
							}
							//i will be 1 greater than start, so rather than i-k+1 we have i-k.
							if(i-k < start){ //go back a bit if we need to; +1 comes from the loop
								i = i-k+1 >= 0 ? i-k : -1; //+1 will come from the loop
								start = i;
							}
						}
					}
				}
			}
		}
		if(first_call && bad_prefix > 0 && (bad_prefix >= this->seq.length() / 2 || bad_prefix >= 2*k)){
#ifndef NDEBUG			
			std::cerr << "bad_prefix: " << bad_prefix << std::endl;
#endif
			CReadData subread = this->substr(0, bad_prefix+1); //2nd argument is length
			std::vector<bool> suberrors = subread.get_errors(trusted, k, minqual, false);
			std::copy(suberrors.begin(), suberrors.end(), this->errors.begin());
		}
		if(first_call && bad_suffix < std::string::npos && bad_suffix < this->seq.length() &&
		(this->seq.length()-bad_suffix > this->seq.length()/2 || this->seq.length()-bad_suffix > 2*k)){
#ifndef NDEBUG
			std::cerr << "bad_suffix: " << bad_suffix << std::endl;
#endif
			CReadData subread = this->substr(bad_suffix, std::string::npos);
			std::vector<bool> suberrors = subread.get_errors(trusted, k, minqual, false);
			std::copy(suberrors.begin(), suberrors.end(), this->errors.begin()+bad_suffix);
		}
		if(std::find(this->errors.begin(), this->errors.end(), true) != this->errors.end()){
			corrected = true;
		}

		this->seq = original_seq;
		return this->errors;
	}

	std::vector<uint8_t> CReadData::recalibrate(const covariateutils::dq_t& dqs, int minqual) const{
		std::vector<int> recalibrated;
		std::copy(this->qual.begin(), this->qual.end(), std::back_inserter(recalibrated));
		int rg = this->get_rg_int();
		for(int i = 0; i < this->seq.length(); ++i){
			uint8_t q = this->qual[i];
			if(q >= minqual){
				recalibrated[i] = dqs.meanq[rg] + dqs.rgdq[rg] + dqs.qscoredq[rg][q] +
					dqs.cycledq[rg][q][this->second][i];
				if(i > 0){
					int first = seq_nt16_int[seq_nt16_table[this->seq[i-1]]];
					int second = seq_nt16_int[seq_nt16_table[this->seq[i]]];
					if(first < 4 && second < 4){
						int8_t dinuc = 15 & ((first << 2) | second); //1111 & (xx00|00xx)
						recalibrated[i] += dqs.dinucdq[rg][q][dinuc];
					}
				}
			}
		}
		std::vector<uint8_t> ret;
		std::transform(recalibrated.begin(), recalibrated.end(), std::back_inserter(ret),
			[](int q)->uint8_t {return q < 0 ? 0 : KBBQ_MAXQ < q ? KBBQ_MAXQ : q;}); //std::clamp in c++17
		return ret;
	}

	CReadData CReadData::substr(size_t pos, size_t count) const{
		CReadData ret = (*this);
		ret.seq = ret.seq.substr(pos, count);
		size_t len = ret.seq.size();
		ret.qual.resize(len);
		ret.skips.resize(len);
		ret.errors.resize(len);
		std::copy(this->qual.cbegin()+pos, this->qual.cbegin()+pos+len, ret.qual.begin());
		std::copy(this->skips.cbegin()+pos, this->skips.cbegin()+pos+len, ret.skips.begin());
		std::copy(this->errors.cbegin()+pos, this->errors.cbegin()+pos+len, ret.errors.begin());
		return ret;
	}


}

