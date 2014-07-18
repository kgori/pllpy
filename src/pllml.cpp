/*
 * pllml.cpp
 *
 *  Created on: Jul 17, 2014
 *      Author: kgori
 */

#include <stdexcept>
#include <fstream>
#include "pllml.h"
extern "C" {
#include <pll/pll.h>
}

pll::pll() {
    attr.rateHetModel = PLL_GAMMA;
    attr.fastScaling = PLL_TRUE;
    attr.saveMemory = PLL_TRUE;
    attr.useRecom = PLL_FALSE;
    attr.numberOfThreads = 1;
}

pll::pll(string alignment_file, string tree, string partition_file, int num_threads) {
	attr.rateHetModel = PLL_GAMMA;
	attr.fastScaling = PLL_TRUE;
	attr.saveMemory = PLL_TRUE;
	attr.useRecom = PLL_FALSE;
	attr.numberOfThreads = num_threads;
	create_instance();
	load_alignment_file(alignment_file);
	if (_is_file(tree)) {
		load_tree_file(tree);
	}
	else if (_is_tree_string(tree)) {
		load_tree_string(tree);
	}
	else {
		cerr << "Didn't understand tree: " << tree << endl;
		throw exception();
	}
	load_partition_file(partition_file);
	initialise_model();
}

pll::~pll() {
    destroy_model();
    alignment_data = nullptr;
    newick = nullptr;
    partitions = nullptr;
    tr = nullptr;
}

void pll::load_alignment_file(string path) {
    alignment_data = pllParseAlignmentFile(PLL_FORMAT_PHYLIP, path.c_str());
    if (!alignment_data) {
        throw exception();
    }
}

void pll::create_instance() {
    tr = pllCreateInstance(&attr);
}

void pll::load_tree_file(string path) {
    newick = pllNewickParseFile(path.c_str());
    if (!newick) {
        cerr << "tree parse error" << endl;
        throw exception();
    }
    if (!pllValidateNewick(newick)) /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */ {
        cerr << "invalid tree" << endl;
        throw exception();
    }
}

void pll::load_partition_file(string path) {
    pllQueue * partitionInfo;
    partitionInfo = pllPartitionParse(path.c_str());
    if (!pllPartitionsValidate(partitionInfo, alignment_data)) {
        cerr << "partitions parse error" << endl;
        throw exception();
    }
    partitions = pllPartitionsCommit(partitionInfo, alignment_data);
    pllAlignmentRemoveDups(alignment_data, partitions);
    pllQueuePartitionsDestroy(&partitionInfo);
}

void pll::load_tree_string(string nwk) {
    newick = pllNewickParseString(nwk.c_str());
    if (!newick) {
        throw exception();
    }
    if (!pllValidateNewick(newick)) /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */ {
        throw exception();
    }
}

void pll::initialise_model() {
    pllTreeInitTopologyNewick(tr, newick, PLL_FALSE);
    if (!pllLoadAlignment(tr, alignment_data, partitions, PLL_DEEP_COPY)) {
        cerr << "Model finalisation error" << endl;
        throw exception();
    }
    pllInitModel(tr, partitions, alignment_data);
    _is_ready = true;
}

void pll::destroy_model() {
    if (alignment_data) pllAlignmentDataDestroy(alignment_data);
    if (newick) pllNewickParseDestroy(&newick);
    if (partitions) pllPartitionsDestroy(tr, &partitions);
    if (tr) pllDestroyInstance(tr);
    _is_ready = false;
}

void pll::optimise(bool estimate_model) {
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    int pll_bool = estimate_model ? PLL_TRUE : PLL_FALSE;
    pllRaxmlSearchAlgorithm(tr, partitions, pll_bool);
}

void pll::optimise_branch_lengths(int num_iter) {
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    pllOptimizeBranchLengths(tr, partitions, num_iter);
}

void pll::optimise_model() {
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    pllOptimizeModelParameters(tr, partitions, tr->likelihoodEpsilon);
}

double pll::get_likelihood() {
    if (!_is_ready) return 0;
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    return tr->likelihood;
}

int pll::get_number_of_threads() {
    return attr.numberOfThreads;
}

int pll::get_number_of_partitions() {
    return partitions->numberOfPartitions;
}

vector<string> pll::get_partition_names() {
    vector<string> names;
    for (size_t i = 0; i < get_number_of_partitions(); ++i) {
        names.push_back(get_partition_name(i));
    }
    return names;
}

string pll::get_partition_name(int partition) {
    _partitions_bounds_check(partition);
    return string(partitions->partitionData[partition]->partitionName);
}

vector<double> pll::get_alphas() {
    vector<double> alphas;
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        alphas.push_back(get_alpha(i));
    }
    return alphas;
}

double pll::get_alpha(int partition) {
    _partitions_bounds_check(partition);
    return pllGetAlpha(partitions, partition);
}

vector<vector<double>> pll::get_frequencies() {
    vector<vector<double>> freqs;
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        freqs.push_back(get_frequencies_vector(i));
    }
    return freqs;
}

vector<double> pll::get_frequencies_vector(int partition) {
    _partitions_bounds_check(partition);
    cout << "optimise freqs: " << partitions->partitionData[partition]->optimizeBaseFrequencies << endl;
    vector<double> freqs_vec;
    if (!_is_ready) {
        cerr << "Model isn't finalised" << endl;
        return freqs_vec;
    }
    int num_states = partitions->partitionData[partition]->states;
    for (int j = 0; j < num_states; ++j) {
        freqs_vec.push_back(partitions->partitionData[partition]->frequencies[j]);
    }
    return freqs_vec;
}

vector<vector<double>> pll::get_rates() {
    vector<vector<double>> rates;
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        rates.push_back(get_rates_vector(i));
    }
    return rates;
}

vector<double> pll::get_rates_vector(int partition) {
    _partitions_bounds_check(partition);
    vector<double> rates_vec;
    if (!_is_ready) {
        cerr << "Model isn't finalised" << endl;
        return rates_vec;
    }
    int num_states = partitions->partitionData[partition]->states;
    int num_rates = (num_states * (num_states - 1)) / 2;
    for (int j = 0; j < num_rates; ++j) {
        rates_vec.push_back(partitions->partitionData[partition]->substRates[j]);
    }
    return rates_vec;
}

double pll::get_epsilon() {
    return tr->likelihoodEpsilon;
}

string pll::get_tree() {
    if (!_is_ready) {
        cerr << "Model isn't finalised" << endl;
        return "";
    }
    pllTreeToNewick(tr->tree_string, tr, partitions,
                    tr->start->back, PLL_TRUE, PLL_TRUE,
                    PLL_TRUE, PLL_FALSE, PLL_FALSE,
                    PLL_SUMMARIZE_LH, 0,0);
    return tr->tree_string;
}

vector<vector<double> > pll::get_empirical_frequencies() {
    if (!_is_ready) {
        cerr << "Model isn't finalised" << endl;
        throw exception();
    }
    double ** ef;
    vector<vector<double>> vec_2d;
    ef = pllBaseFrequenciesGTR(partitions, alignment_data);
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        vector<double> row_vec;
        size_t ns = partitions->partitionData[i]->states;
        for (size_t j = 0; j < ns; ++j) {
            row_vec.push_back(ef[i][j]);
        }
        vec_2d.push_back(row_vec);
    }
    return vec_2d;
}

void pll::set_epsilon(double epsilon) {
    tr->likelihoodEpsilon = epsilon;
}

void pll::set_rates(vector<double> rates,
        int partition, bool optimisable) {
    if (is_dna(partition)) {
        _partitions_bounds_check(partition);
        int num_states = partitions->partitionData[partition]->states;
        size_t num_rates = (num_states * (num_states - 1)) / 2;
        if (rates.size() != num_rates) {
            cerr << "Rates vector is the wrong length. Should be " << num_rates << endl;
            throw exception();
        }
        pllSetFixedSubstitutionMatrix(&(rates[0]), num_rates, partition, partitions, tr);
        set_optimisable_rates(partition, optimisable);
    }
}

void pll::set_alpha(double alpha, int partition, bool optimisable) {
    _partitions_bounds_check(partition);
    pllSetFixedAlpha(alpha, partition, partitions, tr);
    set_optimisable_alpha(partition, optimisable);
}

void pll::set_frequencies(vector<double> freqs, int partition, bool optimisable) {
    if (!_approx_eq(_vector_sum(freqs), 1)) {
        cerr << "Not setting frequencies: Frequencies do not sum to 1" << endl;
        throw exception();
    }
    _partitions_bounds_check(partition);
    size_t num_states = partitions->partitionData[partition]->states;
    if (freqs.size() != num_states) {
        cerr << "Frequencies vector is the wrong length. Should be " << num_states << endl;
        throw exception();
    }
    pllSetFixedBaseFrequencies(&(freqs[0]), num_states, partition, partitions, tr);
    set_optimisable_frequencies(partition, optimisable);
}

void pll::link_alpha_parameters(string linkage) {
    pllLinkAlphaParameters(const_cast<char*>(linkage.c_str()), partitions);
}

void pll::link_frequencies(string linkage) {
    pllLinkFrequencies(const_cast<char*>(linkage.c_str()), partitions);
}

void pll::link_rates(string linkage) {
    pllLinkRates(const_cast<char*>(linkage.c_str()), partitions);
}

bool pll::is_dna(int partition) {
	_partitions_bounds_check(partition);
	return (partitions->partitionData[partition]->dataType == PLL_DNA_DATA);
}

bool pll::is_protein(int partition) {
	_partitions_bounds_check(partition);
		return (partitions->partitionData[partition]->dataType == PLL_AA_DATA);
}

bool pll::_is_file(string filename) {
    ifstream fl(filename.c_str());
    bool result = true;
    if (!fl) {
        result = false;
    }
    fl.close();
    return result;
}

void pll::set_number_of_threads(int threads) {
    attr.numberOfThreads = threads;
}

bool pll::is_optimisable_alpha(int partition) {
    _partitions_bounds_check(partition);
    return (partitions->partitionData[partition]->optimizeAlphaParameter == PLL_TRUE);
}

bool pll::is_optimisable_frequencies(int partition) {
    _partitions_bounds_check(partition);
    return (partitions->partitionData[partition]->optimizeBaseFrequencies == PLL_TRUE);
}

bool pll::is_optimisable_rates(int partition) {
    _partitions_bounds_check(partition);
    return (partitions->partitionData[partition]->optimizeSubstitutionRates == PLL_TRUE);
}

vector<string> pll::get_model_names() {
	vector<string> names;
	size_t np = get_number_of_partitions();
	for (size_t i = 0; i < np; ++i) {
		names.push_back(get_model_name(i));
	}
	return names;
}

string pll::get_model_name(int partition) {
	_partitions_bounds_check(partition);
	string name;
	if (is_dna(partition)) {
		name = "GTR";
	}
	else if (partitions->partitionData[partition]->protModels == PLL_AUTO) {
		name = _model_name(partitions->partitionData[partition]->autoProtModels);
	} else {
		name = _model_name(partitions->partitionData[partition]->protModels);
	}
	return name;
}

bool pll::_is_tree_string(string tree_string) {
    size_t l = tree_string.length();
    return (tree_string[0]=='(' && tree_string[l-1]==';');
}

void pll::set_optimisable_rates(int partition, bool optimisable) {
    if (is_protein(partition)) {
        cerr << "Optimising rates not implemented for protein models" << endl;
        throw exception();
    }
    _partitions_bounds_check(partition);
    int pll_bool = optimisable ? PLL_TRUE : PLL_FALSE;
    if (is_dna(partition)) {
        partitions->partitionData[partition]->optimizeSubstitutionRates = pll_bool;
    //  partitions->dirty = PLL_TRUE;
        _evaluate_likelihood();
    }
}

void pll::set_optimisable_alpha(int partition, bool optimisable) {
    _partitions_bounds_check(partition);
    int pll_bool = optimisable ? PLL_TRUE : PLL_FALSE;
    partitions->partitionData[partition]->optimizeAlphaParameter = pll_bool;
//  partitions->dirty = PLL_TRUE;
    _evaluate_likelihood();
}

void pll::set_optimisable_frequencies(int partition, bool optimisable) {
    _partitions_bounds_check(partition);
    int pll_bool = optimisable ? PLL_TRUE : PLL_FALSE;
    partitions->partitionData[partition]->optimizeBaseFrequencies = pll_bool;
//  partitions->dirty = PLL_TRUE;
    _evaluate_likelihood();
}

void pll::_evaluate_likelihood() {
    if (!_is_ready) {
        cerr << "Model isn't finalised" << endl;
        return;
    }
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
}

double pll::_vector_sum(vector<double> vec) {
    double s = 0;
    for (size_t i = 0; i < vec.size(); ++i) {
        s += vec[i];
    }
    return s;
}

bool pll::_approx_eq(double a, double b) {
    double d = a - b;
    if (d < 0) d = -d;
    return (d < 0.000001);
}

void pll::_partitions_bounds_check(int partition) {
    int max_partitions = get_number_of_partitions();
    if (partition >= max_partitions) {
        cerr << "The model has " << max_partitions << " partitions" << endl;
        throw exception();
    }
}

string pll::_model_name(int model_num) {
    string name;
    switch (model_num) {

        case PLL_DAYHOFF :
        name = "DAYHOFF";
        break;

        case PLL_DCMUT :
        name = "DCMUT";
        break;

        case PLL_JTT :
        name = "JTT";
        break;

        case PLL_MTREV :
        name = "MTREV";
        break;

        case PLL_WAG :
        name = "WAG";
        break;

        case PLL_RTREV :
        name = "RTREV";
        break;

        case PLL_CPREV :
        name = "CPREV";
        break;

        case PLL_VT :
        name = "VT";
        break;

        case PLL_BLOSUM62 :
        name = "BLOSUM62";
        break;

        case PLL_MTMAM :
        name = "MTMAM";
        break;

        case PLL_LG :
        name = "LG";
        break;

        case PLL_MTART :
        name = "MTART";
        break;

        case PLL_MTZOA :
        name = "MTZOA";
        break;

        case PLL_PMB :
        name = "PMB";
        break;

        case PLL_HIVB :
        name = "HIVB";
        break;

        case PLL_HIVW :
        name = "HIVW";
        break;

        case PLL_JTTDCMUT :
        name = "JTTDCMUT";
        break;

        case PLL_FLU :
        name = "FLU";
        break;

        case PLL_AUTO :
        name = "AUTO";
        break;

        case PLL_LG4 :
        name = "LG4";
        break;

        case PLL_GTR :
        name = "GTR";
        break;

        default:
        name = "UNKNOWN";
        break;
    }
    return name;
}
