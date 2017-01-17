/*
pllml.cpp
Copyright (C) 2014  Kevin Gori

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "pllml.h"
extern "C" {
#include <pll/pll.h>
}

// PUBLIC METHODS

pll::pll(std::string alignment_file, std::string partitions, std::string tree, int num_threads, long rns) {
    _init_attr(num_threads, rns);
    _init_instance();
    _init_alignment_file(alignment_file);
    if (_is_file(partitions)) {
        _init_partition_file(partitions);
    }
    else {
        _init_partition_string(partitions);
    }
    if (_is_file(tree)) {
        _init_tree_file(tree);
    }
    else if (_is_tree_string(tree)) {
        _init_tree_string(tree);
    }
    else {
        std::cerr << "Didn't understand tree: " << tree << std::endl;
        throw std::invalid_argument("Didn't understand tree: " + tree);
    }
    _init_model(false);
}

pll::pll(std::string alignment_file, std::string partitions, bool parsimony, int num_threads, long rns) {
    _init_attr(num_threads, rns);
    _init_instance();
    _init_alignment_file(alignment_file);
    if (_is_file(partitions)) {
        _init_partition_file(partitions);
    }
    else {
        _init_partition_string(partitions);
    }
    _init_tree_random();
    _init_model(parsimony);
}

pll::~pll() {
    _destroy_model();
}

void pll::optimise_tree_search(bool estimate_model) {
    _check_model_ready();
    tr->start = tr->nodep[1];
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    int pll_bool = estimate_model ? PLL_TRUE : PLL_FALSE;
    std::cout << "Performing tree search with RAxML search algorithm, "
              << (estimate_model ? "with " : "without ")
              << "model parameter optimisation..." << std::endl;
    pllRaxmlSearchAlgorithm(tr, partitions, pll_bool);
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    std::cout << "Likelihood = " << tr->likelihood << "\nEND." << std::endl;
}

void pll::optimise_model() {
    optimise(true, true, true, false, false, 0, false);
}

void pll::optimise_rates() {
    optimise(false, true, false, false, false, 0, false);
}

void pll::optimise_freqs() {
    optimise(true, false, false, false, false, 0, false);
}

void pll::optimise_alphas() {
    optimise(false, false, true, false, false, 0, false);
}

void pll::optimise_branch_lengths(int num_iter) {
    _check_model_ready();
    tr->start = tr->nodep[1];
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    pllOptimizeBranchLengths(tr, partitions, num_iter);
}

void pll::optimise(bool rates, bool freqs, bool alphas, bool branches,
                   bool topology, int tree_search_interval, bool final_tree_search) {
    if (!rates && !freqs && !alphas && !branches && !topology) return;
    int i = 0;
    _check_model_ready();
    linkageList *alphaList = partitions->alphaList,
                *rateList  = partitions->rateList,
                *freqList  = partitions->freqList;
    double current_likelihood;
    bool tree_optimised = false;
    double modelEpsilon = 0.0001; // same as in modOpt
    tr->start = tr->nodep[1];
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    for (;;) {
        i++;
        current_likelihood = tr->likelihood;
        std::cout << "  iter " << i << " current lnl = " << current_likelihood << std::endl;

        if (rates) {
            std::cout << "    optimising Q-matrix:     " << std::flush;
            pllOptRatesGeneric(tr, partitions, modelEpsilon, rateList);
            pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
            std::cout << tr->likelihood << std::endl;
        }

        if (branches) {
            std::cout << "    optimising branches:     " << std::flush;
            pllOptimizeBranchLengths(tr, partitions, 16);
            pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
            std::cout << tr->likelihood << std::endl;
        }

        if (freqs) {
            std::cout << "    optimising frequencies:  " << std::flush;
            pllOptBaseFreqs(tr, partitions, modelEpsilon, freqList);
            pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
            std::cout << tr->likelihood << std::endl;
        }

        if (branches && freqs) {
            std::cout << "    optimising branches:     " << std::flush;
            pllOptimizeBranchLengths(tr, partitions, 16);
            pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
            std::cout << tr->likelihood << std::endl;
        }

        if (alphas) {
            std::cout << "    optimising alphas:       " << std::flush;
            pllOptAlphasGeneric (tr, partitions, modelEpsilon, alphaList);
            pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
            std::cout << tr->likelihood << std::endl;
        }

        if (branches && alphas) {
            std::cout << "    optimising branches:     " << std::flush;
            pllOptimizeBranchLengths(tr, partitions, 16);
            pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
            std::cout << tr->likelihood << std::endl;
        }

        // Error - likelihood got worse!
        if(tr->likelihood < (current_likelihood - tr->likelihoodEpsilon)) {
            std::cout << tr->likelihood << " " << current_likelihood << std::endl;
            std::cout << "Difference: " << tr->likelihood - current_likelihood << std::endl;
            throw std::logic_error("Likelihood got worse");
        }

        if (topology && (i % tree_search_interval == 0)) {
            std::cout << "Optimising tree topology:    " << std::flush;
            pllRaxmlSearchAlgorithm(tr, partitions, PLL_FALSE);
            pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
            std::cout << tr->likelihood << std::endl;
            tree_optimised = true;
        }

        // Termination step - difference in likelihood between the start and end of this step is negligible
        if (!(fabs(current_likelihood - tr->likelihood) > tr->likelihoodEpsilon) && (tree_optimised || !topology)) {
            if (final_tree_search && topology) {
                std::cout << "Doing final tree search." << std::endl;
                pllRaxmlSearchAlgorithm(tr, partitions, PLL_FALSE);
                if (rates) {
                    pllOptRatesGeneric(tr, partitions, modelEpsilon, rateList);
                    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
                }

                if (branches) {
                    pllOptimizeBranchLengths(tr, partitions, 32);
                    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
                }

                if (freqs) {
                    pllOptBaseFreqs(tr, partitions, modelEpsilon, freqList);
                    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
                }

                if (branches && freqs) {
                    pllOptimizeBranchLengths(tr, partitions, 32);
                    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
                }

                if (alphas) {
                    pllOptAlphasGeneric (tr, partitions, modelEpsilon, alphaList);
                    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
                }

                if (branches && alphas) {
                    pllOptimizeBranchLengths(tr, partitions, 32);
                    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
                }

            }
            std::cout << "current lnl = " << current_likelihood << std::endl
                 << "tr lnl      = " << tr->likelihood << std::endl
                 << "END." << std::endl;
            break;
        }
    }
}

double pll::get_likelihood() {
    _check_model_ready();
    tr->start = tr->nodep[1];
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
    return tr->likelihood;
}

int pll::get_number_of_threads() {
    return attr.numberOfThreads;
}

int pll::get_number_of_partitions() {
    return partitions->numberOfPartitions;
}

std::vector<std::string> pll::get_partition_names() {
    _check_model_ready();
    std::vector<std::string> names;
    for (int i = 0; i < get_number_of_partitions(); ++i) {
        names.push_back(get_partition_name(i));
    }
    return names;
}

std::string pll::get_partition_name(int partition) {
    _check_model_ready();
    _check_partitions_bounds(partition);
    return std::string(partitions->partitionData[partition]->partitionName);
}

std::vector<double> pll::get_alphas() {
    _check_model_ready();
    std::vector<double> alphas;
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        alphas.push_back(get_alpha(i));
    }
    return alphas;
}

double pll::get_alpha(int partition) {
    _check_model_ready();
    _check_partitions_bounds(partition);
    return pllGetAlpha(partitions, partition);
}

std::vector<std::vector<double>> pll::get_frequencies() {
    _check_model_ready();
    std::vector<std::vector<double>> freqs;
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        freqs.push_back(get_frequencies_vector(i));
    }
    return freqs;
}

std::vector<double> pll::get_frequencies_vector(int partition) {
    _check_model_ready();
    _check_partitions_bounds(partition);
    std::vector<double> freqs_vec;
    if (!_model_ready) {
        std::cerr << "Model isn't finalised" << std::endl;
        return freqs_vec;
    }
    int num_states = partitions->partitionData[partition]->states;
    for (int j = 0; j < num_states; ++j) {
        freqs_vec.push_back(partitions->partitionData[partition]->frequencies[j]);
    }
    return freqs_vec;
}

std::vector<std::vector<double>> pll::get_rates() {
    _check_model_ready();
    std::vector<std::vector<double>> rates;
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        rates.push_back(get_rates_vector(i));
    }
    return rates;
}

std::vector<double> pll::get_rates_vector(int partition) {
    _check_model_ready();
    _check_partitions_bounds(partition);
    std::vector<double> rates_vec;
    if (!_model_ready) {
        std::cerr << "Model isn't finalised" << std::endl;
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

std::string pll::get_tree() {
    _check_model_ready();
    pllTreeToNewick(tr->tree_string, tr, partitions,
                    tr->start->back, PLL_TRUE, PLL_TRUE,
                    PLL_FALSE, PLL_FALSE, PLL_FALSE,
                    0, PLL_FALSE, PLL_FALSE);
    std::string tree(tr->tree_string);
    tree.erase(std::remove(tree.begin(), tree.end(), '\n'), tree.end());
    return tree;
}

std::vector<std::vector<double> > pll::get_empirical_frequencies() {
    _check_model_ready();
    double ** ef;
    std::vector<std::vector<double>> vec_2d;
    ef = pllBaseFrequenciesInstance(tr, partitions);
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        std::vector<double> row_vec;
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

void pll::set_rates(std::vector<double> rates,
        int partition, bool optimisable) {
    _check_model_ready();
//    double old_fracchange = get_frac_change();
    if (is_dna(partition)) {
        _check_partitions_bounds(partition);
        int num_states = partitions->partitionData[partition]->states;
        size_t num_rates = (num_states * (num_states - 1)) / 2;
        if (rates.size() != num_rates) {
            std::ostringstream errmsg;
            errmsg << "Rates vector is the wrong length. Should be " << num_rates << std::endl;
            throw std::invalid_argument(errmsg.str());
        }
        pllSetSubstitutionMatrix(&(rates[0]), num_rates, partition, partitions, tr);
        set_optimisable_rates(partition, optimisable);
//        double new_fracchange = get_frac_change();
//        _update_q_matrix_and_brlens(partition, old_fracchange, new_fracchange);
    }
    else {
        std::ostringstream errmsg;
        errmsg << "Can only set rates directly for DNA partitions" << std::endl;
        throw std::domain_error(errmsg.str());
    }
}

void pll::set_alpha(double alpha, int partition, bool optimisable) {
    _check_model_ready();
    _check_partitions_bounds(partition);
    pllSetFixedAlpha(alpha, partition, partitions, tr);
    set_optimisable_alpha(partition, optimisable);
}

void pll::set_frequencies(std::vector<double> freqs, int partition, bool optimisable) {
    _check_model_ready();
    _check_partitions_bounds(partition);
    if (!_approx_eq(_vector_sum(freqs), 1)) {
        std::ostringstream errmsg;
        errmsg << "Not setting frequencies: Frequencies do not sum to 1" << std::endl;
        throw std::invalid_argument(errmsg.str());
    }
    size_t num_states = partitions->partitionData[partition]->states;
    if (freqs.size() != num_states) {
        std::ostringstream errmsg;
        errmsg << "Frequencies vector is the wrong length. Should be " << num_states << std::endl;
        throw std::invalid_argument(errmsg.str());
    }
    set_optimisable_frequencies(partition, true); // frequencies only updated if optimisable flag is true
    pllSetFixedBaseFrequencies(&(freqs[0]), num_states, partition, partitions, tr);
    set_optimisable_frequencies(partition, optimisable);
}

void pll::set_tree(std::string nwk) {
    _init_tree_string(nwk);
    _evaluate_likelihood();
}

void pll::link_alpha_parameters(std::string linkage) {
    _check_model_ready();
    pllLinkAlphaParameters(const_cast<char*>(linkage.c_str()), partitions);
}

void pll::link_frequencies(std::string linkage) {
    _check_model_ready();
    pllLinkFrequencies(const_cast<char*>(linkage.c_str()), partitions);
}

void pll::link_rates(std::string linkage) {
    _check_model_ready();
    pllLinkRates(const_cast<char*>(linkage.c_str()), partitions);
}

bool pll::is_dna(int partition) {
    _check_partitions_bounds(partition);
    return (partitions->partitionData[partition]->dataType == PLL_DNA_DATA);
}

bool pll::is_protein(int partition) {
    _check_partitions_bounds(partition);
        return (partitions->partitionData[partition]->dataType == PLL_AA_DATA);
}

void pll::set_number_of_threads(int threads) {
    attr.numberOfThreads = threads;
}

bool pll::is_optimisable_alpha(int partition) {
    _check_partitions_bounds(partition);
    return (partitions->partitionData[partition]->optimizeAlphaParameter == PLL_TRUE);
}

bool pll::is_optimisable_frequencies(int partition) {
    _check_partitions_bounds(partition);
    return (partitions->partitionData[partition]->optimizeBaseFrequencies == PLL_TRUE);
}

bool pll::is_optimisable_rates(int partition) {
    _check_partitions_bounds(partition);
    return (partitions->partitionData[partition]->optimizeSubstitutionRates == PLL_TRUE);
}

std::vector<std::string> pll::get_model_names() {
    std::vector<std::string> names;
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        names.push_back(get_model_name(i));
    }
    return names;
}

std::string pll::get_model_name(int partition) {
    _check_partitions_bounds(partition);
    std::string name;
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

void pll::set_optimisable_rates(int partition, bool optimisable) {
    if (is_protein(partition)) {
        std::ostringstream errmsg;
        errmsg << "Optimising rates not implemented for protein models" << std::endl;
        throw std::domain_error(errmsg.str());
    }
    _check_partitions_bounds(partition);
    int pll_bool = optimisable ? PLL_TRUE : PLL_FALSE;
    if (is_dna(partition)) {
        partitions->partitionData[partition]->optimizeSubstitutionRates = pll_bool;
        partitions->dirty = PLL_TRUE;
        _evaluate_likelihood();
    }
}

void pll::set_optimisable_alpha(int partition, bool optimisable) {
    _check_partitions_bounds(partition);
    int pll_bool = optimisable ? PLL_TRUE : PLL_FALSE;
    partitions->partitionData[partition]->optimizeAlphaParameter = pll_bool;
    partitions->dirty = PLL_TRUE;
    _evaluate_likelihood();
}

//double pll::get_frac_change() {
//    return tr->fracchange;
//}

void pll::set_optimisable_frequencies(int partition, bool optimisable) {
    _check_partitions_bounds(partition);
    int pll_bool = optimisable ? PLL_TRUE : PLL_FALSE;
    partitions->partitionData[partition]->optimizeBaseFrequencies = pll_bool;
    partitions->dirty = PLL_TRUE;
    _evaluate_likelihood();
}

// PRIVATE METHODS

void pll::_init_attr(int num_threads, long rns) {
    attr.rateHetModel = PLL_GAMMA;
    attr.fastScaling = PLL_TRUE;
    attr.saveMemory = PLL_TRUE;
    attr.useRecom = PLL_FALSE;
    attr.randomNumberSeed = rns;
    attr.numberOfThreads = num_threads;
}

void pll::_init_instance() {
    tr = pllCreateInstance(&attr);
    _instance_ready = true;
}

void pll::_init_alignment_file(std::string path) {
    if (!_is_file(path)) {
        std::ostringstream errmsg;
        errmsg << "Couldn't find the alignment file " << path << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    alignment = std::unique_ptr<pllAlignmentData, AlignmentDeleter>(pllParseAlignmentFile(PLL_FORMAT_PHYLIP, path.c_str()));
    if (!alignment) {
        std::ostringstream errmsg;
        errmsg << "Couldn't parse the alignment at " << path << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    _alignment_ready = true;
}

void pll::_init_partition_file(std::string path) {
    if (!_alignment_ready) {
        std::ostringstream errmsg;
        errmsg << "Must load alignment before partitions" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    pllQueue * partitionInfo;
    partitionInfo = pllPartitionParse(path.c_str());
    if (!pllPartitionsValidate(partitionInfo, alignment.get())) {
        std::ostringstream errmsg;
        errmsg << "partitions parse error" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    partitions = pllPartitionsCommit(partitionInfo, alignment.get());
    pllAlignmentRemoveDups(alignment.get(), partitions);
    pllQueuePartitionsDestroy(&partitionInfo);
    _partitions_ready = true;
}

void pll::_init_partition_string(std::string p_string) {
    if (!_alignment_ready) {
        std::ostringstream errmsg;
        errmsg << "Must load alignment before partitions" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    pllQueue * partitionInfo;
    partitionInfo = pllPartitionParseString(p_string.c_str());
    if (!pllPartitionsValidate(partitionInfo, alignment.get())) {
        std::cerr << "partitions parse error" << std::endl;
        throw std::runtime_error("partitions parse error");
    }
    partitions = pllPartitionsCommit(partitionInfo, alignment.get());
    pllAlignmentRemoveDups(alignment.get(), partitions);
    pllQueuePartitionsDestroy(&partitionInfo);
    _partitions_ready = true;
}

void pll::_init_tree_file(std::string path) {
    if (!_alignment_ready || !_partitions_ready) {
        std::ostringstream errmsg;
        errmsg << "Must load alignment and partitions before tree" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    newick = std::unique_ptr<pllNewickTree, NewickDeleter>(pllNewickParseFile(path.c_str()));
    if (!newick) {
        std::ostringstream errmsg;
        errmsg << "tree parse error" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    if (!pllValidateNewick(newick.get())) /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */ {
        std::ostringstream errmsg;
        errmsg << "invalid tree" << std::endl;
        throw std::invalid_argument(errmsg.str());
    }
    pllTreeInitTopologyNewick(tr, newick.get(), PLL_FALSE);
    newick.reset();
    _tree_ready = true;
}

void pll::_init_tree_string(std::string nwk) {
    if (!_alignment_ready || !_partitions_ready) {
        std::ostringstream errmsg;
        errmsg << "Must load alignment and partitions before tree" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    newick = std::unique_ptr<pllNewickTree, NewickDeleter>(pllNewickParseString(nwk.c_str()));
    if (!newick) {
        std::ostringstream errmsg;
        errmsg << "Newick parse error" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    if (!pllValidateNewick(newick.get())) /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */ {
        std::ostringstream errmsg;
        errmsg << "Invalid newick tree" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    pllTreeInitTopologyNewick(tr, newick.get(), PLL_FALSE);
    newick.reset();
    _tree_ready = true;
}

void pll::_init_tree_random() {
    if (!_alignment_ready || !_partitions_ready) {
        std::ostringstream errmsg;
        errmsg << "Must load alignment and partitions before tree" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    pllTreeInitTopologyRandom(tr, alignment->sequenceCount, alignment->sequenceLabels);
    _tree_ready = true;
}

void pll::_init_model(bool parsimony_tree) {
    if (!_instance_ready || !_alignment_ready || !_partitions_ready || !_tree_ready) {
        std::ostringstream errmsg;
        errmsg << "Must load alignment, tree and partitions before initialising the model" << std::endl;
        throw std::runtime_error(errmsg.str());
    }
    if (!pllLoadAlignment(tr, alignment.get(), partitions)) {
        throw std::runtime_error("Model finalisation error");
    }
    if (parsimony_tree) {
        pllComputeRandomizedStepwiseAdditionParsimonyTree(tr, partitions);
    }
    pllInitModel(tr, partitions);
    alignment.reset();
    _model_ready = true;
}

void pll::_destroy_model() {
    if (partitions) pllPartitionsDestroy(tr, &partitions);
    if (tr) pllDestroyInstance(tr);
}

bool pll::_is_file(std::string filename) {
    std::ifstream fl(filename.c_str());
    bool result = true;
    if (!fl) {
        result = false;
    }
    fl.close();
    return result;
}

bool pll::_is_tree_string(std::string tree_string) {
    size_t l = tree_string.length();
    return (tree_string[0]=='(' && tree_string[l-1]==';');
}

void pll::_evaluate_likelihood() {
    if (!_model_ready) {
        std::cerr << "WARNING: Model isn't finalised." << std::endl;
        return;
    }
    pllEvaluateLikelihood(tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);
}

double pll::_vector_sum(std::vector<double> vec) {
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

void pll::_check_partitions_bounds(int partition) {
    int max_partitions = get_number_of_partitions();
    if (partition >= max_partitions) {
        std::ostringstream errmsg;
        errmsg << "The model has " << max_partitions << " partitions" << std::endl;
        throw std::range_error(errmsg.str());
    }
}

void pll::_check_model_ready() {
    if (!_model_ready) {
        throw std::runtime_error("The model isn't ready for this operation");
    }
}

std::string pll::_model_name(int model_num) {
    std::string name;
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

        case PLL_STMTREV :
        name = "STMTREV";
        break;

        case PLL_LG4M :
        name = "LG4M";
        break;

        case PLL_LG4X :
        name = "LG4X";
        break;

        case PLL_AUTO :
        name = "AUTO";
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

//void pll::_update_q_matrix_and_brlens(int model, double old_fracchange, double new_fracchange) {
//    pllInitReversibleGTR(tr, partitions, model);
//#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
//    pllMasterBarrier (tr, partitions, PLL_THREAD_COPY_RATES);
//#endif
//    _update_all_brlens (old_fracchange, new_fracchange);
//}
//
//void pll::_update_all_brlens(double old_fracchange, double new_fracchange) {
//    nodeptr p;
//
//    p = tr->start;
//    if(!(isTip(p->number, tr->mxtips))) throw std::exception();
//
//    _update_brlens_recursive(p->back, tr->mxtips, old_fracchange, new_fracchange);
//}
//
//void pll::_update_brlens_recursive(nodeptr p, int tips, double old_fracchange, double new_fracchange) {
//    _update_brlen(p, old_fracchange, new_fracchange);
//
//    if (!isTip (p->number, tips)) {
//        _update_brlens_recursive(p->next->back, tips, old_fracchange, new_fracchange);
//        _update_brlens_recursive(p->next->next->back, tips, old_fracchange, new_fracchange);
//    }
//}
//
//void pll::_update_brlen(nodeptr p, double old_fracchange, double new_fracchange) {
//    double z;
//    int j;
//
//    for (j = 0; j < PLL_NUM_BRANCHES; ++ j) {
//        z = exp ((log (p->z[j]) * old_fracchange) / new_fracchange);
//        if (z < PLL_ZMIN) z = PLL_ZMIN;
//        if (z > PLL_ZMAX) z = PLL_ZMAX;
//        p->z[j] = p->back->z[j] = z;
//    }
//}

bool pll::isTip(int number, int maxTips) {
    if(!(number > 0)) {
        throw std::domain_error("Number must be positive.");
    }
    if(number <= maxTips)
        return PLL_TRUE;
    else
        return PLL_FALSE;
}
