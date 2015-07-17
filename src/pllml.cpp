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

#include <stdexcept>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "pllml.h"
extern "C" {
#include <pll/pll.h>
}

// PUBLIC METHODS

pll::pll(string alignment_file, string part_desc, string tree, int num_threads, long rns) {
    auto attr = _init_attr(num_threads, rns);
    alignmentUPtr alignment = _parse_alignment_file(alignment_file);
    queueUPtr queue = _parse_partitions(part_desc, std::move(alignment));
    newickUPtr newick = _parse_tree(tree);
    tr = instanceUPtr(pllCreateInstance(&attr));

    partitions = pllPartitionsCommit(queue.get(), alignment.get());
    pllAlignmentRemoveDups(alignment.get(), partitions);
    pllTreeInitTopologyNewick(tr.get(), newick.get(), PLL_FALSE);
    pllLoadAlignment(tr.get(), alignment.get(), partitions);
    pllInitModel(tr.get(), partitions);
}

pll::pll(string alignment_file, string part_desc, bool parsimony, int num_threads, long rns) {
    auto attr = _init_attr(num_threads, rns);
    alignmentUPtr alignment = _parse_alignment_file(alignment_file);
    queueUPtr queue = _parse_partitions(part_desc, std::move(alignment));
    tr = instanceUPtr(pllCreateInstance(&attr));

    partitions = pllPartitionsCommit(queue.get(), alignment.get());
    pllAlignmentRemoveDups(alignment.get(), partitions);
    pllTreeInitTopologyRandom(tr.get(), alignment->sequenceCount, alignment->sequenceLabels);
    pllLoadAlignment(tr.get(), alignment.get(), partitions);
    if (parsimony) pllComputeRandomizedStepwiseAdditionParsimonyTree(tr.get(), partitions);
    pllInitModel(tr.get(), partitions);
}

pll::~pll() {
    pllPartitionsDestroy(tr.get(), &partitions);
    //std::cout << "Destroyed partitions list" << std::endl;
}

void pll::optimise_tree_search(bool estimate_model) {
    tr->start = tr->nodep[1];
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
    int pll_bool = estimate_model ? PLL_TRUE : PLL_FALSE;
    pllRaxmlSearchAlgorithm(tr.get(), partitions, pll_bool);
}

void pll::optimise_model() {
    optimise(true, true, true, false);
}

void pll::optimise_rates() {
    optimise(true, false, false, false);
}

void pll::optimise_freqs() {
    optimise(false, true, false, false);
}

void pll::optimise_alphas() {
    optimise(false, false, true, false);
}

void pll::optimise_branch_lengths(int num_iter) {
    tr->start = tr->nodep[1];
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
    pllOptimizeBranchLengths(tr.get(), partitions, num_iter);
}

void pll::optimise(bool rates, bool freqs, bool alphas, bool branches) {
    if (!rates && !freqs && !alphas && !branches) return;
    int i = 0;
    double current_likelihood;
    double modelEpsilon = 0.0001; // same as in modOpt
    tr->start = tr->nodep[1];
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
    for (;;) {
        i++;
        current_likelihood = tr->likelihood;
        cerr << "  iter " << i << " current lnl = " << current_likelihood << endl;

        if (rates) {
            pllOptRatesGeneric(tr.get(), partitions, modelEpsilon, partitions->rateList);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            cerr << "    rates:  " << tr->likelihood << endl;
        }

        if (branches) {
            pllOptimizeBranchLengths(tr.get(), partitions, 3);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            cerr << "    brlen1: " << tr->likelihood << endl;
        }

        if (freqs) {
            pllOptBaseFreqs(tr.get(), partitions, modelEpsilon, partitions->freqList);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            cerr << "    freqs:  " << tr->likelihood << endl;
        }

        if (branches) {
            pllOptimizeBranchLengths(tr.get(), partitions, 3);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            cerr << "    brlen2: " << tr->likelihood << endl;
        }

        if (alphas) {
            pllOptAlphasGeneric (tr.get(), partitions, modelEpsilon, partitions->alphaList);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            cerr << "    alphas: " << tr->likelihood << endl;
        }

        if (branches) {
            pllOptimizeBranchLengths(tr.get(), partitions, 3);
            pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
            cerr << "    brlen3: " << tr->likelihood << endl;
        }

        if(tr->likelihood - current_likelihood <= 0) {
            cerr << tr->likelihood << " " << current_likelihood << endl;
            cerr << "Difference: " << tr->likelihood - current_likelihood << endl;
            throw exception();
        }

        if (fabs(current_likelihood - tr->likelihood) <= tr->likelihoodEpsilon){
            cerr << "current lnl = " << current_likelihood << endl
                 << "tr lnl      = " << tr->likelihood << endl
                 << "END" << endl;
            break;
        }
    }
}

double pll::get_likelihood() {
    tr->start = tr->nodep[1];
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
    return tr->likelihood;
}

int pll::get_number_of_partitions() {
    return partitions->numberOfPartitions;
}

vector<string> pll::get_partition_names() {
    vector<string> names;
    for (int i = 0; i < get_number_of_partitions(); ++i) {
        names.push_back(get_partition_name(i));
    }
    return names;
}

string pll::get_partition_name(int partition) {
    _check_partitions_bounds(partition);
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
    _check_partitions_bounds(partition);
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
    _check_partitions_bounds(partition);
    vector<double> freqs_vec;
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
    _check_partitions_bounds(partition);
    vector<double> rates_vec;
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
    pllTreeToNewick(tr->tree_string, tr.get(), partitions,
                    tr->start->back, PLL_TRUE, PLL_TRUE,
                    PLL_FALSE, PLL_FALSE, PLL_FALSE,
                    0, PLL_FALSE, PLL_FALSE);
    string tree(tr->tree_string);
    tree.erase(std::remove(tree.begin(), tree.end(), '\n'), tree.end());
    return tree;
}

vector<vector<double> > pll::get_empirical_frequencies() {
    double ** ef;
    vector<vector<double>> vec_2d;
    ef = pllBaseFrequenciesInstance(tr.get(), partitions);
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        vector<double> row_vec;
        size_t ns = partitions->partitionData[i]->states;
        for (size_t j = 0; j < ns; ++j) {
            row_vec.push_back(ef[i][j]);
        }
        vec_2d.push_back(row_vec);
    }
    pllEmpiricalFrequenciesDestroy (&ef, partitions->numberOfPartitions);
    return vec_2d;
}

void pll::set_epsilon(double epsilon) {
    tr->likelihoodEpsilon = epsilon;
}

void pll::set_rates(vector<double> rates,
        int partition, bool optimisable) {
    double old_fracchange = get_frac_change();
    if (is_dna(partition)) {
        _check_partitions_bounds(partition);
        int num_states = partitions->partitionData[partition]->states;
        size_t num_rates = (num_states * (num_states - 1)) / 2;
        if (rates.size() != num_rates) {
            cerr << "Rates vector is the wrong length. Should be " << num_rates << endl;
            throw exception();
        }
        pllSetSubstitutionMatrix(&(rates[0]), num_rates, partition, partitions, tr.get());
        set_optimisable_rates(partition, optimisable);
        double new_fracchange = get_frac_change();
        _update_q_matrix_and_brlens(partition, old_fracchange, new_fracchange);
    }
}

void pll::set_alpha(double alpha, int partition, bool optimisable) {
    _check_partitions_bounds(partition);
    pllSetFixedAlpha(alpha, partition, partitions, tr.get());
    set_optimisable_alpha(partition, optimisable);
}

void pll::set_frequencies(vector<double> freqs, int partition, bool optimisable) {
    _check_partitions_bounds(partition);
    if (!_approx_eq(_vector_sum(freqs), 1)) {
        cerr << "Not setting frequencies: Frequencies do not sum to 1" << endl;
        throw exception();
    }
    size_t num_states = partitions->partitionData[partition]->states;
    if (freqs.size() != num_states) {
        cerr << "Frequencies vector is the wrong length. Should be " << num_states << endl;
        throw exception();
    }
    set_optimisable_frequencies(partition, true); // frequencies only updated if optimisable flag is true
    pllSetFixedBaseFrequencies(&(freqs[0]), num_states, partition, partitions, tr.get());
    set_optimisable_frequencies(partition, optimisable);
}

void pll::set_tree(string nwk) {
    newickUPtr newick = _parse_tree(nwk);
    pllTreeInitTopologyNewick(tr.get(), newick.get(), PLL_FALSE);
    _evaluate_likelihood();
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
    _check_partitions_bounds(partition);
    return (partitions->partitionData[partition]->dataType == PLL_DNA_DATA);
}

bool pll::is_protein(int partition) {
    _check_partitions_bounds(partition);
        return (partitions->partitionData[partition]->dataType == PLL_AA_DATA);
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

vector<string> pll::get_model_names() {
    vector<string> names;
    size_t np = get_number_of_partitions();
    for (size_t i = 0; i < np; ++i) {
        names.push_back(get_model_name(i));
    }
    return names;
}

string pll::get_model_name(int partition) {
    _check_partitions_bounds(partition);
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

void pll::set_optimisable_rates(int partition, bool optimisable) {
    if (is_protein(partition)) {
        cerr << "Optimising rates not implemented for protein models" << endl;
        throw exception();
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

double pll::get_frac_change() {
    return tr->fracchange;
}

void pll::set_optimisable_frequencies(int partition, bool optimisable) {
    _check_partitions_bounds(partition);
    int pll_bool = optimisable ? PLL_TRUE : PLL_FALSE;
    partitions->partitionData[partition]->optimizeBaseFrequencies = pll_bool;
    partitions->dirty = PLL_TRUE;
    _evaluate_likelihood();
}

// PRIVATE METHODS

pllInstanceAttr pll::_init_attr(int num_threads, long rns) {
    pllInstanceAttr attr;
    attr.rateHetModel = PLL_GAMMA;
    attr.fastScaling = PLL_TRUE;
    attr.saveMemory = PLL_TRUE;
    attr.useRecom = PLL_FALSE;
    attr.randomNumberSeed = rns;
    attr.numberOfThreads = num_threads;
    return attr;
}

alignmentUPtr pll::_parse_alignment_file(string path) {
    if (!_is_file(path)) {
        cerr << "Couldn't find the alignment file " << path << endl;
        throw exception();
    }
    alignmentUPtr alignment;
    alignment = alignmentUPtr(pllParseAlignmentFile(PLL_FORMAT_PHYLIP, path.c_str()));
    if (!alignment) {
        //cout << "Trying to parse as fasta" << endl;
        alignment = alignmentUPtr(pllParseAlignmentFile(PLL_FORMAT_FASTA, path.c_str()));
    }
    if (!alignment) {
        cerr << "Couldn't parse the alignment at " << path << endl;
        throw exception();
    }
    return alignment;
}

queueUPtr pll::_parse_partitions(string partitions, alignmentUPtr&& alignment) {
    queueUPtr q;
    if (_is_file(partitions)) {
        q = queueUPtr(pllPartitionParse(partitions.c_str()));
    }
    else {
        q = queueUPtr(pllPartitionParseString(partitions.c_str()));
    }
    if (!pllPartitionsValidate(q.get(), alignment.get())) {
        cerr << "partitions parse error" << endl;
        throw exception();
    }
    return q;
}

newickUPtr pll::_parse_tree(string tree) {
    newickUPtr newick;
    if (_is_file(tree)) {
        newick = newickUPtr(pllNewickParseFile(tree.c_str()));
    }
    else {
        newick = newickUPtr(pllNewickParseString(tree.c_str()));
    }
    if (!newick) {
        cerr << "tree parse error" << endl;
        throw exception();
    }
    if (!pllValidateNewick(newick.get())) /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */ {
        cerr << "invalid tree" << endl;
        throw exception();
    }
    return newick;
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

bool pll::_is_tree_string(string tree_string) {
    size_t l = tree_string.length();
    return (tree_string[0]=='(' && tree_string[l-1]==';');
}

void pll::_evaluate_likelihood() {
    pllEvaluateLikelihood(tr.get(), partitions, tr->start, PLL_TRUE, PLL_FALSE);
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

void pll::_check_partitions_bounds(int partition) {
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

void pll::_update_q_matrix_and_brlens(int model, double old_fracchange, double new_fracchange) {
    pllInitReversibleGTR(tr.get(), partitions, model);
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
    pllMasterBarrier (tr.get(), partitions, PLL_THREAD_COPY_RATES);
#endif
    _update_all_brlens (old_fracchange, new_fracchange);
}

void pll::_update_all_brlens(double old_fracchange, double new_fracchange) {
    nodeptr p;

    p = tr->start;
    if(!(isTip(p->number, tr->mxtips))) throw exception();

    _update_brlens_recursive(p->back, tr->mxtips, old_fracchange, new_fracchange);
}

void pll::_update_brlens_recursive(nodeptr p, int tips, double old_fracchange, double new_fracchange) {
    _update_brlen(p, old_fracchange, new_fracchange);

    if (!isTip (p->number, tips)) {
        _update_brlens_recursive(p->next->back, tips, old_fracchange, new_fracchange);
        _update_brlens_recursive(p->next->next->back, tips, old_fracchange, new_fracchange);
    }
}

void pll::_update_brlen(nodeptr p, double old_fracchange, double new_fracchange) {
    double z;
    int j;

    for (j = 0; j < PLL_NUM_BRANCHES; ++ j) {
        z = exp ((log (p->z[j]) * old_fracchange) / new_fracchange);
        if (z < PLL_ZMIN) z = PLL_ZMIN;
        if (z > PLL_ZMAX) z = PLL_ZMAX;
        p->z[j] = p->back->z[j] = z;
    }
}

bool pll::isTip(int number, int maxTips) {
    if(number <= 0) throw exception();
    if(number <= maxTips)
        return PLL_TRUE;
    else
        return PLL_FALSE;
}

