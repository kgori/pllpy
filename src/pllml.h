/*
pllml.h
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

#ifndef PLL_WRAPPER_H_
#define PLL_WRAPPER_H_

#include <iostream>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>
extern "C" {
#include <pll/pll.h>
}

using namespace std;

struct AlignmentDeleter {
    void operator()(pllAlignmentData *pAlignment) {
        pllAlignmentDataDestroy(pAlignment);
        //std::cout << "Alignment destroyed" << std::endl;
    }
};

struct QueueDeleter {
    void operator()(pllQueue *pQueue) {
        pllQueuePartitionsDestroy(&pQueue);
        //std::cout << "Queue destroyed" << std::endl;
    }
};

struct NewickDeleter {
    void operator()(pllNewickTree *pNewick) {
        pllNewickParseDestroy(&pNewick);
        //std::cout << "Newick destroyed" << std::endl;
    }
};

struct InstanceDeleter {
    void operator()(pllInstance *pInstance) {
        pllDestroyInstance(pInstance);
        //std::cout << "Instance Destroyed" << std::endl;
    }
};

typedef std::unique_ptr<pllAlignmentData, AlignmentDeleter> alignmentUPtr;
typedef std::unique_ptr<pllNewickTree, NewickDeleter> newickUPtr;
typedef std::unique_ptr<pllQueue, QueueDeleter> queueUPtr;
typedef std::unique_ptr<pllInstance, InstanceDeleter> instanceUPtr;

class pll {
public:
    // Make and destroy
    pll(string alignment_file, string part_desc, string tree, int num_threads = 1, long rns=0xDEADBEEF);
    pll(string alignment_file, string part_desc, bool parsimony, int num_threads = 1, long rns=0xDEADBEEF);
    virtual ~pll();
    pll(pll&& rhs) = delete;               // No copy       - need to deep copy all pointed-at data
    pll& operator=(pll&& rhs) = delete;    // No assignment - (seems like too much work)

    // Run optimisations
    void                   optimise(bool rates, bool freqs, bool alphas, bool branches);
    void                   optimise_alphas();
    void                   optimise_branch_lengths(int num_iter=32);
    void                   optimise_freqs();
    void                   optimise_model();
    void                   optimise_rates();
    void                   optimise_tree_search(bool estimate_model);

    // Getters
    double                 get_likelihood();
    vector<string>         get_partition_names();
    string                 get_partition_name(int partition);
    vector<string>         get_model_names();
    string                 get_model_name(int partition);
    double                 get_epsilon();
    vector<double>         get_alphas();
    double                 get_alpha(int partition);
    vector<vector<double>> get_frequencies();
    vector<double>         get_frequencies_vector(int partition);
    vector<vector<double>> get_rates();
    vector<double>         get_rates_vector(int partition);
    int                    get_number_of_partitions();
    string                 get_tree();
    vector<vector<double>> get_empirical_frequencies();
    double                 get_frac_change();

    // Setters
    void                   set_epsilon(double epsilon);
    void                   set_alpha(double alpha, int partition, bool optimisable);
    void                   set_frequencies(vector<double> freqs, int partition, bool optimisable);
    void                   set_rates(vector<double> rates, int partition, bool optimisable);
    void                   set_optimisable_alpha(int partition, bool optimisable);
    void                   set_optimisable_frequencies(int partition, bool optimisable);
    void                   set_optimisable_rates(int partition, bool optimisable);
    void                   set_tree(string nwk);

    // Partition management
    void                   link_alpha_parameters(string linkage);
    void                   link_frequencies(string linkage);
    void                   link_rates(string linkage);

    // Check settings
    bool                   is_dna(int partition);
    bool                   is_protein(int partition);
    bool                   is_optimisable_alpha(int partition);
    bool                   is_optimisable_frequencies(int partition);
    bool                   is_optimisable_rates(int partition);

private:
    // Model initialisation
    pllInstanceAttr        _init_attr(int num_threads, long rns);
    alignmentUPtr          _parse_alignment_file(string path);
    queueUPtr              _parse_partitions(string partitions, alignmentUPtr&& alignment);
    newickUPtr             _parse_tree(string path);

    // Helper functions
    void                   _evaluate_likelihood();
    double                 _vector_sum(vector<double>);
    bool                   _approx_eq(double a, double b);
    bool                   _is_file(string filename);
    bool                   _is_tree_string(string tree_string);
    string                 _model_name(int model_num);
    void                   _destroy_model();
    void                   _update_q_matrix_and_brlens(int model, double old_fracchange, double new_fracchange);
    void                   _update_all_brlens(double old_fracchange, double new_fracchange);
    void                   _update_brlens_recursive(nodeptr p, int tips, double old_fracchange, double new_fracchange);
    void                   _update_brlen(nodeptr p, double old_fracchange, double new_fracchange);
    bool                   isTip(int number, int maxTips);

    // Error checking
    void                   _check_partitions_bounds(int partition);
    void                   _check_model_ready();

    // Data
    instanceUPtr tr;
    partitionList* partitions;
};

#endif /* PLL_WRAPPER_H_ */
