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
#include <string>
#include <vector>
extern "C" {
#include <pll/pll.h>
}

using namespace std;

class pll {
public:
    // Make and destroy
    pll(string alignment_file, string partitions, string tree, int num_threads = 1, long rns=0xDEADBEEF);
    pll(string alignment_file, string partitions, bool parsimony, int num_threads = 1, long rns=0xDEADBEEF);
    virtual ~pll();

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
    int                    get_number_of_threads();
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
    void                   set_number_of_threads(int threads);

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
    void                   _init_attr(int num_threads, long rns);
    void                   _init_instance();
    void                   _init_alignment_file(string path);
    void                   _init_partition_file(string path);
    void                   _init_partition_string(string p_string);
    void                   _init_tree_file(string path);
    void                   _init_tree_string(string nwk);
    void                   _init_tree_random();
    void                   _init_model(bool parsimony);

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
    pllAlignmentData *alignment  = nullptr;
    pllInstance      *tr         = nullptr;
    pllNewickTree    *newick     = nullptr;
    partitionList    *partitions = nullptr;
    pllInstanceAttr attr;
    string alignment_file;
    string tree_file;
    string partition_file;
    bool _instance_ready   = false;
    bool _model_ready      = false;
    bool _alignment_ready  = false;
    bool _partitions_ready = false;
    bool _tree_ready       = false;
};

#endif /* PLL_WRAPPER_H_ */
