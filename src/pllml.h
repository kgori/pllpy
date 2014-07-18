/*
 * pllml.h
 *
 *  Created on: Jul 17, 2014
 *      Author: kgori
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
	pll();
	pll(string alignment_file, string tree, string partition_file, int num_threads = 1);
	virtual ~pll();

	// Model data IO
	void                   load_alignment_file(string path);
    void                   load_tree_file(string path);
    void                   load_partition_file(string path);
    void                   load_tree_string(string nwk);

    // Model management
    void                   create_instance();
    void                   initialise_model();
    void                   destroy_model();

    // Run optimisations
	void                   optimise(bool estimate_model=true);
	void                   optimise_branch_lengths(int num_iter=32);
	void                   optimise_model();

	// Collect results
	double                 get_likelihood();
	vector<string>         get_partition_names();
	string                 get_partition_name(int partition);
	vector<string>		   get_model_names();
	string				   get_model_name(int partition);
	vector<double>         get_alphas();
	double                 get_alpha(int partition);
	vector<vector<double>> get_frequencies();
	vector<double>         get_frequencies_vector(int partition);
	vector<vector<double>> get_rates();
	vector<double>         get_rates_vector(int partition);
	double                 get_epsilon();
	int                    get_number_of_partitions();
	int                    get_number_of_threads();
	string                 get_tree();
	vector<vector<double>> get_empirical_frequencies();

	// Set parameters
	void                   set_epsilon(double epsilon);
	void                   set_alpha(double alpha, int partition, bool optimisable);
	void                   set_frequencies(vector<double> freqs, int partition, bool optimisable);
	void                   set_rates(vector<double> rates, int partition, bool optimisable);
	void                   set_optimisable_alpha(int partition, bool optimisable);
	void                   set_optimisable_frequencies(int partition, bool optimisable);
	void                   set_optimisable_rates(int partition, bool optimisable);
	void                   set_number_of_threads(int threads);

	// Parameter management
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
	void                   _evaluate_likelihood();
	double                 _vector_sum(vector<double>);
	bool                   _approx_eq(double a, double b);
	void                   _partitions_bounds_check(int partition);
	bool                   _is_file(string filename);
	bool                   _is_tree_string(string tree_string);
	string				   _model_name(int model_num);

	pllAlignmentData * alignment_data;
	pllInstance * tr;
	pllNewickTree * newick;
	partitionList * partitions;
	pllInstanceAttr attr;
	string alignment_file;
	string tree_file;
	string partition_file;
	bool _is_ready = false;
};

#endif /* PLL_WRAPPER_H_ */
