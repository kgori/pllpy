# pllml.pxd
# Copyright (C) 2014  Kevin Gori

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from libcpp.string  cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
cdef extern from "pllml.h":
    cdef cppclass pll:
        # Make
        pll(libcpp_string alignment_file, libcpp_string partitions, libcpp_string tree, int num_threads, long rns) except +
        pll(libcpp_string alignment_file, libcpp_string partitions, bool parsimony, int num_threads, long rns) except +

        # Run optimisations
        void optimise(bool rates, bool freqs, bool alphas, bool branches, bool topology, int tree_search_interval, bool final_tree_search) nogil except +
        void optimise_alphas() nogil except +
        void optimise_branch_lengths(int num_iter) nogil except +
        void optimise_freqs() nogil except +
        void optimise_model() nogil except +
        void optimise_rates() nogil except +
        void optimise_tree_search(bool estimate_model) nogil except +

        # Getters
        double                               get_likelihood() nogil except +
        libcpp_vector[libcpp_string]         get_partition_names() except +
        libcpp_string                        get_partition_name(int partition) except +
        libcpp_vector[libcpp_string]         get_model_names() except +
        libcpp_string                        get_model_name(int partition) except +
        double                               get_epsilon() except +
        libcpp_vector[double]                get_alphas() except +
        double                               get_alpha(int partition) except +
        libcpp_vector[libcpp_vector[double]] get_frequencies() except +
        libcpp_vector[double]                get_frequencies_vector(int partition) except +
        libcpp_vector[libcpp_vector[double]] get_rates() except +
        libcpp_vector[double]                get_rates_vector(int partition) except +
        int                                  get_number_of_partitions() except +
        libcpp_string                        get_tree() except +
        libcpp_vector[libcpp_vector[double]] get_empirical_frequencies() except +
        # double                               get_frac_change() except +

        # Setters
        void set_epsilon(double epsilon) except +
        void set_alpha(double alpha, int partition, bool optimisable) except +
        void set_frequencies(libcpp_vector[double] freqs, int partition, bool optimisable) except +
        void set_rates(libcpp_vector[double] rates, int partition, bool optimisable) except +
        void set_optimisable_alpha(int partition, bool optimisable) except +
        void set_optimisable_frequencies(int partition, bool optimisable) except +
        void set_optimisable_rates(int partition, bool optimisable) except +
        void set_tree(libcpp_string nwk) except +

        # Partition management
        void link_alpha_parameters(libcpp_string linkage) except +
        void link_frequencies(libcpp_string linkage) except +
        void link_rates(libcpp_string linkage) except +

        # Check settings
        bool is_dna(int partition) except +
        bool is_protein(int partition) except +
        bool is_optimisable_alpha(int partition) except +
        bool is_optimisable_frequencies(int partition) except +
        bool is_optimisable_rates(int partition) except +
