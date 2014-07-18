from libcpp.string  cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
cdef extern from "pllml.h":
    cdef cppclass pll:
        # Make
        pll() except +
        pll(libcpp_string alignment_file, libcpp_string tree, libcpp_string partition_file, int num_threads) except +

        # Initialise model
        void load_alignment_file(libcpp_string path) except +
        void load_tree_file(libcpp_string path) except +
        void load_partition_file(libcpp_string path) except +
        void load_tree_string(libcpp_string nwk) except +

        # Model management
        void create_instance()
        void initialise_model() except +
        void destroy_model()

        # Run analysis
        void optimise(bool estimate_model) nogil
        void optimise_branch_lengths(int num_iter) nogil
        void optimise_model() nogil

        # Collect results
        double get_likelihood() nogil
        libcpp_vector[libcpp_string] get_model_names() except +
        libcpp_string get_model_name(int partition) except +
        libcpp_vector[libcpp_string] get_partition_names() except +
        libcpp_string get_partition_name(int partition) except +
        libcpp_vector[double] get_alphas() except +
        double get_alpha(int partition) except +
        libcpp_vector[libcpp_vector[double]] get_frequencies() except +
        libcpp_vector[double] get_frequencies_vector(int partition) except +
        libcpp_vector[libcpp_vector[double]] get_rates() except +
        libcpp_vector[double] get_rates_vector(int partition) except +
        double get_epsilon()
        int get_number_of_partitions()
        int get_number_of_threads()
        libcpp_string get_tree()
        libcpp_vector[libcpp_vector[double]] get_empirical_frequencies() except +

        # Set parameters
        void set_epsilon(double epsilon)
        void set_alpha(double alpha, int partition, bool optimisable) except +
        void set_frequencies(libcpp_vector[double] freqs, int partition, bool optimisable) except +
        void set_rates(libcpp_vector[double] rates, int partition, bool optimisable) except +
        void set_optimisable_alpha(int partition, bool optimisable) except +
        void set_optimisable_frequencies(int partition, bool optimisable) except +
        void set_optimisable_rates(int partition, bool optimisable) except +
        void set_number_of_threads(int threads)

        # Parameter management
        void link_alpha_parameters(libcpp_string linkage)
        void link_frequencies(libcpp_string linkage)
        void link_rates(libcpp_string linkage)

        # Check settings
        bool is_dna(int partition) except +
        bool is_protein(int partition) except +
        bool is_optimisable_alpha(int partition) except +
        bool is_optimisable_frequencies(int partition) except +
        bool is_optimisable_rates(int partition) except +
