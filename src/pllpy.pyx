#cython: c_string_type=str, c_string_encoding=ascii
from  libcpp.string  cimport string as libcpp_string
from  libcpp.set     cimport set as libcpp_set
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair as libcpp_pair
from  libcpp.map     cimport map  as libcpp_map
from  libcpp cimport bool
from  libc.string cimport const_char
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from pllml cimport pll as _pll


cdef class pll:

    cdef _pll *inst

    def __dealloc__(self):
         del self.inst

    def get_model_names(self):
        _r = self.inst.get_model_names()
        cdef list py_result = _r
        return py_result
    
    def set_optimisable_rates(self,  partition ,  optimisable ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
        self.inst.set_optimisable_rates((<int>partition), (<bool>optimisable))
    
    def optimise_freqs(self):
        self.inst.optimise_freqs()
    
    def optimise_alphas(self):
        self.inst.optimise_alphas()
    
    def set_optimisable_alpha(self,  partition ,  optimisable ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
        self.inst.set_optimisable_alpha((<int>partition), (<bool>optimisable))
    
    def get_rates(self):
        _r = self.inst.get_rates()
        cdef list py_result = _r
        return py_result
    
    def get_rates_vector(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        _r = self.inst.get_rates_vector((<int>partition))
        cdef list py_result = _r
        return py_result
    
    def is_optimisable_frequencies(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        cdef bool _r = self.inst.is_optimisable_frequencies((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def optimise_model(self):
        self.inst.optimise_model()
    
    def link_alpha_parameters(self, str linkage ):
        assert isinstance(linkage, str), 'arg linkage wrong type'
        self.inst.link_alpha_parameters((<libcpp_string>linkage))
    
    def set_rates(self, list rates ,  partition ,  optimisable ):
        assert isinstance(rates, list) and all(isinstance(elemt_rec, float) for elemt_rec in rates), 'arg rates wrong type'
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
        cdef libcpp_vector[double] v0 = rates
        self.inst.set_rates(v0, (<int>partition), (<bool>optimisable))
    
    def set_tree(self, str nwk):
        assert isinstance(nwk, str), 'arg nwk wrong type'
        self.inst.set_tree((<libcpp_string>nwk))

    def get_model_name(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        cdef libcpp_string _r = self.inst.get_model_name((<int>partition))
        py_result = <libcpp_string>_r
        return py_result
    
    def optimise_branch_lengths(self,  num_iter ):
        assert isinstance(num_iter, (int, long)), 'arg num_iter wrong type'
        self.inst.optimise_branch_lengths((<int>num_iter))
    
    def is_optimisable_rates(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        cdef bool _r = self.inst.is_optimisable_rates((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def set_epsilon(self, double epsilon ):
        assert isinstance(epsilon, float), 'arg epsilon wrong type'
        self.inst.set_epsilon((<double>epsilon))
    
    def link_rates(self, str linkage ):
        assert isinstance(linkage, str), 'arg linkage wrong type'
        self.inst.link_rates((<libcpp_string>linkage))
    
    def optimise_rates(self):
        self.inst.optimise_rates()
    
    def get_alphas(self):
        _r = self.inst.get_alphas()
        cdef list py_result = _r
        return py_result
    
    def set_frequencies(self, list freqs ,  partition ,  optimisable ):
        assert isinstance(freqs, list) and all(isinstance(elemt_rec, float) for elemt_rec in freqs), 'arg freqs wrong type'
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
        cdef libcpp_vector[double] v0 = freqs
        self.inst.set_frequencies(v0, (<int>partition), (<bool>optimisable))

    def is_protein(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        cdef bool _r = self.inst.is_protein((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def set_optimisable_frequencies(self,  partition,  optimisable):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
        self.inst.set_optimisable_frequencies((<int>partition), (<bool>optimisable))
    
    def is_optimisable_alpha(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        cdef bool _r = self.inst.is_optimisable_alpha((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def get_tree(self):
        cdef libcpp_string _r = self.inst.get_tree()
        py_result = <libcpp_string>_r
        return py_result
    
    def get_frequencies_vector(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        _r = self.inst.get_frequencies_vector((<int>partition))
        cdef list py_result = _r
        return py_result
    
    def get_number_of_partitions(self):
        cdef int _r = self.inst.get_number_of_partitions()
        py_result = <int>_r
        return py_result
    
    def get_epsilon(self):
        cdef double _r = self.inst.get_epsilon()
        py_result = <double>_r
        return py_result
    
    def _init_0(self, str alignment_file , str partitions , str tree ,  num_threads ,  rns ):
        assert isinstance(alignment_file, str), 'arg alignment_file wrong type'
        assert isinstance(partitions, str), 'arg partitions wrong type'
        assert isinstance(tree, str), 'arg tree wrong type'
        assert isinstance(num_threads, (int, long)), 'arg num_threads wrong type'
        assert isinstance(rns, (int, long)), 'arg rns wrong type'
        self.inst = new _pll((<libcpp_string>alignment_file), (<libcpp_string>partitions), (<libcpp_string>tree), (<int>num_threads), (<long int>rns))
    
    def _init_1(self, str alignment_file , str partitions ,  parsimony ,  num_threads ,  rns ):
        assert isinstance(alignment_file, str), 'arg alignment_file wrong type'
        assert isinstance(partitions, str), 'arg partitions wrong type'
        assert isinstance(parsimony, (int, long)), 'arg parsimony wrong type'
        assert isinstance(num_threads, (int, long)), 'arg num_threads wrong type'
        assert isinstance(rns, (int, long)), 'arg rns wrong type'
        self.inst = new _pll((<libcpp_string>alignment_file), (<libcpp_string>partitions), (<bool>parsimony), (<int>num_threads), (<long int>rns))
    
    def __init__(self, *args):
        if (len(args)==5) and (isinstance(args[0], str)) and (isinstance(args[1], str)) and (isinstance(args[2], str)) and (isinstance(args[3], (int, long))) and (isinstance(args[4], (int, long))):
            self._init_0(*args)
        elif (len(args)==5) and (isinstance(args[0], str)) and (isinstance(args[1], str)) and (isinstance(args[2], (int, long))) and (isinstance(args[3], (int, long))) and (isinstance(args[4], (int, long))):
            self._init_1(*args)
        else:
            # raise Exception('can not handle type of %s' % (args,))
            raise Exception("Received wrong argument types:\n"
                            "expected <class 'str'> <class 'str'> <class 'bool|int'> <class 'int'> <class 'int'>\n"
                            "got {} {} {} {} {}".format(*[type(arg) for arg in args]))
    
    def get_alpha(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef double _r = self.inst.get_alpha((<int>partition))
        py_result = <double>_r
        return py_result
    
    def get_partition_names(self):
        _r = self.inst.get_partition_names()
        cdef list py_result = _r
        return py_result
    
    def is_dna(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef bool _r = self.inst.is_dna((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def optimise(self,  rates=True,  freqs=True,  alphas=True,  branches=True,
                 topology=True, tree_search_interval=5, final_tree_search=True):
        """

        :param rates: (bool) Sets optimisation of rates (i.e. Q-matrix parameters) on
        or off. Default=on.
        :param freqs: (bool) Sets optimisation of equilibrium base frequencies on or
        off. Default=on.
        :param alphas: (bool) Sets optimisation of alpha parameter of gamma-distributed
        rate variation on or off. Default=on.
        :param branches: (bool) Sets optimisation of branch lengths on or off. Default=on.
        :param topology: (bool) Sets optimisation of tree topology on or off. Topology
        is updated every few iterations, according to the tree_search_interval. Default=on.
        :param tree_search_interval: (int) Number of model optimisation iterations
        between tree searches. Default=5.
        :param final_tree_search: (bool) Sets whether to do a final optimisation of
        tree topology after the algorithm converges. Also does one final optimisation
        of the activated model parameters. Default=on.
        """
        assert isinstance(rates, (int, long)), 'arg rates wrong type'
        assert isinstance(freqs, (int, long)), 'arg freqs wrong type'
        assert isinstance(alphas, (int, long)), 'arg alphas wrong type'
        assert isinstance(branches, (int, long)), 'arg branches wrong type'
        assert isinstance(topology, (int, long)), 'arg topology wrong type'
        assert isinstance(tree_search_interval, (int, long)), 'arg tree_search_interval wrong type'
        assert isinstance(final_tree_search, (int, long)), 'arg final_tree_search wrong type'
        self.inst.optimise((<bool>rates), (<bool>freqs), (<bool>alphas), (<bool>branches),
                           (<bool>topology), (<int>tree_search_interval), (<bool>final_tree_search))
    
    def set_alpha(self, double alpha,  partition,  optimisable):
        assert isinstance(alpha, float), 'arg alpha wrong type'
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
        self.inst.set_alpha((<double>alpha), (<int>partition), (<bool>optimisable))
    
    def get_empirical_frequencies(self):
        _r = self.inst.get_empirical_frequencies()
        cdef list py_result = _r
        return py_result
    
    def optimise_tree_search(self,  estimate_model ):
        assert isinstance(estimate_model, (int, long)), 'arg estimate_model wrong type'
        self.inst.optimise_tree_search((<bool>estimate_model))
    
    def get_partition_name(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef libcpp_string _r = self.inst.get_partition_name((<int>partition))
        py_result = <libcpp_string>_r
        return py_result
    
    # def get_frac_change(self):
    #     cdef double _r = self.inst.get_frac_change()
    #     py_result = <double>_r
    #     return py_result
    
    def link_frequencies(self, str linkage ):
        assert isinstance(linkage, str), 'arg linkage wrong type'
        self.inst.link_frequencies((<libcpp_string>linkage))
    
    def get_likelihood(self):
        cdef double _r = self.inst.get_likelihood()
        py_result = <double>_r
        return py_result
    
    def get_frequencies(self):
        _r = self.inst.get_frequencies()
        cdef list py_result = _r
        return py_result 
