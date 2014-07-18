#cython: c_string_encoding=ascii  # for cython>=0.19
from  libcpp.string  cimport string as libcpp_string
from  libcpp.set     cimport set as libcpp_set
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair as libcpp_pair
from  libcpp.map     cimport map  as libcpp_map
from  smart_ptr cimport shared_ptr
from  AutowrapRefHolder cimport AutowrapRefHolder
from  libcpp cimport bool
from  libc.string cimport const_char
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from pllml cimport pll as _pll
cdef extern from "autowrap_tools.hpp":
    char * _cast_const_away(char *) 

cdef class pll:

    cdef shared_ptr[_pll] inst

    def __dealloc__(self):
         self.inst.reset()

    
    def get_model_names(self):
        _r = self.inst.get().get_model_names()
        cdef list py_result = _r
        return py_result
    
    def create_instance(self):
        self.inst.get().create_instance()
    
    def destroy_model(self):
        self.inst.get().destroy_model()
    
    def set_optimisable_alpha(self,  partition ,  optimisable ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
    
    
        self.inst.get().set_optimisable_alpha((<int>partition), (<bool>optimisable))
    
    def load_alignment_file(self, bytes path ):
        assert isinstance(path, bytes), 'arg path wrong type'
    
        self.inst.get().load_alignment_file((<libcpp_string>path))
    
    def get_rates(self):
        _r = self.inst.get().get_rates()
        cdef list py_result = _r
        return py_result
    
    def get_rates_vector(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        _r = self.inst.get().get_rates_vector((<int>partition))
        cdef list py_result = _r
        return py_result
    
    def is_optimisable_frequencies(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef bool _r = self.inst.get().is_optimisable_frequencies((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def optimise_model(self):
        self.inst.get().optimise_model()
    
    def load_partition_file(self, bytes path ):
        assert isinstance(path, bytes), 'arg path wrong type'
    
        self.inst.get().load_partition_file((<libcpp_string>path))
    
    def link_alpha_parameters(self, bytes linkage ):
        assert isinstance(linkage, bytes), 'arg linkage wrong type'
    
        self.inst.get().link_alpha_parameters((<libcpp_string>linkage))
    
    def set_optimisable_rates(self,  partition ,  optimisable ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
    
    
        self.inst.get().set_optimisable_rates((<int>partition), (<bool>optimisable))
    
    def get_number_of_threads(self):
        cdef int _r = self.inst.get().get_number_of_threads()
        py_result = <int>_r
        return py_result
    
    def set_rates(self, list rates ,  partition ,  optimisable ):
        assert isinstance(rates, list) and all(isinstance(elemt_rec, float) for elemt_rec in rates), 'arg rates wrong type'
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
        cdef libcpp_vector[double] v0 = rates
    
    
        self.inst.get().set_rates(v0, (<int>partition), (<bool>optimisable))
        
    
    def get_model_name(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef libcpp_string _r = self.inst.get().get_model_name((<int>partition))
        py_result = <libcpp_string>_r
        return py_result
    
    def set_number_of_threads(self,  threads ):
        assert isinstance(threads, (int, long)), 'arg threads wrong type'
    
        self.inst.get().set_number_of_threads((<int>threads))
    
    def optimise_branch_lengths(self,  num_iter ):
        assert isinstance(num_iter, (int, long)), 'arg num_iter wrong type'
    
        self.inst.get().optimise_branch_lengths((<int>num_iter))
    
    def is_optimisable_rates(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef bool _r = self.inst.get().is_optimisable_rates((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def link_rates(self, bytes linkage ):
        assert isinstance(linkage, bytes), 'arg linkage wrong type'
    
        self.inst.get().link_rates((<libcpp_string>linkage))
    
    def initialise_model(self):
        self.inst.get().initialise_model()
    
    def get_alphas(self):
        _r = self.inst.get().get_alphas()
        cdef list py_result = _r
        return py_result
    
    def set_frequencies(self, list freqs ,  partition ,  optimisable ):
        assert isinstance(freqs, list) and all(isinstance(elemt_rec, float) for elemt_rec in freqs), 'arg freqs wrong type'
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
        cdef libcpp_vector[double] v0 = freqs
    
    
        self.inst.get().set_frequencies(v0, (<int>partition), (<bool>optimisable))
        
    
    def is_protein(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef bool _r = self.inst.get().is_protein((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def set_optimisable_frequencies(self,  partition ,  optimisable ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
    
    
        self.inst.get().set_optimisable_frequencies((<int>partition), (<bool>optimisable))
    
    def is_optimisable_alpha(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef bool _r = self.inst.get().is_optimisable_alpha((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def get_tree(self):
        cdef libcpp_string _r = self.inst.get().get_tree()
        py_result = <libcpp_string>_r
        return py_result
    
    def get_frequencies_vector(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        _r = self.inst.get().get_frequencies_vector((<int>partition))
        cdef list py_result = _r
        return py_result
    
    def get_number_of_partitions(self):
        cdef int _r = self.inst.get().get_number_of_partitions()
        py_result = <int>_r
        return py_result
    
    def get_epsilon(self):
        cdef double _r = self.inst.get().get_epsilon()
        py_result = <double>_r
        return py_result
    
    def load_tree_file(self, bytes path ):
        assert isinstance(path, bytes), 'arg path wrong type'
    
        self.inst.get().load_tree_file((<libcpp_string>path))
    
    def _init_0(self):
        self.inst = shared_ptr[_pll](new _pll())
    
    def _init_1(self, bytes alignment_file , bytes tree , bytes partition_file ,  num_threads ):
        assert isinstance(alignment_file, bytes), 'arg alignment_file wrong type'
        assert isinstance(tree, bytes), 'arg tree wrong type'
        assert isinstance(partition_file, bytes), 'arg partition_file wrong type'
        assert isinstance(num_threads, (int, long)), 'arg num_threads wrong type'
    
    
    
    
        self.inst = shared_ptr[_pll](new _pll((<libcpp_string>alignment_file), (<libcpp_string>tree), (<libcpp_string>partition_file), (<int>num_threads)))
    
    def __init__(self, *args):
        if not args:
             self._init_0(*args)
        elif (len(args)==4) and (isinstance(args[0], bytes)) and (isinstance(args[1], bytes)) and (isinstance(args[2], bytes)) and (isinstance(args[3], (int, long))):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def get_alpha(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef double _r = self.inst.get().get_alpha((<int>partition))
        py_result = <double>_r
        return py_result
    
    def get_partition_names(self):
        _r = self.inst.get().get_partition_names()
        cdef list py_result = _r
        return py_result
    
    def load_tree_string(self, bytes nwk ):
        assert isinstance(nwk, bytes), 'arg nwk wrong type'
    
        self.inst.get().load_tree_string((<libcpp_string>nwk))
    
    def optimise(self,  estimate_model ):
        assert isinstance(estimate_model, (int, long)), 'arg estimate_model wrong type'
    
        self.inst.get().optimise((<bool>estimate_model))
    
    def set_alpha(self, double alpha ,  partition ,  optimisable ):
        assert isinstance(alpha, float), 'arg alpha wrong type'
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
    
    
    
        self.inst.get().set_alpha((<double>alpha), (<int>partition), (<bool>optimisable))
    
    def get_empirical_frequencies(self):
        _r = self.inst.get().get_empirical_frequencies()
        cdef list py_result = _r
        return py_result
    
    def set_epsilon(self, double epsilon ):
        assert isinstance(epsilon, float), 'arg epsilon wrong type'
    
        self.inst.get().set_epsilon((<double>epsilon))
    
    def get_partition_name(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef libcpp_string _r = self.inst.get().get_partition_name((<int>partition))
        py_result = <libcpp_string>_r
        return py_result
    
    def is_dna(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef bool _r = self.inst.get().is_dna((<int>partition))
        py_result = <bool>_r
        return py_result
    
    def link_frequencies(self, bytes linkage ):
        assert isinstance(linkage, bytes), 'arg linkage wrong type'
    
        self.inst.get().link_frequencies((<libcpp_string>linkage))
    
    def get_likelihood(self):
        cdef double _r = self.inst.get().get_likelihood()
        py_result = <double>_r
        return py_result
    
    def get_frequencies(self):
        _r = self.inst.get().get_frequencies()
        cdef list py_result = _r
        return py_result 
