#cython: c_string_encoding=ascii  # for cython>=0.19

# pllml.cpp
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
    
    def set_optimisable_rates(self,  partition ,  optimisable ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
    
    
        self.inst.get().set_optimisable_rates((<int>partition), (<bool>optimisable))
    
    def set_optimisable_alpha(self,  partition ,  optimisable ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
        assert isinstance(optimisable, (int, long)), 'arg optimisable wrong type'
    
    
        self.inst.get().set_optimisable_alpha((<int>partition), (<bool>optimisable))
    
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
    
    def link_alpha_parameters(self, bytes linkage ):
        assert isinstance(linkage, bytes), 'arg linkage wrong type'
    
        self.inst.get().link_alpha_parameters((<libcpp_string>linkage))
    
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
    
    def set_parsimony_tree(self):
        self.inst.get().set_parsimony_tree()
    
    def link_rates(self, bytes linkage ):
        assert isinstance(linkage, bytes), 'arg linkage wrong type'
    
        self.inst.get().link_rates((<libcpp_string>linkage))
    
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
    
    def _init_0(self, bytes alignment_file , bytes partition_file , bytes tree ,  num_threads ,  rns ):
        assert isinstance(alignment_file, bytes), 'arg alignment_file wrong type'
        assert isinstance(partition_file, bytes), 'arg partition_file wrong type'
        assert isinstance(tree, bytes), 'arg tree wrong type'
        assert isinstance(num_threads, (int, long)), 'arg num_threads wrong type'
        assert isinstance(rns, (int, long)), 'arg rns wrong type'
    
    
    
    
    
        self.inst = shared_ptr[_pll](new _pll((<libcpp_string>alignment_file), (<libcpp_string>partition_file), (<libcpp_string>tree), (<int>num_threads), (<long int>rns)))
    
    def _init_1(self, bytes alignment_file , bytes partition_file ,  num_threads ,  rns ):
        assert isinstance(alignment_file, bytes), 'arg alignment_file wrong type'
        assert isinstance(partition_file, bytes), 'arg partition_file wrong type'
        assert isinstance(num_threads, (int, long)), 'arg num_threads wrong type'
        assert isinstance(rns, (int, long)), 'arg rns wrong type'
    
    
    
    
        self.inst = shared_ptr[_pll](new _pll((<libcpp_string>alignment_file), (<libcpp_string>partition_file), (<int>num_threads), (<long int>rns)))
    
    def __init__(self, *args):
        if (len(args)==5) and (isinstance(args[0], bytes)) and (isinstance(args[1], bytes)) and (isinstance(args[2], bytes)) and (isinstance(args[3], (int, long))) and (isinstance(args[4], (int, long))):
             self._init_0(*args)
        elif (len(args)==4) and (isinstance(args[0], bytes)) and (isinstance(args[1], bytes)) and (isinstance(args[2], (int, long))) and (isinstance(args[3], (int, long))):
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
    
    def is_dna(self,  partition ):
        assert isinstance(partition, (int, long)), 'arg partition wrong type'
    
        cdef bool _r = self.inst.get().is_dna((<int>partition))
        py_result = <bool>_r
        return py_result
    
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
