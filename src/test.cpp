#include <iostream>
#include <string>
#include "pllml.h"
extern "C" {
#include <pll/pll.h>
}

int main(int argc, char const *argv[])
{
    std::string FILE="/Users/kgori/scratch/basic_pll/GGPS1.phy";
    std::string PART="/Users/kgori/scratch/basic_pll/GGPS1.partitions.txt";
    std::string TREE="/Users/kgori/scratch/basic_pll/GGPS1.tree";

    auto inst = pll(FILE, PART, TREE, 1, 123);
    std::cout << inst.get_likelihood() << std::endl;
    auto ef = inst.get_empirical_frequencies();
    inst.optimise(true, true, true, true);
    inst.optimise_tree_search(true);
    std::cout << inst.get_likelihood() << std::endl;
    std::cout << inst.get_tree() << std::endl;
    return 0;
}
