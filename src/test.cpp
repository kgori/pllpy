#include <iostream>
#include <string>
#include <vector>
#include "pllml.h"
extern "C" {
#include <pll/pll.h>
}

int main(int argc, char const *argv[])
{
    std::string FILE="/homes/kgori/scratch/basic_pll/GGPS1.phy";
    std::string PART="/homes/kgori/scratch/basic_pll/GGPS1.partitions.txt";
    std::string TREE="/homes/kgori/scratch/basic_pll/GGPS1.tree";

    auto inst = make_unique<pll>(FILE, PART, TREE, 1, 123);
    std::cout << inst->get_likelihood() << std::endl;
    std::cout << inst->get_partition_name(0) << std::endl;
    std::cout << inst->get_model_name(0) << std::endl;
    std::cout << inst->get_epsilon() << std::endl;
    std::cout << inst->get_alpha(0) << std::endl;
    std::cout << inst->get_number_of_partitions() << std::endl;
    std::cout << inst->get_tree() << std::endl;
    std::cout << inst->get_frac_change() << std::endl;
    inst->set_tree("((Ptro:0.00147084048849247711,Ppan:0.00106294370976534763):0.00444598321816729036,(Cjac:0.05157798212603657839,(Csab:0.00440006365327790666,(Mmul:0.00652936575529242547,Panu:0.00101047194512476272):0.00381049900890796569):0.01359968787254225639):0.01375788728353777995,Hsap:0.00674547067186638278):0.0;");
    inst->set_epsilon(0.001);
    inst->set_alpha(2, 0, true);
    auto ef = inst->get_empirical_frequencies();

    inst->optimise(true, true, true, true);
    inst->optimise_tree_search(true);
    std::cout << inst->get_likelihood() << std::endl;
    std::cout << inst->get_tree() << std::endl;
    return 0;
}
