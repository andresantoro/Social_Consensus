/**
* \file AnnealedConsensusSolver.hpp
* \brief Header file for class AnnealedConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 27/09/2018
*/

#ifndef ANNEALEDCONSENSUSSOLVER_HPP_
#define ANNEALEDCONSENSUSSOLVER_HPP_

#include <unordered_map>
#include <vector>
#include <random>
#include <utility>
#include <numeric>
#include "ConsensusSolver.hpp"

namespace soc
{//start of namespace soc

/**
* \class AnnealedConsensusSolver AnnealedConsensusSolver.hpp
* \brief Class to solve a continuous opinion formation dynamics on an annealed
* network
*/
class AnnealedConsensusSolver : public ConsensusSolver
{
public:
    AnnealedConsensusSolver(
            std::unordered_map<Node,double>& priority_map,
            std::unordered_map<Node,double>& influence_map,
            std::vector<double>& initial_state_vector,
            double eta = 0.5, int seed = 42);
    AnnealedConsensusSolver(
            std::unordered_map<Node,double>& priority_map,
            std::unordered_map<Node,double>& influence_map,
            double eta = 0.5, int seed = 42);

    void consensus_step();

private:
    std::discrete_distribution<Node> priority_distribution_;
};

}//end of namespace soc

#endif /* ANNEALEDCONSENSUSSOLVER_HPP_ */
