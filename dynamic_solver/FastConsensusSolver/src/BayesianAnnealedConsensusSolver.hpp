/**
* \file BayesianAnnealedConsensusSolver.hpp
* \brief Header file for class BayesianAnnealedConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 27/09/2018
*/

#ifndef BAYESIANANNEALEDCONSENSUSSOLVER_HPP_
#define BAYESIANANNEALEDCONSENSUSSOLVER_HPP_

#include <unordered_map>
#include <vector>
#include <random>
#include <utility>
#include <numeric>
#include "AnnealedConsensusSolver.hpp"

namespace soc
{//start of namespace soc

/**
* \class BayesianAnnealedConsensusSolver BayesianAnnealedConsensusSolver.hpp
* \brief Class to solve a continuous opinion formation dynamics on an annealed
* network
*/
class BayesianAnnealedConsensusSolver : public AnnealedConsensusSolver
{
public:
    BayesianAnnealedConsensusSolver(
            std::unordered_map<Node,double>& priority_map,
            std::unordered_map<Node,double>& influence_map,
            std::vector<double>& initial_state_vector,
            double eta = 0.5, int seed = 42, std::size_t max_cluster = 1,
            bool both_speak = false);
    BayesianAnnealedConsensusSolver(
            std::unordered_map<Node,double>& priority_map,
            std::unordered_map<Node,double>& influence_map,
            double eta = 0.5, int seed = 42, std::size_t max_cluster = 1,
            bool both_speak = false);

    void consensus_step();

private:
    std::vector<double> variance_vector_;
};

}//end of namespace soc

#endif /* BAYESIANANNEALEDCONSENSUSSOLVER_HPP_ */
