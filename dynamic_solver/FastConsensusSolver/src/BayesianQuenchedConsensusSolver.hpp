/**
* \file BayesianQuenchedConsensusSolver.hpp
* \brief Header file for class BayesianQuenchedConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 27/09/2018
*/

#ifndef BAYESIANQUENCHEDCONSENSUSSOLVER_HPP_
#define BAYESIANQUENCHEDCONSENSUSSOLVER_HPP_

#include <unordered_map>
#include <vector>
#include <random>
#include <utility>
#include <numeric>
#include "QuenchedConsensusSolver.hpp"

namespace soc
{//start of namespace soc

/**
* \class BayesianQuenchedConsensusSolver BayesianQuenchedConsensusSolver.hpp
* \brief Class to solve a continuous opinion formation dynamics on an annealed
* network
*/
class BayesianQuenchedConsensusSolver : public QuenchedConsensusSolver
{
public:
    BayesianQuenchedConsensusSolver(
            std::unordered_map<Node,std::vector<Node>>& network_map,
            std::unordered_map<Node,double>& influence_map,
            std::vector<double>& initial_state_vector,
            double eta = 0.5, int seed = 42, std::size_t max_cluster = 1,
            bool both_speak = false);
    BayesianQuenchedConsensusSolver(
            std::unordered_map<Node,std::vector<Node>>& network_map,
            std::unordered_map<Node,double>& influence_map,
            double eta = 0.5, int seed = 42, std::size_t max_cluster = 1,
            bool both_speak = false);

    void consensus_step();

private:
    std::vector<double> variance_vector_;
};

}//end of namespace soc

#endif /* BAYESIANQUENCHEDCONSENSUSSOLVER_HPP_ */
