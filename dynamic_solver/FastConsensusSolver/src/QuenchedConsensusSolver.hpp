/**
* \file QuenchedConsensusSolver.hpp
* \brief Header file for class QuenchedConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 27/09/2018
*/

#ifndef QUENCHEDCONSENSUSSOLVER_HPP_
#define QUENCHEDCONSENSUSSOLVER_HPP_

#include <unordered_map>
#include <vector>
#include <random>
#include <utility>
#include <numeric>
#include "ConsensusSolver.hpp"

namespace soc
{//start of namespace soc

/**
* \class QuenchedConsensusSolver QuenchedConsensusSolver.hpp
* \brief Class to solve a continuous opinion formation dynamics on a quenched
* network
*/
class QuenchedConsensusSolver : public ConsensusSolver
{
public:
    QuenchedConsensusSolver(
            std::unordered_map<Node,std::vector<Node>>& network_map,
            std::unordered_map<Node,double>& influence_map,
            std::vector<double>& initial_state_vector,
            double eta = 0.5, int seed = 42);
    QuenchedConsensusSolver(
            std::unordered_map<Node,std::vector<Node>>& network_map,
            std::unordered_map<Node,double>& influence_map,
            double eta = 0.5, int seed = 42);

    void consensus_step();

private:
    std::unordered_map<Node,std::vector<Node>> network_map_;

};

}//end of namespace soc

#endif /* QUENCHEDCONSENSUSSOLVER_HPP_ */
