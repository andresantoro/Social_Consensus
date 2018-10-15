/**
* \file QuenchedConsensusSolver.cpp
* \brief Methods for the class QuenchedConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 27/09/2018
*/

#include "QuenchedConsensusSolver.hpp"

using namespace std;

namespace soc
{//start of namespace net

/*---------------------------
 *      Constructor
 *---------------------------*/

/**
* \brief Constructor of the class QuenchedConsensusSolver
* \param[in] network_map unordered_map of <Node, vector<Node>> for structure
* \param[in] influence_map unordered_map of <Node, double> to fix the influence
* \param[in] initial_state_vector vector<double> to represent initial opinion
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
QuenchedConsensusSolver::QuenchedConsensusSolver(
    unordered_map<Node, vector<Node> >& network_map,
    unordered_map<Node, double>& influence_map,
    vector<double>& initial_state_vector, double eta, int seed,
    size_t max_cluster) :
        ConsensusSolver(influence_map, initial_state_vector, eta, seed,
                max_cluster), network_map_(network_map)
{
}

/**
* \brief Constructor of the class QuenchedConsensusSolver
* \param[in] network_map unordered_map of <Node, vector<Node>> for structure
* \param[in] influence_map unordered_map of <Node, double> to fix the influence
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
QuenchedConsensusSolver::QuenchedConsensusSolver(
    unordered_map<Node, vector<Node> >& network_map,
    unordered_map<Node, double>& influence_map,
    double eta, int seed, size_t max_cluster) :
        ConsensusSolver(influence_map, eta, seed, max_cluster),
        network_map_(network_map)
{
}

/*---------------------------
 *     Mutators
 *---------------------------*/
/**
* \brief Make a step of the consensus dynamics
*/
void QuenchedConsensusSolver::consensus_step()
{
    Node listener = random_int(size_, gen_);
    Node speaker = network_map_.at(listener)[
            random_int(network_map_.at(listener).size(), gen_)];
    //determine variation
    double variation = ((state_vector_[speaker]-state_vector_[listener])*eta_
                       *(1+influence_map_[speaker]-influence_map_[listener]));
    state_vector_[listener] += variation;
    //store the variation step
    history_vector_.push_back(pair<Node,double>(listener, variation));
}

}//end of namespace soc
