
/**
* \file AnnealedConsensusSolver.cpp
* \brief Methods for the class AnnealedConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 27/09/2018
*/

#include "AnnealedConsensusSolver.hpp"

using namespace std;

namespace soc
{//start of namespace net

/*---------------------------
 *      Constructor
 *---------------------------*/

/**
* \brief Constructor of the class AnnealedConsensusSolver
* \param[in] priority_map unordered_map of <Node,double> for speaker priority
* \param[in] influence_map unordered_map of <Node, double> to fix the influence
* \param[in] initial_state_vector vector<double> to represent initial opinion
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
AnnealedConsensusSolver::AnnealedConsensusSolver(
    unordered_map<Node,double>& priority_map,
    unordered_map<Node, double>& influence_map,
    vector<double>& initial_state_vector, double eta, int seed,
    size_t max_cluster) :
        ConsensusSolver(influence_map, initial_state_vector, eta, seed,
                max_cluster), priority_distribution_()
{
    //initialize priority distribution
    vector<double> temp(size_, 0.);
    for (Node i = 0; i < size_; i++)
    {
        temp[i] = priority_map[i];
    }
    priority_distribution_ = discrete_distribution<Node>(temp.begin(),
            temp.end());
}

/**
* \brief Constructor of the class AnnealedConsensusSolver
* \param[in] priority_map unordered_map of <Node,double> for speaker priority
* \param[in] influence_map unordered_map of <Node, double> to fix the influence
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
AnnealedConsensusSolver::AnnealedConsensusSolver(
    unordered_map<Node,double>& priority_map,
    unordered_map<Node, double>& influence_map,
    double eta, int seed, size_t max_cluster) :
        ConsensusSolver(influence_map, eta, seed, max_cluster),
        priority_distribution_()
{
    //initialize priority distribution
    vector<double> temp(size_, 0.);
    for (Node i = 0; i < size_; i++)
    {
        temp[i] = priority_map[i];
    }
    priority_distribution_ = discrete_distribution<Node>(temp.begin(),
            temp.end());
}

/*---------------------------
 *     Mutators
 *---------------------------*/
/**
* \brief Make a step of the consensus dynamics
*/
void AnnealedConsensusSolver::consensus_step()
{
    Node listener = random_int(size_, gen_);
    Node speaker = priority_distribution_(gen_);
    while (speaker == listener)
    {
        speaker = priority_distribution_(gen_);
    }
    //determine variation
    double variation = ((state_vector_[speaker]-state_vector_[listener])*eta_
                       *(1+influence_map_[speaker]-influence_map_[listener]));
    state_vector_[listener] += variation;
    //store the variation step
    history_vector_.push_back(pair<Node,double>(listener, variation));
}

}//end of namespace soc
