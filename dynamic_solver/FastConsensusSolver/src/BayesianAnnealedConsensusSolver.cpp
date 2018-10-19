/**
* \file BayesianAnnealedConsensusSolver.cpp
* \brief Methods for the class BayesianAnnealedConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 27/09/2018
*/

#include "BayesianAnnealedConsensusSolver.hpp"
#include <iostream>

using namespace std;

namespace soc
{//start of namespace net

/*---------------------------
 *      Constructor
 *---------------------------*/

/**
* \brief Constructor of the class BayesianAnnealedConsensusSolver
* \param[in] priority_map unordered_map of <Node,double> for speaker priority
* \param[in] influence_map unordered_map of <Node, double> to fix the influence
* \param[in] initial_state_vector vector<double> to represent initial opinion
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
BayesianAnnealedConsensusSolver::BayesianAnnealedConsensusSolver(
    unordered_map<Node,double>& priority_map,
    unordered_map<Node, double>& influence_map,
    vector<double>& initial_state_vector, double eta, int seed,
    size_t max_cluster, bool both_speak) :
        AnnealedConsensusSolver(priority_map, influence_map,
                initial_state_vector, eta, seed,
                max_cluster, both_speak),
        variance_vector_(priority_map.size(), 1.)
{
}

/**
* \brief Constructor of the class BayesianAnnealedConsensusSolver
* \param[in] priority_map unordered_map of <Node,double> for speaker priority
* \param[in] influence_map unordered_map of <Node, double> to fix the influence
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
BayesianAnnealedConsensusSolver::BayesianAnnealedConsensusSolver(
    unordered_map<Node,double>& priority_map,
    unordered_map<Node, double>& influence_map,
    double eta, int seed, size_t max_cluster, bool both_speak) :
        AnnealedConsensusSolver(priority_map, influence_map, eta, seed,
                max_cluster, both_speak),
        variance_vector_(priority_map.size(), 1.)
{
}

/*---------------------------
 *     Mutators
 *---------------------------*/
/**
* \brief Make a step of the consensus dynamics
*/
void BayesianAnnealedConsensusSolver::consensus_step()
{
    Node listener = random_int(size_, gen_);
    Node speaker = priority_distribution_(gen_);
    while (speaker == listener)
    {
        speaker = priority_distribution_(gen_);
    }
    //initial state
    double raw_variation = state_vector_[speaker]-state_vector_[listener];
    double initial_speaker_variance = variance_vector_[speaker];
    double initial_listener_variance = variance_vector_[listener];

    //get "scaled" kalman gain
    double k = (eta_*(1 + influence_map_[speaker] - influence_map_[listener])*
        initial_listener_variance/(initial_speaker_variance +
                initial_listener_variance));


    //determine variation on state and update variance
    double variation_listener = (raw_variation*k);
    state_vector_[listener] += variation_listener;
    variance_vector_[listener] = ((1-k)*(1-k)*initial_listener_variance +
            k*k*initial_speaker_variance);


    //store the variation step
    history_vector_.push_back(pair<Node,double>(listener, variation_listener));

    //change the state of the speaker also
    if (both_speak_)
    {
        //get "scaled" kalman gain
        double k = (eta_*(1 + influence_map_[listener] -
                    influence_map_[speaker])*initial_speaker_variance/
                (initial_speaker_variance + initial_listener_variance));

        //determine variation
        double variation_speaker = (-raw_variation*k);
        state_vector_[speaker] += variation_speaker;
        variance_vector_[speaker] = ((1-k)*(1-k)*initial_speaker_variance +
            k*k*initial_listener_variance);

        //store the variation step
        history_vector_.push_back(pair<Node,double>(speaker,
                    variation_speaker));
    }
}

}//end of namespace soc
