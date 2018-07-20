/**
* \file ConsensusSolver.cpp
* \brief Methods for the class ConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 03/02/2018
*/

#include <ConsensusSolver.h>
#include <limits>

using namespace std;

namespace soc
{//start of namespace net

//utility function
double standard_deviation(vector<double>& v); 
unsigned int random_int(size_t size, RNGType gen);

/*---------------------------
 *      Constructor
 *---------------------------*/

/**
* \brief Constructor of the class ConsensusSolver
* \param[in] network_map unordered_map of <Node,vector<Node>> for the structure 
* \param[in] influence_map unordered_map of <Node, double> to fix the influence 
* \param[in] initial_state_vector vector<double> to represent initial opinion
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
ConsensusSolver::ConsensusSolver(unordered_map<Node,vector<Node>>& network_map,
    unordered_map<Node, double>& influence_map, 
    vector<double>& initial_state_vector, double eta, int seed) : 
        network_map_(network_map), influence_map_(influence_map),
        initial_state_vector_(initial_state_vector), eta_(eta), gen_(seed),
        history_vector_(), state_vector_(initial_state_vector) 
{
}

/**
* \brief Constructor of the class ConsensusSolver
* \param[in] network_map unordered_map of <Node,vector<Node>> for the structure 
* \param[in] influence_map unordered_map of <Node, double> to fix the influence 
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
ConsensusSolver::ConsensusSolver(unordered_map<Node,vector<Node>>& network_map,
    unordered_map<Node, double>& influence_map, double eta, int seed) : 
        network_map_(network_map), influence_map_(influence_map),
        initial_state_vector_(), eta_(eta), gen_(seed),
        history_vector_(), state_vector_() 
{
    //set the initial state at random
    for (Node node = 0; node < network_map_.size(); node++)
    {
        double opinion = generate_canonical<double, 
            numeric_limits<double>::digits>(gen_);
        initial_state_vector_.push_back(opinion);
        state_vector_.push_back(opinion);
    }
}

/*---------------------------
 *     Mutators 
 *---------------------------*/
/**
* \brief Reset the solver, but keep the initial state
*/
void ConsensusSolver::reset()
{
    state_vector_.clear();
    history_vector_.clear();
    //set the at the initial state 
    for (Node node = 0; node < network_map_.size(); node++)
    {
        state_vector_.push_back(initial_state_vector_[node]);
    }
}

/**
* \brief Reset the solver, and set the initial state at random 
*/
void ConsensusSolver::reset_all()
{
    initial_state_vector_.clear();
    state_vector_.clear();
    history_vector_.clear();
    //set the initial state at random
    for (Node node = 0; node < network_map_.size(); node++)
    {
        double opinion = generate_canonical<double, 
            numeric_limits<double>::digits>(gen_);
        initial_state_vector_.push_back(opinion);
        state_vector_.push_back(opinion);
    }
}

/**
* \brief Make a step of the consensus dynamics
*/
void ConsensusSolver::consensus_step()
{
    Node listener = random_int(network_map_.size(), gen_);
    Node speaker = network_map_.at(listener)[
            random_int(network_map_.at(listener).size(), gen_)];
    //determine variation
    double variation = ((state_vector_[speaker]-state_vector_[listener])*eta_
                       *(1+influence_map_[speaker]-influence_map_[listener]));
    state_vector_[listener] += variation;    
    //store the variation step
    history_vector_.push_back(pair<Node,double>(listener, variation));
}

/**
* \brief Make consensus step until tolerance is reached
* \param[in] tol double representing the tolerance on the standard deviation 
* of opinion
*/
void ConsensusSolver::reach_consensus(double tol)
{
    while (standard_deviation(state_vector_) > tol)
    {
        consensus_step();
    }
}

/*---------------------------
 *    Utility functions 
 *---------------------------*/

double standard_deviation(vector<double>& v)
{
    double mean = std::accumulate(v.begin(), v.end(), 0.0)/v.size();
    double variance = (std::inner_product(v.begin(), v.end(), v.begin(), 
                       0.0)/v.size() - mean*mean);
    return sqrt(variance);
}

unsigned int random_int(size_t size, RNGType& gen)
{
    return floor(generate_canonical<double,
        numeric_limits<double>::digits>(gen)*size);
}


}//end of namespace soc
