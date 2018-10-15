/**
* \file ConsensusSolver.cpp
* \brief Methods for the class ConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 03/02/2018
*/

#include "ConsensusSolver.hpp"

using namespace std;

namespace soc
{//start of namespace net

/*---------------------------
 *      Constructor
 *---------------------------*/

/**
* \brief Constructor of the class ConsensusSolver
* \param[in] influence_map unordered_map of <Node, double> to fix the influence
* \param[in] initial_state_vector vector<double> to represent initial opinion
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
ConsensusSolver::ConsensusSolver(
    unordered_map<Node, double>& influence_map,
    vector<double>& initial_state_vector, double eta, int seed,
    size_t max_cluster) :
        size_(influence_map.size()), influence_map_(influence_map),
        initial_state_vector_(initial_state_vector), eta_(eta), gen_(seed),
        history_vector_(), state_vector_(initial_state_vector),
        cluster_vector_(max_cluster), cluster_mean_vector_(max_cluster),
        cluster_std_vector_(max_cluster), cluster_label_vector_(size_, 0)
{
}

/**
* \brief Constructor of the class ConsensusSolver
* \param[in] influence_map unordered_map of <Node, double> to fix the influence
* \param[in] eta double parameter for the dynamics
* \param[in] seed int seed for the random number generator
*/
ConsensusSolver::ConsensusSolver(
    unordered_map<Node, double>& influence_map, double eta, int seed,
    size_t max_cluster) :
        size_(influence_map.size()), influence_map_(influence_map),
        initial_state_vector_(), eta_(eta), gen_(seed),
        history_vector_(), state_vector_(), cluster_vector_(max_cluster),
        cluster_mean_vector_(max_cluster),
        cluster_std_vector_(max_cluster),
        cluster_label_vector_(size_, 0)
{
    //set the initial state at random
    for (Node node = 0; node < size_; node++)
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
    for (Node node = 0; node < size_; node++)
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
    for (Node node = 0; node < size_; node++)
    {
        double opinion = generate_canonical<double,
            numeric_limits<double>::digits>(gen_);
        initial_state_vector_.push_back(opinion);
        state_vector_.push_back(opinion);
    }
}


/**
* \brief Make consensus step until tolerance is reached
* \param[in] tol double representing the tolerance on the standard deviation
* of opinion
*/
void ConsensusSolver::reach_consensus(double tol)
{
    //let us perform the same number of step as the number of agents
    //it is unecessary to look at the standard deviation at each step
    for (int i = 0; i < size_ ; i++)
    {
        consensus_step();
    }

    //refine the state until convergence
    if (cluster_vector_.size() > 1)
    {
        //initialize clusters mean randomly according to data
        initialize_cluster_mean(cluster_vector_.size());
        //initialize the assignements to the clusters
        for (size_t node = 0; node < size_; node++)
        {
            assign_to_cluster(node);
        }
        //get an appropriate separation in clusters
        converge_cluster();

        //iterate until convergence
        while (*max_element(cluster_std_vector_.begin(),
                    cluster_std_vector_.end()) > tol)
        {
            consensus_step();
            converge_cluster();
        }
    }
    else
    {
        //there is only one cluster
        while (standard_deviation(state_vector_) > tol)
        {
            consensus_step();
        }
    }
}

/*----------------------------
 *    Utility methods
 *----------------------------*/
//Initialize the cluster mean at random using the data
void ConsensusSolver::initialize_cluster_mean(size_t max_cluster)
{
    for (int cluster; cluster < max_cluster; cluster ++)
    {
        cluster_mean_vector_[cluster] = state_vector_[random_int(size_, gen_)];
    }
}

//Assign the node to the appropriate cluster; return true if it changed
bool ConsensusSolver::assign_to_cluster(size_t node)
{
    //search for the appropriate cluster with minimal deviation
    size_t appropriate_cluster = 0;
    double min_distance = abs(state_vector_[node] - cluster_mean_vector_[0]);
    for (size_t cluster = 1; cluster < cluster_mean_vector_.size(); cluster++)
    {
        double distance = abs(state_vector_[node]
                - cluster_mean_vector_[cluster]);
        if (distance < min_distance)
        {
            appropriate_cluster = cluster;
            min_distance = distance;
        }
    }

    bool variation = appropriate_cluster != cluster_label_vector_[node];

    if (variation)
    {
        //remove node from previous cluster and assign to new one
        cluster_vector_[cluster_label_vector_[node]].erase(node);
        cluster_label_vector_[node] = appropriate_cluster;
        cluster_vector_[appropriate_cluster].insert(node);
    }
    return variation;
}

//update the mean value for a particular cluster
void ConsensusSolver::update_cluster_mean(size_t cluster)
{
    if(cluster_vector_[cluster].size() > 0)
    {
        cluster_mean_vector_[cluster] = average(cluster_vector_[cluster]);
    }
    else
    {
        //this implies the cluster mean is off
        //choose another mean randomly
        cluster_mean_vector_[cluster] = state_vector_[random_int(size_, gen_)];
    }
}

//uses the k-mean algorithm to converge the clusters and their mean
void ConsensusSolver::converge_cluster()
{
    bool variation = true;
    while(variation)
    {
        //update mean of the clusters
        for (size_t cluster = 0; cluster < cluster_mean_vector_.size();
                cluster++)
        {
            update_cluster_mean(cluster);
        }

        //assignement step
        variation = false;
        for (size_t node = 0; node < size_; node++)
        {
            variation |= assign_to_cluster(node);
        }
    }

    //calculate std for the clusters
    for (size_t cluster = 0; cluster < cluster_std_vector_.size();
            cluster++)
    {
        cluster_std_vector_[cluster] = standard_deviation(
                cluster_vector_[cluster]);
    }
}


/*---------------------------
 *    External functions
 *---------------------------*/
unsigned int random_int(std::size_t size, RNGType& gen)
{
    return floor(std::generate_canonical<double,
        std::numeric_limits<double>::digits>(gen)*size);
}


}//end of namespace soc
