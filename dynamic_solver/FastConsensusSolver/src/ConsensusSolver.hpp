/**
* \file ConsensusSolver.hpp
* \brief Header file for class ConsensusSolver
* \author Guillaume St-Onge
* \version 1.0
* \date 11/07/2018
*/

#ifndef CONSENSUSSOLVER_HPP_
#define CONSENSUSSOLVER_HPP_
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <random>
#include <utility>
#include <numeric>
#include <limits>
#include <algorithm>

namespace soc
{//start of namespace soc

typedef unsigned int Node;
typedef std::mt19937 RNGType;

/**
* \class ConsensusSolver ConsensusSolver.hpp
* \brief Abstract class to solve a continuous opinion formation dynamics
*/
class ConsensusSolver
{
public:
    //Constructor
    ConsensusSolver(std::unordered_map<Node,double>& influence_map,
                    std::vector<double>& initial_state_vector,
                    double eta = 0.5, int seed = 42,
                    std::size_t max_cluster = 1);
    ConsensusSolver(std::unordered_map<Node,double>& influence_map,
                    double eta = 0.5, int seed = 42,
                    std::size_t max_cluster = 1);

    //Destructor
    virtual ~ConsensusSolver() = default;

    //Accessors
    std::vector<double> get_initial_state_vector() const
        {return initial_state_vector_;}
    std::vector<double> get_state_vector() const
        {return state_vector_;}
    double get_mean() const
        {return (std::accumulate(state_vector_.begin(), state_vector_.end(),
            0.0)/state_vector_.size());}
    int get_time() const
        {return history_vector_.size();}
    std::vector<std::pair<Node, double>> get_history_vector() const
        {return history_vector_;}

    //Mutators
    void reset();
    void reset_all();
    virtual void consensus_step() = 0; //to be overloaded
    void reach_consensus(double tol);

protected:
    //utility methods
    void initialize_cluster_mean(std::size_t max_cluster);
    bool assign_to_cluster(std::size_t node);
    void update_cluster_mean(std::size_t cluster);
    void converge_cluster();

    //members
    std::size_t size_;
    std::unordered_map<Node,double> influence_map_;
    std::vector<double> initial_state_vector_;
    std::vector<double> state_vector_;
    std::vector<std::pair<Node, double>> history_vector_;
    std::vector<std::unordered_set<std::size_t>> cluster_vector_;
    std::vector<double> cluster_mean_vector_;
    std::vector<double> cluster_std_vector_;
    std::vector<size_t> cluster_label_vector_;
    RNGType gen_;
    double eta_;

};

/*--------------------------
 * External functions
 *--------------------------*/

//compute standard deviation for a container
template< typename T>
double standard_deviation(T& v)
{
    double mean = std::accumulate(v.begin(), v.end(), 0.0)/v.size();
    double variance = (std::inner_product(v.begin(), v.end(), v.begin(),
                       0.0)/v.size() - mean*mean);
    return sqrt(variance);
}

//compute average for a container
template< typename T>
double average(T& v)
{
    return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

//get a random integer upper bounded by size
unsigned int random_int(std::size_t size, RNGType& gen);


}//end of namespace soc

#endif /* CONSENSUSSOLVER_HPP_ */
