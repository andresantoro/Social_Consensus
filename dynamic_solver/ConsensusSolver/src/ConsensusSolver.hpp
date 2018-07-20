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
#include <vector>
#include <random>
#include <utility>
#include <numeric>

namespace soc
{//start of namespace soc

typedef unsigned int Node; 
typedef std::mt19937 RNGType;

/**
* \class ConsensusSolver ConsensusSolver.hpp
* \brief Class to solve a continuous opinion formation dynamics
*/
class ConsensusSolver 
{
public:
    //Constructor
    ConsensusSolver(std::unordered_map<Node,vector<Node>>& network_map,
                    std::unordered_map<Node,double>& influence_map,
                    std::vector<double>& initial_state_vector,
                    double eta = 0.5, int seed = 42);
    ConsensusSolver(std::unordered_map<Node,vector<Node>>& network_map,
                    std::unordered_map<Node,double>& influence_map,
                    double eta = 0.5, int seed = 42);

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
    void consensus_step();
    void reach_consensus(tol);
    
private:
    std::unordered_map<Node,vector<Node>> network_map_;
    std::unsigned_map<Node,double>& influence_map_;
    std::vector<double> initial_state_vector_;
    std::vector<double> state_vector_;
    std::vector<std::pair<Node, double>> history_vector_;
    RNGType gen_;
    double eta_;
};


}//end of namespace soc 

#endif /* CONSENSUSSOLVER_HPP_ */
