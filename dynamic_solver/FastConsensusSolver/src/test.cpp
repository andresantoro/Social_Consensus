#include "AnnealedConsensusSolver.hpp"
#include <iostream>

using namespace std;
using namespace soc;

int main(int argc, const char *argv[])
{
    //unordered_map<Node, vector<Node> > network_map;
    //vector<Node> v0,v1;
    //v0.push_back(1);
    //v1.push_back(0);
    //network_map[0] = v0;
    //network_map[1] = v1;
    unordered_map<Node, double> influence_map;
    influence_map[0] = 1;
    influence_map[1] = 1;
    AnnealedConsensusSolver solver(influence_map, influence_map, 0.5,  10);


    solver.reach_consensus(0.01);
    cout << solver.get_time() << endl;


    return 0;
}
