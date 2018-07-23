#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <ConsensusSolver.h>

using namespace std;
using namespace soc;

namespace py = pybind11;

PYBIND11_MODULE(FastConsensusSolver, m)
{
    py::class_<ConsensusSolver>(m, "ConsensusSolver")
        .def(py::init<unordered_map<Node,vector<Node>>&,
                    unordered_map<Node,double>&, vector<double>&,
                    double, int>())
        .def(py::init<unordered_map<Node,vector<Node>>&,
                    unordered_map<Node,double>&,double, int>())
        .def("get_initial_state_vector", 
            &ConsensusSolver::get_initial_state_vector)
        .def("get_state_vector", &ConsensusSolver::get_state_vector)
        .def("get_mean", &ConsensusSolver::get_mean)
        .def("get_time", &ConsensusSolver::get_time)
        .def("get_history_vector", &ConsensusSolver::get_history_vector)
        .def("reset", &ConsensusSolver::reset)
        .def("reset_all", &ConsensusSolver::reset_all)
        .def("consensus_step", &ConsensusSolver::consensus_step)
        .def("reach_consensus", &ConsensusSolver::reach_consensus);
}
