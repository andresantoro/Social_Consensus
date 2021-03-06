#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "QuenchedConsensusSolver.hpp"
#include "AnnealedConsensusSolver.hpp"
#include "BayesianAnnealedConsensusSolver.hpp"
#include "BayesianQuenchedConsensusSolver.hpp"

using namespace std;
using namespace soc;

namespace py = pybind11;

PYBIND11_MODULE(FastConsensusSolver, m)
{
    //Quenched network solver
    py::class_<QuenchedConsensusSolver>(m,
            "QuenchedSolver")
        .def(py::init<unordered_map<Node,vector<Node> >&,
                    unordered_map<Node,double>&, vector<double>&, double,
                    int, size_t, bool>(),
            py::arg("network_map"), py::arg("influence_map"),
            py::arg("initial_state"), py::arg("eta")=0.5,
            py::arg("seed")=42, py::arg("max_cluster")=1,
            py::arg("both_speak")=false)
        .def(py::init<unordered_map<Node,vector<Node> >&,
                    unordered_map<Node,double>&, double, int, size_t, bool>(),
            py::arg("network_map"), py::arg("influence_map"),
            py::arg("eta")=0.5, py::arg("seed")=42, py::arg("max_cluster")=1,
            py::arg("both_speak")=false)
        .def("get_initial_state_vector",
            &QuenchedConsensusSolver::get_initial_state_vector)
        .def("get_state_vector", &QuenchedConsensusSolver::get_state_vector)
        .def("get_mean", &QuenchedConsensusSolver::get_mean)
        .def("get_time", &QuenchedConsensusSolver::get_time)
        .def("get_history_vector", &QuenchedConsensusSolver::get_history_vector)
        .def("reset", &QuenchedConsensusSolver::reset)
        .def("reset_all", &QuenchedConsensusSolver::reset_all)
        .def("consensus_step", &QuenchedConsensusSolver::consensus_step)
        .def("reach_consensus", &QuenchedConsensusSolver::reach_consensus);

    //Annealed network solver
    py::class_<AnnealedConsensusSolver>(m,
            "AnnealedSolver")
        .def(py::init<unordered_map<Node,double>&,
                    unordered_map<Node,double>&, vector<double>&, double,
                    int, size_t, bool>(),
            py::arg("priority_map"), py::arg("influence_map"),
            py::arg("initial_state"), py::arg("eta")=0.5,
            py::arg("seed")=42, py::arg("max_cluster")=1,
            py::arg("both_speak")=false)
        .def(py::init<unordered_map<Node,double>&,
                    unordered_map<Node,double>&, double, int, size_t, bool>(),
            py::arg("priority_map"), py::arg("influence_map"),
            py::arg("eta")=0.5, py::arg("seed")=42, py::arg("max_cluster")=1,
            py::arg("both_speak")=false)
        .def("get_initial_state_vector",
            &AnnealedConsensusSolver::get_initial_state_vector)
        .def("get_state_vector", &AnnealedConsensusSolver::get_state_vector)
        .def("get_mean", &AnnealedConsensusSolver::get_mean)
        .def("get_time", &AnnealedConsensusSolver::get_time)
        .def("get_history_vector", &AnnealedConsensusSolver::get_history_vector)
        .def("reset", &AnnealedConsensusSolver::reset)
        .def("reset_all", &AnnealedConsensusSolver::reset_all)
        .def("consensus_step", &AnnealedConsensusSolver::consensus_step)
        .def("reach_consensus", &AnnealedConsensusSolver::reach_consensus);

    //Bayesian annealed network solver
    py::class_<BayesianAnnealedConsensusSolver>(m,
            "BayesianAnnealedSolver")
        .def(py::init<unordered_map<Node,double>&,
                    unordered_map<Node,double>&, vector<double>&, double,
                    int, size_t, bool>(),
            py::arg("priority_map"), py::arg("influence_map"),
            py::arg("initial_state"), py::arg("eta")=0.5,
            py::arg("seed")=42, py::arg("max_cluster")=1,
            py::arg("both_speak")=false)
        .def(py::init<unordered_map<Node,double>&,
                    unordered_map<Node,double>&, double, int, size_t, bool>(),
            py::arg("priority_map"), py::arg("influence_map"),
            py::arg("eta")=0.5, py::arg("seed")=42, py::arg("max_cluster")=1,
            py::arg("both_speak")=false)
        .def("get_initial_state_vector",
            &BayesianAnnealedConsensusSolver::get_initial_state_vector)
        .def("get_state_vector", &BayesianAnnealedConsensusSolver::get_state_vector)
        .def("get_mean", &BayesianAnnealedConsensusSolver::get_mean)
        .def("get_time", &BayesianAnnealedConsensusSolver::get_time)
        .def("get_history_vector", &BayesianAnnealedConsensusSolver::get_history_vector)
        .def("reset", &BayesianAnnealedConsensusSolver::reset)
        .def("reset_all", &BayesianAnnealedConsensusSolver::reset_all)
        .def("consensus_step", &BayesianAnnealedConsensusSolver::consensus_step)
        .def("reach_consensus", &BayesianAnnealedConsensusSolver::reach_consensus);

    //Bayesian quenched network solver
    py::class_<BayesianQuenchedConsensusSolver>(m,
            "BayesianQuenchedSolver")
        .def(py::init<unordered_map<Node,vector<Node> >&,
                    unordered_map<Node,double>&, vector<double>&, double,
                    int, size_t, bool>(),
            py::arg("network_map"), py::arg("influence_map"),
            py::arg("initial_state"), py::arg("eta")=0.5,
            py::arg("seed")=42, py::arg("max_cluster")=1,
            py::arg("both_speak")=false)
        .def(py::init<unordered_map<Node,vector<Node> >&,
                    unordered_map<Node,double>&, double, int, size_t, bool>(),
            py::arg("network_map"), py::arg("influence_map"),
            py::arg("eta")=0.5, py::arg("seed")=42, py::arg("max_cluster")=1,
            py::arg("both_speak")=false)
        .def("get_initial_state_vector",
            &BayesianQuenchedConsensusSolver::get_initial_state_vector)
        .def("get_state_vector", &BayesianQuenchedConsensusSolver::get_state_vector)
        .def("get_mean", &BayesianQuenchedConsensusSolver::get_mean)
        .def("get_time", &BayesianQuenchedConsensusSolver::get_time)
        .def("get_history_vector", &BayesianQuenchedConsensusSolver::get_history_vector)
        .def("reset", &BayesianQuenchedConsensusSolver::reset)
        .def("reset_all", &BayesianQuenchedConsensusSolver::reset_all)
        .def("consensus_step", &BayesianQuenchedConsensusSolver::consensus_step)
        .def("reach_consensus", &BayesianQuenchedConsensusSolver::reach_consensus);
}
