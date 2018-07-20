#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <ConsensusSolver.hpp>

using namespace std;
using namespace soc;

namespace py = pybind11;

PYBIND11_MODULE(spreading_CR, m)
{
    py::class_<ConsensusSolver>(m, "ConsensusSolver")
        .def(py::init<>())
        .def("get_time_vector", &ConsensusSolver::get_time_vector)
        .def("get_Inode_number_vector", 
            &ConsensusSolver::get_Inode_number_vector)
        .def("get_Rnode_number_vector", 
            &ConsensusSolver::get_Rnode_number_vector)
        .def("is_absorbed", &ConsensusSolver::is_absorbed)
        .def("initialize", (void (ConsensusSolver::*)(double, unsigned int))
            &ConsensusSolver::initialize)
        .def("initialize", (void (ConsensusSolver::*)(vector<NodeLabel>&,
            unsigned int)) &ConsensusSolver::initialize)
        .def("reset", &ConsensusSolver::reset)
        .def("next_state", &ConsensusSolver::next_state)
        .def("evolve", &ConsensusSolver::evolve);
}
