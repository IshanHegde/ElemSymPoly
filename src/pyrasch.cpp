#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

#include <pybind11/eigen.h>

#include "Rsm.cpp"

namespace py = pybind11;






PYBIND11_MODULE(pyrasch,m)
    {
        m.doc() = "pyrasch";
        
        py::class_<RSM>(m,"RSM")
        .def(py::init<const Eigen::MatrixXd &>())
        .def("PROX",&RSM::PROX)
        .def("JMLE",&RSM::JMLE)
        .def("estimate_thresholds",&RSM::estimate_thresholds)
        .def("estimate_counts",&RSM::estimate_counts)
        .def_readwrite("expected_value",&RSM::expected_value,py::return_value_policy::reference_internal)
        .def_readwrite("variance",&RSM::variance,py::return_value_policy::reference_internal)
        .def_readwrite("data_probability",&RSM::data_probability,py::return_value_policy::reference_internal)
        .def_readwrite("RA_Thresholds",&RSM::RA_Thresholds,py::return_value_policy::reference_internal)
        .def_readwrite("observed_counts",&RSM::observed_counts,py::return_value_policy::reference_internal)
        .def_readwrite("ability",&RSM::ability,py::return_value_policy::reference_internal)
        .def_readwrite("difficulty",&RSM::difficulty,py::return_value_policy::reference_internal);

    }
