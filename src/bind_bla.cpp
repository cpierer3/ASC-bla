#include <pybind11/pybind11.h>
#include <sstream>

#include <matrixview.hpp>

#include "vector.hpp"
#include "matrixview.hpp"
#include <pybind11/numpy.h>

using namespace ASC_bla;
namespace py = pybind11;

PYBIND11_MODULE(bla, m) {
  m.doc() = "Basic linear algebra module"; // optional module docstring

  py::class_<Vector<double>>(m, "Vector")
      .def(py::init<size_t>(), py::arg("size"), "create vector of given size")
      .def("__len__", &Vector<double>::size, "return size of vector")

      .def("__setitem__",
           [](Vector<double> &self, int i, double v) {
             if (i < 0)
               i += self.size();
             if (i < 0 || i >= self.size())
               throw py::index_error("vector index out of range");
             self(i) = v;
           })
      .def("__getitem__", [](Vector<double> &self, int i) { return self(i); })

      .def("__setitem__",
           [](Vector<double> &self, py::slice inds, double val) {
             size_t start, stop, step, n;
             if (!inds.compute(self.size(), &start, &stop, &step, &n))
               throw py::error_already_set();
             self.range(start, stop).slice(0, step) = val;
           })

      .def("__add__",
           [](Vector<double> &self, Vector<double> &other) {
             return Vector<double>(self + other);
           })

      .def("__rmul__", [](Vector<double> &self,
                          double scal) { return Vector<double>(scal * self); })

      .def("__str__",
           [](const Vector<double> &self) {
             std::stringstream str;
             str << self;
             return str.str();
           })

      .def(py::pickle(
          [](Vector<double> &self) { // __getstate__
            /* return a tuple that fully encodes the state of the object */
            return py::make_tuple(self.size(),
                                  py::bytes((char *)(void *)&self(0),
                                            self.size() * sizeof(double)));
          },
          [](py::tuple t) { // __setstate__
            if (t.size() != 2)
              throw std::runtime_error("should be a 2-tuple!");

          Vector<double> v(t[0].cast<size_t>());
          py::bytes mem = t[1].cast<py::bytes>();
          std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()), v.size()*sizeof(double));
          return v;
        }))
    ;

     py::class_<Matrix<double, ASC_bla::RowMajor>>(m, "Matrix")
      .def(py::init<size_t, size_t>(), py::arg("rows"), py::arg("cols"),
           "Construct a Matrix with specified dimensions\n\n"
           "Args:\n"
           "    rows: Number of rows in the matrix\n"
           "    cols: Number of columns in the matrix")
      .def("__getitem__",
           [](Matrix<double, RowMajor> &self, std::tuple<int, int> ind) {
             return self(std::get<0>(ind), std::get<1>(ind));
           })
      .def("__setitem__",
           [](Matrix<double, RowMajor> &self, std::tuple<int, int> ind,
              double v) { self(std::get<0>(ind), std::get<1>(ind)) = v; })
      .def_property_readonly("shape",
                             [](const Matrix<double, RowMajor> &self) {
                               return std::tuple(self.rows(), self.cols());
                             })
      .def("__str__", [](const Matrix<double, ASC_bla::RowMajor> &self) {
        std::stringstream str;
        str << self;
        return str.str();
      })
      .def("__mul__", [](const Matrix<double, ASC_bla::RowMajor> &self,
                          const Matrix<double, ASC_bla::RowMajor> &other) { return Matrix<double, ASC_bla::RowMajor>(self* other); })

      .def("row", [](Matrix<double, RowMajor> &self, int i) {
      return self.row(i);
      }, py::return_value_policy::reference_internal)

      .def("col", [](Matrix<double, RowMajor> &self, int j) {
      return self.col(j);
      }, py::return_value_policy::reference_internal)
      
      .def("transpose", [](const Matrix<double, RowMajor> &self) {
      return self.transpose();
      })
      .def("inverse", [](const Matrix<double, RowMajor> &self) {
      return self.inverse();
      })
      .def("to_numpy", [](Matrix<double, RowMajor> &self) {
      return py::array_t<double>(
        {self.rows(), self.cols()},
        {sizeof(double) * self.cols(), sizeof(double)},
        &self(0,0)
    );
})
      ;
}
