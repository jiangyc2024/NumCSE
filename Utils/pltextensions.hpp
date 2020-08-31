/**
 * @brief Include this file to access some extensions to matplotlibcpp. It
 * mimics the original library by injection of its functions into the usual
 * namespace. Use this as a substitute for (yet) missing pyplot functions in
 * matplotlibcpp, to add overloads/template versions of existing functions or to
 * alter default behaviour.
 */

#pragma once

#include "matplotlibcpp.h"

namespace matplotlibcpp {

typedef std::map<std::string, std::string> Keywords_t;

namespace detail {

PyObject *safe_import(PyObject &pymod, const std::string &function_name) {
  auto func = PyObject_GetAttrString(&pymod, "triplot");
  if (!func || !PyFunction_Check(func)) {
    throw std::runtime_error("Error finding function " + function_name);
  }
  return func;
}

struct _extension {
  
  // functions go here
  PyObject *_triplot; 

  static _extension &get() {
    static _extension ext;
    return ext;
  }

 private:
  _extension() {
    _interpreter::get();  // Trigger construction of the vanilla interpreter
    auto pyplotname = PyString_FromString("matplotlib.pyplot");
    if (pyplotname) {
      Py_DECREF(pyplotname);
      auto pymod = PyImport_Import(pyplotname);
      if (pymod) {
        _triplot = safe_import(*pymod, "triplot");
      } else
        std::runtime_error("Error loading pyplot");
    } else
      std::runtime_error("Error creating string");
  }

  ~_extension() { Py_Finalize(); }
};

}  // end namespace detail

PyObject *kwargs_from_keywords(const Keywords_t &keywords = {}) {
  auto kwargs = PyDict_New();
  for (std::map<std::string, std::string>::const_iterator it = keywords.begin();
       it != keywords.end(); ++it) {
    PyDict_SetItemString(kwargs, it->first.c_str(),
                         PyUnicode_FromString(it->second.c_str()));
  }

  return kwargs;
}

// Extend functions so one can plot 2d meshes
template <typename Vector, typename Matrix>
bool triplot(const Vector &x, const Vector &y, const Matrix &T,
             const Keywords_t &keywords = {}) {
  assert(x.size() == y.size());

  auto plot_args = PyTuple_New(3);
  PyTuple_SetItem(plot_args, 0, get_array(x));
  PyTuple_SetItem(plot_args, 1, get_array(y));
  PyTuple_SetItem(plot_args, 2, get_2darray(T));

  auto kwargs = kwargs_from_keywords(keywords);
  auto res =
      PyObject_Call(detail::_extension::get()._triplot, plot_args, kwargs);

  Py_DECREF(plot_args);
  Py_DECREF(kwargs);

  if (res) Py_DECREF(res);
  return res;
}

// Overload the existing text() function to allow keyword arguments
template <typename Numeric>
void text(Numeric x, Numeric y, const std::string &s,
          const Keywords_t &keywords = {}) {
  assert(x.size() == y.size());

  auto args = PyTuple_New(3);
  PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));
  PyTuple_SetItem(args, 1, PyFloat_FromDouble(y));
  PyTuple_SetItem(args, 2, PyString_FromString(s.c_str()));

  auto kwargs = kwargs_from_keywords(keywords);

  auto res = PyObject_Call(detail::_interpreter::get().s_python_function_text,
                           args, kwargs);

  Py_DECREF(args);
  Py_DECREF(kwargs);

  if (res)
    Py_DECREF(res);
  else
    throw std::runtime_error("Call to text() failed.");
}

}  // namespace matplotlibcpp
