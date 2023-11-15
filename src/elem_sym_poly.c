//
// Created by ishan on 11/12/23.
//


#include <elementary_symmetric_polynomial.h>
#include <stdlib.h>
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#if PY_MAJOR_VERSION < 3
#error "Requires Python 3"
#include "stopcompilation"
#endif


#define DEFAULT_PRECISION 53

static PyObject * ElemSymPolyError;

// This assumes it is passed a list of python floats ( C doubles ) of length power of 2
static PyObject * py_compute_elem_sym_poly(PyObject * Py_UNUSED(self), PyObject *args){

	PyObject * elements_lst;
	int N;
	int actual_size;
	int precision = DEFAULT_PRECISION;

	if (!PyArg_ParseTuple(args, "O!ii", &PyList_Type, &elements_lst, &actual_size, &precision)){
		return NULL;
	}

	N = (int)PyList_Size(elements_lst);

	double * elements = malloc(sizeof(double) * N);

	for (int i = 0;i< N; i++){
		PyObject* item = PyList_GetItem(elements_lst, i);

		if (!PyFloat_Check(item) || PyErr_Occurred()){
			PyErr_SetString(ElemSymPolyError, "List must contain only floats");
			free(elements);
			return NULL;
		}
		else{
			elements[i] = PyFloat_AsDouble(item);
		}
	}

	elementary_symmetric_state_t  state = init_elementary_symmetric_state(N, precision);

	update_elementary_symmetric_state(state, elements, N);

	double * poly = compute_elementary_symmetric_polynomials(state);

	free_elementary_symmetric_state(state);
	mpfr_mp_memory_cleanup();

	PyObject * poly_lst = PyList_New(actual_size+1);

	for (int i = 0;i< actual_size+1; i++){
		PyObject * item = PyFloat_FromDouble(poly[i]);

		if (item == NULL){
			PyErr_SetString(ElemSymPolyError, "Internal error");
			free(elements);
			free(poly);
			return NULL;
		}

		PyList_SetItem(poly_lst, i, item);
	}

	free(poly);
	free(elements);

	return poly_lst;

};

static PyMethodDef ElemSymPolyMethods[] = {
		{"compute_elem_sym_poly", py_compute_elem_sym_poly, METH_VARARGS, "Compute the elementary symmetric polynomials"},
		{NULL, NULL, 0, NULL}
};

static struct PyModuleDef ElemSymPolyDef = {
		PyModuleDef_HEAD_INIT,
		"elem_sym_poly",
		NULL,
		-1,
		ElemSymPolyMethods,
};

PyMODINIT_FUNC PyInit_pyElemSymPoly(void){


	PyObject * module = PyModule_Create(&ElemSymPolyDef);


	if (module == NULL){
		return NULL;
	}

	ElemSymPolyError = PyErr_NewException("elem_sym_poly.error", NULL, NULL);
	Py_INCREF(ElemSymPolyError);

	if (PyModule_AddObject(module, "error", ElemSymPolyError ) < 0){
		Py_DECREF(ElemSymPolyError);
		Py_CLEAR(ElemSymPolyError);
		Py_DECREF(module);
		return NULL;
	}

	return module;
}

#undef DEFAULT_PRECISION