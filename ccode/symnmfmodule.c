#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symnmf.h"

static int index1(int pointIndex, int elementIndex, int d)
{
	return d * pointIndex + elementIndex;
}

PyObject* CreateReturnValue(double* points, int n, int d)
{
	int i, elem, len, len2, j;
	PyObject *temp, *point, *ret = PyList_New(0);
	if (NULL == ret)
		return NULL;
	for (i = 0; i < n; i++)
	{
		temp = PyList_New(0);
		if (NULL == temp)
			goto error;
		for (elem = 0; elem < d; elem++)
		{
			point = Py_BuildValue("d", points[index1(i, elem, d)]);
			if (NULL == point)
			{
				for (j = 0; j < elem; j++)
				{
					free(PyList_GetItem(temp, j));
				}
				goto error;
			}
			PyList_Append(temp, point);
		}
		PyList_Append(ret, temp);
	}
	return ret;

error:
	len = PyList_Size(ret);
	for (i = 0; i < len; i++)
	{
		PyObject *point = PyList_GetItem(ret, i);
		len2 = PyList_Size(point);
		for (j = 0; j < len2; j++)
		{
			free(PyList_GetItem(point, j));
		}
		free(point);
	}
	free(ret);
	return NULL;
}

//n - number of points
//d - dimension of each point
//unpacks input into points 

void ParseInput(int n, int d, double *points, PyObject *input)
{
    int i, j;
    PyObject *temp;
    for (i = 0; i < n; i++)
	{
		temp = PyList_GetItem(input, i);
		for (j = 0; j < d; j++)
		{
			points[index1(i, j, d)] = PyFloat_AsDouble(PyList_GetItem(temp, j));
		}
	}
}

PyObject *SymWrapper(PyObject* self, PyObject* args)
{
    int n, d, status;
    double *points, *res;
    PyObject *object;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &object))
	{
		return NULL; 
	}
    points = calloc(n * d, sizeof(double));
	if (NULL == points)
		return NULL;
    ParseInput(n, d, points, object);
    status = 0;
    res = sym(n, d, points, &status);
    free(points);
    if (0 != status)
    {
        return NULL;
    }
    else 
    {
        object = CreateReturnValue(res, n, n);
        free(res);
        return object;
    }
}

PyObject *DDGWrapper(PyObject* self, PyObject* args)
{
    int n, d, status;
    double *points, *res;
    PyObject *object;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &object))
	{
		return NULL; 
	}
    points = calloc(n * d, sizeof(double));
	if (NULL == points)
		return NULL;
    ParseInput(n, d, points, object);
    status = 0;
    res = ddg(n, d, points, &status);
    free(points);
    if (0 != status)
    {
        return NULL;
    }
    else 
    {
        object = CreateReturnValue(res, n, n);
        free(res);
        return object;
    }
}

PyObject *NormWrapper(PyObject* self, PyObject* args)
{
    int n, d, status;
    double *points, *res;
    PyObject *object;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &object))
	{
		return NULL; 
	}
    points = calloc(n * d, sizeof(double));
	if (NULL == points)
		return NULL;
    ParseInput(n, d, points, object);
    status = 0;
    res = norm(n, d, points, &status);
    free(points);
    if (0 != status)
    {
        return NULL;
    }
    else 
    {
        object = CreateReturnValue(res, n, n);
        free(res);
        return object;
    }
}

PyObject* SymNMFWrapper(PyObject* self, PyObject* args)
{
    int n, k, status;
    double *w, *h, *res;
    PyObject *object1, *object2;
    if (!PyArg_ParseTuple(args, "iiOO", &n, &k, &object1, &object2))
	{
		return NULL; 
	}
    w = calloc(n * n, sizeof(double));
	if (NULL == w)
		return NULL;
    h = calloc(n * k, sizeof(double));
	if (NULL == h)
    {
        free(w);
		return NULL;
    }
    ParseInput(n, n, w, object1);
    ParseInput(n, k, h, object2);
    status = 0;
    res = symnmf(n, k, w, h, &status);
    free(h);
	free(w);
    if (0 != status)
    {
        return NULL;
    }
    else 
    {
        object1 = CreateReturnValue(res, n, k);
        free(res);
        return object1;
    }
}

static PyMethodDef methods[] = {
	{
		"sym",
		(PyCFunction)SymWrapper,
		METH_VARARGS,
		PyDoc_STR("Creates the symmetric A matrix.")
	}, 
    {
		"ddg",
		(PyCFunction)DDGWrapper,
		METH_VARARGS,
		PyDoc_STR("Creates the diagonal D matrix.")
	}, 
    {
		"norm",
		(PyCFunction)NormWrapper,
		METH_VARARGS,
		PyDoc_STR("Creates the normalized w matrix.")
	}, 
    {
		"symnmf",
		(PyCFunction)SymNMFWrapper,
		METH_VARARGS,
		PyDoc_STR("Provided with the w norm matrix and an initialized h-matrix, converges into the h-classifier matrix.")
	}, 
	{
		NULL, NULL, 0, NULL
	}
};

static struct PyModuleDef moduleDef =
{
	PyModuleDef_HEAD_INIT, "symnmf_c_api", NULL, -1, methods
};

PyMODINIT_FUNC PyInit_symnmf_c_api(void)
{
	PyObject* m = PyModule_Create(&moduleDef);
	if (!m) {
		return NULL;
	}
	return m;
}