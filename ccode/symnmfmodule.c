#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symnmf.h"

#ifndef SYMNMFModule
#define SYMNMFModule

static int index1(int pointIndex, int elementIndex, int d)
{
	return d * pointIndex + elementIndex;
}

static int index2(int pointIndex, int d)
{
	return d * pointIndex;
}

PyObject* createReturnValue(double* points, int k, int d)
{
	int i, elem, len, len2, j;
	PyObject *temp, *point, *ret = PyList_New(0);
	if (NULL == ret)
		return NULL;
	for (i = 0; i < k; i++)
	{
		temp = PyList_New(0);
		if (NULL == temp)
			goto error;
		for (elem = 0; elem < d; elem++)
		{
			point = Py_BuildValue("lf", points[index1(i, elem, d)]);
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
/*
n - number of points
d - dimension of each point
unpacks input into points 
*/
void ParseInput(int n, int d, double *points, PyObject *input)
{
    int i;
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
    
}

PyObject *DDGWrapper(PyObject* self, PyObject* args)
{
    
}

PyObject *NormWrapper(PyObject* self, PyObject* args)
{
    
}

PyObject* SymNMFWrapper(PyObject* self, PyObject* args)
{
    
}

#endif
