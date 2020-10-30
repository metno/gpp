%module gridpp
%include "typemaps.i"
%{
#define SWIG_FILE_WITH_INIT

#define INVALID_DIMENSIONS_ERROR(N, DTYPE) SWIG_exception(SWIG_TypeError, "Could not convert input to " #N "D array of type '" #DTYPE "'");
#if 0
    #include <iostream>
    #define PRINT_DEBUG(message) std::cout << message << std::endl;
#else
    #define PRINT_DEBUG(message) ;
#endif
%}
#if defined(SWIGPYTHON)
%include "numpy.i"
#endif
%init %{
#if defined(SWIGPYTHON)
    import_array();
#endif
    gridpp::initialize_omp();
%}

%include exception.i
/*
      SWIG_MemoryError
      SWIG_IOError
      SWIG_RuntimeError
      SWIG_IndexError
      SWIG_TypeError
      SWIG_DivisionByZero
      SWIG_OverflowError
      SWIG_SyntaxError
      SWIG_ValueError
      SWIG_SystemError
*/
%exception {
    try {
        $action
    }
    catch (std::invalid_argument &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_ValueError, s.c_str());
    }
    catch (gridpp::not_implemented_exception &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_RuntimeError, s.c_str());
    }
    catch (std::exception &e) {
        std::string s(e.what());
        SWIG_exception(SWIG_RuntimeError, s.c_str());
    }
    catch (...) {
        SWIG_exception(SWIG_RuntimeError, "Unknown exception");
    }
}
%include "std_vector.i"
%include "std_string.i"
namespace std {
%template(IntVector) vector<int>;
%template(FloatVector) vector<float>;
%template(FloatVector3) vector<vector<vector<float> > >;
%template(IntVector2) vector<vector<int> >;
%template(DoubleVector) vector<double>;
%template(FloatVector2) vector<vector<float> >;
%template(DoubleVector2) vector<vector<double> >;
}

%apply (float* IN_ARRAY1, int DIM1) {(float* v, int n)}
%apply (double* IN_ARRAY1, int DIM1) {(double* v, int n)}

#if defined(SWIGPYTHON)
%define %np_vector_typemaps(DTYPE, NPY_DTYPE)
/*
 * 1D vectors
 */
%typemap(in) std::vector<DTYPE> (std::vector<DTYPE>*ptr=NULL, PyArrayObject* py_array=NULL, PyObject* py_obj=NULL, PyObject* py_obj0=NULL){
    PRINT_DEBUG("Typemap(in) std::vector<DTYPE>");
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        if(num_dims != 1)
            SWIG_exception(SWIG_RuntimeError, "Unknown exception");

        if(array_type($input) == NPY_DTYPE) {
            py_obj = PyArray_FROMANY($input, array_type($input), 1, 1, NPY_ARRAY_DEFAULT);
        }
        else {
            py_obj0 = PyArray_FROMANY($input, array_type($input), 1, 1, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_DTYPE), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int size = array_size($input, 0);
        $1 = std::vector<DTYPE>(arg, arg+size);
    }
    else {
        ptr = new std::vector<DTYPE>();
        int test = swig::asptr($input, &ptr);
        if(ptr == NULL || !SWIG_IsOK(test)) {
            INVALID_DIMENSIONS_ERROR(1, DTYPE);
        }
        $1 = *ptr;
    }
}

%typemap(in) const std::vector<DTYPE> & (std::vector<DTYPE>*ptr=NULL, std::vector<DTYPE> temp, PyArrayObject* py_array=NULL, PyObject* py_obj=NULL, PyObject* py_obj0=NULL){
    PRINT_DEBUG("Typemap(in) const std::vector<DTYPE> &");
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        if(num_dims != 1)
            SWIG_exception(SWIG_RuntimeError, "Vector must be 1 dimensional");
        if(array_type($input) == NPY_DTYPE) {
            py_obj = PyArray_FROMANY($input, array_type($input), 1, 1, NPY_ARRAY_DEFAULT);
        }
        else {
            py_obj0 = PyArray_FROMANY($input, array_type($input), 1, 1, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_DTYPE), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int size = array_size($input, 0);
        temp = std::vector<DTYPE>(arg, arg+size);
        $1 = &temp;
    }
    else {
        ptr = new std::vector<DTYPE>();
        int test = swig::asptr($input, &ptr);
        if(ptr == NULL || !SWIG_IsOK(test)) {
            INVALID_DIMENSIONS_ERROR(1, DTYPE);
        }
        $1 = ptr;
    }
}

/* Same as the const version above */
%typemap(in) std::vector<DTYPE> & (std::vector<DTYPE>*ptr=NULL, std::vector<DTYPE> temp, PyArrayObject* py_array=NULL, PyObject* py_obj=NULL, PyObject* py_obj0=NULL){
    PRINT_DEBUG("Typemap(in) const std::vector<DTYPE> &");
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        if(num_dims != 1)
            SWIG_exception(SWIG_RuntimeError, "Vector must be 1 dimensional");
        if(array_type($input) == NPY_DTYPE) {
            py_obj = PyArray_FROMANY($input, array_type($input), 1, 1, NPY_ARRAY_DEFAULT);
        }
        else {
            py_obj0 = PyArray_FROMANY($input, array_type($input), 1, 1, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_DTYPE), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int size = array_size($input, 0);
        temp = std::vector<DTYPE>(arg, arg+size);
        $1 = &temp;
    }
    else {
        ptr = new std::vector<DTYPE>();
        int test = swig::asptr($input, &ptr);
        if(ptr == NULL || !SWIG_IsOK(test)) {
            INVALID_DIMENSIONS_ERROR(1, DTYPE);
        }
        $1 = ptr;
    }
}

/* When an input is defined as an output, then don't read the variable */
%typemap(in, numinputs=0) std::vector<DTYPE> & OUTPUT (std::vector<DTYPE> temp, PyObject* py_obj=NULL){
    PRINT_DEBUG("Typemap(in) std::vector<DTYPE> & OUTPUT");
    $1 = &temp;
}

/* Inputs that are moved to outputs must be appended */
%typemap(argout) std::vector<DTYPE>& OUTPUT (PyObject* py_obj=NULL) {
    PRINT_DEBUG("Typemap(argout) std::vector<DTYPE>& OUTPUT");
    const std::vector<DTYPE>& temp = *$1;
    int s = temp.size();
    npy_intp dims[1] = {s};
    py_obj = PyArray_ZEROS(1, dims, NPY_DTYPE, 0);
    for(long i = 0; i < s; i++) {
        DTYPE* ref = (DTYPE*) PyArray_GETPTR1((PyArrayObject*) py_obj, i);
        ref[0] = temp[i];
    }
    %append_output(py_obj);
}

%typemap(freearg) std::vector<DTYPE>, const std::vector<DTYPE>&, std::vector<std::vector<DTYPE> >, const std::vector<std::vector<DTYPE> >&, std::vector<std::vector<std::vector<DTYPE> > >, const std::vector<std::vector<std::vector<DTYPE> > >& {
    if(py_obj0$argnum != NULL)
        Py_DECREF(py_obj0$argnum);

    if(py_obj$argnum != NULL)
        Py_DECREF(py_obj$argnum);

    if(ptr$argnum != NULL) {
        delete ptr$argnum;
    }
}

/*
 * We get a segfault if we enable this. Seems python should handle the objects memory
%typemap(freearg) std::vector<DTYPE> & OUTPUT, std::vector<std::vector<DTYPE> >& OUTPUT, std::vector<std::vector<std::vector<DTYPE> > >& OUTPUT {
    if(py_obj$argnum != NULL) {
        Py_DECREF(py_obj$argnum);
    }
}
*/

%typemap(out) std::vector<DTYPE> {
    PRINT_DEBUG("Typemap(out) std::vector<DTYPE>");
    npy_intp dims[1] = {$1.size()};
    $result = PyArray_ZEROS(1, dims, NPY_DTYPE, 0);
    for(long i = 0; i < $1.size(); i++) {
        DTYPE* ref = (DTYPE*) PyArray_GETPTR1((PyArrayObject*) $result, i);
        ref[0] = $1[i];
    }
}

%typecheck(SWIG_TYPECHECK_INTEGER) std::vector<DTYPE>, const std::vector<DTYPE> & {
    PRINT_DEBUG("typecheck std::vector<DTYPE>");
    if(is_array($input)) {
        $1 = array_numdims($input) == 1 ? 1 : 0;
    }
    else {
        bool isok = swig::check<std::vector<DTYPE> >($input);
        $1 = isok;
    }
}

/*
 * 2D vectors
 */
%typemap(in) std::vector<std::vector<DTYPE> > (std::vector<std::vector<DTYPE> >*ptr=NULL, PyArrayObject* py_array=NULL, PyObject* py_obj=NULL, PyObject* py_obj0=NULL){
    PRINT_DEBUG("Typemap(in) std::vector<std::vector<DTYPE> >");
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        if(num_dims != 2)
            SWIG_exception(SWIG_RuntimeError, "Vector must be 2 dimensional");
        if(array_type($input) == NPY_DTYPE) {
            py_obj = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
        }
        else {
            py_obj0 = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_DTYPE), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        std::vector<std::vector<DTYPE> > temp = std::vector<std::vector<DTYPE> >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i] = std::vector<DTYPE>(arg + i*s1, arg + i*s1 + s1);
        }
        $1 = temp;
    }
    else {
        ptr = new std::vector<std::vector<DTYPE> >();
        int test = swig::asptr($input, &ptr);
        if(ptr == NULL || !SWIG_IsOK(test)) {
            INVALID_DIMENSIONS_ERROR(2, DTYPE);
        }
        $1 = *ptr;
    }
}

%typemap(in) const std::vector<std::vector<DTYPE> > & (std::vector<std::vector<DTYPE> >*ptr=NULL, std::vector<std::vector<DTYPE> > temp, PyArrayObject* py_array=NULL, PyObject* py_obj=NULL, PyObject* py_obj0=NULL){
    PRINT_DEBUG("Typemap(in) const std::vector<std::vector<DTYPE> > &");
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        if(num_dims != 2)
            SWIG_exception(SWIG_RuntimeError, "Vector must be 2 dimensional");
        py_obj;
        if(array_type($input) == NPY_DTYPE) {
            py_obj = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
        }
        else {
            py_obj0 = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_DTYPE), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        assert(py_array != NULL);
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        temp = std::vector<std::vector<DTYPE> >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i] = std::vector<DTYPE>(arg + i*s1, arg + i*s1 + s1);
        }
        $1 = &temp;
    }
    else {
        ptr = new std::vector<std::vector<DTYPE> >();
        swig::asptr($input, &ptr);
        int test = swig::asptr($input, &ptr);
        if(ptr == NULL || !SWIG_IsOK(test)) {
            INVALID_DIMENSIONS_ERROR(2, DTYPE);
        }
        $1 = ptr;
    }
}
/* Same as the const version above */
%typemap(in) std::vector<std::vector<DTYPE> > & (std::vector<std::vector<DTYPE> >*ptr=NULL, std::vector<std::vector<DTYPE> > temp, PyArrayObject* py_array=NULL, PyObject* py_obj=NULL, PyObject* py_obj0=NULL){
    PRINT_DEBUG("Typemap(in) std::vector<std::vector<DTYPE> > &");
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        if(num_dims != 2)
            SWIG_exception(SWIG_RuntimeError, "Vector must be 2 dimensional");
        py_obj;
        if(array_type($input) == NPY_DTYPE) {
            py_obj = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
        }
        else {
            py_obj0 = PyArray_FROMANY($input, array_type($input), 2, 2, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_DTYPE), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        assert(py_array != NULL);
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        temp = std::vector<std::vector<DTYPE> >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i] = std::vector<DTYPE>(arg + i*s1, arg + i*s1 + s1);
        }
        $1 = &temp;
    }
    else {
        ptr = new std::vector<std::vector<DTYPE> >();
        swig::asptr($input, &ptr);
        int test = swig::asptr($input, &ptr);
        if(ptr == NULL || !SWIG_IsOK(test)) {
            INVALID_DIMENSIONS_ERROR(2, DTYPE);
        }
        $1 = ptr;
    }
}

/* When an input is defined as an output, then don't read the variable */
%typemap(in, numinputs=0) std::vector<std::vector<DTYPE> > & OUTPUT (std::vector<std::vector<DTYPE> > temp){
    PRINT_DEBUG("Typemap(in) std::vector<std::vector<DTYPE> > & OUTPUT");
    $1 = &temp;
}

/* Inputs that are moved to outputs must be appended */
%typemap(argout) std::vector<std::vector<DTYPE> >& OUTPUT (PyObject* py_obj=NULL){
    PRINT_DEBUG("Typemap(argout) std::vector<std::vector<DTYPE> >& OUTPUT");
    const std::vector<std::vector<DTYPE> >& temp = *$1;
    int s0 = temp.size();
    int s1 = 0;
    if(s0 != 0)
        s1 = temp[0].size();
    npy_intp dims[2] = {s0, s1};
    py_obj = PyArray_ZEROS(2, dims, NPY_DTYPE, 0);
    for(long i = 0; i < s0; i++) {
        for(long j = 0; j < s1; j++) {
            DTYPE* ref = (DTYPE*) PyArray_GETPTR2((PyArrayObject*) py_obj, i, j);
            ref[0] = temp[i][j];
        }
    }
    %append_output(py_obj);
}

%typemap(out) std::vector<std::vector<DTYPE> > {
    PRINT_DEBUG("Typemap(out) std::vector<std::vector<DTYPE> >");
    const std::vector<std::vector<DTYPE> >& temp = $1;
    int s0 = temp.size();
    int s1 = 0;
    if(s0 != 0)
        s1 = temp[0].size();
    npy_intp dims[2] = {s0, s1};
    $result = PyArray_ZEROS(2, dims, NPY_DTYPE, 0);
    for(long i = 0; i < s0; i++) {
        for(long j = 0; j < s1; j++) {
            DTYPE* ref = (DTYPE*) PyArray_GETPTR2((PyArrayObject*) $result, i, j);
            ref[0] = temp[i][j];
        }
    }
}

%typecheck(SWIG_TYPECHECK_INTEGER) std::vector<std::vector<DTYPE> >, const std::vector<std::vector<DTYPE> > & {
    PRINT_DEBUG("typecheck std::vector<std::vector<DTYPE> >");
    if(is_array($input))
        $1 = array_numdims($input) == 2 ? 1 : 0;
    else {
        bool isok = swig::check<std::vector<std::vector<DTYPE> > >($input);
        $1 = isok;
    }
}

/*
 * 3D vectors
 */
%typemap(in) std::vector<std::vector<std::vector<DTYPE> > > (std::vector<std::vector<std::vector<DTYPE> > >*ptr=NULL, PyArrayObject* py_array=NULL, PyObject* py_obj=NULL, PyObject* py_obj0=NULL){
    PRINT_DEBUG("Typemap(in) std::vector<std::vector<std::vector<DTYPE> > >");
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        if(num_dims != 3)
            SWIG_exception(SWIG_RuntimeError, "Vector must be 3 dimensional");
        py_obj;
        if(array_type($input) == NPY_DTYPE) {
            py_obj = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
        }
        else {
            py_obj0 = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_DTYPE), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        int s2 = array_size($input, 2);
        std::vector<std::vector<std::vector<DTYPE> > > temp = std::vector<std::vector<std::vector<DTYPE> > >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i].resize(s1);
            for(int j = 0; j < s1; j++) {
                temp[i][j] = std::vector<DTYPE>(arg + i * s1 * s2 + j * s2, arg + i * s1 * s2 + (j + 1) * s2);
            }
        }
        $1 = temp;
    }
    else {
        ptr = new std::vector<std::vector<std::vector<DTYPE> > >();
        int test = swig::asptr($input, &ptr);
        if(ptr == NULL || !SWIG_IsOK(test)) {
            INVALID_DIMENSIONS_ERROR(3, DTYPE);
        }
        $1 = *ptr;
    }
}

%typemap(in) const std::vector<std::vector<std::vector<DTYPE> > > & (std::vector<std::vector<std::vector<DTYPE> > >*ptr=NULL, std::vector<std::vector<std::vector<DTYPE> > > temp, PyArrayObject* py_array=NULL, PyObject* py_obj=NULL, PyObject* py_obj0=NULL){
    PRINT_DEBUG("Typemap(in) const std::vector<std::vector<std::vector<DTYPE> > > &");
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        if(num_dims != 3)
            SWIG_exception(SWIG_RuntimeError, "Vector must be 3 dimensional");
        py_obj;
        if(array_type($input) == NPY_DTYPE) {
            py_obj = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
        }
        else {
            py_obj0 = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_DTYPE), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        assert(py_array != NULL);
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        int s2 = array_size($input, 2);
        temp = std::vector<std::vector<std::vector<DTYPE> > >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i].resize(s1);
            for(int j = 0; j < s1; j++) {
                temp[i][j] = std::vector<DTYPE>(arg + i * s1 * s2 + j * s2, arg + i * s1 * s2 + (j + 1) * s2);
            }
        }
        $1 = &temp;
    }
    else {
        ptr = new std::vector<std::vector<std::vector<DTYPE> > >();
        int test = swig::asptr($input, &ptr);
        if(ptr == NULL || !SWIG_IsOK(test)) {
            INVALID_DIMENSIONS_ERROR(3, DTYPE);
        }
        $1 = ptr;
    }
}

/* Same as the const version above */
%typemap(in) std::vector<std::vector<std::vector<DTYPE> > > & (std::vector<std::vector<std::vector<DTYPE> > >*ptr=NULL, std::vector<std::vector<std::vector<DTYPE> > > temp, PyArrayObject* py_array=NULL, PyObject* py_obj=NULL, PyObject* py_obj0=NULL){
    PRINT_DEBUG("Typemap(in) const std::vector<std::vector<std::vector<DTYPE> > > &");
    if(is_array($input)) {
        int num_dims = array_numdims($input);
        if(num_dims != 3)
            SWIG_exception(SWIG_RuntimeError, "Vector must be 3 dimensional");
        py_obj;
        if(array_type($input) == NPY_DTYPE) {
            py_obj = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
        }
        else {
            py_obj0 = PyArray_FROMANY($input, array_type($input), 3, 3, NPY_ARRAY_DEFAULT);
            assert(py_obj0 != NULL);
            py_obj = PyArray_CastToType((PyArrayObject*) py_obj0, PyArray_DescrFromType(NPY_DTYPE), 0);
        }
        assert(py_obj != NULL);
        py_array = (PyArrayObject*) py_obj;
        assert(py_array != NULL);
        DTYPE *arg = (DTYPE*) array_data(py_array);
        int s0 = array_size($input, 0);
        int s1 = array_size($input, 1);
        int s2 = array_size($input, 2);
        temp = std::vector<std::vector<std::vector<DTYPE> > >(s0);
        for(int i = 0; i < s0; i++) {
            temp[i].resize(s1);
            for(int j = 0; j < s1; j++) {
                temp[i][j] = std::vector<DTYPE>(arg + i * s1 * s2 + j * s2, arg + i * s1 * s2 + (j + 1) * s2);
            }
        }
        $1 = &temp;
    }
    else {
        ptr = new std::vector<std::vector<std::vector<DTYPE> > >();
        int test = swig::asptr($input, &ptr);
        if(ptr == NULL || !SWIG_IsOK(test)) {
            INVALID_DIMENSIONS_ERROR(3, DTYPE);
        }
        $1 = ptr;
    }
}

/* When an input is defined as an output, then don't read the variable */
%typemap(in, numinputs=0) std::vector<std::vector<std::vector<DTYPE> > > & OUTPUT (std::vector<std::vector<std::vector<DTYPE> > > temp){
    PRINT_DEBUG("Typemap(in) std::vector<std::vector<std::vector<DTYPE> > > & OUTPUT");
    $1 = &temp;
}

/* Inputs that are moved to outputs must be appended */
%typemap(argout) std::vector<std::vector<std::vector<DTYPE> > > OUTPUT (PyObject* py_obj=NULL){
    PRINT_DEBUG("Typemap(argout) std::vector<std::vector<std::vector<DTYPE> > >");
    const std::vector<std::vector<std::vector<DTYPE> > >& temp = $1;
    int s0 = temp.size();
    int s1 = 0;
    if(s0 != 0)
        s1 = temp[0].size();
    int s2 = 0;
    if(s0 != 0 && s1 != 0)
        s2 = temp[0][0].size();
    npy_intp dims[3] = {s0, s1, s2};
    py_obj = PyArray_ZEROS(3, dims, NPY_DTYPE, 0);
    for(long i = 0; i < s0; i++) {
        for(long j = 0; j < s1; j++) {
            for(long k = 0; k < s2; k++) {
                DTYPE* ref = (DTYPE*) PyArray_GETPTR3((PyArrayObject*) py_obj, i, j, k);
                ref[0] = temp[i][j][k];
            }
        }
    }
    %append_output(py_obj);
}

%typemap(out) std::vector<std::vector<std::vector<DTYPE> > > {
    PRINT_DEBUG("Typemap(out) std::vector<std::vector<std::vector<DTYPE> > >");
    const std::vector<std::vector<std::vector<DTYPE> > >& temp = $1;
    int s0 = temp.size();
    int s1 = 0;
    if(s0 != 0)
        s1 = temp[0].size();
    int s2 = 0;
    if(s0 != 0 && s1 != 0)
        s2 = temp[0][0].size();
    npy_intp dims[3] = {s0, s1, s2};
    $result = PyArray_ZEROS(3, dims, NPY_DTYPE, 0);
    for(long i = 0; i < s0; i++) {
        for(long j = 0; j < s1; j++) {
            for(long k = 0; k < s2; k++) {
                DTYPE* ref = (DTYPE*) PyArray_GETPTR3((PyArrayObject*) $result, i, j, k);
                ref[0] = temp[i][j][k];
            }
        }
    }
}

%typecheck(SWIG_TYPECHECK_INTEGER) std::vector<std::vector<std::vector<DTYPE> > >, const std::vector<std::vector<std::vector<DTYPE> > > & {
    PRINT_DEBUG( "typecheck std::vector<std::vector<std::vector<DTYPE> > >");
    if(is_array($input))
        $1 = array_numdims($input) == 3 ? 1 : 0;
    else {
        bool isok = swig::check<std::vector<std::vector<std::vector<DTYPE> > > >($input);
        $1 = isok;
    }
}


//%typecheck(SWIG_TYPECHECK_INTEGER) std::vector<std::vector<DTYPE> > {
//    std::cout << "HERE" << "TYPECHECK 2D" << "DTYPE" << " " << "NPY_DTYPE" << std::endl;
//    $1 = array_numdims($input) == 2 ? 1 : 0;
//}


%enddef
%np_vector_typemaps(int, NPY_INT)
%np_vector_typemaps(long, NPY_LONG)
%np_vector_typemaps(float, NPY_FLOAT)
%np_vector_typemaps(double, NPY_DOUBLE)
#endif

%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& output };
%apply std::vector<std::vector<int> >& OUTPUT { std::vector<std::vector<int> >& count };
%apply std::vector<std::vector<int> >& OUTPUT { std::vector<std::vector<int> >& y_coord };
%apply std::vector<std::vector<int> >& OUTPUT { std::vector<std::vector<int> >& x_coord };
%apply std::vector<float>& OUTPUT { std::vector<float>& distances };
%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& distances };
%apply int& OUTPUT { int& X1_out };
%apply int& OUTPUT { int& Y1_out };
%apply int& OUTPUT { int& X2_out };
%apply int& OUTPUT { int& Y2_out };

%{
#include "gridpp.h"
%}

%include "gridpp.h"
