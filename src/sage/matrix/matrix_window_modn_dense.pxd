from matrix_window cimport MatrixWindow

cdef extern from "sage/ext/multi_modular.h":
    ctypedef unsigned long mod_int

cdef class MatrixWindow_modn_dense(MatrixWindow):
    pass