#pragma once

class matrix {
  private:
    double *_data;
    size_t _rows;
    size_t _cols;
  public:
    size_t rows() const { return _rows; }
    size_t cols() const { return _cols; }
    size_t size() const { return _rows*_cols; }
    double* data() const { return _data; }

    // Getters
    double operator()(const size_t i) const { return _data[i]; }
    double operator()(const size_t i, const size_t j) const { return _data[i*_cols+j]; }
    // Setters
    double& operator()(const size_t i) { return _data[i]; }
    double& operator()(const size_t i, const size_t j) { return _data[i*_cols+j]; }

    // Constructors
    matrix() {
        _data = nullptr;
        _rows = _cols = 0;
        //std::cout << "Default constructor" << std::endl;
    }

    matrix(size_t _rows, size_t _cols = 1) {
        this->_rows = _rows;
        this->_cols = _cols;
        this->_data = new double[_rows*_cols];
        //std::cout << "Base constructor" << std::endl;
    }

    matrix(const matrix &M) {
        this->_rows = M.rows();
        this->_cols = M.cols();
        this->_data = new double[_rows*_cols];
        std::copy(M.data(), M.data()+M.size(), this->_data);
        //std::cout << "Copy constructor" << std::endl;
    }

    // Destructor
    ~matrix() {
        if ( _data != nullptr )
            delete [] _data;
        //std::cout << "Destructor" << std::endl;
    }

    // Assignment operator
    matrix& operator=(const matrix& M) {
        std::cout << "Assignment" << std::endl;
       // if different size or not allocated
        if ( _rows*_cols != M.size() ) {
            std::cout << "Assignment.alloc" << std::endl;
            delete [] _data;
            this->_data = new double[M.size()];
        }
        this->_rows = M.rows();
        this->_cols = M.cols();
        std::copy(M.data(), M.data()+M.size(), this->_data);
        return *this;
    }

    // Matrix-Matrix Multiplication
    matrix operator*(const matrix &M) {
        matrix R(this->rows(), M.cols());
        for (size_t i = 0; i < this->rows(); i++) {
            for (size_t j = 0; j < M.cols(); j++) {
                double sum = 0.0;
                for (size_t k = 0; k < this->cols(); k++) {
                    sum += (*this)(i,k) * M(k,j);
                }
                R(i,j) = sum;
            }
        }
        return R;
    }

    // Matrix Transpose
    matrix transpose() const {
        matrix T(_cols, _rows);
        for (size_t i = 0; i < T.rows() ; i++) {
            for (size_t j = 0; j < T.cols(); j++) {
                T(i,j) = (*this)(j,i);
            }
        }
        return T;
    }

    // Set all elements to a value
    void fill(const double val) {
        for (size_t i = 0; i < _rows*_cols; i++) {
            _data[i] = val;
        }
    }
};

// Utility functions
void print(const matrix &M) {
    for (size_t i = 0; i < M.rows(); i++) {
        std::cout << std::endl;
        for (size_t j = 0; j < M.cols(); j++) {
            std::cout << "\t" << M(i,j);
        }
    }
    std::cout << std::endl;
}


