#ifndef MATRIXBASE_H
#define MATRIXBASE_H
#endif

#include<iostream>
#include<string>

using namespace std;

template <typename T>
class MatrixBase{
	
	public:
		// constructors
		MatrixBase();
		MatrixBase(size_t num_rows, size_t num_cols);
		MatrixBase(size_t num_rows, size_t num_cols, string type);
		MatrixBase(T **matrix, size_t num_rows, size_t num_cols);
		~MatrixBase();

		// copy constructor
		MatrixBase(const MatrixBase<T>& matrix_to_copy);
		
		// member functions
		void PrintMatrix();
		MatrixBase Transpose();
		// operator overloading
		MatrixBase operator= (const MatrixBase<T>& matrix_to_copy);
		// assign elements at row idx and col idx
		T& operator() (size_t r_idx, size_t c_idx);
		// read elements at row idx and col idx
		T operator() (size_t r_idx, size_t c_idx) const;
		// multiply two matrices
		MatrixBase operator* (const MatrixBase& matrix_to_multipy);
		// add two matrices
		MatrixBase operator+ (const MatrixBase& matrix_to_add);
		// subtract two matrices
		MatrixBase operator- (const MatrixBase& matrix_to_subtract);
		// add a scalar
		MatrixBase operator+ (T x);
		// subtract a scalar
		MatrixBase operator- (T x);
		// multiply a scalar
		MatrixBase operator* (T x);
		// divide by a scalar
		MatrixBase operator/ (T x);
		// getter and setters
		size_t get_ncols();
		size_t get_nrows();
		T get_element(size_t idx_r, size_t idx_c);
		void set_element(size_t idx_r, size_t idx_c, T val);
		bool is_empty();

	protected:
		size_t nrows_, ncols_;
		T **matrix_;


};