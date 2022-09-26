#include "matrix_base_class.h"

template <typename T>
MatrixBase<T>::MatrixBase(){
	nrows_ = 0;
	ncols_ = 0;
	matrix_ = nullptr;
}

template <typename T>
MatrixBase<T>::MatrixBase(size_t num_rows, size_t num_cols){
	nrows_ = num_rows;
	ncols_ = num_cols;
	matrix_ = new T*[nrows_];
	// assign memory for the columns
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r] = new T[ncols_];
	}
	// initialize MatrixBase with zeros
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			matrix_[idx_r][idx_c] = 0;
		}
	}

}

template <typename T>
MatrixBase<T>::MatrixBase(size_t num_rows, size_t num_cols, string type){
	nrows_ = num_rows;
	ncols_ = num_cols;
	matrix_ = new T*[nrows_];
	// assign memory for the columns
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r] = new T[ncols_];
	}

	if (type == "eye"){
		if (num_rows != num_cols){
			throw invalid_argument("To create identity Matrix number of rows and columns should be equal\n");
		}
			// initialize MatrixBase with zeros and then assign 1 to diagonal
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
				if (idx_r == idx_c)
					matrix_[idx_r][idx_c] = 1;
				else
					matrix_[idx_r][idx_c] = 0;
			}
		}
	}

}

template <typename T>
MatrixBase<T>::MatrixBase(T **matrix, size_t num_rows, size_t num_cols){
	nrows_ = num_rows;
	ncols_ = num_cols;
	matrix_ = new T*[nrows_];

	// assign memory for the columns
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r] = new T[ncols_];
	}

	// initialize MatrixBase with zeros and then assign 1 to diagonal
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
				matrix_[idx_r][idx_c] = matrix[idx_r][idx_c];
		}
	}
}

template <typename T>
MatrixBase<T>::MatrixBase(const MatrixBase<T>& matrix_to_copy){
	nrows_ = matrix_to_copy.nrows_;
	ncols_ = matrix_to_copy.ncols_;
	matrix_ = new T*[nrows_];
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r] = new T[ncols_];
	}
	// initialize MatrixBase with zeros
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			matrix_[idx_r][idx_c] = matrix_to_copy.matrix_[idx_r][idx_c];
		}
	}
}

template <typename T>
MatrixBase<T>::~MatrixBase(){
	if (ncols_ > 0){
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			delete[] matrix_[idx_r];
		}
	}

	if(nrows_ > 0){
		delete[] matrix_;
	}

}

// member functions
template <typename T>
void MatrixBase<T>::PrintMatrix(){
	printf("*****************************\n");
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			printf("%g\t", matrix_[idx_r][idx_c]);
		}
		printf("\n");
	}
	printf("*****************************\n");
}

template <typename T>
MatrixBase<T> MatrixBase<T>::Transpose(){
	MatrixBase<T> transpose_matrix(ncols_, nrows_);
	for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			transpose_matrix(idx_c, idx_r) = matrix_[idx_r][idx_c];
		}
	}
	return transpose_matrix;
}

// operator overloading
template <typename T>
MatrixBase<T> MatrixBase<T>::operator= (const MatrixBase<T>& matrix_to_copy){
	if( &matrix_to_copy == this){
		return *this;
	}

	if (ncols_ > 0){
		for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
			delete[] matrix_[idx_r];
		}
	}

	if(nrows_ > 0){
		delete[] matrix_;
	}

	nrows_ = matrix_to_copy.nrows_;
	ncols_ = matrix_to_copy.ncols_;
	matrix_ = new T*[nrows_];
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		matrix_[idx_r] = new T[ncols_];
	}
	// initialize Matrix with the matrix_to_copy values
	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			matrix_[idx_r][idx_c] = matrix_to_copy.matrix_[idx_r][idx_c];
		}
	}
	return *this;
}

template <typename T>
T& MatrixBase<T>::operator()(size_t r_idx, size_t c_idx){
	if (r_idx >= nrows_ || c_idx >= ncols_)
	{
		throw out_of_range("Matrx subscript out of bounds\n");
	}

	return matrix_[r_idx][c_idx];
}

template <typename T>
T MatrixBase<T>::operator() (size_t r_idx, size_t c_idx) const {
	if (r_idx >= nrows_ || c_idx >= ncols_)
	{
		throw out_of_range("Matrx subscript out of bounds\n");
	}

	return matrix_[r_idx][c_idx];
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator* (const MatrixBase& matrix_to_multipy){
	if (ncols_ != matrix_to_multipy.nrows_)
		throw invalid_argument("MatrixBase dimensions invalid for multiplication");
	MatrixBase<T> return_MatrixBase(nrows_, matrix_to_multipy.ncols_);
	for (size_t idx_c = 0; idx_c < ncols_; idx_c++) {
    	for (size_t idx_r = 0; idx_r < nrows_; idx_r++) {
        	for (size_t idx2_c = 0; idx2_c < matrix_to_multipy.ncols_; idx2_c++) {
            	return_MatrixBase.matrix_[idx_r][idx2_c] += matrix_[idx_r][idx_c] * matrix_to_multipy.matrix_[idx_c][idx2_c];
        	}
    	}
	}

	return return_MatrixBase;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator+ (const MatrixBase& matrix_to_add){
	// check dimensions
	if (nrows_ != matrix_to_add.nrows_ || ncols_ != matrix_to_add.ncols_)
		throw invalid_argument("MatrixBase dimensions invalid for addition");

	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] + matrix_to_add.matrix_[idx_r][idx_c];
	}

	return return_MatrixBase;

}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator- (const MatrixBase& matrix_to_subtract){
	// check dimensions
	if (nrows_ != matrix_to_subtract.nrows_ || ncols_ != matrix_to_subtract.ncols_)
		throw invalid_argument("MatrixBase dimensions invalid for subtraction");

	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] - matrix_to_subtract.matrix_[idx_r][idx_c];
	}

	return return_MatrixBase;

}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator+ (T x){
	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] + x;
	}

	return return_MatrixBase;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator- (T x){
	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] - x;
	}

	return return_MatrixBase;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator* (T x){
	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] * x;
	}

	return return_MatrixBase;
}

template <typename T>
MatrixBase<T> MatrixBase<T>::operator/ (T x){
	MatrixBase<T> return_MatrixBase(nrows_, ncols_);

	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			return_MatrixBase(idx_r, idx_c) = matrix_[idx_r][idx_c] / x;
	}

	return return_MatrixBase;
}

// getters and setters
template<typename T>
size_t MatrixBase<T>::get_ncols(){
	return ncols_;
}

template<typename T>
size_t MatrixBase<T>::get_nrows(){
	return nrows_;
}

template<typename T>
T MatrixBase<T>::get_element(size_t idx_r, size_t idx_c){
	return matrix_[idx_r][idx_c];
}

template<typename T>
void MatrixBase<T>::set_element(size_t idx_r, size_t idx_c, T val){
	matrix_[idx_r][idx_c] = val;
}

template<typename T>
bool MatrixBase<T>::is_empty(){
	if (matrix_ == NULL){
		return true;
	}else{
		return false;
	}
}

// Explicit template instantiation
template class MatrixBase<float>;
template class MatrixBase<double>;