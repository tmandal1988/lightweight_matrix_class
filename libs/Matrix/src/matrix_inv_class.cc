#include "matrix_inv_class.h"

template <typename T>
MatrixInv<T> MatrixInv<T>::Inverse(){
	size_t nrows_ = this->nrows_;
	size_t ncols_ = this->ncols_;
	T **matrix_ = this->matrix_;
	if (nrows_ != ncols_)
		throw invalid_argument("Matrix is not square\n");

	T orig_matrix[ncols_][nrows_];
	for (size_t idx_r = 0; idx_r < nrows_; idx_r++){
		for (size_t idx_c = 0; idx_c < ncols_; idx_c++)
			orig_matrix[idx_r][idx_c] = matrix_[idx_r][idx_c];
	}

	MatrixInv<T> inverse_matrix(nrows_, ncols_, "eye");

	for(size_t idx_r = 0; idx_r < nrows_; idx_r++){
		T diag_val = orig_matrix[idx_r][idx_r];
		if ( diag_val == 0){
			for( size_t idx_non_zero = idx_r + 1; idx_non_zero < nrows_; idx_r ++){
				if ( orig_matrix[idx_non_zero][idx_r] != 0){
					// swap rows
					for (size_t idx_c = 0; idx_c < ncols_; ++idx_c)
					{
						T temp_val_orig = orig_matrix[idx_r][idx_c];
						orig_matrix[idx_r][idx_c] = orig_matrix[idx_non_zero][idx_c];
						orig_matrix[idx_non_zero][idx_c] = temp_val_orig;

						temp_val_orig = inverse_matrix(idx_r, idx_c);
						inverse_matrix(idx_r, idx_c) = inverse_matrix(idx_non_zero, idx_c);
						inverse_matrix(idx_non_zero, idx_c) = temp_val_orig;
					}
					diag_val = orig_matrix[idx_r][idx_r];
					break;
				}
			}
		}

		// divide all the columns for idx_r row by diag_val
		for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
			orig_matrix[idx_r][idx_c] = orig_matrix[idx_r][idx_c]/diag_val;
			inverse_matrix(idx_r, idx_c) = inverse_matrix(idx_r, idx_c)/diag_val;
		}
		for(size_t idx2_r = 0; idx2_r < nrows_; idx2_r++){
			if (idx2_r != idx_r){
				T value_to_multiply = orig_matrix[idx2_r][idx_r];
				for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
					orig_matrix[idx2_r][idx_c] = orig_matrix[idx2_r][idx_c] - orig_matrix[idx_r][idx_c]*value_to_multiply;
					inverse_matrix(idx2_r, idx_c) = inverse_matrix(idx2_r, idx_c) - inverse_matrix(idx_r, idx_c)*value_to_multiply;

				}
			}
		}
	}

	return inverse_matrix;
}

template <typename T>
MatrixInv<T> MatrixInv<T>::InverseUsingQr(){
	size_t nrows_ = this->nrows_;
	size_t ncols_ = this->ncols_;
	MatrixInv<T> q(nrows_, ncols_, "eye");
	MatrixInv<T> r(this->matrix_, nrows_, ncols_);
	T v1, v2;
	if (nrows_ != ncols_)
		throw invalid_argument("Matrix is not square\n");

	for(size_t idx_c = 0; idx_c < ncols_; idx_c++){
		for(size_t idx_r = idx_c + 1; idx_r < nrows_; idx_r++){
			if( abs(r(idx_c, idx_c)) < 1e-20 ){
				for(size_t idx_cc = idx_c; idx_cc < ncols_; idx_cc++){
					v1 = r(idx_r, idx_cc);
					r(idx_r, idx_cc) = r(idx_c, idx_cc);
					r(idx_c, idx_cc) = v1;
				}

					for(size_t idx_qc = 0; idx_qc < nrows_; idx_qc++){
						v2 = q(idx_r, idx_qc);
						q(idx_r, idx_qc) = q(idx_c, idx_qc);
						q(idx_c, idx_qc) = v2;
					}
			}

			T temp = r(idx_r, idx_c)/r(idx_c, idx_c);
			for(size_t idx_qcc = 0; idx_qcc < ncols_; idx_qcc++){
				q(idx_r, idx_qcc) -= temp*q(idx_c, idx_qcc);
			}
			
			for (size_t idx_rcc = idx_c; idx_rcc < ncols_; idx_rcc++){		
				r(idx_r, idx_rcc) -= temp*r(idx_c, idx_rcc);
			}
		}
	}

	return r.BackSubstitution()*q;
}

template<typename T>
MatrixInv<T> MatrixInv<T>::BackSubstitution(){
	size_t nrows_ = this->nrows_;
	size_t ncols_ = this->ncols_;
	MatrixInv<T> rb(nrows_, ncols_);

	if (nrows_ != ncols_)
		throw invalid_argument("Matrix is not square\n");

	for(int idx_r = nrows_; idx_r-- > 0;){

		rb(idx_r, idx_r) = 1/(this->matrix_[idx_r][idx_r]);

		for(int idx_c = idx_r; idx_c-- > 0;){

			for(size_t idx_k = idx_c; idx_k <= idx_r - 1; idx_k++){

				rb(idx_c, idx_r) = rb(idx_c, idx_r) -  this->matrix_[idx_c][idx_k + 1]*rb(idx_k+1, idx_r);
			}
			rb(idx_c, idx_r) = rb(idx_c, idx_r)/this->matrix_[idx_c][idx_c];
		}
	}

	return rb;
}

// Explicit template instantiation
template class MatrixInv<float>;
template class MatrixInv<double>;