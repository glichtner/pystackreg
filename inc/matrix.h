#include <vector>

#ifndef MATRIX_H_
#define MATRIX_H_

template<typename T>
class matrix {
    std::vector<T>m;

    public:
    	matrix() {
    	}

        matrix(unsigned int rows, unsigned int cols) {
        	resize(rows, cols);
        }

        void resize(unsigned int rows, unsigned int cols) {
            this->rows = rows;
            this->cols = cols;
        	m.resize(rows*cols);
        }

		T& operator()(unsigned int const i, unsigned int const j) noexcept
		{
			return m[i * cols + j];
		}

		T* getRowPtr(unsigned int const i)
		{
			return &(m[i * cols]);
		}

		T* begin() {
			return &m[0];
		}

		T* end() {
			return &m[rows * cols];
		}

		T* getPtr(unsigned int const i, unsigned int const j) {
			return &(m[i * cols + j]);
		}

		unsigned int nrows() {
			return rows;
		}

		unsigned int ncols() {
			return cols;
		}

    private:
		unsigned int cols = 0;
		unsigned int rows = 0;

};
#endif
