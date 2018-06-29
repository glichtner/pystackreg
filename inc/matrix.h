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

		T& operator()(unsigned int const i, unsigned int const j)
		{
			return m[i * cols + j];
		}

		T* getRowPtr(unsigned int const i)
		{
			return &m[0] + i * cols;
		}

		T* begin() {
			return &m[0];
		}

		T* end() {
			return &m[0] + rows * cols;
		}

		T* getPtr(unsigned int const i, unsigned int const j) {
			return &m[0] + i * cols + j;
		}

		unsigned int nrows() {
			return rows;
		}

		unsigned int ncols() {
			return cols;
		}

    private:
		unsigned int cols;
		unsigned int rows;

};
#endif
