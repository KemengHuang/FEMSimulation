#include "fem_math.h"

MatrixXd vec_double(MatrixXd F) {
    const int cols = F.cols();
    const int rows = F.rows();
    const int nums = cols * rows;
    MatrixXd result(nums, 1);
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            result(i * rows + j, 0) = F(j, i);
        }
    }
    return result;
}

MatrixXf vec_float(MatrixXf F) {
    const int cols = F.cols();
    const int rows = F.rows();
    const int nums = cols * rows;
    MatrixXf result(nums, 1);
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            result(i * rows + j, 0) = F(j, i);
        }
    }
    return result;
}