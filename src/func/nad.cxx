// nad.cxx
// Non-affine displacement calculation.
module;

#include <vector>
#include <stdexcept>
#include <thread>
#include <iostream>
#include <cmath>

import binning;
import types;

export module nad;

export{

// Function to perform matrix multiplication
auto matmul(const std::vector<std::vector<float>>& A,
                const std::vector<std::vector<float>>& B) -> std::vector<std::vector<float>> {
    int rowsA = A.size();
    int colsA = A[0].size();
    int rowsB = B.size();
    int colsB = B[0].size();

    // Check if dimensions match for multiplication
    if (colsA != rowsB) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
    }

    // Initialize result matrix with zeros
    std::vector<std::vector<float>> result(rowsA, std::vector<float>(colsB, 0.0f));

    // Perform multiplication
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsB; ++j) {
            for (int k = 0; k < colsA; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}


// Compute D2min on a span of neighbor vectors.
void thread_dtmin(int start, int end){
    for (int id=start;id<end;++id){
        // std::vector<std::vector<std::array<float, 3>>> delta1(delta[id].begin(), delta[id].begin()+neigh_list[id].n);
        // = std::vector<<>> vec2(vec.begin(), vec.begin() + n)
        // j=
        // v=
        // w=
        int a=0;
    }
}

//
void compute_dtmin(){


}
}

