//
// Created by au642261 on 8/3/25.
//

#ifndef LINALGROUTINES_H
#define LINALGROUTINES_H

#include <type_traits>
#include "Types.h"

// Solve Ax=C assuming A is tridiagonal using the Thomas algorithm
template <typename vecType, typename matType>
void tridiagSolver(vecType& resVec, const matType& ATriag, const vecType& C) {
    // Compile time check to see that vecType is actually compatible with matType...
    using mat_val_t = typename matType::value_type::value_type;
    using vec_val_t = typename vecType::value_type;
    using result_t  = std::common_type_t<mat_val_t, vec_val_t>;
    static_assert(std::is_convertible_v<result_t, typename vecType::value_type>,
        "resVec must be able to store the result of matType × vecType");

    // Prepare temp storage
    matType tempMat = ATriag;
    vecType tempVec = C;
    int N = C.size();

    // First sweep forward
    for (int i=1; i<N; i++) {
        tempMat[1][i] = tempMat[1][i] - tempMat[0][i]/tempMat[1][i-1] * tempMat[2][i-1];  // b'_n
        tempVec[i] = tempVec[i] - tempMat[0][i]/tempMat[1][i-1] * tempVec[i-1];  // d'_n
    }

    // Second sweep backwards
    resVec[N-1] = tempVec[N-1] / tempMat[1][N-1];
    for (int i=N-2; i>=0; --i) {
        resVec[i] = (tempVec[i] - tempMat[2][i]*resVec[i+1]) / tempMat[1][i];
    }
}

// Tridiag matrix vector mult - implement this better or just use LAPACK...
template <typename vecType, typename matType>
void tridiagMatVecMult(vecType& resVec, const matType& triMat, const vecType& vec) {
    using vec_val_t = typename vecType::value_type;
    using mat_val_t = typename matType::value_type::value_type;
    using result_t  = std::common_type_t<mat_val_t, vec_val_t>;
    static_assert(std::is_convertible_v<result_t, vec_val_t>,
        "resVec must be able to store the result of matType × vecType");

    int N = static_cast<int>(vec.size());
    resVec[0] = triMat[1][0]*vec[0] + triMat[2][0]*vec[1];
    resVec[N-1] = triMat[0][N-1]    *vec[N-2] + triMat[1][N-1]*vec[N-1];
    for (int i=1; i<N-1; i++) {
        vec_val_t sum_i(0.);
        for (int j=0; j<3; j++) {
            sum_i += triMat[j][i] * vec[i+j-1];
        }
        resVec[i] = sum_i;
    }
}


dMat dTridiagDiagMult(const dMat& triMat, const dVec& diagMat) {
    int N = static_cast<int>(diagMat.size());
    dMat res(3, dVec(N, 0.));

    res[1][0] = triMat[1][0] * diagMat[0];
    res[1][N-1] = triMat[1][N-1] * diagMat[N-1];
    res[2][0] = triMat[2][0] * diagMat[1];
    res[0][N-1] = triMat[0][N-1] * diagMat[N-2];
    for (int i=1; i<N-1; i++) {
        res[1][i] = triMat[1][i] * diagMat[i];
        res[0][i] = triMat[0][i] * diagMat[i-1];
        res[2][i] = triMat[2][i] * diagMat[i+1];
    }
    return res;
}

#endif //LINALGROUTINES_H
