// This is a template for implementing the Davidson algorithm
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "davidson_class.h"

int Davidson_class::diagonalize_with_davidson(const size_t nroots) {
    int error = 0;

    // set parameters
    size_t m_dim = matrix.n_rows;
    size_t kmin = nroots + 1;
    size_t kmax = nroots + 10;
    size_t max_iter = nroots * 40;
    double rtol = 0.0000001;

    // Allocate memory for basis vectors, v, and sigma vectors, c.
    arma::mat v = arma::zeros<arma::mat>(m_dim, kmax);
    arma::mat c = arma::zeros<arma::mat>(m_dim, kmax);

    // Allocate memory for residual and new vector
    arma::vec resid = arma::zeros<arma::vec>(m_dim);
    arma::vec newv = arma::zeros<arma::vec>(m_dim);

    // Create instances for Krylov subspace
    arma::mat kmat;
    arma::mat kevec;
    arma::vec keval;

    // Create instance of rnorm
    double rnorm;

    size_t citer = 1;
    size_t croot = 0;
    size_t kdim = kmin;

    // Create initial guess.
    create_initial_guess(v, kdim);

    while (citer < max_iter && croot <= nroots) {

        // reset kdim
        kdim = kmin;

        c = matrix * v;

        // build the kmatrix, v.Ac, and diagonalize.
        build_kmat(v.cols(0, kdim - 1), c.cols(0, kdim - 1), kmat);

        arma::eig_sym(keval, kevec, kmat);

        while (kdim < kmax && citer < max_iter) {

            // Compute residual vector, and its norm
            get_residual(keval, kevec, v.cols(0, kdim - 1), c.cols(0, kdim - 1),
                         croot, kdim, resid);

            rnorm = arma::norm(resid);

            if (rnorm < rtol)
                break;

            print_iteration_information(keval, citer, rnorm);

            // Form new vector. Orthonormalize and add to space
            make_newvector(resid, keval, croot, newv);
            orthonormalize_vector(v.cols(0, kdim - 1), kdim, newv);
            v.col(kdim) = newv;

            kdim++;

            check_basis_orthonormality(v.cols(0, kdim - 1));

            c = matrix * v;

            // Build new kmatrix, v.Ac, and diagonalize.
            build_kmat(v.cols(0, kdim - 1), c.cols(0, kdim - 1), kmat);
            arma::eig_sym(keval, kevec, kmat);

            citer++;
        }

        if (rnorm > rtol) {
            // truncate basis
            truncate_basis(kevec, kmin, kmax, v);
            check_basis_orthonormality(v.cols(0, kmin - 1));
            // v.print();
        }

        else {
            truncate_basis(kevec, kmin, kmax, v);
            check_basis_orthonormality(v.cols(0, kmin - 1));

            // increment croot to optimize next root
            croot++;
        }
    }

    // Set eigval and eigvec to results
    eigval = keval.subvec(0, nroots - 1);
    eigvec = v.cols(0, nroots - 1);

    std::cout << "Final Residual: " << rnorm << std::endl;

    if (rnorm >= rtol)
        error++;  // Calculation did not converge

    return error;
}

arma::vec Davidson_class::get_eigval() {
    return (eigval);
}

arma::mat Davidson_class::get_eigvec() {
    return (eigvec);
}

int Davidson_class::create_initial_guess(arma::mat& v, const size_t kdim) {
    double zero_conv = 1e-12;
    arma::arma_rng::set_seed(1);

    arma::vec guess_vec_one(matrix.n_cols, arma::fill::randu);
    v.col(0) = guess_vec_one;

    v.col(0) = v.col(0) / norm(v.col(0));

    arma::vec guess_vec_two(matrix.n_cols, arma::fill::randu);
    orthonormalize_vector(v.col(0), kdim, guess_vec_two);

    v.col(1) = guess_vec_two;

    for (int i = 2; i < kdim; i++) {
        arma::vec guess_vec_next(matrix.n_cols, arma::fill::randu);
        orthonormalize_vector(v.cols(0, i - 1), kdim, guess_vec_next);
        v.col(i) = guess_vec_next;
    }

    // std::cout << "Initial Guess Vectors" << std::endl;
    // v.cols(0,kdim-1).print();

    // Check orthonormality
    check_basis_orthonormality(v.cols(0, kdim - 1));

    return 0;
}

int Davidson_class::build_kmat(const arma::mat& v, const arma::mat& c, arma::mat& kmat) {
    kmat = v.t() * c;
    return 0;
}

int Davidson_class::get_residual(const arma::vec& keval, const arma::mat& kevec,
                 const arma::mat& v, const arma::mat& c, const int root,
                 const size_t kdim, arma::vec& resid) {
    resid.zeros();
    for (int i = 0; i < kdim; i++) {
        resid +=
            kevec(i, root) * c.col(i) - kevec(i, root) * keval(root) * v.col(i);
    }

    return 0;
}

int Davidson_class::print_iteration_information(const arma::vec& keval, const size_t citer,
                                const double rnorm) {
    std::cout << "Iteration " << citer << " eigenvalues: " << std::endl;
    keval.print();
    std::cout << "Residual: " << rnorm << std::endl;
    return 0;
}

int Davidson_class::make_newvector(const arma::vec& resid,
                   const arma::vec& keval, const int root, arma::vec& newv) {
    for (int I = 0; I < matrix.n_rows; I++) {
        newv(I) = resid(I) / (keval(root) - matrix(I, I) + 1e-7);
    }
    return 0;
}

int Davidson_class::orthonormalize_vector(const arma::mat& v, const size_t kdim,
                          arma::vec& newv) {
    double zero_conv = 1e-12;

    // Check if the matrix is orthonormal

    check_basis_orthonormality(v);

    // Orthonormalize vector
    arma::vec orthnorm_vec = newv;

    for (int i = 0; i < v.n_cols; i++) {
        orthnorm_vec -=
            (arma::dot(newv, v.col(i)) / (arma::dot(v.col(i), v.col(i)))) *
            v.col(i);
    }

    newv = orthnorm_vec / norm(orthnorm_vec);

    // check if new vector is orthonormal to basis

    double newv_norm = norm(newv);
    if (abs(newv_norm - 1) > zero_conv) {
        std::cout << "The orthnonormalized vector is not normal" << std::endl;
        std::cout << "The norm is " << newv_norm << std::endl;
        exit(1);
    }

    arma::vec newv_dot_prod(v.n_cols, arma::fill::ones);
    for (int i = 0; i < v.n_cols; i++) {
        newv_dot_prod(i) = arma::dot(newv, v.col(i));
    }

    if (newv_dot_prod.is_zero(zero_conv) == false) {
        std::cout << "The orthonormalized vector is not orthogonal to the basis"
                  << std::endl;
        std::cout << "The vector has dot products of " << newv_dot_prod
                  << " with the basis" << std::endl;
        exit(1);
    }

    return 0;
}

int Davidson_class::truncate_basis(const arma::mat& kevec, const size_t kmin, const size_t kmax,
                   arma::mat& v) {
    // Truncation during the optimization of the root
    if (kevec.n_cols == v.n_cols) {
        arma::mat prev_v = v;
        v = prev_v * kevec;

        v.cols(kmin, kmax - 1).zeros();
    }

    // Truncation to create the guess vectors for the optimization of the next
    // root
    //(since the optimization of each root may finish before kdim reaches kmax)
    else {
        arma::mat prev_v = v.cols(0, kevec.n_cols - 1);
        v.cols(0, kevec.n_cols - 1) = prev_v * kevec;

        v.cols(kmin, kmax - 1).zeros();
    }

    return 0;
}

int Davidson_class::check_basis_orthonormality(const arma::mat& basis) {
    double zero_conv = 1e-12;
    for (int i = 0; i < basis.n_cols; i++) {
        if (abs(norm(basis.col(i)) - 1) > zero_conv) {
            std::cout << "The input basis is not normal" << std::endl;
            std::cout << "vector " << i << " has a norm of "
                      << norm(basis.col(i)) << std::endl;
            exit(1);
        }
    }

    for (int i = 0; i < basis.n_cols; i++) {
        for (int j = 0; j < i; j++) {
            double dot_prod = arma::dot(basis.col(i), basis.col(j));
            if (abs(dot_prod) > zero_conv) {
                std::cout << "The input basis is not orthogonal" << std::endl;
                std::cout << "The vectors " << i << " and " << j
                          << " have a dot product of " << dot_prod << std::endl;
                exit(1);
            }
        }
    }
    return 0;
}