#include <stdio.h>
#include <stdlib.h>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "davidson_class.h"

template <class Matrix>
void print_matrix(Matrix matrix) {
    matrix.print(std::cout);
}

// provide explicit instantiations of the template function for
// every matrix type you use somewhere in your program.
template void print_matrix<arma::mat>(arma::mat matrix);
template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);

int main() {
    // Create an example matrix

    arma::mat example_matrix(100, 100, arma::fill::zeros);

    srand(1);

    for (int i = 0; i < 100; i++) {
        int random = rand() % 20;
        example_matrix(i, i) = random;
        for (int n = 0; n < 100; n++) {
            int random_2 = rand() % 20;
            if (random_2 == 3) {
                double random_3 = rand() % 200;
                example_matrix(i, n) = example_matrix(n, i) = random_3 / 100;
            }
        }
    }

    // getting roots by calling davidson just once
    int nroot = 5;

    Davidson_class Davidson_object(example_matrix);

    int error = Davidson_object.diagonalize_with_davidson(nroot);

    arma::vec eval; 
    arma::mat evec;
    eval = Davidson_object.get_eigval();
    evec = Davidson_object.get_eigvec();
    eval.print();

    if (error != 0) {
        std::cout << "Error code: " << error << std::endl;
        exit(error);
    }

    std::cout << "Lowest " << nroot << " eigenvalues are:" << std::endl;
    eval.print();

    std::cout << "Lowest " << nroot << "eigenvectors are:" << std::endl;
    evec.print();

    // Diagonalize with armadillo to check answer
    arma::vec arma_eigvals;
    arma::mat arma_eigvecs;

    arma::eig_sym(arma_eigvals, arma_eigvecs, example_matrix);

    std::cout << "Armadillo lowest eigenvalues: " << std::endl;
    arma_eigvals.subvec(0, nroot - 1).print();

    std::cout << "Armadillo lowest eigenvectors: " << std::endl;
    arma_eigvecs.cols(0, nroot - 1).print();

    // print difference in eigenvalues and eigenvectors
    std::cout << "Difference between Davidson and Armadillo eigenvalues: "
              << std::endl;
    (abs(arma_eigvals.subvec(0, nroot - 1)) - abs(eval)).print();

    std::cout << "Difference between Davidson and Armadillo eigenvectors: "
              << std::endl;
    (abs(arma_eigvecs.cols(0, nroot - 1)) - abs(evec)).print();

    return 0;
}