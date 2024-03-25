#ifndef davidson_template_h
#define davidson_template_h

#include <armadillo>

//Davidson Class
class Davidson_class {
    arma::mat matrix;  //The matrix
    arma::vec eigval;  //eigenvalues
    arma::mat eigvec;  //eigenvectors

    /*
    * create_initial_guess: Creates the initial guess vectors
    * Input:
    *  matrix = symmetric matrix
    *  kdim = dimensions of subspace (number of initial guess vectors)
    * Output:
    *  v = Guess vectors
    * Returns:
    *  No return
    */
    int create_initial_guess(arma::mat& v, const size_t kdim);

    /*
    * build_kmat: Build subspace matrix
    * Input:
    *  v = Guess vectors
    *  c = sigmoidal vectors
    * Output:
    *  kmat = subspace matrix
    * Returns:
    *  No return
    */
    int build_kmat(const arma::mat& v, const arma::mat& c, arma::mat& kmat);

    /*
    * get_residual: Get residual vector
    * Input:
    *  keval = subspace eigenvalues
    *  kevec = subspace eigenvectors
    *  v = Guess vectors
    *  c = sigmoidal vectors
    *  root = optimized root
    * Output:
    *  resid = residual vector
    * Returns:
    *  No return
    */
    int get_residual(const arma::vec& keval, const arma::mat& kevec,
                    const arma::mat& v, const arma::mat& c, const int root,
                    const size_t kdim, arma::vec& resid);
    /*
    * print_iteration_information: Print iteration information
    * Input:
    *  keval = subspace eigenvalues
    *  citer = iteration number
    *  rnorm = norm of residual vector
    * Output:
    *  No output
    * Returns:
    *  No return
    */
    int print_iteration_information(const arma::vec& keval, const size_t citer,
                                    const double rnorm);

    /*
    * make_newvector: Make correction vector to be added to subspace
    * Input:
    *  matrix = symmetric matrix
    *  resid = residual vector
    *  keval = subspace eigenvalues
    * Output:
    *  newv = correction vector to be added to subspace
    * Returns:
    *  No return
    */
    int make_newvector(const arma::vec& resid,
                    const arma::vec& keval, const int root, arma::vec& newv);

    /*
    * orthonormalize_vector: orthonormalize a vector to an orthogonal basis
    * Input:
    *  v = orthonormal basis
    *  kdim = number of vectors in the basis
    * Output:
    *  newv = orthonormalized vector
    * Returns:
    *  No return
    */
    int orthonormalize_vector(const arma::mat& v, const size_t kdim,
                            arma::vec& newv);

    /*
    * truncate_basis: Get linear combination of guess vectors to truncate basis
    * Input:
    *  kevec = subspace eigenvectors
    *  kmin = smallest dimensions for subspace
    *  kmax = larges dimensions for subspace
    * Output:
    *  v = new truncated guess vectors
    * Returns:
    *  exits program with message if orthonormalized vector is found to be non-orthonormal
    */
    int truncate_basis(const arma::mat& kevec, const size_t kmin, const size_t kmax,
                    arma::mat& v);
    /*
    * check_basis_orthonormality: Checks orthonormality of basis
    * Input:
    *  basis = basis
    * Output:
    *  No output
    * Returns:
    *  exits program with message if basis is found to be non-orthonormal
    */
    int check_basis_orthonormality(const arma::mat& basis);

    public:
    Davidson_class(arma::mat m){
        matrix = m;
    }
    /*
    * diagonalize_with_davidson: Diagonalize a symmetric matrix using the Davidson
    * algorithm.
    * Input:
    *  matrix = symmetric matrix
    *  nroots = number of eigenvalues to find
    * Output:
    *  eigval = eigenvalues
    *  eigvec = eigenvectors
    * Returns:
    *  error = 0 success
    */
    int diagonalize_with_davidson(const size_t nroots);
    /*
    * get_eigval: get eigenvalues
    * Input:
    *  none
    * Output:
    *  none
    * Returns:
    *  armadillo vector of eigenvalues
    */
    arma::vec get_eigval();
    /*
    * get_eigvec: get eigenvectors
    * Input:
    *  none
    * Output:
    *  none
    * Returns:
    *  armadillo matrix of eigenvectors
    */
    arma::mat get_eigvec();
};


#endif
