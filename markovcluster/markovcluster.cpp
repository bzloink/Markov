 /* Markov clustering algorithm
    by Andrew M. Launder
    
    See README for proper code usage. */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

unsigned long int Dim(char* filename) {
 // Determine number of lines in input file (i.e., the order of the square transition matrix).
    std::ifstream data(filename);
    
    unsigned long int dim = 0;
    
    if (!data) {
        dim = 1;
        return dim;
    }
    
    std::string line;
    
    while (std::getline(data, line)) {
        ++dim;
    }
    
    return dim;
}

std::vector<double> LineSplit(std::ifstream& data) {
 // Splits a string and stores it in a vector of doubles.
    std::string line, item;
    std::getline(data, line);
    std::istringstream linestream(line);
    std::vector<std::string>* linesplitstr = new std::vector<std::string>;
    
    while (linestream >> item) {
        linesplitstr->push_back(item);
    }
    
    unsigned long int linelength = linesplitstr->size();
    std::vector<double> linesplitdouble(linelength);
    
    unsigned long int i;
    for (i = 0; i < linelength; ++i) {
        linesplitdouble[i] = std::stod((*linesplitstr)[i]);
    }
    delete linesplitstr;
    
    return linesplitdouble;
}

std::vector<std::vector<double>> ParseFile(char* filename, unsigned long int dim) {
 // Parses input into transition matrix object.
    std::ifstream data(filename);
    
    if (!data) {
        std::vector<std::vector<double>> errorreturn(1);
        errorreturn[0].push_back(0);
        return errorreturn;
    }
    
    std::vector<double>* valueline = new std::vector<double>(dim);
    std::vector<std::vector<double>> values(dim);
    
    unsigned long int i;
    for (i = 0; i < dim; ++i) {
        *valueline = LineSplit(data);
        values[i] = *valueline;
    }
    delete valueline;
    
    return values;
}

void PrintMatrix(std::vector<std::vector<double>> matrix, unsigned long int dim) {
 // Prints a square matrix (vector of vectors of doubles).
    std::cout << std::endl;
    
    unsigned long int i, j;
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int MatrixCheck(std::vector<std::vector<double>> matrix, unsigned long int dim) {
 /* Returns 2 if input file has been specified incorrectly.
    Returns 1 if matrix is not square.
    Returns 0 otherwise.                                    */
    int matrixcheck = 0;
    
    if (dim == 1) {
        matrixcheck = 2;
        return matrixcheck;
    }
    
    unsigned long int i;
    for (i = 0; i < dim; ++i) {
        if (matrix[i].size() != dim) {
            matrixcheck = 1;
        }
    }
    
    return matrixcheck;
}

std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>> initialmatrix, unsigned long int dim) {
 // Transposes a square matrix (vector of vectors of doubles).
    decltype(initialmatrix) finalmatrix(dim);
    
    unsigned long int i, j;
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            finalmatrix[i].push_back(initialmatrix[j][i]);
        }
    }
    
    return finalmatrix;
}

std::vector<std::vector<double>> DivideSum(std::vector<std::vector<double>> initialmatrix, unsigned long int dim) {
 // Makes a matrix left stochastic.
    decltype(initialmatrix) finalmatrix(dim);
    
    std::vector<double>* column = new std::vector<double>(dim);
    
    unsigned long int i, j;
    for (i = 0; i < dim; ++i) {
        (*column)[i] = 0;
        for (j = 0; j < dim; ++j) {
            if (initialmatrix[j][i] == 0) {
                continue;
            }
            (*column)[i] += initialmatrix[j][i];
        }
    }
    
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            if ((*column)[j] == 0) {
                finalmatrix[i].push_back(0);
                continue;
            }
            finalmatrix[i].push_back(initialmatrix[i][j] / (*column)[j]);
        }
    }
    
    return finalmatrix;
}

std::vector<std::vector<double>> Expand(std::vector<std::vector<double>> initialmatrix, unsigned long int dim) {
 // Squares a matrix via a standard matrix product.
    decltype(initialmatrix) finalmatrix(dim);
    double element = 0;
    
    unsigned long int i, j, k;
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            for (k = 0; k < dim; ++k) {
                if (initialmatrix[i][k] == 0 && initialmatrix[k][j] == 0) {
                    continue;
                }
                element += initialmatrix[i][k] * initialmatrix[k][j];
            }
            finalmatrix[i].push_back(element);
            element = 0;
        }
    }
    
    return finalmatrix;
}

std::vector<std::vector<double>> Inflate(std::vector<std::vector<double>> initialmatrix, unsigned long int dim, double power) {
 // Returns the power of a matrix via a Hadamard matrix product.
    decltype(initialmatrix) finalmatrix(dim);
    
    unsigned long int i, j;
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            if (initialmatrix[i][j] == 0) {
                finalmatrix[i].push_back(0);
                continue;
            }
            finalmatrix[i].push_back(pow(initialmatrix[i][j], power));
        }
    }
    
    return finalmatrix;
}

std::vector<std::vector<double>> MakeZero(std::vector<std::vector<double>> initialmatrix, unsigned long int dim, double zeroconv) {
 // Sets any element less than a given threshold to 0.
    decltype(initialmatrix) finalmatrix(dim);
    
    unsigned long int i, j;
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            if (initialmatrix[i][j] < zeroconv) {
                finalmatrix[i].push_back(0);
                continue;
            }
            finalmatrix[i].push_back(initialmatrix[i][j]);
        }
    }
    
    return finalmatrix;
}

int ConvCheck(std::vector<std::vector<double>> newmatrix, std::vector<std::vector<double>> oldmatrix, double conv, unsigned long int dim) {
 /* Returns 1 if RMS of transition matrix elements not yet converged.
    Returns 0 otherwise.                                              */
    int convcheck = 1;
    double RMS = 0;
    
    unsigned long int i, j;
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            if (newmatrix[i][j] == 0 && oldmatrix[i][j] == 0) {
                continue;
            }
            RMS += pow(newmatrix[i][j] - oldmatrix[i][j], 2);
        }
    }
    RMS = pow(RMS, 0.5);
    
    if (RMS < conv) {
        convcheck = 0;
    }
    
    return convcheck;
}

int main(int argc, char** argv) {
 // Markov clustering algorithm.
    double conv, inflation, zeroconv;
    int convexp, maxiter, zeroconvexp;
    char* filename;
    
    if (argc <= 2) {
        convexp = -6;
        maxiter = 100;
        zeroconvexp = -15;
        conv = pow(10, convexp);
        inflation = 2;
        zeroconv = pow(10, zeroconvexp);
        filename = argv[1];
    }
    else {
        if (!argv[1]) {
            std::cerr << "No convergence criterion given!" << std::endl;
            return EXIT_FAILURE;
        }
        std::istringstream convstr(argv[1]);
        convstr >> convexp;
        if (convexp < 0) {
            std::cerr << "Convergence criterion must be greater than 0." << std::endl;
            return EXIT_FAILURE;
        }
        convexp = convexp * -1;
        conv = pow(10, convexp);
        
        if (!argv[2]) {
            std::cerr << "No max iterations limit given!" << std::endl;
            return EXIT_FAILURE;
        }
        std::istringstream maxiterstr(argv[2]);
        maxiterstr >> maxiter;
        if (maxiter <= 0 || maxiter >= 1000) {
            std::cerr << "Max iterations limit must be greater than 0 and less than 1000." << std::endl;
            return EXIT_FAILURE;
        }
        
        if (!argv[3]) {
            std::cerr << "No inflation coefficient given!" << std::endl;
            return EXIT_FAILURE;
        }
        std::istringstream inflationstr(argv[3]);
        inflationstr >> inflation;
        if (inflation <= 1) {
            std::cerr << "Inflation coefficient must be greater than 1." << std::endl;
            return EXIT_FAILURE;
        }
        
        if (!argv[4]) {
            std::cerr << "No zero threshold given!" << std::endl;
            return EXIT_FAILURE;
        }
        std::istringstream zeroexpstr(argv[4]);
        zeroexpstr >> zeroconvexp;
        if (zeroconvexp <= 0) {
            std::cerr << "Zero threshold must be greater than 0." << std::endl;
            return EXIT_FAILURE;
        }
        zeroconvexp = zeroconvexp * -1;
        zeroconv = pow(10, zeroconvexp);
        
        filename = argv[5];
    }
    
    if (!filename) {
        std::cerr << "No input file given!" << std::endl;
        return EXIT_FAILURE;
    }
    
    unsigned long int dim = Dim(filename);
    std::vector<std::vector<double>>* initialtransmat = new std::vector<std::vector<double>>;
    *initialtransmat = ParseFile(filename, dim);
    
    int matrixcheck = MatrixCheck(*initialtransmat, dim);
    if (matrixcheck == 2) {
        std::cerr << "Have you correctly specified your input file?" << std::endl;
        return EXIT_FAILURE;
    }
    if (matrixcheck == 1) {
        std::cerr << "Check your transition matrix. Is it square?" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::vector<std::vector<double>>* oldtransmat = new std::vector<std::vector<double>>;
    *oldtransmat = Transpose(*initialtransmat, dim);
    
    std::cout << std::endl;
    std::cout << "Inital transition matrix:" << std::endl;
    PrintMatrix(*initialtransmat, dim);
    delete initialtransmat;
    
    *oldtransmat = DivideSum(*oldtransmat, dim);
    
    *oldtransmat = Expand(*oldtransmat, dim);
    
    *oldtransmat = Inflate(*oldtransmat, dim, inflation);
    *oldtransmat = DivideSum(*oldtransmat, dim);
    
    *oldtransmat = MakeZero(*oldtransmat, dim, zeroconv);
    
    std::vector<std::vector<double>>* newtransmat = new std::vector<std::vector<double>>;
    *newtransmat = *oldtransmat;
    double convcheck = 1;
    int iter = 1;
    while (convcheck == 1 && iter <= maxiter) {
        *newtransmat = Expand(*newtransmat, dim);
        
        *newtransmat = Inflate(*newtransmat, dim, inflation);
        *newtransmat = DivideSum(*newtransmat, dim);
        
        *newtransmat = MakeZero(*newtransmat, dim, zeroconv);
        
        convcheck = ConvCheck(*newtransmat, *oldtransmat, conv, dim);
        ++iter;
        *oldtransmat = *newtransmat;
    }
    delete oldtransmat;
    
    if (convcheck == 1) {
        std::cerr << "Max iterations reached!" << std::endl;
        std::cerr << std::endl;
        std::cerr << "No convergence!" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Final transition matrix:" << std::endl;
        PrintMatrix(*newtransmat, dim);
        std::cerr << "Better luck next time!" << std::endl;
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }
    
    std::cout << "Converged in " << iter << " steps! Congratulations!" << std::endl;
    std::cout << std::endl;
    std::cout << "Final transition matrix:" << std::endl;
    PrintMatrix(*newtransmat, dim);
    delete newtransmat;
    std::cout << "Program exiting successfully." << std::endl;
    std::cout << std::endl;
    
    return 0;
}
