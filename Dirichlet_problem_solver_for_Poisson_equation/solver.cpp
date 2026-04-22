#include "inmost.h"
#include <stdio.h>

#include <iostream>
#include <fstream>

// u = sin(3x)cos(4y) - точное решение
double u(double x, double y) { return sin(3 * x) * cos(4 * y); } 
// f(x) = 25sin(3x)cos(4y) - производная, посчитанная аналитически
double f(double x, double y) { return 25 * sin(3 * x) * cos(4 * y);}

// заполнение матрицы и правой части для системы, умниженной на h^2
void fill_system(INMOST::Sparse::Matrix &A, INMOST::Sparse::Vector &b, size_t N_inner, double h) {
    auto idx = [N_inner](size_t i, size_t j) {
        return (i-1) + (j-1) * N_inner;
    };

    for (size_t j = 1; j <= N_inner; ++j){
        for (size_t i = 1; i <= N_inner; ++i){
            double x = i * h;
            double y = j * h;
            size_t row = idx(i, j);
            A[row][row] = 4.0;
            b[row] = f(x, y) * h * h;

            if (i > 1) {
                A[row][idx(i - 1, j)] = -1.0;
            } else {
                b[row] += u(0.0, y);
            }

            if (i < N_inner) {
                A[row][idx(i + 1, j)] = -1.0;
            } else {
                b[row] += u(1.0, y);
            }

            if (j > 1) {
                A[row][idx(i, j - 1)] = -1.0;
            } else {
                b[row] += u(x, 0.0);
            }

            if (j < N_inner) {
                A[row][idx(i, j + 1)] = -1.0;
            } else {
                b[row] += u(x, 1.0);
            }
        }
    }
}

double C_norm(INMOST::Sparse::Vector &solution, size_t N_inner, double h){
    double max_err = 0.0;
    auto idx = [N_inner](size_t i, size_t j) {
        return (i-1) + (j-1) * N_inner;
    };

    for (size_t j = 1; j <= N_inner; ++j){
        for (size_t i = 1; i <= N_inner; ++i){
            size_t row = idx(i, j);
            double x = i * h;
            double y = j * h;

            double err = std::abs(solution[row] - u(x,y));
            if (err > max_err) {
                max_err = err;
            }
        }
    }
    return max_err;
}

double L2_norm(INMOST::Sparse::Vector &solution, size_t N_inner, double h){
    double sum_err = 0.0;
    auto idx = [N_inner](size_t i, size_t j){
        return ((i - 1) + (j-1) * N_inner);
    };

    for (size_t j = 1; j <= N_inner; ++j){
        for (size_t i = 1; i <= N_inner; ++i){
            size_t row = idx(i, j);
            double x = i * h;
            double y = j * h;

            sum_err += (solution[row] - u(x,y)) * (solution[row] - u(x,y));
        }
    }
    return std::sqrt(sum_err * h * h);
}

int main(int argc, char *argv[]){

    unsigned n = 500; // разбер матриццы
    if (argc > 1){
        n = std::atoi(argv[1]);
    }

    double h = 1.0 / n;
    unsigned n_inner = n - 1;
    unsigned n_system = n_inner * n_inner;

    INMOST::Solver::Initialize(&argc, &argv);
    INMOST::Solver S(INMOST::Solver::INNER_ILU2);

    S.SetParameter("absolute_tolerance", "1e-14");
    S.SetParameter("relative_tolerance", "1e-11");

    INMOST::Sparse::Matrix A;
    INMOST::Sparse::Vector b;
    INMOST::Sparse::Vector solution;

    A.SetInterval(0, n_system);
    b.SetInterval(0, n_system);
    solution.SetInterval(0, n_system);

    fill_system(A, b, n_inner, h);

    S.SetMatrix(A);

    bool solved = S.Solve(b, solution);
    std::cout << "Iterations: " << S.Iterations() << std::endl;
    std::cout << "Preconditioner time: " << S.PreconditionerTime() << std::endl;
    std::cout << "Iterations time: " << S.IterationsTime() << std::endl;

    if (!solved){
        std::cout << "Linear solver failure!" << std::endl;
        std::cout << "Reason: " << S.ReturnReason() << std::endl;
    } else {
        double err_C = C_norm(solution, n_inner, h);
        double err_L2 = L2_norm(solution, n_inner, h);

        std::cout << "L2-norm error: " << err_L2 << std::endl;
        std::cout << "C-norm error: " << err_C << std::endl;

        std::ofstream outfile("results.csv", std::ios::app);
        if (outfile.is_open()) {
            outfile << n << "," << h << "," << err_C << "," << err_L2 << "," 
                    << S.PreconditionerTime() << "," << S.IterationsTime() << "\n";
            outfile.close();
        } else {
            std::cout << "Не удалось открыть файл results.csv для записи!" << std::endl;
        }
    }
    
    INMOST::Solver::Finalize();
    return 0;
}
