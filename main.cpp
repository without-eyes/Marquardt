#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

double f(vector<double> x) {
    double x1 = x[0];
    double x2 = x[1];
    return 5 * pow(x1, 2) + 2 * x1 * x2 + 2 * pow(x2, 2) - x1 + x2;
}

vector<double> gradient(vector<double> x) {
    double x1 = x[0];
    double x2 = x[1];
    vector<double> grad(2);
    grad[0] = 10 * x1 + 2 * x2 - 1;
    grad[1] = 2 * x1 + 4 * x2 + 1;
    return grad;
}

vector<vector<double>> hessian() {
    vector<vector<double>> hess(2, vector<double>(2));
    hess[0][0] = 10;
    hess[0][1] = 2;
    hess[1][0] = 2;
    hess[1][1] = 4;
    return hess;
}

vector<double> marquardt(vector<double> x0, double epsilon) {
    double lambda = 10000;
    vector<double> x = x0;
    vector<double> grad = gradient(x);
    int iter = 0;

    cout << "Iteration №" << iter
         << ":\nx" << iter << " = [" << fixed << setprecision(4) << x[0] << ", " << x[1] << "]^T\n" <<
         "a" << iter << " = " << lambda <<
         "f(x) = " << f(x) <<
         endl << endl;

    while (sqrt(grad[0]*grad[0] + grad[1]*grad[1]) > epsilon) {
        iter++;
        vector<vector<double>> hess = hessian();
        for (int i = 0; i < 2; ++i) {
            hess[i][i] += lambda;
        }
        vector<double> delta_x(2);
        delta_x[0] = -1 * (hess[0][0]*grad[0] + hess[0][1]*grad[1]) / (hess[0][0]*hess[1][1] - hess[0][1]*hess[1][0]);
        delta_x[1] = -1 * (hess[1][0]*grad[0] + hess[1][1]*grad[1]) / (hess[0][0]*hess[1][1] - hess[0][1]*hess[1][0]);

        vector<double> new_x(2);
        new_x[0] = x[0] + delta_x[0];
        new_x[1] = x[1] + delta_x[1];

        cout << "Iteration №" << iter << ":" << endl <<
             "x" << iter-1 << " = [" << x[0] << ", " << x[1] << "]^T" << endl <<
             "a" << iter-1 << " = " << lambda << endl <<
             "delf" " = [" << grad[0] << ", " << grad[1] << "]" << endl <<
             "||delf||" " = " << sqrt(grad[0]*grad[0] + grad[1]*grad[1]) << endl <<
             "x" << iter << " = [" << new_x[0] << ", " << new_x[1] << "]" << endl <<
             "f" << iter << " = " << f(new_x) << endl;

        if (f(new_x) < f(x)) {
            lambda /= 2;
            x = new_x;
            grad = gradient(x);
        } else {
            lambda *= 2;
        }
        cout << "f" << iter << " < f" << iter-1 << " -> " << "c" << iter << " = " << lambda << endl << endl;
    }

    return x;
}

int main() {
    vector<double> initial_guess = {5, 2};
    double epsilon = 1e-1;

    cout << "Initial guess: x = (" << fixed << setprecision(3) << initial_guess[0] << ", " << initial_guess[1] << "), f(x) = " << f(initial_guess) << endl;

    vector<double> solution = marquardt(initial_guess, epsilon);

    cout << "Minimum point: x = (" << fixed << setprecision(3) << solution[0] << ", " << solution[1] << "), Minimum value: " << f(solution) << endl;

    return 0;
}
