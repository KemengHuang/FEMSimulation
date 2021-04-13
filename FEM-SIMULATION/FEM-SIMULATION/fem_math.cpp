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

double f(double x, double a, double b, double c, double d) {
	double f = a * x * x * x + b * x * x + c * x + d;
	return f;
}

double df(double x, double a, double b, double c) {
	double df = 3 * a * x * x + 2 * b * x + c;
	return df;
}

std::vector<double> NewtonSolverForCubicEquation(const double& a, const double& b, const double& c, const double& d)
{
	double EPS = 1e-6;
	double DX = 0;
	std::vector<double> results;
	double specialPoint = -b / a / 3;
	double pos[2];
	int solves = 1;
	double delta = 4 * b * b - 12 * a * c;
	if (delta > 0) {
		pos[0] = (sqrt(delta) - 2 * b) / 6 / a;
		pos[1] = (-sqrt(delta) - 2 * b) / 6 / a;
		double v1 = f(pos[0], a, b, c, d);
		double v2 = f(pos[1], a, b, c, d);
		if (abs(v1) < 1e-10) {
			v1 = 0;
		}
		if (abs(v2) < 1e-10) {
			v2 = 0;
		}
		double sign = v1 * v2;
		DX = (pos[0] - pos[1]);
		if (sign <= 0) {
			solves = 3;
		}
		else if (sign > 0) {
			if ((a < 0 && f(pos[0], a, b, c, d) > 0) || (a > 0 && f(pos[0], a, b, c, d) < 0)) {
				DX = -DX;
			}
		}
	}
	else if (delta == 0) {
		if (abs(f(specialPoint, a, b, c, d)) < 1e-10) {
			for (int i = 0; i < 3; i++) {
				double tempReuslt = specialPoint;
				results.push_back(tempReuslt);
			}
			return results;
		}
		if (a > 0) {
			if (f(specialPoint, a, b, c, d) > 0) {
				DX = 1;
			}
			else if (f(specialPoint, a, b, c, d) < 0) {
				DX = -1;
			}
		}
		else if (a < 0) {
			if (f(specialPoint, a, b, c, d) > 0) {
				DX = -1;
			}
			else if (f(specialPoint, a, b, c, d) < 0) {
				DX = 1;
			}
		}

	}

	double start = specialPoint - DX;
	double x0 = start;
	double result[3];

	for (int i = 0; i < solves; i++) {
		double x1 = 0;
		int itCount = 0;
		do
		{
			if (itCount)
				x0 = x1;

			x1 = x0 - ((f(x0, a, b, c, d)) / (df(x0, a, b, c)));
			itCount++;

		} while (abs(x1 - x0) > EPS);
		results.push_back(x1);
		start = start + DX;
		x0 = start;
	}

	return results;
}