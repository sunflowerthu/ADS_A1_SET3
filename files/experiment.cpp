#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <random>
#include <fstream>
#include <cmath>

bool isInside(double x, double y,
    double x1, double y1, double r1,
    double x2, double y2, double r2,
    double x3, double y3, double r3) {

    return (x - x1)*(x - x1) + (y - y1)*(y - y1) <= r1*r1
        && (x - x2)*(x - x2) + (y - y2)*(y - y2) <= r2*r2
        && (x - x3)*(x - x3) + (y - y3)*(y - y3) <= r3*r3;
}

double computeExactArea() {
    return 0.25 * std::numbers::pi + 1.25 * asin(0.8) - 1.0;
}

// Широкая область
double wideMonteCarlo(int numPoints) {
    double x1 = 1.0, y1 = 1.0, r1 = 1.0;
    double x2 = 1.5, y2 = 2.0, r2 = sqrt(5.0) / 2.0;
    double x3 = 2.0, y3 = 1.5, r3 = sqrt(5.0) / 2.0;

    double a2 = std::max(x1 + r1, std::max(x2 + r2, x3 + r3));
    double a1 = std::min(x1 - r1, std::min(x2 - r2, x3 - r3));
    double b2 = std::max(y1 + r1, std::max(y2 + r2, y3 + r3));
    double b1 = std::min(y1 - r1, std::min(y2 - r2, y3 - r3));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_x(a1, a2);
    std::uniform_real_distribution<double> dist_y(b1, b2);

    int inside = 0;
    for (int i = 0; i < numPoints; i++) {
        double x = dist_x(gen);
        double y = dist_y(gen);

        if (isInside(x, y, x1, y1, r1, x2, y2, r2, x3, y3, r3)) {
            inside++;
        }
    }

    double rect_area = (a2 - a1) * (b2 - b1);
    return rect_area * static_cast<double>(inside) / static_cast<double>(numPoints);
}

// Узкая область
double narrowMonteCarlo(int numPoints) {
    double x1 = 1.0, y1 = 1.0, r1 = 1.0;
    double x2 = 1.5, y2 = 2.0, r2 = sqrt(5.0) / 2.0;
    double x3 = 2.0, y3 = 1.5, r3 = sqrt(5.0) / 2.0;

    // Из рисунка видно, что пересечение находится примерно в [1.0, 2.0] x [1.0, 2.0]
    double a1 = 0.8, a2 = 2.2;  // x-границы
    double b1 = 0.8, b2 = 2.2;  // y-границы

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_x(a1, a2);
    std::uniform_real_distribution<double> dist_y(b1, b2);

    int inside = 0;
    for (int i = 0; i < numPoints; i++) {
        double x = dist_x(gen);
        double y = dist_y(gen);

        if (isInside(x, y, x1, y1, r1, x2, y2, r2, x3, y3, r3)) {
            inside++;
        }
    }

    double rect_area = (a2 - a1) * (b2 - b1);
    return rect_area * static_cast<double>(inside) / static_cast<double>(numPoints);
}

int main() {
    // Вычисляем точную площадь
    double exact_area = computeExactArea();

    // Записываем данные экспериментов
    std::ofstream dataFile("experiment_data.csv");
    dataFile << "rectangle_type,num_points,computed_area,relative_error\n";

    for (int n = 100; n <= 100000; n += 500) {
        // Широкая область
        double wide_area = wideMonteCarlo(n);
        double wide_error = std::abs(wide_area - exact_area) / exact_area;

        // Узкая область
        double narrow_area = narrowMonteCarlo(n);
        double narrow_error = std::abs(narrow_area - exact_area) / exact_area;

        dataFile << "wide," << n << "," << wide_area << "," << wide_error << "\n";
        dataFile << "narrow," << n << "," << narrow_area << "," << narrow_error << "\n";
    }

    dataFile.close();

    return 0;
}
