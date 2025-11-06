#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <random>

bool isInside(double x, double y,
    double x1, double y1, double r1,
    double x2, double y2, double r2,
    double x3, double y3, double r3) {

    return (x - x1)*(x - x1) + (y - y1)*(y - y1) <= r1*r1
    && (x - x2)*(x - x2) + (y - y2)*(y - y2) <= r2*r2
    && (x - x3)*(x - x3) + (y - y3)*(y - y3) <= r3*r3;
}

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());

    double x1, x2, x3, y1, y2, y3, r1, r2, r3;
    std::cin >> x1 >> y1 >> r1 >> x2 >> y2 >> r2 >> x3 >> y3 >> r3;;

    double a2 = std::max(x1 + r1, std::max(x2 + r2, x3 + r3));
    double a1 = std::max(x1 - r1, std::max(x2 - r2, x3 - r3));
    double b2 = std::max(y1 + r1, std::max(y2 + r2, y3 + r3));
    double b1 = std::max(y1 - r1, std::max(y2 - r2, y3 - r3));

    std::uniform_real_distribution<double> dist1(a1, a2);
    std::uniform_real_distribution<double> dist2(b1, b2);
    int inside = 0;
    int n = 1000000;
    for (int i = 0; i < n; i++) {
        double x = dist1(gen);
        double y = dist2(gen);

        if (isInside(x, y, x1, y1, r1, x2, y2, r2, x3, y3, r3)) {
            inside++;
        }
    }

    double s_rect = (a2 - a1) * (b2 - b1);
    double a = static_cast<double>(inside) / static_cast<double>(n); ;

    double result = s_rect * a;

    std::cout << result;
}
