#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <chrono>

using namespace std;
struct LookupElement {
    double k;
    double ln_k;
};

std::vector<LookupElement> createLookupTable() {
    std::vector<LookupElement> table;
    int numEntries = 52;
    table.push_back({pow(2, 8), std::log(pow(2, 8))});
    table.push_back({pow(2, 4), std::log(pow(2, 4))});
    table.push_back({pow(2, 2), std::log(pow(2, 2))});
    table.push_back({pow(2, 1), std::log(pow(2, 1))});

    for (int i = 4; i < numEntries; ++i) {
        double k = pow(2, -(i - 3)) + 1;
        table.push_back({k, std::log(k)});
    }
    return table;
}

double findLargestK(double x, const std::vector<LookupElement>& lookupTable, bool useDivision) {
    for (const auto& element : lookupTable) {
        if (useDivision && x / element.k < 1) {
            return element.k;
        } else if (!useDivision && x * element.k < 1) {
            return element.k;
        }
    }
    return 1;
}

struct InvariantElement {
    double xj;
    double yj;
    double k;
    double ln_k;
    double kxj;
    double ylnk;
};

std::vector<InvariantElement> createInvariantTable(double x, const std::vector<LookupElement>& lookupTable) {
    std::vector<InvariantElement> invariantTable;
    double xj = x / findLargestK(x, lookupTable, true);
    double yj = std::log(findLargestK(x, lookupTable, true));
    double prev_yj = yj;
    for (const auto& element : lookupTable) {
        double k = findLargestK(xj, lookupTable, false);
        double ln_k = std::log(k);
        double kxj = k * xj;
        double ylnk = yj - ln_k;
        invariantTable.push_back({xj, yj, k, ln_k, kxj, ylnk});
        xj = kxj;
        yj -= ln_k;
        if (std::abs(prev_yj - yj) < 1e-2) {
            break;
        }
        prev_yj = yj;
    }

    return invariantTable;
}

double calculateLn(double x, const std::vector<InvariantElement>& invariantTable) {
    double residuo = abs(1 - invariantTable.back().kxj);
    double yj1 = invariantTable.back().yj - invariantTable.back().ln_k;
    return yj1 - residuo;
}
int main() {
    double x = 100;
    std::vector<LookupElement> lookupTable = createLookupTable();
    std::vector<InvariantElement> invariantTable = createInvariantTable(x, lookupTable);
    auto start = std::chrono::high_resolution_clock::now();
    double ln_x = calculateLn(x, invariantTable);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    /*for (auto a : lookupTable){
        std::cout << a.k << " - " << a.ln_k << std::endl;
    }*/
    for (auto a : invariantTable){
        std::cout << a.xj << " - " << a.yj << " - " << a.k << " - " << a.ln_k << " - " << a.kxj << " - " << a.ylnk << std::endl;
    }
    return 0;
}
