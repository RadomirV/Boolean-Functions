#include <chrono>
#include <thread>
#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <unordered_set>

#include "BF.h"

using namespace bf;

typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
typedef int int32_t;
// typedef unsigned long uint64_t;

int main()
{
    srand(time(NULL));
    const int numberOfVariables = 5;
    unsigned int GoodPairsCount = 0;
    unsigned int BadPairsCount = 0;
    int iterationCount = 0;
    auto limits = BF::generateBorderBalancesFunctions();
    // auto limits = std::make_pair(BF("0000000011111111"), BF("1111111100000000"));
    // auto limits = std::make_pair(BF("00001111"), BF("11110000"));
    // auto limits = std::make_pair(BF("0011"), BF("1100"));
    BF func = limits.first;
    std::vector<unsigned int> good_freq;
    std::vector<unsigned int> bad_freq;
    good_freq.resize(300);
    bad_freq.resize(300);
    while (func != limits.second)
    {
        auto wht_coef = func.WHTransform();

        auto badPairs = func.pairsToWorsen();      // Pairs which decrease the nonlinearity
        auto improvePairs = func.pairsToImprove(); // Pairs which improve nonlinearity
        GoodPairsCount += improvePairs.size();
        BadPairsCount += badPairs.size();
        good_freq[improvePairs.size()]++;
        bad_freq[badPairs.size()]++;

        if (iterationCount % 100000 == 0)
        {
            uint32_t transfer = 8;
            uint32_t inRow = 0;
            func.printPerBit();
            std::cout << "\t\t\t\t\t\tIterN = " << iterationCount << "\nGood = " << improvePairs.size() << " Bad = " << badPairs.size() << " nonlinearity = " << func.nonlinearity(wht_coef) << "\n";
            for (int i = 0; i < 300; ++i)
            {
                if (good_freq[i] != 0)
                {
                    std::cout << "[" << i << "]= " << good_freq[i] << " ";
                    inRow++;
                    if (inRow % transfer == 0)
                        std::cout << "\n";
                }
            }
            inRow = 0;
            std::cout << "\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" << std::endl;
            for (int i = 0; i < 300; ++i)
            {
                if (bad_freq[i] != 0)
                {
                    std::cout << "[" << i << "]= " << bad_freq[i] << " ";
                    inRow++;
                    if (inRow % transfer == 0)
                        std::cout << "\n";
                }
            }
            std::cout << std::endl;
        }
        if (improvePairs.size() == 0 && badPairs.size() == 0)
        {
            func.printPerBit();
            std::cout << "NO PAIRS AT ALL GOOD AND BAD NO PAIRS AT ALL GOOD AND BAD NO PAIRS AT ALL GOOD AND BAD\n";
            std::cout << "\t\t\t\t\t\tIterN = " << iterationCount << "\nGood = " << improvePairs.size() << " Bad = " << badPairs.size() << " nonlinearity = " << func.nonlinearity(wht_coef) << "\n";
        }
        if (func.nonlinearity(wht_coef) > 12)
        {
            std::cout << "Good = " << improvePairs.size() << " Bad = " << badPairs.size() << " nonlinearity = " << func.nonlinearity(wht_coef) << "\t\t\t\t\t\tIterN = " << iterationCount << "\n";
        }

        func.nextBalanced();
        iterationCount++;
    }
    auto wht_coef = func.WHTransform();

    auto badPairs = func.pairsToWorsen();      // Pairs which decrease the nonlinearity
    auto improvePairs = func.pairsToImprove(); // Pairs which improve nonlinearity

    GoodPairsCount += improvePairs.size();
    BadPairsCount += badPairs.size();
    func.print();
    func.printPerBit();
    std::cout << "\t\t\t\t\t\tIterN = " << iterationCount << "\nGood = " << improvePairs.size() << " Bad = " << badPairs.size() << "\n";

    std::cout << "GoodPairsCount= " << GoodPairsCount << "\n";
    std::cout << "BadPairsCount= " << BadPairsCount << "\n";
    uint32_t transfer = 8;
    uint32_t inRow = 0;
    good_freq[improvePairs.size()]++;
    bad_freq[badPairs.size()]++;
    std::cout << "Distribution of the powers of sets of improving pairs" << std::endl;
    for (int i = 0; i < 300; ++i)
    {
        if (good_freq[i] != 0)
        {
            std::cout << "[" << i << "]= " << good_freq[i] << " ";
            inRow++;
            if (inRow % transfer == 0)
                std::cout << "\n";
        }
    }
    inRow = 0;
    std::cout << "distribution of the powers of sets of worsening pairs" << std::endl;
    for (int i = 0; i < 300; ++i)
    {
        if (bad_freq[i] != 0)
        {
            std::cout << "[" << i << "]= " << bad_freq[i] << " ";
            inRow++;
            if (inRow % transfer == 0)
                std::cout << "\n";
        }
    }

    return 0;
}
