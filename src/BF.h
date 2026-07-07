#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

class BF
{
private:
    std::vector<unsigned int> vec_;
    int n_;

public:
    BF();
    BF(const std::string &);
    BF(int n_, int par = 0);
    BF(const BF &);
    ~BF() {};

    BF operator=(const BF &);
    bool operator==(const BF &) const;
    bool operator!=(const BF &) const;
    bool operator[](unsigned int) const;

    int getNumberOfVariables() const;
    static bool getValue(const BF &, unsigned int);
    void setValue(unsigned int, bool);
    std::vector<unsigned int> values();
    void print(bool onlyVector = true) const;
    void printPerBit() const;

    static unsigned int weight(const BF &);
    static unsigned int degree(const BF &);
    static unsigned int correlationImmunity(const BF &, std::vector<int> &WHT);
    static unsigned int nonlinearity(const BF &, std::vector<int> &WHT);
    static int cnF(const BF &, std::vector<int> &);
    static unsigned int propagationCriteria(const BF &, std::vector<int> &);

    static void printANF(const BF &);
    static void printBAA(const BF &, std::vector<int> &WHT);

    static BF mobiusTransform(const BF &);
    static std::vector<int> WHTransform(const BF &);
    static std::vector<int> autoCor(const BF &, std::vector<int>);

    BF generateAffine(uint32_t maskOfVariables, bool addOne) const;
    std::pair<uint32_t, uint32_t> generatePair() const;
    std::pair<uint32_t, uint32_t> generateImprovePair(uint32_t) const;
    static std::pair<BF, BF> generateBorderBalancesFunctions();
    static BF genBalancedFunc(int numberOfVariables);
    void nextBalanced();

    static std::pair<std::vector<unsigned int>, std::vector<unsigned int>> goodPairsVec(const BF &, std::vector<int> &);
    static std::pair<std::vector<unsigned int>, std::vector<unsigned int>> optGoodPairsVec(const BF &, std::vector<int> &);

    static void fillWSets(std::vector<uint32_t> &, std::vector<uint32_t> &, std::vector<uint32_t> &, std::vector<uint32_t> &, const std::vector<int> &);
    static BF swapOnSets(BF &, uint32_t, uint32_t);

    std::vector<std::pair<uint32_t, uint32_t>> pairsToWorsen() const;
    std::vector<std::pair<uint32_t, uint32_t>> pairsToImprove() const;
    std::vector<std::pair<uint32_t, uint32_t>> pairsToImproveStraight() const;

    bool isNeutralOrImprovePair(const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::pair<uint32_t, uint32_t> &) const;
    bool isImprovePair(const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::pair<uint32_t, uint32_t> &) const;
    bool isImprovePairEasy(const std::vector<int> &wht_coef, const std::pair<uint32_t, uint32_t> &pair) const;

    void nonlinearityImprove(uint64_t, uint64_t &, uint64_t &);
    bool testPairsToImproveFunctions(std::vector<std::pair<uint32_t, uint32_t>> &, std::vector<std::pair<uint32_t, uint32_t>> &) const;
};
