#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace bf
{

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
    bool getValue(unsigned int) const;
    void setValue(unsigned int, bool);
    std::vector<unsigned int> values();
    void print(bool onlyVector = true) const;
    void printPerBit() const;

    unsigned int weight() const;
    unsigned int degree() const;
    unsigned int correlationImmunity(std::vector<int> &WHT) const;
    unsigned int nonlinearity(std::vector<int> &WHT) const;
    int cnF(std::vector<int> &) const;
    unsigned int propagationCriteria(std::vector<int> &) const;

    void printANF() const;
    void printBAA(std::vector<int> &WHT) const;

    BF mobiusTransform() const;
    std::vector<int> WHTransform() const;
    std::vector<int> autoCor(std::vector<int>) const;

    BF generateAffine(uint32_t maskOfVariables, bool addOne) const;
    std::pair<uint32_t, uint32_t> generatePair() const;
    std::pair<uint32_t, uint32_t> generateImprovePair(uint32_t) const;
    static std::pair<BF, BF> generateBorderBalancesFunctions();
    static BF genBalancedFunc(int numberOfVariables);
    void nextBalanced();

    std::pair<std::vector<unsigned int>, std::vector<unsigned int>> goodPairsVec(std::vector<int> &) const;
    std::pair<std::vector<unsigned int>, std::vector<unsigned int>> optGoodPairsVec(std::vector<int> &) const;

    static void fillWSets(std::vector<uint32_t> &, std::vector<uint32_t> &, std::vector<uint32_t> &, std::vector<uint32_t> &, const std::vector<int> &);
    BF swapOnSets(uint32_t, uint32_t) const;

    std::vector<std::pair<uint32_t, uint32_t>> pairsToWorsen() const;
    std::vector<std::pair<uint32_t, uint32_t>> pairsToImprove() const;
    std::vector<std::pair<uint32_t, uint32_t>> pairsToImproveStraight() const;

    bool isNeutralOrImprovePair(const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::pair<uint32_t, uint32_t> &) const;
    bool isImprovePair(const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::pair<uint32_t, uint32_t> &) const;
    bool isImprovePairEasy(const std::vector<int> &wht_coef, const std::pair<uint32_t, uint32_t> &pair) const;

    void nonlinearityImprove(uint64_t, uint64_t &, uint64_t &);
    bool testPairsToImproveFunctions(std::vector<std::pair<uint32_t, uint32_t>> &, std::vector<std::pair<uint32_t, uint32_t>> &) const;
};

} // namespace bf
