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

    static unsigned int weight(const BF &);
    static unsigned int degree(const BF &);
    static unsigned int COR(const BF &, std::vector<int> &WHT);
    static unsigned int Nonlinearity(const BF &, std::vector<int> &WHT);
    static int CN_f(const BF &, std::vector<int> &);
    static unsigned int PC(const BF &, std::vector<int> &);

    static void print_ANF(const BF &);
    static void print_BAA(const BF &, std::vector<int> &WHT);
    void print(bool onlyVector = true) const;
    void print_per_bit() const;

    static BF mobius_transform(const BF &);
    static std::vector<int> WH_transform(const BF &);
    static std::vector<int> Auto_Cor(const BF &, std::vector<int>);

    static std::pair<std::vector<unsigned int>, std::vector<unsigned int>> good_pairsVec(const BF &, std::vector<int> &);
    static std::pair<std::vector<unsigned int>, std::vector<unsigned int>> OptGoodPairsVec(const BF &, std::vector<int> &);

    static void FillWSets(std::vector<uint32_t> &, std::vector<uint32_t> &, std::vector<uint32_t> &, std::vector<uint32_t> &, const std::vector<int> &);
    static BF SwapOnSets(BF &, uint32_t, uint32_t);

    BF generateAffine(uint32_t maskOfVariables, bool addOne) const;
    std::pair<uint32_t, uint32_t> generatePair() const;
    std::pair<uint32_t, uint32_t> generateImprovePair(uint32_t) const;
    static std::pair<BF, BF> generateBorderBalancesFunctions();
    static BF GenBalancedFunc(int numberOfVariables);
    void nextBalanced();

    std::vector<std::pair<uint32_t, uint32_t>> PairsToWorsen() const;
    std::vector<std::pair<uint32_t, uint32_t>> pairsToImprove() const;
    std::vector<std::pair<uint32_t, uint32_t>> pairsToImproveStraight() const;

    bool isNeutralOrImprovePair(const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::pair<uint32_t, uint32_t> &) const;
    bool isImprovePair(const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::vector<uint32_t> &, const std::pair<uint32_t, uint32_t> &) const;
    bool isImprovePairEasy(const std::vector<int> &wht_coef, const std::pair<uint32_t, uint32_t> &pair) const;

    void nonlinearityImprove(uint64_t, uint64_t &, uint64_t &);
    bool testPairsToImproveFunctions(std::vector<std::pair<uint32_t, uint32_t>> &, std::vector<std::pair<uint32_t, uint32_t>> &) const;
};

int CountBoolMatrixRankAndStepTransform(std::vector<unsigned int> &Matrix, int columns = 32);
int CountBoolMatrixRankAndStepTransform(std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns = 32);
bool IsMatrixSpecialStepTransform(const std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns = 32);
std::vector<uint32_t> getSolutionsOfSystem(const std::vector<std::pair<uint32_t, uint8_t>> &PairMatrix, int n_);

std::vector<std::pair<uint32_t, uint32_t>> ConvertToPairs(const std::pair<std::vector<unsigned int>, std::vector<unsigned int>> &GPVec);
uint64_t convertTo_uint64(uint32_t num1, uint32_t num2);
std::pair<uint32_t, uint32_t> convertTo_pair_uint32_t(uint64_t num);

uint32_t count_first_zeros(uint32_t x);
int getElderBitPos(unsigned int num, int start = 31);
bool bit_value(unsigned int num, unsigned char bit);
void SetBit(uint32_t &num, uint8_t bit, bool is_one);
void print_monom(unsigned int monom);

unsigned int uint_weight(unsigned int x);
bool weight_mod(unsigned int x);

unsigned int next_combination(unsigned int prev);
