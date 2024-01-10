#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <bitset>
#include <algorithm>
#include <set>
#pragma once

class BF
{
private:
    std::vector<unsigned int> vec;
    int n;

public:
    BF();
    BF(const std::string &);
    BF(int n, int par = 0);
    BF(const BF &);
    ~BF() {}
    BF operator=(const BF &);
    static unsigned int weight(const BF &);
    int get_number_of_variables();
    static BF mobius_transform(const BF &);
    static void print_ANF(const BF &);
    static unsigned int degree(const BF &);
    static unsigned int COR(const BF &, std::vector<int> &WHT);
    static unsigned int Nonlinearity(const BF &, std::vector<int> &WHT);
    static void print_BAA(const BF &, std::vector<int> &WHT);
    inline static bool get_val(const BF &, unsigned int);
    void print();
    void print_per_bit();
    bool operator==(const BF &) const;
    bool operator!=(const BF &) const;
    bool operator[](unsigned int);
    void set_val(unsigned int, bool);
    static std::vector<int> WH_transform(const BF &);
    static std::vector<int> Auto_Cor(const BF &, std::vector<int>);
    static int CN_f(const BF &, std::vector<int> &);
    static unsigned int PC(const BF &, std::vector<int> &);
    static std::pair<std::vector<unsigned int>, std::vector<unsigned int>> good_pairsVec(BF &, std::vector<int> &);
    static std::pair<std::vector<unsigned int>, std::vector<unsigned int>> OptGoodPairsVec(BF &, std::vector<int> &);
    static void FillWSets(std::vector<uint32_t> &, std::vector<uint32_t> &, std::vector<uint32_t> &, std::vector<uint32_t> &, std::vector<int> &);
    static BF SwapOnSets(BF &, uint32_t, uint32_t);
    static BF GenBalancedFunc(int numberOfVariables);
    std::vector<unsigned int> get_vec()
    {
        return this->vec;
    }

    std::vector<std::pair<uint32_t, uint32_t>> PairsToImprove();
    std::vector<std::pair<uint32_t, uint32_t>> PairsToImproveStraight();
    bool TestPairsToImproveFunctions(std::vector<std::pair<uint32_t, uint32_t>> &, std::vector<std::pair<uint32_t, uint32_t>> &);
};

int CountBoolMatrixRankAndStepTransform(std::vector<unsigned int> &Matrix, int columns = 32);
int CountBoolMatrixRankAndStepTransform(std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns = 32);
std::vector<std::pair<uint32_t, uint32_t>> ConvertToPairs(std::pair<std::vector<unsigned int>, std::vector<unsigned int>> &GPVec);
