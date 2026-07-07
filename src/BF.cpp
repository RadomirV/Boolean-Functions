#include "BF.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unordered_set>

uint32_t count_first_zeros(uint32_t x)
{
    if (x == 0)
        return 0;

    uint32_t n_ = 1;
    if (x << 16 == 0)
    {
        n_ += 16;
        x >>= 16;
    }
    if (x << 24 == 0)
    {
        n_ += 8;
        x >>= 8;
    }
    if (x << 28 == 0)
    {
        n_ += 4;
        x >>= 4;
    }
    if (x << 30 == 0)
    {
        n_ += 2;
        x >>= 2;
    }
    uint32_t tmp = x & 1;
    uint32_t res = n_ - tmp;
    return res;
}

bool bit_value(unsigned int num, unsigned char bit)
{
    return ((uint32_t(1) << (bit)) & num) != 0;
}
void SetBit(uint32_t &num, uint8_t bit, bool is_one)
{
    {
        if (bit >= 32)
        {
            return;
        }

        const uint32_t mask = 1u << static_cast<unsigned>(bit);
        if (is_one)
            num |= mask;
        else
            num &= ~mask;
    }
}
int getElderBitPos(unsigned int num, int start)
{
    int i = start;
    while (i >= 0)
    {
        if (bit_value(num, i))
            break;
        i--;
    }
    return i;
}
void print_monom(unsigned int monom)
{
    if (monom == 0)
    {
        std::cout << 1;
        return;
    }
    for (unsigned int i = 0; i < 32; i++)
        if (bit_value(monom, i))
            std::cout << "x" << i;
}
uint64_t convertTo_uint64(uint32_t num1, uint32_t num2)
{
    return static_cast<uint64_t>((static_cast<uint64_t>(num1) << 32) | static_cast<uint64_t>(num2));
}
std::pair<uint32_t, uint32_t> convertTo_pair_uint32_t(uint64_t num)
{
    uint32_t num0 = static_cast<uint32_t>(num >> 32);
    uint32_t num1 = static_cast<uint32_t>(num);
    return std::make_pair(num0, num1);
}
unsigned int uint_weight(unsigned int x)
{
    x -= ((x >> 1) & 0x55555555);
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    x = (x + (x >> 4)) & 0x0f0f0f0f;
    x += (x >> 8);
    x += (x >> 16);
    x = x & 0x3f;
    return x;
}
bool weight_mod(unsigned int x)
{
    x = x ^ (x >> 1);
    x = x ^ (x >> 2);
    x = x ^ (x >> 4);
    x = x ^ (x >> 8);
    x = x ^ (x >> 16);
    return (x & 0x00000001) == 1;
}

unsigned int next_combination(unsigned int prev)
{
    unsigned int b = (prev + 1) & prev;
    unsigned int c = uint_weight((b - 1) ^ prev) - 2;
    return (((((prev + 1) ^ prev) << 1) + 1) << c) ^ b;
}
int CountBoolMatrixRankAndStepTransform(std::vector<unsigned int> &Matrix, int columns)
{
    int Rank = 0;

    columns--;
    for (int i = 0; i < Matrix.size() && (columns >= 0); i++)
    {
        auto max = std::max_element(Matrix.begin() + i, Matrix.end());
        if (*max)
        {
            std::swap(Matrix[i], *max);

            while (bit_value(Matrix[i], columns) == 0)
            {
                columns--;
            }

            for (int j = i + 1; j < Matrix.size(); j++) // annigilate all 1 under max in column
                if (bit_value(Matrix[j], columns))
                    Matrix[j] ^= Matrix[i];
            columns--;
        }
        else
            break;
    }
    while (Matrix.back() == 0)
        Matrix.pop_back();

    Rank = Matrix.size();
    return Rank;
}

std::vector<std::pair<uint32_t, uint32_t>> ConvertToPairs(const std::pair<std::vector<unsigned int>, std::vector<unsigned int>> &GPVec)
{ // из 2х множеств построить декартово произведение(массив пар)
    std::vector<std::pair<uint32_t, uint32_t>> pairs;

    for (uint32_t i = 0; i < GPVec.first.size(); i++)
    {
        for (uint32_t j = 0; j < GPVec.second.size(); j++)
        {
            pairs.push_back(std::make_pair(GPVec.first[i], GPVec.second[j]));
        }
    }
    return pairs;
}
int CountBoolMatrixRankAndStepTransform(std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns)
{
    int Rank = 0;
    columns--;
    int numberOfVariables = columns;

    for (int i = 0; i < PairMatrix.size() && (columns >= 0); i++)
    {
        auto max = std::max_element(PairMatrix.begin() + i, PairMatrix.end(), [](std::pair<uint32_t, uint8_t> pair_1, std::pair<uint32_t, uint8_t> pair_2)
                                    { return pair_1.first < pair_2.first; });
        if (max->first)
        {
            std::swap(PairMatrix[i], *max);
            while (bit_value(PairMatrix[i].first, columns) == 0)
            {
                columns--;
            }
            for (int j = i + 1; j < PairMatrix.size(); j++) // annigilate all 1 under max
                if (bit_value(PairMatrix[j].first, columns))
                {
                    PairMatrix[j].first ^= PairMatrix[i].first;
                    PairMatrix[j].second ^= PairMatrix[i].second;
                }

            columns--;
        }
        else
            break;
    }

    while (!PairMatrix.empty() && PairMatrix.back().first == 0 && PairMatrix.back().second == 0)
        PairMatrix.pop_back();

    Rank = PairMatrix.size();

    // Here after step matrix make special step matrix
    // uint32_t mainVarVec = 0;
    for (int k = PairMatrix.size() - 1; k >= 0; k--)
    {
        int elderBitPos = getElderBitPos(PairMatrix[k].first, numberOfVariables);
        // SetBit(mainVarVec, elderBitPos, 1);
        if (elderBitPos == -1)
        {
            PairMatrix = std::vector<std::pair<unsigned int, uint8_t>>();
            return -1;
        }

        for (int j = k - 1; j >= 0; j--)
        {
            if (bit_value(PairMatrix[j].first, elderBitPos))
            {
                PairMatrix[j].first ^= PairMatrix[k].first;
                PairMatrix[j].second ^= PairMatrix[k].second;
            }
        }
    }

    return Rank;
}

bool IsMatrixSpecialStepTransform(const std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns)
{
    if (PairMatrix.empty() || columns < 0)
        return false;
    for (auto it = PairMatrix.begin(); it != PairMatrix.end(); it++)
    {
        int elderPos = getElderBitPos(it->first, columns);
        for (auto it2 = it + 1; it2 != PairMatrix.end(); it2++)
        {
            if (bit_value(it2->first, elderPos))
                return false;
        }
    }

    for (auto it = PairMatrix.end() - 1;; it--)
    {
        if (it == PairMatrix.begin())
            break;
        int elderPos = getElderBitPos(it->first, columns);

        for (auto it2 = it - 1;; it2--)
        {
            if (bit_value(it2->first, elderPos))
                return false;
            if (it2 == PairMatrix.begin())
                break;
        }
    }

    return true;
}
std::vector<uint32_t> getSolutionsOfSystem(const std::vector<std::pair<uint32_t, uint8_t>> &PairMatrix, int n)
{
    if (!IsMatrixSpecialStepTransform(PairMatrix, n))
    {
        return std::vector<uint32_t>();
    }
    uint32_t mainVec = 0, dependendVec = 0;
    for (int i = 0; i < PairMatrix.size(); i++)
    {
        int elderBitPos = getElderBitPos(PairMatrix[i].first, n - 1);
        SetBit(mainVec, elderBitPos, 1);
    }
    dependendVec = ~mainVec;
    uint32_t mask = 0;
    mask = ~mask;
    mask >>= 32 - n;
    dependendVec &= mask;
    uint32_t currPreviousToDependentVec = dependendVec;

    std::vector<uint32_t> solutions;

    while (true)
    {
        uint32_t currentSolution = 0;
        currentSolution |= currPreviousToDependentVec;
        for (int i = 0; i < PairMatrix.size(); i++)
        {
            int elderBitPos = getElderBitPos(PairMatrix[i].first, n - 1);
            uint32_t r = PairMatrix[i].first & currPreviousToDependentVec;

            uint32_t bitToSet = ((uint32_t)weight_mod(r) + PairMatrix[i].second) % 2;

            bool isXone = bitToSet;

            SetBit(currentSolution, elderBitPos, isXone);
        }
        solutions.push_back(currentSolution);
        if (currPreviousToDependentVec == 0)
            break;
        currPreviousToDependentVec = (currPreviousToDependentVec - 1) & dependendVec;
    }

    return solutions;
}

BF::BF()
{
    unsigned int n_ = 1;
    std::vector<unsigned int> val;
    val.push_back(0);
}
BF::BF(const std::string &str)
{

    // the length of string must be a power of n_
    unsigned int len = str.length();

    try
    {
        bool res = len & (len - 1);
        if (res == true)
        {
            throw "incorrect string length for BF::BF(std::string &) constructor";
        }
    }
    catch (std::string msg)
    {
        std::cerr << msg << '\n';
    }

    try
    {
        bool is_string_correct = true;
        int i = 0;
        while ((str[i] == 48 || str[i] == 49) && i < len)
            i++;
        if (i != len)
            throw "incorrect symbols in string";
    }
    catch (std::string msg)
    {
        std::cerr << msg << '\n';
    }

    std::vector<unsigned int> val;
    this->n_ = count_first_zeros(len);
    int vector_size = ((len - 1) >> 5) + 1;

    for (int i = 0; i < vector_size; i++)
    {
        unsigned int mask = 1;
        unsigned int element = 0;
        for (int j = 0; (j < 32) && (j < len); j++)
        {
            if (str[(i << 5) + j] == '1')
                element |= mask;
            mask <<= 1;
        }
        val.push_back(element);
    }

    this->vec_ = val;
}
BF::BF(int n_, int par) // par=0 -> identically 0; par=1 -> identically 1; par =(any other value)random func
{
    try
    {
        if (n_ <= 0)
            throw "invalid number of variables";
    }
    catch (std::string msg)
    {
        std::cerr << msg << '\n';
    }

    unsigned int len = uint32_t(1) << n_;
    unsigned int vector_size = len >> 5;
    if (n_ <= 5)
        vector_size = 1;
    std::vector<unsigned int> val;
    switch (par)
    {
    case 0:
        for (int i = 0; i < vector_size; ++i)
            val.push_back(0);
        break;
    case 1:
        if (len < 32)
        {
            unsigned int num = 0;
            num = ~num;
            num >>= (32 - len);
            val.push_back(num);
            break;
        }
        for (int i = 0; i < vector_size; ++i)
            val.push_back(std::numeric_limits<unsigned int>::max());
        break;
    default:
        if (len < 32)
        {
            unsigned int num = rand() - rand();
            num >>= (32 - len);
            val.push_back(num);
            break;
        }
        for (int i = 0; i < vector_size; ++i)
            val.push_back(rand() - rand());
        break;
    }

    this->n_ = n_;
    this->vec_ = val;
}
BF::BF(const BF &func)
{
    this->n_ = func.n_;
    this->vec_ = func.vec_;
}
unsigned int BF::weight(const BF &a)
{
    unsigned int size = a.vec_.size();
    uint32_t total_weight = 0;

    for (int i = 0; i < size; ++i)
    {
        total_weight += uint_weight(a.vec_[i]);
    }

    return total_weight;
}
void BF::print(bool onlyVector) const
{
    if (!onlyVector)
        std::cout << "variables= " << this->n_ << std::endl;
    unsigned int len = (uint32_t)1 << this->n_;
    unsigned int mask = 1;
    unsigned int size = this->vec_.size();

    for (int j = 0; j < size; j++)
    {
        mask = 1;
        for (int i = 0; (i < 32) && (i < len); i++)
        {
            if (this->vec_[j] & mask)
                std::cout << "1";
            else
                std::cout << "0";
            mask <<= 1;
        }
        if (j % 4 == 3)
            std::cout << '\n';
        else
            std::cout << ' ';
    }
}
BF BF::operator=(const BF &func)
{
    this->n_ = func.n_;
    this->vec_ = func.vec_;
    return *this;
}
int BF::getNumberOfVariables() const
{
    return this->n_;
}
BF BF::mobius_transform(const BF &func)
{
    unsigned int vec_size = func.vec_.size();
    std::vector<unsigned int> val = func.vec_;

    if (func.n_ >= 1)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 1) & 0xaaaaaaaa);
    if (func.n_ >= 2)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 2) & 0xcccccccc);
    if (func.n_ >= 3)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 4) & 0xf0f0f0f0);
    if (func.n_ >= 4)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 8) & 0xff00ff00);
    if (func.n_ >= 5)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 16) & 0xffff0000);

    if (func.n_ > 5)
    {
        for (unsigned int i = 0; i < func.n_ - 5; ++i)
        {
            unsigned int half = 1 << i;

            for (unsigned int j = 0; j < vec_size; j = j + 2 * half)
                for (unsigned int k = j + half; k < (j + 2 * half); k++)
                    val[k] = val[k] ^ val[k - half];
        }
    }
    BF res;
    res.vec_ = val;
    res.n_ = func.n_;
    return res;
}
void BF::print_per_bit() const
{
    for (int i = 0; i < this->vec_.size(); ++i)
        std::cout << std::bitset<32>(this->vec_[i]) << std::endl;
}
bool BF::operator==(const BF &func) const
{
    if (this->n_ == func.n_)
    {
        unsigned int size = func.vec_.size();
        int i = 0;
        for (; i < size; i++)
        {
            if (this->vec_[i] != func.vec_[i])
                return false;
        }
        return true;
    }
    std::cerr << "the number of variables in operator == is different" << '\n';
    return false;
}
bool BF::operator!=(const BF &func) const
{
    if (this->n_ == func.n_)
    {
        unsigned int size = func.vec_.size();
        int i = 0;
        for (; i < size; i++)
        {
            if (this->vec_[i] != func.vec_[i])
                return true;
        }
        return false;
    }
    std::cerr << "the number of variables in operator != is different" << '\n';
    return false;
}
void BF::print_ANF(const BF &func)
{
    unsigned int curr = 0, size = func.vec_.size();
    int count_zero = std::count_if(func.vec_.begin(), func.vec_.end(), [](int i)
                                   { return i == 0; });
    if (count_zero == size)
    {
        std::cout << "0";
        return;
    }
    for (unsigned int i = 0; i < size; i++)
        for (unsigned int j = 0; j < 32; j++)
        {
            curr = i << 5;
            curr += j;
            if (bit_value(func.vec_[i], j))
            {
                print_monom(curr);
                std::cout << " + ";
            }
        }
}
unsigned int BF::degree(const BF &func)
{

    unsigned int max = 0, size = func.vec_.size();

    unsigned int weight = 0, curr_val = 0;
    for (int i = size - 1; i >= 0; i--)
    {

        for (int j = 31; j >= 0; j--)
        {
            curr_val = i << 5;
            curr_val += j;
            if (bit_value(func.vec_[i], j))
            {
                weight = uint_weight(curr_val);
                if (weight > max)
                    max = weight;
            }
        }
    }
    return max;
}
std::vector<int> BF::WH_transform(const BF &func)
{
    unsigned int size = func.vec_.size();

    std::vector<int> char_vec;
    // make characteristic vector
    for (unsigned int i = 0; i < size; ++i)
    {
        for (unsigned int j = 0; j < 32; j++)
        {
            if (bit_value(func.vec_[i], j))
                char_vec.push_back(-1);
            else
                char_vec.push_back(1);
        }
    }
    if (func.n_ < 5)
    {
        for (int i = 0; i < 32 - (1 << func.n_); i++)
            char_vec.pop_back();
    }

    unsigned int char_vec_size = char_vec.size();
    int sum = 0, dif = 0;

    for (unsigned int i = 0; i < func.n_; ++i)
    {
        unsigned int half = 1 << i;

        for (unsigned int j = 0; j < char_vec_size; j = j + 2 * half)
        {
            for (unsigned int k = j; k < j + half; k++)
            {
                sum = char_vec[k] + char_vec[k + half];
                dif = char_vec[k] - char_vec[k + half];
                char_vec[k] = sum;
                char_vec[k + half] = dif;
            }
        }
    }
    return char_vec;
}
unsigned int BF::COR(const BF &func, std::vector<int> &WHT_coef)
{

    auto gen_limits = [](unsigned int n_, unsigned int k) -> std::pair<unsigned int, unsigned int>
    {
        return std::make_pair((uint32_t)((1 << k) - 1) << (n_ - k), 0xFFFFFFFF >> (32 - k));
    };

    if (WHT_coef.empty())
        WHT_coef = WH_transform(func);

    int i = 1;

    for (; i <= func.n_; i++)
    {
        auto limits = gen_limits(func.n_, i);
        unsigned int current = limits.first;

        if (WHT_coef[current] != 0)
            return i - 1;

        while (current != limits.second)
        {
            current = next_combination(current);
            if (WHT_coef[current] != 0)
                return i - 1;
        }
    }
    return func.n_;
}
unsigned int BF::Nonlinearity(const BF &func, std::vector<int> &WHT)
{
    if (WHT.empty())
        WHT = WH_transform(func);
    auto min_max = std::minmax_element(WHT.begin(), WHT.end());
    int max = 0;
    if (abs(*min_max.first) > abs(*min_max.second))
        max = abs(*min_max.first);
    else
        max = abs(*min_max.second);

    return (1 << (func.n_ - 1)) - (max >> 1);
}
void BF::print_BAA(const BF &func, std::vector<int> &WHT)
{
    if (WHT.empty())
        WHT = WH_transform(func);

    auto min_max = std::minmax_element(WHT.begin(), WHT.end());
    unsigned int max_pos = 0;
    if (abs(*min_max.first) > abs(*min_max.second))
        max_pos = min_max.first - WHT.begin();
    else
        max_pos = min_max.second - WHT.begin();

    bool is_max_positive = WHT[max_pos] > 0 ? false : true;

    for (unsigned int i = 0; i < 31; i++)
        if (bit_value(max_pos, i))
            std::cout << "x" << i << "+";
    if (is_max_positive)
        std::cout << "1";
}
bool BF::getValue(const BF &func, unsigned int x_set)
{
    unsigned int ix = x_set >> 5;
    unsigned int bit = x_set & 0x0000001f;
    return bit_value(func.vec_[ix], bit);
}
bool BF::operator[](unsigned int x_set) const
{
    unsigned int ix = x_set >> 5;
    unsigned int bit = x_set & 0x0000001f;
    return bit_value(this->vec_[ix], bit);
}
void BF::setValue(unsigned int x_set, bool val)
{

    unsigned int ix = x_set >> 5;
    unsigned int bit = x_set & 0x0000001f;
    unsigned int mask = 1;

    this->vec_[ix] = (this->vec_[ix] & ~(1u << bit)) | (static_cast<unsigned int>(val) << bit);
}
std::vector<unsigned int> BF::values()
{
    return this->vec_;
}
std::vector<int> BF::Auto_Cor(const BF &func, std::vector<int> wht_coef)
{
    if (wht_coef.empty())
        wht_coef = BF::WH_transform(func);
    unsigned int size = wht_coef.size();

    std::for_each(wht_coef.begin(), wht_coef.end(), [](int &n_)
                  { n_ *= n_; });

    int sum = 0, dif = 0;

    for (unsigned int i = 0; i < func.n_; ++i)
    {
        unsigned int half = 1 << i;

        for (unsigned int j = 0; j < size; j = j + 2 * half)
        {
            for (unsigned int k = j; k < j + half; k++)
            {
                sum = wht_coef[k] + wht_coef[k + half];
                dif = wht_coef[k] - wht_coef[k + half];
                wht_coef[k] = sum;
                wht_coef[k + half] = dif;
            }
        }
    }
    unsigned int n_variables = func.n_;
    std::for_each(wht_coef.begin(), wht_coef.end(), [&n_variables](int &n_)
                  { n_ >>= n_variables; });

    return wht_coef;
}
int BF::CN_f(const BF &func, std::vector<int> &auto_cor)
{
    if (auto_cor.empty())
        auto_cor = BF::Auto_Cor(func, std::vector<int>());

    auto min_max = std::minmax_element(auto_cor.begin() + 1, auto_cor.end());
    int max = 0;
    if (abs(*min_max.first) > abs(*min_max.second))
        max = abs(*min_max.first);
    else
        max = abs(*min_max.second);

    return ((1 << (func.n_ - 2)) - (max >> 2));
}
unsigned int BF::PC(const BF &func, std::vector<int> &auto_cor)
{
    if (auto_cor.empty())
        auto_cor = BF::Auto_Cor(func, std::vector<int>());

    auto gen_limits = [](unsigned int n_, unsigned int k) -> std::pair<unsigned int, unsigned int>
    {
        return std::make_pair((uint32_t)((1 << k) - 1) << (n_ - k), 0xFFFFFFFF >> (32 - k));
    };

    int i = 1;

    for (; i <= func.n_; i++)
    {
        auto limits = gen_limits(func.n_, i);
        unsigned int current = limits.first;

        if (auto_cor[current] != 0)
            return i - 1;

        while (current != limits.second)
        {
            current = next_combination(current);
            if (auto_cor[current] != 0)
                return i - 1;
        }
    }
    return func.n_;
}

std::pair<std::vector<unsigned int>, std::vector<unsigned int>> BF::OptGoodPairsVec(const BF &func, std::vector<int> &wht_coef)
{
    auto min_max = std::minmax_element(wht_coef.begin(), wht_coef.end());
    int max;
    if (abs(*min_max.first) > abs(*min_max.second))
        max = abs(*min_max.first);
    else
        max = abs(*min_max.second);
    std::vector<uint32_t> W_1_pos, W_1_neg;
    unsigned int current = 0;
    for (int i = 0; i < wht_coef.size(); i++)
    {
        if (wht_coef[i] == max)
            W_1_pos.push_back(i);
        if (wht_coef[i] == max * (-1))
            W_1_neg.push_back(i);
    }

    std::vector<std::pair<uint32_t, uint8_t>> PairMatrix;
    for (int i = 0; i < W_1_pos.size(); i++)
        PairMatrix.push_back(std::make_pair(W_1_pos[i], 0));
    for (int i = 0; i < W_1_neg.size(); i++)
        PairMatrix.push_back(std::make_pair(W_1_neg[i], 1));

    int Rank = CountBoolMatrixRankAndStepTransform(PairMatrix, func.n_);

    if (Rank == -1)
    {
        return std::pair<std::vector<unsigned int>, std::vector<unsigned int>>();
    }

    auto solutionsFirst = getSolutionsOfSystem(PairMatrix, func.n_);

    std::vector<uint32_t> finalSolutionsFirst;

    for (int i = 0; i < solutionsFirst.size(); i++)
        if (func[solutionsFirst[i]] == 0)
            finalSolutionsFirst.push_back(solutionsFirst[i]);

    PairMatrix.clear();

    // now we found the first element of good pairs
    for (int i = 0; i < W_1_pos.size(); i++)
        PairMatrix.push_back(std::make_pair(W_1_pos[i], 1));
    for (int i = 0; i < W_1_neg.size(); i++)
        PairMatrix.push_back(std::make_pair(W_1_neg[i], 0));

    Rank = CountBoolMatrixRankAndStepTransform(PairMatrix, func.n_);

    if (Rank == -1)
    {
        return std::pair<std::vector<unsigned int>, std::vector<unsigned int>>();
    }

    auto solutionsSecond = getSolutionsOfSystem(PairMatrix, func.n_);
    std::vector<uint32_t> finalSolutionsSecond;

    for (int i = 0; i < solutionsSecond.size(); i++)
        if (func[solutionsSecond[i]] == 1)
            finalSolutionsSecond.push_back(solutionsSecond[i]);
    std::pair<std::vector<unsigned int>, std::vector<unsigned int>> ResultGoodPairs;
    ResultGoodPairs = std::make_pair(finalSolutionsFirst, finalSolutionsSecond);

    return ResultGoodPairs;
}
std::pair<BF, BF> BF::generateBorderBalancesFunctions()
{
    return std::make_pair(BF(std::bitset<32>(0xFFFFFFFF >> (16)).to_string()), BF(std::bitset<32>(((1 << 16) - 1) << (16)).to_string()));
}
void BF::nextBalanced()
{
    this->vec_[0] = next_combination(this->vec_[0]);
}
void BF::FillWSets(std::vector<uint32_t> &W_1_pos, std::vector<uint32_t> &W_1_neg, std::vector<uint32_t> &W_3_pos, std::vector<uint32_t> &W_3_neg, const std::vector<int> &wht_coef)
{
    auto min_max = std::minmax_element(wht_coef.begin(), wht_coef.end());
    int max;
    if (abs(*min_max.first) > abs(*min_max.second))
        max = abs(*min_max.first);
    else
        max = abs(*min_max.second);
    for (int i = 0; i < wht_coef.size(); i++)
    {
        if (wht_coef[i] == max - 4)
        {
            W_3_pos.push_back(i);
        }
        if (wht_coef[i] == (max * (-1)) + 4)
        {
            W_3_neg.push_back(i);
        }
        if (wht_coef[i] == max)
        {
            W_1_pos.push_back(i);
        }
        if (wht_coef[i] == (max * (-1)))
        {
            W_1_neg.push_back(i);
        }
    }
}
BF BF::SwapOnSets(BF &func, uint32_t set_x1, uint32_t set_x2)
{
    if (set_x1 >= ((uint32_t)1 << func.n_) || set_x2 >= ((uint32_t)1 << func.n_))
    {
        throw "pair is out of range";
    }
    auto result = func;
    result.setValue(set_x1, func[set_x2]);
    result.setValue(set_x2, func[set_x1]);
    return result;
}
BF BF::GenBalancedFunc(int numberOfVariables) // Генерируем функцию а потом подправляем ее до уравновешенной
{
    BF func(numberOfVariables, 10);

    int32_t diff = ((uint32_t)1 << (numberOfVariables - 1)) - BF::weight(func); // разница чего больше единиц или нулей

    while (diff < 0) // если больше единиц, то в случайных местах где f(x)=1 заменить на 0
    {
        uint32_t mask = 0;
        mask = ~mask;
        mask >>= 32 - numberOfVariables;
        uint32_t x_set = rand() & mask;
        if (func[x_set] == true)
        {
            func.setValue(x_set, false);
            diff++;
        }
    }
    while (diff > 0) // если больше нулей, то в случайных местах где f(x)=0 заменить на 1
    {
        uint32_t mask = 0;
        mask = ~mask;
        mask >>= 32 - numberOfVariables;
        uint32_t x_set = rand() & mask;
        if (func[x_set] == false)
        {
            func.setValue(x_set, true);
            diff--;
        }
    }
    return func;
}
std::pair<std::vector<unsigned int>, std::vector<unsigned int>> BF::good_pairsVec(const BF &func, std::vector<int> &wht_coef)
{
    const unsigned int n_ = func.n_;
    auto min_max = std::minmax_element(wht_coef.begin(), wht_coef.end());
    int max;
    if (abs(*min_max.first) > abs(*min_max.second))
        max = abs(*min_max.first);
    else
        max = abs(*min_max.second);
    std::vector<uint32_t> W_1_pos, W_1_neg;
    unsigned int current = 0;
    for (int i = 0; i < wht_coef.size(); i++)
    {
        if (wht_coef[i] == max)
            W_1_pos.push_back(i);
        if (wht_coef[i] == max * (-1))
            W_1_neg.push_back(i);
    }
    std::vector<unsigned int> A_00, B_01, A_11, B_10;
    for (unsigned int curr = 0; curr < 1 << n_; curr++)
    {
        if (func[curr] == false)
        {
            int a = 0;
            for (; a < W_1_pos.size(); a++)
                if (weight_mod(W_1_pos[a] & curr) == true)
                    break;
            if (a == W_1_pos.size() && W_1_pos.size())
                A_00.push_back(curr);

            int b = 0;
            for (; b < W_1_neg.size(); b++)
                if (weight_mod(W_1_neg[b] & curr) == false)
                    break;
            if (b == W_1_neg.size() && W_1_neg.size())
                B_01.push_back(curr);
        }
        else
        {
            int a = 0;
            for (; a < W_1_pos.size(); a++)
                if (weight_mod(W_1_pos[a] & curr) == false)
                    break;
            if (a == W_1_pos.size() && W_1_pos.size())
                A_11.push_back(curr);

            int b = 0;
            for (; b < W_1_neg.size(); b++)
                if (weight_mod(W_1_neg[b] & curr) == true)
                    break;
            if (b == W_1_neg.size() && W_1_neg.size())
                B_10.push_back(curr);
        }
    }

    std::vector<unsigned int> AB_0, AB_1;
    std::set_intersection(A_00.begin(), A_00.end(), B_01.begin(), B_01.end(), std::insert_iterator<std::vector<unsigned int>>(AB_0, AB_0.begin()));
    std::set_intersection(A_11.begin(), A_11.end(), B_10.begin(), B_10.end(), std::insert_iterator<std::vector<unsigned int>>(AB_1, AB_1.begin()));
    if (W_1_pos.size() == 0)
    {
        AB_0 = B_01;
        AB_1 = B_10;
    }
    if (W_1_neg.size() == 0)
    {
        AB_0 = A_00;
        AB_1 = A_11;
    }

    return std::make_pair(AB_0, AB_1);
}

std::vector<std::pair<uint32_t, uint32_t>> BF::PairsToWorsen() const
{
    auto wht_coef = BF::WH_transform(*this);
    std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
    BF::FillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);

    std::vector<std::vector<uint32_t>> A_01, A_10, B_00, B_11;
    A_01.resize(W_1_pos.size());
    A_10.resize(W_1_pos.size());
    B_00.resize(W_1_neg.size());
    B_11.resize(W_1_neg.size());
    for (unsigned int curr = 0; curr < 1 << n_; curr++)
    {
        if (this->operator[](curr) == false)
        {
            for (int i = 0; i < W_1_pos.size(); ++i)
            {
                if (weight_mod(curr & W_1_pos[i]) == true)
                    A_01[i].push_back(curr);
            }
            for (int i = 0; i < W_1_neg.size(); ++i)
            {
                if (weight_mod(curr & W_1_neg[i]) == false)
                    B_00[i].push_back(curr);
            }
        }
        else
        {
            for (int i = 0; i < W_1_pos.size(); ++i)
            {
                if (weight_mod(curr & W_1_pos[i]) == false)
                    A_10[i].push_back(curr);
            }
            for (int i = 0; i < W_1_neg.size(); ++i)
            {
                if (weight_mod(curr & W_1_neg[i]) == true)
                    B_11[i].push_back(curr);
            }
        }
    }

    std::unordered_set<uint64_t> Allvalues;
    for (int i = 0; i < W_1_pos.size(); ++i)
    {
        if (!A_01[i].empty() && !A_10[i].empty())
        {
            for (int j = 0; j < A_01[i].size(); ++j)
            {
                for (int k = 0; k < A_10[i].size(); ++k)
                {
                    Allvalues.insert(convertTo_uint64(A_01[i][j], A_10[i][k]));
                }
            }
        }
    }
    for (int i = 0; i < W_1_neg.size(); ++i)
    {
        if (!B_00[i].empty() && !B_11[i].empty())
        {
            for (int j = 0; j < B_00[i].size(); ++j)
            {
                for (int k = 0; k < B_11[i].size(); ++k)
                {
                    Allvalues.insert(convertTo_uint64(B_00[i][j], B_11[i][k]));
                }
            }
        }
    }

    std::vector<std::pair<uint32_t, uint32_t>> result;
    for (auto &iter : Allvalues)
    {
        result.push_back(convertTo_pair_uint32_t(iter));
    }

    return result;
}
std::vector<std::pair<uint32_t, uint32_t>> BF::pairsToImprove() const
{ // построить множесто пар улучшающих нелинейность
    auto wht_coef = BF::WH_transform(*this);
    auto GoodPairs = BF::OptGoodPairsVec(*this, wht_coef); // Вернуть пару множеств(2 множества из теоремы слева от разности) с помощью решения СЛУ

    std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
    BF::FillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);

    std::vector<std::pair<uint32_t, uint32_t>> pairs = ConvertToPairs(GoodPairs); // конвертируем в декартово произведение
    std::vector<std::pair<uint32_t, uint32_t>> res;
    for (uint64_t i = 0; i < pairs.size(); i++) // проверяю удовлетворяет ли пара условиям для улучшение нелинейность
    {
        int flag = 0;

        for (int j = 0; j < W_3_pos.size(); j++)
        {
            bool correct = weight_mod(W_3_pos[j] & pairs[i].first) == false || weight_mod(W_3_pos[j] & pairs[i].second) == true; // (с,x1)=0 ⋁ (с,x2)=1
            if (correct)                                                                                                         //
                flag++;
            else
                break;
        }
        if (flag != W_3_pos.size()) // Если пара не удовлетворила какому-то элементу из W_3^+ то смотрим следующую пару
            continue;

        flag = 0;

        for (int j = 0; j < W_3_neg.size(); j++)
        {
            bool correct = weight_mod(W_3_neg[j] & pairs[i].first) == true || weight_mod(W_3_neg[j] & pairs[i].second) == false; //(d,x1)=1 ⋁ (d,x2)=0
            if (correct)
                flag++;
            else
                break;
        }
        if (flag != W_3_neg.size()) // Если пара не удовлетворила какому-то элементу из W_3^- то смотрим следующую пару
            continue;

        res.push_back(pairs[i]); // Пара удовлетворила всем условиям -> она улучшает нелинейность
    }

    return res;
}
std::vector<std::pair<uint32_t, uint32_t>> BF::pairsToImproveStraight() const // отличие только в том как вычисялются множества слева от разности
{
    auto wht_coef = BF::WH_transform(*this);
    auto GoodPairs = BF::good_pairsVec(*this, wht_coef);

    std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
    BF::FillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);

    std::vector<std::pair<uint32_t, uint32_t>> pairs = ConvertToPairs(GoodPairs), res;
    for (uint64_t i = 0; i < pairs.size(); i++)
    {
        int flag = 0;

        for (int j = 0; j < W_3_pos.size(); j++)
        {
            bool correct = weight_mod(W_3_pos[j] & pairs[i].first) == false || weight_mod(W_3_pos[j] & pairs[i].second) == true;
            if (correct)
                flag++;
            else
                break;
        }
        if (flag != W_3_pos.size())
            continue;
        ;

        flag = 0;

        for (int j = 0; j < W_3_neg.size(); j++)
        {
            bool correct = weight_mod(W_3_neg[j] & pairs[i].first) == true || weight_mod(W_3_neg[j] & pairs[i].second) == false;
            if (correct)
                flag++;
            else
                break;
        }
        if (flag != W_3_neg.size())
            continue;

        res.push_back(pairs[i]);
    }

    return res;
}
bool BF::testPairsToImproveFunctions(std::vector<std::pair<uint32_t, uint32_t>> &Optimized, std::vector<std::pair<uint32_t, uint32_t>> &Straight) const
{
    if (Optimized.size() != Straight.size())
        return false;
    for (auto &j : Optimized)
    {
        auto it = std::find_if(begin(Straight), end(Straight), [j](std::pair<uint32_t, uint32_t> my_pair)
                               { if (j.first == my_pair.first && j.second == my_pair.second) return true; else return false; });
        if (it == std::end(Straight))
            return false;
    }
    return true;
}
bool BF::isNeutralOrImprovePair(const std::vector<uint32_t> &W_1_pos, const std::vector<uint32_t> &W_1_neg, const std::pair<uint32_t, uint32_t> &pair) const
{
    for (auto &i : W_1_pos)
    {
        bool isWorsening = weight_mod(i & pair.first) == true && weight_mod(i & pair.second) == false; // (w,x1)==1 && (w,x2)==0
        if (isWorsening)
            return false;
    }
    for (auto &i : W_1_neg)
    {
        bool isWorsening = weight_mod(i & pair.first) == false && weight_mod(i & pair.second) == true; // (w,x1)==0 && (w,x2)==1
        if (isWorsening)
            return false;
    }

    return true;
}

bool BF::isImprovePair(const std::vector<uint32_t> &W_1_pos, const std::vector<uint32_t> &W_1_neg, const std::vector<uint32_t> &W_3_pos, const std::vector<uint32_t> &W_3_neg, const std::pair<uint32_t, uint32_t> &pair) const
{

    for (auto &i : W_1_pos)
    {
        bool correct = weight_mod(i & pair.first) == false && weight_mod(i & pair.second) == true;
        if (!correct)
            return false;
    }

    for (auto &i : W_1_neg)
    {
        bool correct = weight_mod(i & pair.first) == true && weight_mod(i & pair.second) == false;
        if (!correct)
            return false;
    }

    for (int j = 0; j < W_3_pos.size(); j++)
    {
        bool correct = weight_mod(W_3_pos[j] & pair.first) == false || weight_mod(W_3_pos[j] & pair.second) == true; // (с,x1)=0 ⋁ (с,x2)=1
        if (!correct)                                                                                                //
            return false;
    }

    for (int j = 0; j < W_3_neg.size(); j++)
    {
        bool correct = weight_mod(W_3_neg[j] & pair.first) == true || weight_mod(W_3_neg[j] & pair.second) == false; //(d,x1)=1 ⋁ (d,x2)=0
        if (!correct)
            return false;
    }

    return true; // Пара удовлетворила всем условиям -> она улучшает нелинейность
}
bool BF::isImprovePairEasy(const std::vector<int> &wht_coef, const std::pair<uint32_t, uint32_t> &pair) const
{
    std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
    BF::FillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);
    for (int j = 0; j < W_3_pos.size(); j++)
    {
        bool correct = weight_mod(W_3_pos[j] & pair.first) == false || weight_mod(W_3_pos[j] & pair.second) == true; // (с,x1)=0 ⋁ (с,x2)=1
        if (!correct)                                                                                                //
            return false;
    }

    for (int j = 0; j < W_3_neg.size(); j++)
    {
        bool correct = weight_mod(W_3_neg[j] & pair.first) == true || weight_mod(W_3_neg[j] & pair.second) == false; //(d,x1)=1 ⋁ (d,x2)=0
        if (!correct)
            return false;
    }

    return false;
}
BF BF::generateAffine(uint32_t maskOfVariables, bool addOne) const
{
    BF affineFunc(this->n_);
    uint32_t mask = 1;
    for (int i = 0; i <= this->n_; ++i)
    {
        if (mask & maskOfVariables)
        {
            affineFunc.setValue(mask, true);
        }
        mask <<= 1;
    }
    affineFunc.setValue(0, addOne);
    affineFunc = BF::mobius_transform(affineFunc);
    return affineFunc;
}
std::pair<uint32_t, uint32_t> BF::generatePair() const
{
    std::pair<uint32_t, uint32_t> pair = std::make_pair(rand() % (1 << this->n_), rand() % (1 << this->n_));
    while (this->operator[](pair.first) == this->operator[](pair.second))
    {
        pair.second = rand() % (1 << this->n_);
    }
    if (this->operator[](pair.first) == true)
    {
        std::swap(pair.first, pair.second);
    }
    return pair;
}
std::pair<uint32_t, uint32_t> BF::generateImprovePair(uint32_t attemps) const
{
    std::pair<uint32_t, uint32_t> improvePair;
    auto wht_coef = BF::WH_transform(*this);
    auto goodSet = BF::OptGoodPairsVec(*this, wht_coef);
    for (int i = 0; i <= attemps; ++i)
    {
        improvePair.first = goodSet.first[rand() % goodSet.first.size()];
        improvePair.second = goodSet.second[rand() % goodSet.second.size()];
        if (BF::isImprovePairEasy(wht_coef, improvePair))
        {
            return improvePair;
        }
    }
    return std::pair<uint32_t, uint32_t>();
}
void BF::nonlinearityImprove(uint64_t linDiff, uint64_t &neutralAttemps, uint64_t &improveAttemps)
{
    if (linDiff % 2 != 0)
    {
        return;
    }
    auto wht_coef = BF::WH_transform(*this);

    auto tmp = std::roundf((1 << (this->n_ - 1)) - std::pow(2, this->n_ / 2 - 1));
    int Nf_boundaries = this->n_ % 2 == 0 ? tmp - 2 : tmp;

    std::pair<uint32_t, uint32_t> neutralPair;
    std::pair<uint32_t, uint32_t> pair;
    for (int j = 0; j <= improveAttemps && linDiff > 0; ++j)
    {

        std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
        wht_coef = BF::WH_transform(*this);
        BF::FillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);
        auto NF = BF::Nonlinearity(*this, wht_coef);
        if (BF::Nonlinearity(*this, wht_coef) >= Nf_boundaries)
            break;

        bool needToBuildImproveSet = true;
        std::vector<std::pair<uint32_t, uint32_t>> improvePairs;

        for (int i = 0; i < 100000; ++i)
        {
            pair = this->generatePair();
            if (this->isImprovePair(W_1_pos, W_1_neg, W_3_pos, W_3_neg, pair))
            {
                *this = BF::SwapOnSets(*this, pair.first, pair.second);
                linDiff -= 2;
                j = 0;
                needToBuildImproveSet = false;
                break;
            }
        }
        // std::cout << " W_1_pos = " << W_1_pos.size()
        //           << " W_1_neg = " << W_1_neg.size()
        //           << " W_3_pos = " << W_3_pos.size()
        //           << " W_3_neg = " << W_3_neg.size()
        //           << std::endl;
        if (!needToBuildImproveSet)
        {
            continue;
        }

        improvePairs = this->pairsToImprove();

        if (improvePairs.size() != 0)
        {
            pair = improvePairs[rand() % improvePairs.size()];
            *this = BF::SwapOnSets(*this, pair.first, pair.second);
            linDiff -= 2;
            j = 0;
        }
        else
        {
            int i = 0;
            bool foundNeutral = false;
            while (i <= neutralAttemps)
            {
                neutralPair = this->generatePair();
                if (this->isNeutralOrImprovePair(W_1_pos, W_1_neg, neutralPair))
                {
                    foundNeutral = true;
                    break;
                }
                ++i;
            }
            if (foundNeutral)
            {
                *this = BF::SwapOnSets(*this, neutralPair.first, neutralPair.second);
            }
            else
            {
                neutralAttemps = 0;
                return;
            }
        }
    }

    improveAttemps = 0;
}
