#include "BF.h"

int count_first_zeros(unsigned int x)
{
    if (x == 0)
        return 0;

    unsigned int n = 1;
    if (x << 16 == 0)
    {
        n += 16;
        x >>= 16;
    }
    if (x << 24 == 0)
    {
        n += 8;
        x >>= 8;
    }
    if (x << 28 == 0)
    {
        n += 4;
        x >>= 4;
    }
    if (x << 30 == 0)
    {
        n += 2;
        x >>= 2;
    }
    unsigned int tmp = x & 1;
    unsigned int res = n - tmp;
    return res;
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

inline unsigned int next_combination(unsigned int prev)
{
    unsigned int b = (prev + 1) & prev;
    unsigned int c = uint_weight((b - 1) ^ prev) - 2;
    return (((((prev + 1) ^ prev) << 1) + 1) << c) ^ b;
}

BF::BF()
{
    unsigned int n = 1;
    std::vector<unsigned int> val;
    val.push_back(0);
}
BF::BF(const std::string &str)
{

    // the length of string must be a power of n
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
        std::cout << msg << '\n';
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
        std::cout << msg << '\n';
    }

    std::vector<unsigned int> val;
    this->n = count_first_zeros(len);
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

    this->vec = val;
}
BF::BF(int n, int par) // par=0 -> identically 0; par=1 -> identically 1; par =(any other value)random func
{
    try
    {
        if (n <= 0)
            throw "invalid number of variables";
    }
    catch (std::string msg)
    {
        std::cout << msg << '\n';
    }

    unsigned int len = uint32_t(1) << n;
    unsigned int vector_size = len >> 5;
    if (n <= 5)
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

    this->n = n;
    this->vec = val;
}
BF::BF(const BF &func)
{
    this->n = func.n;
    this->vec = func.vec;
}

unsigned int BF::weight(const BF &a)
{
    unsigned int size = a.vec.size();
    uint32_t total_weight = 0;

    for (int i = 0; i < size; ++i)
    {
        total_weight += uint_weight(a.vec[i]);
    }

    return total_weight;
}

void BF::print()
{
    std::cout << "variables= " << this->n << std::endl;
    unsigned int len = (uint32_t)1 << this->n;
    unsigned int mask = 1;
    unsigned int size = this->vec.size();

    for (int j = 0; j < size; j++)
    {
        mask = 1;
        for (int i = 0; (i < 32) && (i < len); i++)
        {
            if (this->vec[j] & mask)
                std::cout << "1";
            else
                std::cout << "0";
            mask <<= 1;
        }
        std::cout << '\n';
    }

    // for (int i = 0; i < this->vec.size(); ++i)
    //     std::cout << std::bitset<32>(this->vec[i]) << std::endl;
}

BF BF::operator=(const BF &func)
{
    this->n = func.n;
    this->vec = func.vec;
    return *this;
}

int BF::get_number_of_variables()
{
    return this->n;
}

BF BF::mobius_transform(const BF &func)
{
    unsigned int vec_size = func.vec.size();
    std::vector<unsigned int> val = func.vec;

    if (func.n >= 1)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 1) & 0xaaaaaaaa);
    if (func.n >= 2)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 2) & 0xcccccccc);
    if (func.n >= 3)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 4) & 0xf0f0f0f0);
    if (func.n >= 4)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 8) & 0xff00ff00);
    if (func.n >= 5)
        for (int i = 0; i < vec_size; ++i)
            val[i] = val[i] ^ ((val[i] << 16) & 0xffff0000);

    if (func.n > 5)
    {
        for (unsigned int i = 0; i < func.n - 5; ++i)
        {
            unsigned int half = 1 << i;

            for (unsigned int j = 0; j < vec_size; j = j + 2 * half)
                for (unsigned int k = j + half; k < (j + 2 * half); k++)
                    val[k] = val[k] ^ val[k - half];
        }
    }
    BF res;
    res.vec = val;
    res.n = func.n;
    return res;
}

void BF::print_per_bit()
{
    for (int i = 0; i < this->vec.size(); ++i)
        std::cout << std::bitset<32>(this->vec[i]) << std::endl;
}

bool BF::operator==(const BF &func) const
{
    if (this->n == func.n)
    {
        unsigned int size = func.vec.size();
        int i = 0;
        for (; i < size; i++)
        {
            if (this->vec[i] != func.vec[i])
                return false;
        }
        return true;
    }
    std::cout << "the number of variables in operator == is different" << '\n';
    return false;
}

bool BF::operator!=(const BF &func) const
{
    if (this->n == func.n)
    {
        unsigned int size = func.vec.size();
        int i = 0;
        for (; i < size; i++)
        {
            if (this->vec[i] != func.vec[i])
                return true;
        }
        return false;
    }
    std::cout << "the number of variables in operator != is different" << '\n';
    return false;
}

inline bool bit_value(unsigned int num, unsigned char bit)
{
    return (uint32_t(1) << (bit)) & num;
}

inline void SetBit(uint32_t &num, uint8_t bit, bool is_one)
{
    if (is_one)
        num |= (uint32_t)1 << bit;
    else
    {
        uint32_t mask = (uint32_t)1 << bit;
        mask = ~mask;
        num &= mask;
    }
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

void BF::print_ANF(const BF &func)
{
    unsigned int curr = 0, size = func.vec.size();
    int count_zero = std::count_if(func.vec.begin(), func.vec.end(), [](int i)
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
            if (bit_value(func.vec[i], j))
            {
                print_monom(curr);
                std::cout << " + ";
            }
        }
}

unsigned int BF::degree(const BF &func)
{

    unsigned int max = 0, size = func.vec.size();

    unsigned int weight = 0, curr_val = (1 << func.n) - 1;
    for (int i = size - 1; i >= 0; i--)
    {

        for (int j = 31; j >= 0; j--)
        {
            curr_val = i << 5;
            curr_val += j;
            if (bit_value(func.vec[i], j))
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
    unsigned int size = func.vec.size();

    std::vector<int> char_vec;
    // make characteristic vector
    for (unsigned int i = 0; i < size; ++i)
    {
        for (unsigned int j = 0; j < 32; j++)
        {
            if (bit_value(func.vec[i], j))
                char_vec.push_back(-1);
            else
                char_vec.push_back(1);
        }
    }
    if (func.n < 5)
    {
        for (int i = 0; i < 32 - (1 << func.n); i++)
            char_vec.pop_back();
    }

    unsigned int char_vec_size = char_vec.size();
    int sum = 0, dif = 0;

    for (unsigned int i = 0; i < func.n; ++i)
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
    // 🤢🤢🤢🤢🤢
    return char_vec;
}

unsigned int BF::COR(const BF &func, std::vector<int> &WHT_coef)
{

    auto gen_limits = [](unsigned int n, unsigned int k) -> std::pair<unsigned int, unsigned int>
    {
        return std::make_pair((uint32_t)((1 << k) - 1) << (n - k), 0xFFFFFFFF >> (32 - k));
    };

    if (WHT_coef.empty())
        WHT_coef = WH_transform(func);

    int i = 1;

    for (; i <= func.n; i++)
    {
        auto limits = gen_limits(func.n, i);
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
    return func.n;
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

    return (1 << (func.n - 1)) - (max >> 1);
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

bool BF::get_val(const BF &func, unsigned int x_set)
{
    unsigned int ix = x_set >> 5;
    unsigned int bit = x_set & 0x0000001f;
    return bit_value(func.vec[ix], bit);
}

bool BF::operator[](unsigned int x_set)
{
    unsigned int ix = x_set >> 5;
    unsigned int bit = x_set & 0x0000001f;
    return bit_value(this->vec[ix], bit);
}
void BF::set_val(unsigned int x_set, bool val)
{

    unsigned int ix = x_set >> 5;
    unsigned int bit = x_set & 0x0000001f;
    unsigned int mask = 1;
    mask = mask << bit;
    mask = ~mask;
    if (val)
        this->vec[ix] |= (unsigned int)1 << bit;
    else
        this->vec[ix] &= mask;
}

std::vector<int> BF::Auto_Cor(const BF &func, std::vector<int> wht_coef)
{
    if (wht_coef.empty())
        wht_coef = BF::WH_transform(func);
    unsigned int size = wht_coef.size();
    // for (int i = 0; i < size; i++)
    //    wht_coef[i] *= wht_coef[i];

    std::for_each(wht_coef.begin(), wht_coef.end(), [](int &n)
                  { n *= n; });

    int sum = 0, dif = 0;

    for (unsigned int i = 0; i < func.n; ++i)
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
    unsigned int n_variables = func.n;
    std::for_each(wht_coef.begin(), wht_coef.end(), [&n_variables](int &n)
                  { n >>= n_variables; });

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

    // std::cout << "MAX_autocor= " << max << "\n\n";

    return ((1 << (func.n - 2)) - (max >> 2));
}

unsigned int BF::PC(const BF &func, std::vector<int> &auto_cor)
{
    if (auto_cor.empty())
        auto_cor = BF::Auto_Cor(func, std::vector<int>());

    auto gen_limits = [](unsigned int n, unsigned int k) -> std::pair<unsigned int, unsigned int>
    {
        return std::make_pair((uint32_t)((1 << k) - 1) << (n - k), 0xFFFFFFFF >> (32 - k));
    };

    int i = 1;

    for (; i <= func.n; i++)
    {
        auto limits = gen_limits(func.n, i);
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
    return func.n;
}

int CountBoolMatrixRankAndStepTransform(std::vector<unsigned int> &Matrix, int columns)
{
    int Rank = 0;

    columns--;
    for (int i = 0; i < Matrix.size() && (columns >= 0); i++)
    {
        auto max = std::max_element(Matrix.begin() + i, Matrix.end());
        // std::cout << "i=" << i << " max= " << std::bitset<32>(*max) << '\n';
        if (*max)
        {
            std::swap(*(Matrix.begin() + i), *(max));

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

    // std::for_each(Matrix.begin(), Matrix.end(), [](unsigned int num)
    //               { std::cout << std::bitset<32>(num) << '\n'; });

    while (*(Matrix.end() - 1) == 0)
        Matrix.pop_back();

    Rank = Matrix.size();
    return Rank;
}

int getElderBitPos(unsigned int num, int start = 31)
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

int CountBoolMatrixRankAndStepTransform(std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns)
{
    int Rank = 0;
    int i = 0;
    columns--;
    int numberOfVariables = columns;

    for (; i < PairMatrix.size() && (columns >= 0); i++)
    {
        auto max = std::max_element(PairMatrix.begin() + i, PairMatrix.end(), [](std::pair<uint32_t, uint8_t> pair_1, std::pair<uint32_t, uint8_t> pair_2)
                                    { return pair_1.first < pair_2.first; });
        if (max->first)
        {
            std::swap(*(PairMatrix.begin() + i), *(max));
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

    while ((PairMatrix.end() - 1)->first == 0 && (PairMatrix.end() - 1)->second == 0)
        PairMatrix.pop_back();

    // std::for_each(PairMatrix.begin(), PairMatrix.end(), [](std::pair<uint32_t, uint8_t> pair)
    //               { std::cout << std::bitset<32>(pair.first) << "|" << (uint16_t)pair.second << '\n'; });

    Rank = PairMatrix.size();

    // Here after step matrix make special step matrix
    // uint32_t mainVarVec = 0;
    for (int k = PairMatrix.size() - 1; k >= 0; k--)
    {
        int elderBitPos = getElderBitPos(PairMatrix[k].first, numberOfVariables);
        // SetBit(mainVarVec, elderBitPos, 1);
        if (elderBitPos == -1)
        {
            //std::cout << "No solution\n";
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
std::vector<uint32_t> getSolutionsOfSystem(const std::vector<std::pair<uint32_t, uint8_t>> &PairMatrix, int n)
{
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

std::pair<std::vector<unsigned int>, std::vector<unsigned int>> BF::OptGoodPairsVec(BF &func, std::vector<int> &wht_coef)
{
    const unsigned int n = func.n;
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

    int Rank = CountBoolMatrixRankAndStepTransform(PairMatrix, func.n);

    if (Rank == -1)
    {
        //std::cout << "no solution for equations\n";
        return std::pair<std::vector<unsigned int>, std::vector<unsigned int>>();
    }

    auto solutionsFirst = getSolutionsOfSystem(PairMatrix, func.n);

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

    Rank = CountBoolMatrixRankAndStepTransform(PairMatrix, func.n);

    if (Rank == -1)
    {
        //std::cout << "no solution for equations\n";
        return std::pair<std::vector<unsigned int>, std::vector<unsigned int>>();
    }

    auto solutionsSecond = getSolutionsOfSystem(PairMatrix, func.n);
    std::vector<uint32_t> finalSolutionsSecond;

    for (int i = 0; i < solutionsSecond.size(); i++)
        if (func[solutionsSecond[i]] == 1)
            finalSolutionsSecond.push_back(solutionsSecond[i]);
    std::pair<std::vector<unsigned int>, std::vector<unsigned int>> ResultGoodPairs;
    ResultGoodPairs=std::make_pair(finalSolutionsFirst,finalSolutionsSecond);
    

    return ResultGoodPairs;
}

void BF::FillWSets(std::vector<uint32_t> &W_1_pos, std::vector<uint32_t> &W_1_neg, std::vector<uint32_t> &W_3_pos, std::vector<uint32_t> &W_3_neg, std::vector<int> &wht_coef)
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
            continue;
        }
        if (wht_coef[i] == (max * (-1)) + 4)
        {
            W_3_neg.push_back(i);
            continue;
        }
        if (wht_coef[i] == max)
        {
            W_1_pos.push_back(i);
            continue;
        }
        if (wht_coef[i] == max * (-1))
        {
            W_1_neg.push_back(i);
            continue;
        }
    }
}

BF BF::SwapOnSets(BF &func, uint32_t set_x1, uint32_t set_x2)
{
    auto result = func;
    result.set_val(set_x1, func[set_x2]);
    result.set_val(set_x2, func[set_x1]);
    return result;
}

BF BF::GenBalancedFunc(int numberOfVariables)
{
    BF func(numberOfVariables, 10);

    int32_t diff = ((uint32_t)1 << (numberOfVariables - 1)) - BF::weight(func);

    while (diff < 0)
    {
        uint32_t mask = 0;
        mask = ~mask;
        mask >>= 32 - numberOfVariables;
        uint32_t x_set = rand() & mask;
        if (func[x_set] == true)
        {
            func.set_val(x_set, false);
            diff++;
        }
    }
    while (diff > 0)
    {
        uint32_t mask = 0;
        mask = ~mask;
        mask >>= 32 - numberOfVariables;
        uint32_t x_set = rand() & mask;
        if (func[x_set] == false)
        {
            func.set_val(x_set, true);
            diff--;
        }
    }
    return func;
}

std::pair<std::vector<unsigned int>, std::vector<unsigned int>> BF::good_pairsVec(BF &func, std::vector<int> &wht_coef)
{
    const unsigned int n = func.n;
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
    for (unsigned int curr = 0; curr < 1 << n; curr++)
    {
        if (func[curr] == false)
        {
            int a = 0;
            for (; a < W_1_pos.size(); a++)
                if (weight_mod(W_1_pos[a] & curr) == true)
                    break;
            if (a == W_1_pos.size())
                A_00.push_back(curr);

            int b = 0;
            for (; b < W_1_neg.size(); b++)
                if (weight_mod(W_1_neg[b] & curr) == false)
                    break;
            if (b == W_1_neg.size())
                B_01.push_back(curr);
        }
        else
        {
            int a = 0;
            for (; a < W_1_pos.size(); a++)
                if (weight_mod(W_1_pos[a] & curr) == false)
                    break;
            if (a == W_1_pos.size())
                A_11.push_back(curr);

            int b = 0;
            for (; b < W_1_neg.size(); b++)
                if (weight_mod(W_1_neg[b] & curr) == true)
                    break;
            if (b == W_1_neg.size())
                B_10.push_back(curr);
        }
    }

    std::vector<unsigned int> AB_0, AB_1;
    std::set_intersection(A_00.begin(), A_00.end(), B_01.begin(), B_01.end(), std::insert_iterator<std::vector<unsigned int>>(AB_0, AB_0.begin()));
    std::set_intersection(A_11.begin(), A_11.end(), B_10.begin(), B_10.end(), std::insert_iterator<std::vector<unsigned int>>(AB_1, AB_1.begin()));

    return std::make_pair(AB_0, AB_1);
}