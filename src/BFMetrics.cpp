#include "BF.h"
#include "Utils.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace bf
{

unsigned int BF::weight() const
{
    unsigned int size = vec_.size();
    uint32_t total_weight = 0;

    for (int i = 0; i < size; ++i)
    {
        total_weight += uintWeight(vec_[i]);
    }

    return total_weight;
}

void BF::printANF() const
{
    unsigned int curr = 0, size = vec_.size();
    int count_zero = std::count_if(vec_.begin(), vec_.end(), [](int i)
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
            if (bitValue(vec_[i], j))
            {
                printMonom(curr);
                std::cout << " + ";
            }
        }
}

unsigned int BF::degree() const
{

    unsigned int max = 0, size = vec_.size();

    unsigned int weight = 0, curr_val = 0;
    for (int i = size - 1; i >= 0; i--)
    {

        for (int j = 31; j >= 0; j--)
        {
            curr_val = i << 5;
            curr_val += j;
            if (bitValue(vec_[i], j))
            {
                weight = uintWeight(curr_val);
                if (weight > max)
                    max = weight;
            }
        }
    }
    return max;
}

unsigned int BF::correlationImmunity(std::vector<int> &WHT_coef) const
{

    auto gen_limits = [](unsigned int n_, unsigned int k) -> std::pair<unsigned int, unsigned int>
    {
        return std::make_pair((uint32_t)((1 << k) - 1) << (n_ - k), 0xFFFFFFFF >> (32 - k));
    };

    if (WHT_coef.empty())
        WHT_coef = WHTransform();

    int i = 1;

    for (; i <= n_; i++)
    {
        auto limits = gen_limits(n_, i);
        unsigned int current = limits.first;

        if (WHT_coef[current] != 0)
            return i - 1;

        while (current != limits.second)
        {
            current = nextCombination(current);
            if (WHT_coef[current] != 0)
                return i - 1;
        }
    }
    return n_;
}

unsigned int BF::nonlinearity(std::vector<int> &WHT) const
{
    if (WHT.empty())
        WHT = WHTransform();
    auto min_max = std::minmax_element(WHT.begin(), WHT.end());
    int max = 0;
    if (abs(*min_max.first) > abs(*min_max.second))
        max = abs(*min_max.first);
    else
        max = abs(*min_max.second);

    return (1 << (n_ - 1)) - (max >> 1);
}

void BF::printBAA(std::vector<int> &WHT) const
{
    if (WHT.empty())
        WHT = WHTransform();

    auto min_max = std::minmax_element(WHT.begin(), WHT.end());
    unsigned int max_pos = 0;
    if (abs(*min_max.first) > abs(*min_max.second))
        max_pos = min_max.first - WHT.begin();
    else
        max_pos = min_max.second - WHT.begin();

    bool is_max_positive = WHT[max_pos] > 0 ? false : true;

    for (unsigned int i = 0; i < 31; i++)
        if (bitValue(max_pos, i))
            std::cout << "x" << i << "+";
    if (is_max_positive)
        std::cout << "1";
}

int BF::cnF(std::vector<int> &auto_cor) const
{
    if (auto_cor.empty())
        auto_cor = autoCor(std::vector<int>());

    auto min_max = std::minmax_element(auto_cor.begin() + 1, auto_cor.end());
    int max = 0;
    if (abs(*min_max.first) > abs(*min_max.second))
        max = abs(*min_max.first);
    else
        max = abs(*min_max.second);

    return ((1 << (n_ - 2)) - (max >> 2));
}

unsigned int BF::propagationCriteria(std::vector<int> &auto_cor) const
{
    if (auto_cor.empty())
        auto_cor = autoCor(std::vector<int>());

    auto gen_limits = [](unsigned int n_, unsigned int k) -> std::pair<unsigned int, unsigned int>
    {
        return std::make_pair((uint32_t)((1 << k) - 1) << (n_ - k), 0xFFFFFFFF >> (32 - k));
    };

    int i = 1;

    for (; i <= n_; i++)
    {
        auto limits = gen_limits(n_, i);
        unsigned int current = limits.first;

        if (auto_cor[current] != 0)
            return i - 1;

        while (current != limits.second)
        {
            current = nextCombination(current);
            if (auto_cor[current] != 0)
                return i - 1;
        }
    }
    return n_;
}

} // namespace bf
