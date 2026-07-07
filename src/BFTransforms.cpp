#include "BF.h"
#include "Utils.h"

#include <algorithm>

BF BF::mobiusTransform(const BF &func)
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

std::vector<int> BF::WHTransform(const BF &func)
{
    unsigned int size = func.vec_.size();

    std::vector<int> char_vec;
    // make characteristic vector
    for (unsigned int i = 0; i < size; ++i)
    {
        for (unsigned int j = 0; j < 32; j++)
        {
            if (bitValue(func.vec_[i], j))
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

std::vector<int> BF::autoCor(const BF &func, std::vector<int> wht_coef)
{
    if (wht_coef.empty())
        wht_coef = BF::WHTransform(func);
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

