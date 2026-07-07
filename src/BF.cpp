#include "BF.h"
#include "Utils.h"

#include <bitset>
#include <cstdlib>
#include <iostream>
#include <limits>

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
    this->n_ = countFirstZeros(len);
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

void BF::printPerBit() const
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

bool BF::getValue(const BF &func, unsigned int x_set)
{
    unsigned int ix = x_set >> 5;
    unsigned int bit = x_set & 0x0000001f;
    return bitValue(func.vec_[ix], bit);
}

bool BF::operator[](unsigned int x_set) const
{
    unsigned int ix = x_set >> 5;
    unsigned int bit = x_set & 0x0000001f;
    return bitValue(this->vec_[ix], bit);
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

