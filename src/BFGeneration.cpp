#include "BF.h"
#include "Utils.h"

#include <bitset>
#include <cstdlib>

std::pair<BF, BF> BF::generateBorderBalancesFunctions()
{
    return std::make_pair(BF(std::bitset<32>(0xFFFFFFFF >> (16)).to_string()), BF(std::bitset<32>(((1 << 16) - 1) << (16)).to_string()));
}

void BF::nextBalanced()
{
    this->vec_[0] = nextCombination(this->vec_[0]);
}

BF BF::genBalancedFunc(int numberOfVariables) // Генерируем функцию а потом подправляем ее до уравновешенной
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
    affineFunc = BF::mobiusTransform(affineFunc);
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
    auto wht_coef = BF::WHTransform(*this);
    auto goodSet = BF::optGoodPairsVec(*this, wht_coef);
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
