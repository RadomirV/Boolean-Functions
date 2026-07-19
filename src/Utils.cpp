#include "Utils.h"

#include <algorithm>
#include <iostream>


namespace bf
{

uint32_t countFirstZeros(uint32_t x)
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

bool bitValue(unsigned int num, unsigned char bit)
{
    return ((uint32_t(1) << (bit)) & num) != 0;
}
void setBit(uint32_t &num, uint8_t bit, bool is_one)
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
        if (bitValue(num, i))
            break;
        i--;
    }
    return i;
}
void printMonom(unsigned int monom)
{
    if (monom == 0)
    {
        std::cout << 1;
        return;
    }
    for (unsigned int i = 0; i < 32; i++)
        if (bitValue(monom, i))
            std::cout << "x" << i;
}
uint64_t convertToUint64(uint32_t num1, uint32_t num2)
{
    return static_cast<uint64_t>((static_cast<uint64_t>(num1) << 32) | static_cast<uint64_t>(num2));
}
std::pair<uint32_t, uint32_t> convertToPairUint32T(uint64_t num)
{
    uint32_t num0 = static_cast<uint32_t>(num >> 32);
    uint32_t num1 = static_cast<uint32_t>(num);
    return std::make_pair(num0, num1);
}
unsigned int uintWeight(unsigned int x)
{
    x -= ((x >> 1) & 0x55555555);
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    x = (x + (x >> 4)) & 0x0f0f0f0f;
    x += (x >> 8);
    x += (x >> 16);
    x = x & 0x3f;
    return x;
}
bool weightMod(unsigned int x)
{
    x = x ^ (x >> 1);
    x = x ^ (x >> 2);
    x = x ^ (x >> 4);
    x = x ^ (x >> 8);
    x = x ^ (x >> 16);
    return (x & 0x00000001) == 1;
}

unsigned int nextCombination(unsigned int prev)
{
    unsigned int b = (prev + 1) & prev;
    unsigned int c = uintWeight((b - 1) ^ prev) - 2;
    return (((((prev + 1) ^ prev) << 1) + 1) << c) ^ b;
}
int countBoolMatrixRankAndStepTransform(std::vector<unsigned int> &Matrix, int columns)
{
    int Rank = 0;

    columns--;
    for (int i = 0; i < Matrix.size() && (columns >= 0); i++)
    {
        auto max = std::max_element(Matrix.begin() + i, Matrix.end());
        if (*max)
        {
            std::swap(Matrix[i], *max);

            while (bitValue(Matrix[i], columns) == 0)
            {
                columns--;
            }

            for (int j = i + 1; j < Matrix.size(); j++) // annigilate all 1 under max in column
                if (bitValue(Matrix[j], columns))
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

std::vector<std::pair<uint32_t, uint32_t>> convertToPairs(const std::pair<std::vector<unsigned int>, std::vector<unsigned int>> &GPVec)
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
int countBoolMatrixRankAndStepTransform(std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns)
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
            while (bitValue(PairMatrix[i].first, columns) == 0)
            {
                columns--;
            }
            for (int j = i + 1; j < PairMatrix.size(); j++) // annigilate all 1 under max
                if (bitValue(PairMatrix[j].first, columns))
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
        // setBit(mainVarVec, elderBitPos, 1);
        if (elderBitPos == -1)
        {
            PairMatrix = std::vector<std::pair<unsigned int, uint8_t>>();
            return -1;
        }

        for (int j = k - 1; j >= 0; j--)
        {
            if (bitValue(PairMatrix[j].first, elderBitPos))
            {
                PairMatrix[j].first ^= PairMatrix[k].first;
                PairMatrix[j].second ^= PairMatrix[k].second;
            }
        }
    }

    return Rank;
}

bool isMatrixSpecialStepTransform(const std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns)
{
    if (PairMatrix.empty() || columns < 0)
        return false;
    for (auto it = PairMatrix.begin(); it != PairMatrix.end(); it++)
    {
        int elderPos = getElderBitPos(it->first, columns);
        for (auto it2 = it + 1; it2 != PairMatrix.end(); it2++)
        {
            if (bitValue(it2->first, elderPos))
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
            if (bitValue(it2->first, elderPos))
                return false;
            if (it2 == PairMatrix.begin())
                break;
        }
    }

    return true;
}
std::vector<uint32_t> getSolutionsOfSystem(const std::vector<std::pair<uint32_t, uint8_t>> &PairMatrix, int n)
{
    if (!isMatrixSpecialStepTransform(PairMatrix, n))
    {
        return std::vector<uint32_t>();
    }
    uint32_t mainVec = 0, dependendVec = 0;
    for (int i = 0; i < PairMatrix.size(); i++)
    {
        int elderBitPos = getElderBitPos(PairMatrix[i].first, n - 1);
        setBit(mainVec, elderBitPos, 1);
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

            uint32_t bitToSet = ((uint32_t)weightMod(r) + PairMatrix[i].second) % 2;

            bool isXone = bitToSet;

            setBit(currentSolution, elderBitPos, isXone);
        }
        solutions.push_back(currentSolution);
        if (currPreviousToDependentVec == 0)
            break;
        currPreviousToDependentVec = (currPreviousToDependentVec - 1) & dependendVec;
    }

    return solutions;
}

} // namespace bf
