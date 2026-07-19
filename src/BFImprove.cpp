#include "BF.h"
#include "Utils.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <unordered_set>

namespace bf
{

std::pair<std::vector<unsigned int>, std::vector<unsigned int>> BF::optGoodPairsVec(std::vector<int> &wht_coef) const
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

    int Rank = countBoolMatrixRankAndStepTransform(PairMatrix, n_);

    if (Rank == -1)
    {
        return std::pair<std::vector<unsigned int>, std::vector<unsigned int>>();
    }

    auto solutionsFirst = getSolutionsOfSystem(PairMatrix, n_);

    std::vector<uint32_t> finalSolutionsFirst;

    for (int i = 0; i < solutionsFirst.size(); i++)
        if ((*this)[solutionsFirst[i]] == 0)
            finalSolutionsFirst.push_back(solutionsFirst[i]);

    PairMatrix.clear();

    // now we found the first element of good pairs
    for (int i = 0; i < W_1_pos.size(); i++)
        PairMatrix.push_back(std::make_pair(W_1_pos[i], 1));
    for (int i = 0; i < W_1_neg.size(); i++)
        PairMatrix.push_back(std::make_pair(W_1_neg[i], 0));

    Rank = countBoolMatrixRankAndStepTransform(PairMatrix, n_);

    if (Rank == -1)
    {
        return std::pair<std::vector<unsigned int>, std::vector<unsigned int>>();
    }

    auto solutionsSecond = getSolutionsOfSystem(PairMatrix, n_);
    std::vector<uint32_t> finalSolutionsSecond;

    for (int i = 0; i < solutionsSecond.size(); i++)
        if ((*this)[solutionsSecond[i]] == 1)
            finalSolutionsSecond.push_back(solutionsSecond[i]);
    std::pair<std::vector<unsigned int>, std::vector<unsigned int>> ResultGoodPairs;
    ResultGoodPairs = std::make_pair(finalSolutionsFirst, finalSolutionsSecond);

    return ResultGoodPairs;
}

void BF::fillWSets(std::vector<uint32_t> &W_1_pos, std::vector<uint32_t> &W_1_neg, std::vector<uint32_t> &W_3_pos, std::vector<uint32_t> &W_3_neg, const std::vector<int> &wht_coef)
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

BF BF::swapOnSets(uint32_t set_x1, uint32_t set_x2) const
{
    if (set_x1 >= ((uint32_t)1 << n_) || set_x2 >= ((uint32_t)1 << n_))
    {
        throw "pair is out of range";
    }
    auto result = *this;
    result.setValue(set_x1, (*this)[set_x2]);
    result.setValue(set_x2, (*this)[set_x1]);
    return result;
}

std::pair<std::vector<unsigned int>, std::vector<unsigned int>> BF::goodPairsVec(std::vector<int> &wht_coef) const
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
    std::vector<unsigned int> A_00, B_01, A_11, B_10;
    for (unsigned int curr = 0; curr < 1 << n_; curr++)
    {
        if ((*this)[curr] == false)
        {
            int a = 0;
            for (; a < W_1_pos.size(); a++)
                if (weightMod(W_1_pos[a] & curr) == true)
                    break;
            if (a == W_1_pos.size() && W_1_pos.size())
                A_00.push_back(curr);

            int b = 0;
            for (; b < W_1_neg.size(); b++)
                if (weightMod(W_1_neg[b] & curr) == false)
                    break;
            if (b == W_1_neg.size() && W_1_neg.size())
                B_01.push_back(curr);
        }
        else
        {
            int a = 0;
            for (; a < W_1_pos.size(); a++)
                if (weightMod(W_1_pos[a] & curr) == false)
                    break;
            if (a == W_1_pos.size() && W_1_pos.size())
                A_11.push_back(curr);

            int b = 0;
            for (; b < W_1_neg.size(); b++)
                if (weightMod(W_1_neg[b] & curr) == true)
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


std::vector<std::pair<uint32_t, uint32_t>> BF::pairsToWorsen() const
{
    auto wht_coef = WHTransform();
    std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
    BF::fillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);

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
                if (weightMod(curr & W_1_pos[i]) == true)
                    A_01[i].push_back(curr);
            }
            for (int i = 0; i < W_1_neg.size(); ++i)
            {
                if (weightMod(curr & W_1_neg[i]) == false)
                    B_00[i].push_back(curr);
            }
        }
        else
        {
            for (int i = 0; i < W_1_pos.size(); ++i)
            {
                if (weightMod(curr & W_1_pos[i]) == false)
                    A_10[i].push_back(curr);
            }
            for (int i = 0; i < W_1_neg.size(); ++i)
            {
                if (weightMod(curr & W_1_neg[i]) == true)
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
                    Allvalues.insert(convertToUint64(A_01[i][j], A_10[i][k]));
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
                    Allvalues.insert(convertToUint64(B_00[i][j], B_11[i][k]));
                }
            }
        }
    }

    std::vector<std::pair<uint32_t, uint32_t>> result;
    for (auto &iter : Allvalues)
    {
        result.push_back(convertToPairUint32T(iter));
    }

    return result;
}

std::vector<std::pair<uint32_t, uint32_t>> BF::pairsToImprove() const
{ // построить множесто пар улучшающих нелинейность
    auto wht_coef = WHTransform();
    auto GoodPairs = optGoodPairsVec(wht_coef); // Вернуть пару множеств(2 множества из теоремы слева от разности) с помощью решения СЛУ

    std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
    BF::fillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);

    std::vector<std::pair<uint32_t, uint32_t>> pairs = convertToPairs(GoodPairs); // конвертируем в декартово произведение
    std::vector<std::pair<uint32_t, uint32_t>> res;
    for (uint64_t i = 0; i < pairs.size(); i++) // проверяю удовлетворяет ли пара условиям для улучшение нелинейность
    {
        int flag = 0;

        for (int j = 0; j < W_3_pos.size(); j++)
        {
            bool correct = weightMod(W_3_pos[j] & pairs[i].first) == false || weightMod(W_3_pos[j] & pairs[i].second) == true; // (с,x1)=0 ⋁ (с,x2)=1
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
            bool correct = weightMod(W_3_neg[j] & pairs[i].first) == true || weightMod(W_3_neg[j] & pairs[i].second) == false; //(d,x1)=1 ⋁ (d,x2)=0
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
    auto wht_coef = WHTransform();
    auto GoodPairs = goodPairsVec(wht_coef);

    std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
    BF::fillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);

    std::vector<std::pair<uint32_t, uint32_t>> pairs = convertToPairs(GoodPairs), res;
    for (uint64_t i = 0; i < pairs.size(); i++)
    {
        int flag = 0;

        for (int j = 0; j < W_3_pos.size(); j++)
        {
            bool correct = weightMod(W_3_pos[j] & pairs[i].first) == false || weightMod(W_3_pos[j] & pairs[i].second) == true;
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
            bool correct = weightMod(W_3_neg[j] & pairs[i].first) == true || weightMod(W_3_neg[j] & pairs[i].second) == false;
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
        bool isWorsening = weightMod(i & pair.first) == true && weightMod(i & pair.second) == false; // (w,x1)==1 && (w,x2)==0
        if (isWorsening)
            return false;
    }
    for (auto &i : W_1_neg)
    {
        bool isWorsening = weightMod(i & pair.first) == false && weightMod(i & pair.second) == true; // (w,x1)==0 && (w,x2)==1
        if (isWorsening)
            return false;
    }

    return true;
}


bool BF::isImprovePair(const std::vector<uint32_t> &W_1_pos, const std::vector<uint32_t> &W_1_neg, const std::vector<uint32_t> &W_3_pos, const std::vector<uint32_t> &W_3_neg, const std::pair<uint32_t, uint32_t> &pair) const
{

    for (auto &i : W_1_pos)
    {
        bool correct = weightMod(i & pair.first) == false && weightMod(i & pair.second) == true;
        if (!correct)
            return false;
    }

    for (auto &i : W_1_neg)
    {
        bool correct = weightMod(i & pair.first) == true && weightMod(i & pair.second) == false;
        if (!correct)
            return false;
    }

    for (int j = 0; j < W_3_pos.size(); j++)
    {
        bool correct = weightMod(W_3_pos[j] & pair.first) == false || weightMod(W_3_pos[j] & pair.second) == true; // (с,x1)=0 ⋁ (с,x2)=1
        if (!correct)                                                                                                //
            return false;
    }

    for (int j = 0; j < W_3_neg.size(); j++)
    {
        bool correct = weightMod(W_3_neg[j] & pair.first) == true || weightMod(W_3_neg[j] & pair.second) == false; //(d,x1)=1 ⋁ (d,x2)=0
        if (!correct)
            return false;
    }

    return true; // Пара удовлетворила всем условиям -> она улучшает нелинейность
}

bool BF::isImprovePairEasy(const std::vector<int> &wht_coef, const std::pair<uint32_t, uint32_t> &pair) const
{
    std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
    BF::fillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);
    for (int j = 0; j < W_3_pos.size(); j++)
    {
        bool correct = weightMod(W_3_pos[j] & pair.first) == false || weightMod(W_3_pos[j] & pair.second) == true; // (с,x1)=0 ⋁ (с,x2)=1
        if (!correct)                                                                                                //
            return false;
    }

    for (int j = 0; j < W_3_neg.size(); j++)
    {
        bool correct = weightMod(W_3_neg[j] & pair.first) == true || weightMod(W_3_neg[j] & pair.second) == false; //(d,x1)=1 ⋁ (d,x2)=0
        if (!correct)
            return false;
    }

    return false;
}

void BF::nonlinearityImprove(uint64_t linDiff, uint64_t &neutralAttemps, uint64_t &improveAttemps)
{
    if (linDiff % 2 != 0)
    {
        return;
    }
    auto wht_coef = WHTransform();

    auto tmp = std::roundf((1 << (this->n_ - 1)) - std::pow(2, this->n_ / 2 - 1));
    int Nf_boundaries = this->n_ % 2 == 0 ? tmp - 2 : tmp;

    std::pair<uint32_t, uint32_t> neutralPair;
    std::pair<uint32_t, uint32_t> pair;
    for (int j = 0; j <= improveAttemps && linDiff > 0; ++j)
    {

        std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
        wht_coef = WHTransform();
        BF::fillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);
        auto NF = nonlinearity(wht_coef);
        if (nonlinearity(wht_coef) >= Nf_boundaries)
            break;

        bool needToBuildImproveSet = true;
        std::vector<std::pair<uint32_t, uint32_t>> improvePairs;

        for (int i = 0; i < 100000; ++i)
        {
            pair = this->generatePair();
            if (this->isImprovePair(W_1_pos, W_1_neg, W_3_pos, W_3_neg, pair))
            {
                *this = swapOnSets(pair.first, pair.second);
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
            *this = swapOnSets(pair.first, pair.second);
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
                *this = swapOnSets(neutralPair.first, neutralPair.second);
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

} // namespace bf
