#include <gtest/gtest.h>

#include <unordered_set>

#include "BF.h"

TEST(BFUtilsTest, BitOperations)
{
    EXPECT_EQ(count_first_zeros(0), 0u);
    EXPECT_EQ(count_first_zeros(1), 0u);
    EXPECT_EQ(count_first_zeros(8), 3u);
    EXPECT_EQ(count_first_zeros(1u << 16), 16u);
    EXPECT_EQ(count_first_zeros(0xFFFF), 0u);

    EXPECT_EQ(getElderBitPos(0xFFFFu, 12), 12);
    EXPECT_EQ(getElderBitPos(0xFFFFu, 25), 15);
    EXPECT_EQ(getElderBitPos(1u << 16, 16), 16);
    EXPECT_EQ(getElderBitPos(0u, 31), -1);
    EXPECT_EQ(getElderBitPos(0b100100u, 4), 2);

    EXPECT_FALSE(bit_value(0b1010, 0));
    EXPECT_TRUE(bit_value(0b1010, 1));
    EXPECT_FALSE(bit_value(0b1010, 2));
    EXPECT_TRUE(bit_value(0b1010, 3));

    uint32_t testNum = 0;
    SetBit(testNum, 3, true);
    EXPECT_EQ(testNum, 0b1000u);
    SetBit(testNum, 1, true);
    EXPECT_EQ(testNum, 0b1010u);
    SetBit(testNum, 3, false);
    EXPECT_EQ(testNum, 0b0010u);

    testing::internal::CaptureStdout();
    print_monom(0b110u);
    std::string out = testing::internal::GetCapturedStdout();
    EXPECT_EQ(out, "x1x2");

    testing::internal::CaptureStdout();
    print_monom(0b1101101u);
    out = testing::internal::GetCapturedStdout();
    EXPECT_EQ(out, "x0x2x3x5x6");
}

TEST(BFUtilsTest, NextCombinationKeepsBitCount)
{
    // Valid iteration range used by COR/PC: from 11000 to 00011 (n=5, k=2).
    EXPECT_EQ(next_combination(0b11000u), 0b10100u);
    EXPECT_EQ(next_combination(0b10100u), 0b10010u);
    EXPECT_EQ(next_combination(0b10010u), 0b10001u);
    EXPECT_EQ(next_combination(0b10001u), 0b01100u);
}

TEST(BFUtilsTest, WeightOperations)
{
    EXPECT_EQ(uint_weight(0u), 0u);
    EXPECT_EQ(uint_weight(1u), 1u);
    EXPECT_EQ(uint_weight(0b10101010u), 4u);
    EXPECT_EQ(uint_weight(0xFFFFFFFFu), 32u);
    EXPECT_EQ(uint_weight(0x80000000u), 1u);

    EXPECT_FALSE(weight_mod(0u));        // even parity
    EXPECT_TRUE(weight_mod(1u));         // odd parity
    EXPECT_FALSE(weight_mod(0b1010u));   // 2 ones -> even parity
    EXPECT_TRUE(weight_mod(0b1011u));    // 3 ones -> odd parity
    EXPECT_FALSE(weight_mod(0xFFFFFFu)); // 32 ones -> even parity
}
TEST(BFUtilsTest, Convertions)
{
    EXPECT_EQ(convertTo_uint64(0x12345678u, 0x9ABCDEF0u), 0x123456789ABCDEF0ull);
    auto pair = convertTo_pair_uint32_t(0x123456789ABCDEF0ull);
    EXPECT_EQ(pair.first, 0x12345678u);
    EXPECT_EQ(pair.second, 0x9ABCDEF0u);

    std::pair<std::vector<unsigned int>, std::vector<unsigned int>> pairVec = {{0x1u, 0x2u, 0x3u}, {0x4u, 0x5u}};
    auto convertedPairs = ConvertToPairs(pairVec);
    EXPECT_EQ(convertedPairs.size(), pairVec.first.size() * pairVec.second.size());
    EXPECT_EQ(convertedPairs[0], std::make_pair(0x1u, 0x4u));
    EXPECT_EQ(convertedPairs[1], std::make_pair(0x1u, 0x5u));
    EXPECT_EQ(convertedPairs.back(), std::make_pair(0x3u, 0x5u));
}

TEST(BFUtilsTest, MatrixRankAndStepTransform)
{
    {
        std::vector<unsigned int> matrix = {0b110u, 0b101u, 0b011u};
        std::vector<unsigned int> transformed = {0b110u, 0b011u};
        int rank = CountBoolMatrixRankAndStepTransform(matrix, 3);
        EXPECT_EQ(rank, 2);
        EXPECT_EQ(matrix, transformed);
        EXPECT_EQ(matrix.size(), 2u);
    }

    {
        std::vector<unsigned int> matrix = {0b110u, 0b110u, 0b110u};
        int rank = CountBoolMatrixRankAndStepTransform(matrix, 3);
        EXPECT_EQ(rank, 1);
        EXPECT_EQ(matrix.size(), 1u);
        EXPECT_EQ(matrix[0], 0b110u);
    }

    {
        std::vector<std::pair<unsigned int, uint8_t>> inconsistent = {
            {0b11110u, 0u},
            {0b11110u, 0u},
            {0b00110u, 1u},
            {0b01011u, 1u},
            {0b10011u, 1u}};
        int rank = CountBoolMatrixRankAndStepTransform(inconsistent, 5);
        EXPECT_EQ(rank, -1);
        ASSERT_EQ(inconsistent.size(), 0u);
        EXPECT_TRUE(inconsistent.empty());
    }

    {
        std::vector<std::pair<unsigned int, uint8_t>> PairMatrix = {
            {0b10000u, 0u},
            {0b10011u, 1u},
            {0b11010u, 0u},
            {0b10011u, 1u}};
        std::vector<std::pair<unsigned int, uint8_t>> SStepMatrix = {
            {0b10000u, 0u},
            {0b01001u, 1u},
            {0b00011u, 1u}};
        int rank = CountBoolMatrixRankAndStepTransform(PairMatrix, 5);
        EXPECT_EQ(rank, 3);
        EXPECT_TRUE(PairMatrix == SStepMatrix);
    }
    {
        std::vector<std::pair<unsigned int, uint8_t>> Matrix1 = {
            {0b11110u, 0u},
            {0b00110u, 1u},
            {0b01011u, 1u},
            {0b10011u, 1u}};
        EXPECT_FALSE(IsMatrixSpecialStepTransform(Matrix1, 5));

        std::vector<std::pair<unsigned int, uint8_t>> Matrix2 = {
            {0b10000u, 0u},
            {0b01001u, 1u},
            {0b00011u, 1u}};
        EXPECT_TRUE(IsMatrixSpecialStepTransform(Matrix2, 5));

        std::vector<std::pair<unsigned int, uint8_t>> Matrix3 = {
            {0b10000u, 0u},
            {0b01011u, 1u},
            {0b00001u, 1u}};

        EXPECT_FALSE(IsMatrixSpecialStepTransform(Matrix3, 5));

        std::vector<std::pair<unsigned int, uint8_t>> Matrix4 = {
            {0b1100u, 0u},
            {0b0010u, 1u},
            {0b0100u, 1u}};
        EXPECT_FALSE(IsMatrixSpecialStepTransform(Matrix4, 4));
    }
}
TEST(BFUtilsTest, GetSolutionsOfSystem)
{
    {
        std::vector<std::pair<unsigned int, uint8_t>> PairMatrix = {
            {0b10000u, 0u},
            {0b11010u, 0u},
            {0b10011u, 1u}};
        auto rank = CountBoolMatrixRankAndStepTransform(PairMatrix, 5);
        EXPECT_EQ(rank, 3);
        std::vector<uint32_t> solutions = getSolutionsOfSystem(PairMatrix, 5);
        EXPECT_EQ(solutions.size(), 4);
        EXPECT_TRUE(std::find(solutions.begin(), solutions.end(), 0b01010u) != solutions.end());
        EXPECT_TRUE(std::find(solutions.begin(), solutions.end(), 0b00001u) != solutions.end());
        EXPECT_TRUE(std::find(solutions.begin(), solutions.end(), 0b01110u) != solutions.end());
        EXPECT_TRUE(std::find(solutions.begin(), solutions.end(), 0b00101u) != solutions.end());
    }

    {
        std::vector<std::pair<unsigned int, uint8_t>> PairMatrix = {
            {0b11110u, 0u},
            {0b00110u, 1u},
            {0b01011u, 1u},
            {0b10011u, 1u}};
        auto rank = CountBoolMatrixRankAndStepTransform(PairMatrix, 5);
        auto solutions = getSolutionsOfSystem(PairMatrix, 5);
        EXPECT_EQ(solutions.size(), 0);
    }

    {
        std::vector<std::pair<unsigned int, uint8_t>> PairMatrix = {
            {0b11100u, 0u},
            {0b00100u, 0u},
            {0b00001u, 1u}};
        auto rank = CountBoolMatrixRankAndStepTransform(PairMatrix, 5);
        EXPECT_EQ(rank, 3);
        auto solutions = getSolutionsOfSystem(PairMatrix, 5);
        EXPECT_EQ(solutions.size(), 4);

        auto areEqualByContent = [](const std::vector<uint32_t> &a, const std::vector<uint32_t> &b) -> bool
        {
            return std::unordered_multiset<uint32_t>(a.begin(), a.end()) ==
                   std::unordered_multiset<uint32_t>(b.begin(), b.end());
        };
        std::vector<uint32_t> expectedSolutions = {
            0b00001u,
            0b00011u,
            0b11001u,
            0b11011u};
        EXPECT_TRUE(areEqualByContent(solutions, expectedSolutions));
    }
}

TEST(BFCoreTest, ConstructFromStringAndWeight)
{
    BF f("0001");
    EXPECT_EQ(f.getNumberOfVariables(), 2);
    EXPECT_EQ(BF::weight(f), 1u);
    EXPECT_FALSE(BF::getValue(f, 0));
    EXPECT_TRUE(BF::getValue(f, 3));
}

TEST(BFCoreTest, SetAndGetValue)
{
    BF f(3, 0);
    EXPECT_FALSE(f[5]);
    f.setValue(5, true);
    EXPECT_TRUE(f[5]);
    f.setValue(5, false);
    EXPECT_FALSE(f[5]);
}

TEST(BFCoreTest, MobiusTransformIsSelfInverse)
{
    for (int i = 0; i <= 5; ++i)
    {
        BF f(10, 2);
        BF transformed = BF::mobius_transform(f);
        BF restored = BF::mobius_transform(transformed);
        EXPECT_EQ(restored, f);
    }
}

TEST(BFCoreTest, WalshTransformForZeroFunction)
{
    BF f(3, 0);
    std::vector<int> wht = BF::WH_transform(f);
    ASSERT_EQ(wht.size(), 8u);
    EXPECT_EQ(wht[0], 8);
    for (size_t i = 1; i < wht.size(); ++i)
    {
        EXPECT_EQ(wht[i], 0);
    }
}

static std::vector<int> WalshTransformByDefinition(const BF &f, int n)
{
    const unsigned int size = 1u << n;
    std::vector<int> result(size, 0);

    for (unsigned int a = 0; a < size; ++a)
    {
        int sum = 0;
        for (unsigned int x = 0; x < size; ++x)
        {
            const bool fx = BF::getValue(f, x);
            const bool ax = weight_mod(a & x);
            sum += (fx ^ ax) ? -1 : 1;
        }
        result[a] = sum;
    }
    return result;
}

TEST(BFCoreTest, WalshTransformMatchesDefinition)
{

    for (uint32_t i = 0; i <= 10; i++)
    {
        const int n = 8;
        BF f(n, 2);

        const std::vector<int> expected = WalshTransformByDefinition(f, n);
        const std::vector<int> actual = BF::WH_transform(f);

        ASSERT_EQ(actual.size(), expected.size());
        EXPECT_EQ(actual, expected);
    }
}

TEST(BFCoreTest, WalshTransformParsevalIdentity)
{
    for (uint32_t i = 0; i <= 10; i++)
    {
        const int n = 8;
        BF f(n, 2);
        const unsigned int size = 1u << n;

        const std::vector<int> wht = BF::WH_transform(f);
        long long energy = 0;
        for (int v : wht)
            energy += static_cast<long long>(v) * v;

        EXPECT_EQ(energy, static_cast<long long>(size) * size);
    }
}

TEST(BFCoreTest, AffineFunctionHasZeroNonlinearity)
{
    BF f("01010101");
    std::vector<int> wht = BF::WH_transform(f);
    EXPECT_EQ(BF::Nonlinearity(f, wht), 0u);
}

TEST(BFCoreTest, DegreeForMonomialX0X1)
{
    BF f("00010001");
    BF anf = BF::mobius_transform(f);
    EXPECT_EQ(BF::degree(anf), 2u);
}
