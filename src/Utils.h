#pragma once

#include <cstdint>
#include <utility>
#include <vector>

int countBoolMatrixRankAndStepTransform(std::vector<unsigned int> &Matrix, int columns = 32);
int countBoolMatrixRankAndStepTransform(std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns = 32);
bool isMatrixSpecialStepTransform(const std::vector<std::pair<unsigned int, uint8_t>> &PairMatrix, int columns = 32);
std::vector<uint32_t> getSolutionsOfSystem(const std::vector<std::pair<uint32_t, uint8_t>> &PairMatrix, int n_);

std::vector<std::pair<uint32_t, uint32_t>> convertToPairs(const std::pair<std::vector<unsigned int>, std::vector<unsigned int>> &GPVec);
uint64_t convertToUint64(uint32_t num1, uint32_t num2);
std::pair<uint32_t, uint32_t> convertToPairUint32T(uint64_t num);

uint32_t countFirstZeros(uint32_t x);
int getElderBitPos(unsigned int num, int start = 31);
bool bitValue(unsigned int num, unsigned char bit);
void setBit(uint32_t &num, uint8_t bit, bool is_one);
void printMonom(unsigned int monom);

unsigned int uintWeight(unsigned int x);
bool weightMod(unsigned int x);

unsigned int nextCombination(unsigned int prev);
