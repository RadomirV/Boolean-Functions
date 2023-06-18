#include <iostream>
#include "BF.cpp"
#include <chrono>
#include <set>

int main()
{
    srand(time(NULL));
    int numOfVariables = 21;
    int W_set_SUM = 3;
    std::cout << "Enter number of variables:";
    std::cin >> numOfVariables;
    std::cout << "Enter sum of W sets:";
    std::cin >> W_set_SUM;
    std::cout << std::endl;

    // func.print();
    std::vector<int64_t> statistic_stupid, statistic_Opt;
    std::cout << "Time in nanoseconds\n";
    uint64_t attempsToGen = 0;

    for (int i = 0; i < 1000;)
    {
        auto start = std::chrono::steady_clock::now();

        auto start_gen = std::chrono::steady_clock::now();
        BF func = BF::GenBalancedFunc(numOfVariables);
        auto wht_coef = BF::WH_transform(func);
        auto end_gen = std::chrono::steady_clock::now();
        std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
        BF::FillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht_coef);
        if (W_1_pos.size() + W_1_neg.size() < W_set_SUM)
        {
            attempsToGen++;
            continue;
        }

        auto start_Stupid = std::chrono::steady_clock::now();
        auto good_pairs = BF::good_pairsVec(func, wht_coef);
        auto end_Stupid = std::chrono::steady_clock::now();
        if (good_pairs.first.empty() || good_pairs.second.empty())
            continue;

        // std::cout << "number of pairsSTUPID= " << good_pairs.first.size() * good_pairs.second.size() << "\n";

        auto start_Opt = std::chrono::steady_clock::now();
        auto opt_good_pairs = BF::OptGoodPairsVec(func, wht_coef);
        auto end_Opt = std::chrono::steady_clock::now();

        auto end = std::chrono::steady_clock::now();

        auto StupidTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_Stupid - start_Stupid).count();
        auto OptTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_Opt - start_Opt).count();
        auto TimeForIteration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        auto TimeForGenWHT = std::chrono::duration_cast<std::chrono::nanoseconds>(end_gen - start_gen).count();

        std::cout << "Brute= " << StupidTime << " Opt= " << OptTime << " Power of W sets= " << W_1_pos.size() + W_1_neg.size();
        // std::cout << " pairsSTUPID= " << good_pairs.first.size() * good_pairs.second.size();
        std::cout << " pairs= " << opt_good_pairs.first.size() * opt_good_pairs.second.size();
        std::cout << " IterationTime= " << TimeForIteration << " time wht= " << TimeForGenWHT << " Attemps to gen W= " << attempsToGen << '\n';

        statistic_stupid.push_back(StupidTime);
        statistic_Opt.push_back(OptTime);
        attempsToGen = 0;
        i++;
    }

    return 0;
}
