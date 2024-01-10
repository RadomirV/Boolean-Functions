#include "BF.cpp"
#include <chrono>

int main()
{
    srand(time(NULL));
    const int numberOfVariables = 5;
    int printInRow = 6;
    BF func = BF::GenBalancedFunc(numberOfVariables); // сгенерировать уравновешенную функцию с заданным числом переменных
    // BF func("10000010011101111100101001001101");
    while (true)
    {

        std::cout << '\n';
        func.print();
        std::cout << '\n';

        auto wht = BF::WH_transform(func);
        /*
       for (int i = 0; i < wht.size(); i++)
            std::cout << wht[i] << " ";
        std::cout << '\n';
        */

        std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
        BF::FillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht); // заполняем множесва W соответствующими элементами

        auto resPairs = func.PairsToImprove();              // Множество пар улучшающих полученых решением СЛУ
        auto straightPairs = func.PairsToImproveStraight(); // Множество пар улучшающих полученых прямым алгоритмом
        int transfer = 0;
        if (resPairs.empty())
        {
            std::cout << "NO PAIRS\n";
        }
        for (auto &i : resPairs)
        {
            transfer++;
            std::cout << "(" << std::bitset<numberOfVariables>(i.first) << "[" << i.first << "], " << std::bitset<numberOfVariables>(i.second) << "[" << i.second << "]) ";
            if (transfer % printInRow == 0)
                std::cout << std::endl;
        }
        auto testFunc = func.TestPairsToImproveFunctions(resPairs, straightPairs); // проверить сов
        if (!testFunc)
        {
            std::cout << "ERROR IN FUNCTION TESTING\n";
            return 0;
        }

        std::cout << "\nW3_pos_size= " << W_3_pos.size() << " \nW3_neg_size= " << W_3_neg.size() << " \nNf= " << BF::Nonlinearity(func, wht) << std::endl;
        uint32_t x1 = 0, x2 = 0;

        std::cout << "Enter pairs to change:\n x1= ";
        std::cin >> x1;
        std::cout << " x2= ";
        std::cin >> x2;
        try
        {
            func = BF::SwapOnSets(func, x1, x2); // обменять значения на наборах x1, x2
        }
        catch (const char *err)
        {
            std::cout << err << std::endl;
            return 0;
        }
    }

    return 0;
}

/*
int main()
{
    srand(time(NULL));
    int numOfVariables = 6;
    int W_set_SUM = 2;
    std::cout << "Enter number of variables:";
    std::cin >> numOfVariables;
    std::cout << "Enter sum of W sets:";
    std::cin >> W_set_SUM;
    std::cout << std::endl;

    // func.print();
    std::vector<int64_t> statistic_Straight, statistic_Opt;
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

        auto start_Straight = std::chrono::steady_clock::now();
        auto good_pairs = BF::good_pairsVec(func, wht_coef);
        auto end_Straight = std::chrono::steady_clock::now();
        if (good_pairs.first.empty() || good_pairs.second.empty())
            continue;

        // std::cout << "number of pairsStraight= " << good_pairs.first.size() * good_pairs.second.size() << "\n";

        auto start_Opt = std::chrono::steady_clock::now();
        auto opt_good_pairs = BF::OptGoodPairsVec(func, wht_coef);
        auto end_Opt = std::chrono::steady_clock::now();

        auto end = std::chrono::steady_clock::now();

        auto StraightTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_Straight - start_Straight).count();
        auto OptTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_Opt - start_Opt).count();
        auto TimeForIteration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        auto TimeForGenWHT = std::chrono::duration_cast<std::chrono::nanoseconds>(end_gen - start_gen).count();

        std::cout << "Brute= " << StraightTime << " Opt= " << OptTime << " Power of W sets= " << W_1_pos.size() + W_1_neg.size();
        // std::cout << " pairsStraight= " << good_pairs.first.size() * good_pairs.second.size();
        std::cout << " pairs= " << opt_good_pairs.first.size() * opt_good_pairs.second.size();
        std::cout << " IterationTime= " << TimeForIteration << " time wht= " << TimeForGenWHT << " Attemps to gen W= " << attempsToGen << '\n';

        statistic_Straight.push_back(StraightTime);
        statistic_Opt.push_back(OptTime);
        attempsToGen = 0;
        i++;
    }

    return 0;
}
*/
