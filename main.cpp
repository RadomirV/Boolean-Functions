#include "BF.cpp"
#include <chrono>
#include <thread>
typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
typedef int int32_t;
// typedef unsigned long uint64_t;

int main()
{
    srand(time(NULL));
    const int numberOfVariables = 8;
    BF func(numberOfVariables);
    uint64_t lindiff;
    uint64_t neutralAtt;
    uint64_t improveAtt;
    std::map<uint32_t, std::map<uint32_t, uint32_t>> stat;
    std::pair<uint32_t, uint32_t> exitStat;
    std::pair<double, double> exitNf(0, 0);
    int attemps = 10000;
    unsigned int Nf = 0;
    auto wht_coef = BF::WH_transform(func);
    // for n = 8 optimal     80<=Nf<=108
    // for n = 9 optimal    182 <=Nf<= 202
    // for n = 10 optimal    384 <=Nf<=412
    for (int i = 0; i < attemps; ++i)
    {
        func = BF::GenBalancedFunc(numberOfVariables);
        wht_coef = BF::WH_transform(func);
        Nf = BF::Nonlinearity(func, wht_coef);
        while ((Nf < 80) || (Nf > 108))
        {
            func = BF::GenBalancedFunc(numberOfVariables);
            wht_coef = BF::WH_transform(func);
            Nf = BF::Nonlinearity(func, wht_coef);
        }
        //std::cout << "NF = " << Nf << std::endl;
        lindiff = 1000;
        neutralAtt = 16 * numberOfVariables;
        improveAtt = 4 * numberOfVariables;

        func.nonlinearityImprove(lindiff, neutralAtt, improveAtt);
        wht_coef = BF::WH_transform(func);
        unsigned int Nf_new = BF::Nonlinearity(func, wht_coef);

        stat[Nf][Nf_new]++;
        if (neutralAtt == 0)
        {
            exitStat.first++;
            exitNf.first += Nf_new;
        }
        if (improveAtt == 0)
        {
            exitStat.second++;
            exitNf.second += Nf_new;
        }
        // std::this_thread::sleep_for(std::chrono::milliseconds(100));
        std::cout << "I= " << i << " ";
    }
    int transfer = 0;
    int printInRow = 12;
    std::cout << "\nexit of neutral = " << exitStat.first << " exit of improve = " << exitStat.second << std::endl;
    exitNf.first = exitNf.first / exitStat.first;
    exitNf.second = exitNf.second / exitStat.second;
    std::cout << "exit NF of neutral = " << exitNf.first << " exit Nf of improve = " << exitNf.second << std::endl;

    for (auto it = stat.begin(); it != stat.end(); ++it)
    {
        transfer = 0;
        std::cout << " Original NF= " << it->first << "|||";
        for (auto IN_iter = it->second.begin(); IN_iter != it->second.end(); ++IN_iter)
        {
            transfer++;
            std::cout << " " << IN_iter->first << "[" << IN_iter->second << "]";
            if (transfer % printInRow == 0)
                std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}

int main_nnn()
{
    srand(time(NULL));
    const int numberOfVariables = 10;
    BF func(numberOfVariables);
    int iter = 1 << numberOfVariables;
    std::map<uint32_t, uint32_t> stat;
    uint64_t lindiff;
    uint64_t neutralAtt;
    uint64_t improveAtt;
    std::pair<uint32_t, uint32_t> exitStat;
    std::pair<double, double> exitNf(0, 0);
    int isSecondHalf = 0;
    int neuAtt = 8 * numberOfVariables;
    int impAtt = 5;
    for (int i = 1; i < iter; ++i)
    {

        lindiff = 1000;
        neutralAtt = neuAtt;
        improveAtt = impAtt;

        func = func.generateAffine(i, false);
        func.nonlinearityImprove(lindiff, neutralAtt, improveAtt);
        auto wht_coef = BF::WH_transform(func);
        unsigned int Nf = BF::Nonlinearity(func, wht_coef);
        stat[Nf]++;
        if (neutralAtt == 0)
        {
            exitStat.first++;
            exitNf.first += Nf;
        }
        if (improveAtt == 0)
        {
            exitStat.second++;
            exitNf.second += Nf;
        }

        lindiff = 1000;
        neutralAtt = neuAtt;
        improveAtt = impAtt;

        func = func.generateAffine(i, true);
        func.nonlinearityImprove(lindiff, neutralAtt, improveAtt);
        wht_coef = BF::WH_transform(func);
        Nf = BF::Nonlinearity(func, wht_coef);
        stat[Nf]++;
        if (neutralAtt == 0)
        {
            exitStat.first++;
            exitNf.first += Nf;
        }
        if (improveAtt == 0)
        {
            exitStat.second++;
            exitNf.second += Nf;
        }

        // if (i % 100 == 0)
        std::cout << "\t\t\t\titer = " << i << "\n";
    }
    exitNf.first = exitNf.first / exitStat.first;
    exitNf.second = exitNf.second / exitStat.second;

    std::cout << "exit of neutral = " << exitStat.first << " exit of improve = " << exitStat.second << std::endl;
    std::cout << "exit NF of neutral = " << exitNf.first << " exit Nf of improve = " << exitNf.second << std::endl;
    for (auto it = stat.begin(); it != stat.end(); ++it)
    {
        std::cout << "[Nf = " << it->first << "] = " << it->second << "\n";
    }

    return 0;
}

int tryt()
{
    srand(time(NULL));
    const int numberOfVariables = 26;
    BF func = BF::GenBalancedFunc(numberOfVariables);
    // func.print();
    // std::cout << "\n";
    // func.print_per_bit();
    // std::cout << "\n";
    auto wht = BF::WH_transform(func);
    // for (auto &i : wht)
    // {
    //     std::cout << i << " ";
    // }
    // std::cout << std::endl;

    std::vector<uint32_t> W_1_pos, W_1_neg, W_3_pos, W_3_neg;
    BF::FillWSets(W_1_pos, W_1_neg, W_3_pos, W_3_neg, wht);
    std::cout << "Nf= " << BF::Nonlinearity(func, wht) << std::endl;
    // func.nonlinearityImprove(10000, 100000, 10000);
    //  func.print();
    wht = BF::WH_transform(func);
    std::cout << "Nf= " << BF::Nonlinearity(func, wht) << std::endl;

    return 0;
}

int AllBalanced()
{
    srand(time(NULL));
    const int numberOfVariables = 5;
    unsigned int GoodPairsCount = 0;
    unsigned int BadPairsCount = 0;
    int iterationCount = 0;
    auto limits = BF::generateBorderBalancesFunctions();
    // auto limits = std::make_pair(BF("0000000011111111"), BF("1111111100000000"));
    // auto limits = std::make_pair(BF("00001111"), BF("11110000"));
    //    auto limits = std::make_pair(BF("0011"), BF("1100"));
    BF func = limits.first;
    std::vector<unsigned int> good_freq;
    std::vector<unsigned int> bad_freq;
    good_freq.resize(300);
    bad_freq.resize(300);
    while (func != limits.second)
    {
        auto wht_coef = BF::WH_transform(func);

        auto badPairs = func.PairsToWorsen();      // Pairs which decrease the nonlinearity
        auto improvePairs = func.PairsToImprove(); // Pairs which improve nonlinearity
        GoodPairsCount += improvePairs.size();
        BadPairsCount += badPairs.size();
        good_freq[improvePairs.size()]++;
        bad_freq[badPairs.size()]++;
        /*
        func.print();
        std::cout << "\t\t\t\t\t\tIterN = " << iterationCount << "\nGood = " << improvePairs.size() << " Bad = " << badPairs.size() << " Nonlinearity = " << BF::Nonlinearity(func, wht_coef) << "\n";
        for (auto &i : wht_coef)
        {
            std::cout << i << " ";
        }
        std::cout << "FFFFFFFFFFFFFFFF\n";
        */
        /*
        if (improvePairs.size() == 0 && badPairs.size() == 0)
        {
            std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
        }
        for (auto &i : improvePairs)
        {
            std::cout << "(" << std::bitset<numberOfVariables>(i.first) << "[" << i.first << "], " << std::bitset<numberOfVariables>(i.second) << "[" << i.second << "]) ";
        }
        std::cout<<"------------------------------------------------------------------------\n";
        for (auto &i : badPairs)
        {
            std::cout << "(" << std::bitset<numberOfVariables>(i.first) << "[" << i.first << "], " << std::bitset<numberOfVariables>(i.second) << "[" << i.second << "]) ";
        }
        */

        if (iterationCount % 100000 == 0)
        {
            uint32_t transfer = 8;
            uint32_t inRow = 0;
            func.print_per_bit();
            std::cout << "\t\t\t\t\t\tIterN = " << iterationCount << "\nGood = " << improvePairs.size() << " Bad = " << badPairs.size() << " Nonlinearity = " << BF::Nonlinearity(func, wht_coef) << "\n";
            for (int i = 0; i < 300; ++i)
            {
                if (good_freq[i] != 0)
                {
                    std::cout << "[" << i << "]= " << good_freq[i] << " ";
                    inRow++;
                    if (inRow % transfer == 0)
                        std::cout << "\n";
                }
            }
            inRow = 0;
            std::cout << "\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" << std::endl;
            for (int i = 0; i < 300; ++i)
            {
                if (bad_freq[i] != 0)
                {
                    std::cout << "[" << i << "]= " << bad_freq[i] << " ";
                    inRow++;
                    if (inRow % transfer == 0)
                        std::cout << "\n";
                }
            }
            std::cout << std::endl;
        }
        if (improvePairs.size() == 0 && badPairs.size() == 0)
        {
            func.print_per_bit();
            std::cout << "NO PAIRS AT ALL GOOD AND BAD NO PAIRS AT ALL GOOD AND BAD NO PAIRS AT ALL GOOD AND BAD\n";
            std::cout << "\t\t\t\t\t\tIterN = " << iterationCount << "\nGood = " << improvePairs.size() << " Bad = " << badPairs.size() << " Nonlinearity = " << BF::Nonlinearity(func, wht_coef) << "\n";
        }
        if (BF::Nonlinearity(func, wht_coef) > 12)
        {
            std::cout << "Good = " << improvePairs.size() << " Bad = " << badPairs.size() << " Nonlinearity = " << BF::Nonlinearity(func, wht_coef) << "\t\t\t\t\t\tIterN = " << iterationCount << "\n";
        }

        func.nextBalanced();
        iterationCount++;
    }
    auto wht_coef = BF::WH_transform(func);

    auto badPairs = func.PairsToWorsen();      // Pairs which decrease the nonlinearity
    auto improvePairs = func.PairsToImprove(); // Pairs which improve nonlinearity

    GoodPairsCount += improvePairs.size();
    BadPairsCount += badPairs.size();
    func.print();
    func.print_per_bit();
    std::cout << "\t\t\t\t\t\tIterN = " << iterationCount << "\nGood = " << improvePairs.size() << " Bad = " << badPairs.size() << "\n";

    std::cout << "GoodPairsCount= " << GoodPairsCount << "\n";
    std::cout << "BadPairsCount= " << BadPairsCount << "\n";
    uint32_t transfer = 8;
    uint32_t inRow = 0;
    good_freq[improvePairs.size()]++;
    bad_freq[badPairs.size()]++;
    std::cout << "Distribution of the powers of sets of improving pairs" << std::endl;
    for (int i = 0; i < 300; ++i)
    {
        if (good_freq[i] != 0)
        {
            std::cout << "[" << i << "]= " << good_freq[i] << " ";
            inRow++;
            if (inRow % transfer == 0)
                std::cout << "\n";
        }
    }
    inRow = 0;
    std::cout << "distribution of the powers of sets of worsening pairs" << std::endl;
    for (int i = 0; i < 300; ++i)
    {
        if (bad_freq[i] != 0)
        {
            std::cout << "[" << i << "]= " << bad_freq[i] << " ";
            inRow++;
            if (inRow % transfer == 0)
                std::cout << "\n";
        }
    }

    return 0;
}
int mainmain()
{
    srand(time(NULL));
    const int numberOfVariables = 5;
    int printInRow = 6;
    // BF func = BF::GenBalancedFunc(numberOfVariables); // сгенерировать уравновешенную функцию с заданным числом переменных
    std::string str;
    // std::cin >> str;
    // BF func(str);
    // BF func = BF::GenBalancedFunc(numberOfVariables);
    BF func("1010");
    while (true)
    {

        std::cout << '\n';
        func.print();
        std::cout << '\n';
        auto wht = BF::WH_transform(func);
        for (auto &i : wht)
        {
            std::cout << i << " ";
        }
        std::cout << std::endl;

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

        auto badPairs = func.PairsToWorsen();
        std::cout << "\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";

        for (auto &i : badPairs)
        {
            transfer++;
            std::cout << "(" << std::bitset<numberOfVariables>(i.first) << "[" << i.first << "], " << std::bitset<numberOfVariables>(i.second) << "[" << i.second << "]) ";
            if (transfer % printInRow == 0)
                std::cout << std::endl;
        }

        auto testFunc = func.TestPairsToImproveFunctions(resPairs, straightPairs); // проверить совпадение
        if (!testFunc)
        {
            std::cout << "ERROR IN FUNCTION TESTING\n";
            return 0;
        }

        std::cout << "\nW3_pos_size= " << W_3_pos.size() << " \nW3_neg_size= " << W_3_neg.size() << " \nNf= " << BF::Nonlinearity(func, wht) << std::endl;
        for (auto &i : wht)
        {
            std::cout << i << " ";
        }
        std::cout << std::endl;

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
