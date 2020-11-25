#include <iostream>
#include <random>
#include <chrono>
#include <memory>
const int N = 1 << 28;
int main()
{
    const int seed = 0;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dis;
    std::normal_distribution<double> normal;
    std::unique_ptr<bool[]> array = std::make_unique<bool[]>(N);

    for (int i = 0; i < N; i++)
    {
        if (dis(gen) > 0.5)
            array[i] = true;
        else
            array[i] = false;
    }

    int sum = 0;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++)
    {
        if (array[i])
            sum++;
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " microsecond" << std::endl;
//    std::cout << sum << std::endl;




         sum = 0;
     t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++)
    {
            sum+=array[i];
    }
     t2 = std::chrono::high_resolution_clock::now();

    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " microsecond" << std::endl;
//    std::cout << sum << std::endl;
        std::cout << N << std::endl;
}
