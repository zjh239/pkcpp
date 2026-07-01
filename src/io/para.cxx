// para.cxx
module;

#include <iostream>
#include <vector>
#include <thread>

export module para;
export {

// Function that computes the sum of a subarray
void partial_sum(const std::vector<int>& arr, size_t start, size_t end, long long& result) {
    long long sum = 0;
    for (size_t i = start; i < end; ++i) {
        sum += arr[i];
    }
    result = sum;
}

int test() {
    const size_t N = 1'000'000;
    std::vector<int> data(N, 1);  // Initialize with 1s

    const size_t num_threads = std::thread::hardware_concurrency();
    std::cout << "Number of threads supported by CPU: " << num_threads << std::endl;

    // We will use 4 threads
    // const size_t num_threads = 4;
    std::vector<std::thread> threads(num_threads);
    std::vector<long long> partial_results(num_threads, 0);

    size_t chunk_size = N / num_threads;

    // Launch threads
    for (size_t t = 0; t < num_threads; ++t) {
        // std::cout<<"This is thread "<<t<<std::endl;
        size_t start = t * chunk_size;
        size_t end = (t == num_threads - 1) ? N : start + chunk_size;
        threads[t] = std::thread(partial_sum, std::cref(data), start, end, std::ref(partial_results[t]));

    }

    // Join threads
    for (size_t t = 0; t < num_threads; ++t) {
        threads[t].join();
        // std::cout<<"Partial sum result is "<<partial_results[t]<<std::endl;
    }

    // Combine partial results
    long long total_sum = 0;
    for (const auto& sum : partial_results) {
        total_sum += sum;
    }

    std::cout << "Total sum = " << total_sum << std::endl;

    return 0;
}
}
