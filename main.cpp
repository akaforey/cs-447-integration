#include <pthread.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>

using namespace std;

/*
compile with:
    g++ -pthread -o integrate main.cpp
run with:
    ./integrate a b n_samples n_threads
output:
    result, time taken
*/

// TODO Time calculation
// TODO Optimize for speed


struct parameters{
    long thread_id;
    int N; // Number of samples
    int T; // Number of threads
    pthread_mutex_t mutex; // Mutex to protect the shared result variable
    // Gaussian quadrature weights and points
    vector<double> x;
    vector<double> w;
    double* result;
};


// Function to be integrated
double f(double x) {
    return sin(x)/x;
}

// double average(const std::vector<double>& values) {
//     double sum = 0.0;
//     double error = 0.0;
//     for (double value : values) {
//         double y = value - error;
//         double temp = sum + y;
//         error = (temp - sum) - y;
//         sum = temp;
//     }
//     return sum / values.size();
// }

// double pairwiseSum(const std::vector<double>& values) {
//     // Recursively add pairs of numbers together
//     if (values.size() == 1) {
//         return values[0];
//     }
//     if (values.size() == 2) {
//         return values[0] + values[1];
//     }
//     std::vector<double> newValues;
//     for (int i = 0; i < values.size() - 1; i += 2) {
//         newValues.push_back(values[i] + values[i + 1]);
//     }
//     if (values.size() % 2 != 0) {
//         newValues.push_back(values.back());
//     }
//     return pairwiseSum(newValues);
// }

void* integrate(void* params) {
    struct parameters* vars = (struct parameters*) params;
    double thread_result = 0;

    // Calculate the integral for this thread
    for (long i = vars->thread_id; i < vars->N; i += vars->T) {
        thread_result += vars->w[i] * f(vars->x[i]);
        // vars->result[i] = vars->w[i] * f(vars->x[i]);
    }

    // // Add this thread's result to the global result
    pthread_mutex_lock(&(vars->mutex));
    *(vars->result) += thread_result;
    pthread_mutex_unlock(&(vars->mutex));

    return 0;
}

int main(int num_args, char** args) {


    if (num_args != 5){
        cout << "Wrong number of arguments" << endl;
        cout << num_args << endl;
        return 0;
    }

    auto start = std::chrono::high_resolution_clock::now();

    // Set the integral limits
    double a = stod(args[1]);
    double b = stod(args[2]);
    int num_samples = stoi(args[3]);
    int num_threads = stoi(args[4]);

    // vector<double> result;
    // result.resize(num_samples);
    double result;     // Result of the integral calculation
    // double x[num_samples];
    vector<double> x;
    x.resize(num_samples);
    // double w[num_samples];
    vector<double> w;
    w.resize(num_samples);

    // Initialize the mutex
    pthread_mutex_t mutex; // Mutex to protect the shared result variable
    pthread_mutex_init(&mutex, NULL);


    // Initialize the Gaussian quadrature weights and points
    // TODO: Is this the correct method? wiki looks more complicated
    for (int i = 0; i < num_samples; i++) {
        x[i] = a + (b - a) * (i + 0.5) / num_samples;
        w[i] = (b - a) / num_samples;
    }

    // Create the threads
    pthread_t threads[num_threads];
    parameters param_list[num_threads];
    for (long i = 0; i < num_threads; i++) {
        param_list[i] = {i, num_samples, num_threads, mutex, x, w, &result};
        pthread_create(&threads[i], NULL, integrate, (void*) &(param_list[i]));
    }

    // Wait for all threads to finish
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // double final_result = pairwiseSum(result)/num_samples;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Duration: " << duration.count() << " milliseconds" << std::endl;

    // Print the result
    cout << "Result: " << fixed << setprecision(18) << result << endl;

    pthread_mutex_destroy(&mutex);

    return 0;
}
