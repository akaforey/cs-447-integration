#include <pthread.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>
#include <immintrin.h>
#include <stdlib.h>

using namespace std;

/*
compile with:
    g++ -pthread -o integrate main.cpp
run with:
    ./integrate a b n_samples n_threads
output:
    result, time taken


integral_0^1 sin(x)/x dx | 18 digits 
= 0.946083070367183015
*/

// TODO Optimize for speed


struct parameters{
    long thread_id;
    int N; // Number of samples
    int T; // Number of threads
    pthread_mutex_t mutex; // Mutex to protect the shared result variable
    // Gaussian quadrature weights and points
    double x;
    double w;
    double* global_res;
    double a;
    double b;
};


// Function to be integrated
double f(double x) {
    return sin(x)/x;
}

// double sum(vector<double>& values){
//     double sum = 0;
//     for (ulong i=0; i < values.size(); i++){
//         sum += values[i];
//     }
//     return sum;
// }



// double pairwiseSum(vector<double>& values) {
//     // Recursively add pairs of numbers together
//     if (values.size() == 1) {
//         return values[0];
//     }
//     if (values.size() == 2) {
//         return values[0] + values[1];
//     }
//     vector<double> newValues;
//     for (ulong i = 0; i < values.size() - 1; i += 2) {
//         newValues.push_back(values[i] + values[i + 1]);
//     }
//     if (values.size() % 2 != 0) {
//         newValues.push_back(values.back());
//     }
//     return pairwiseSum(newValues);
// }

// double ksum(const vector<double>& values) {
//     double sum = 0.0;
//     double error = 0.0;
//     for (double value : values) {
//         double y = value - error;
//         double temp = sum + y;
//         error = (temp - sum) - y;
//         sum = temp;
//     }
//     return sum;
// }

double ksum(const vector<double>& values) {
    __m256d sum = _mm256_setzero_pd();
    __m256d error = _mm256_setzero_pd();
    for (ulong i = 0; i < values.size(); i += 4) {
        __m256d value = _mm256_loadu_pd(&values[i]);
        __m256d y = _mm256_sub_pd(value, error);
        __m256d temp = _mm256_add_pd(sum, y);
        error = _mm256_sub_pd(_mm256_sub_pd(temp, sum), y);
        sum = temp;
    }
    __m256d sum2 = _mm256_hadd_pd(sum, sum);
    __m128d sum3 = _mm_add_pd(_mm256_extractf128_pd(sum2, 1), _mm256_castpd256_pd128(sum2));
    return _mm_cvtsd_f64(sum3);
}

void* integrate(void* params) {
    struct parameters* vars = (struct parameters*) params;
    // srand((uint) vars->thread_id + 123);
    // uint seed = (uint) vars->thread_id + 123;
    // cout << ((double)rand_r(&seed)/RAND_MAX)*(vars->b - vars->a) + vars->a << endl;
    __m256d w = _mm256_set1_pd(vars->w);
    __m256d x = _mm256_set1_pd(vars->x);
    __m256d h = _mm256_set1_pd( 0.5);
    __m256d sum = _mm256_setzero_pd();
    for (int i = 4*(vars->thread_id); i < vars->N; i += 4*(vars->T)) {
        __m256d xi = _mm256_set_pd(i + 3, i + 2, i + 1, i);
        __m256d xih = _mm256_add_pd(xi, h);
        __m256d xihx = _mm256_mul_pd(xih, x);
        __m256d fi = _mm256_set_pd(f(xihx[3]), f(xihx[2]), f(xihx[1]), f(xihx[0]));
        __m256d temp = _mm256_mul_pd(w, fi);
        sum = _mm256_add_pd(sum, temp);
    }
    __m256d sum2 = _mm256_hadd_pd(sum, sum);
    __m128d sum3 = _mm_add_pd(_mm256_extractf128_pd(sum2, 1), _mm256_castpd256_pd128(sum2));
    double thread_result = _mm_cvtsd_f64(sum3);
    pthread_mutex_lock(&(vars->mutex));
    *(vars->global_res) += thread_result;
    pthread_mutex_unlock(&(vars->mutex));

    return 0;
}

// void* integrate(void* params) {
//     struct parameters* vars = (struct parameters*) params;
//     vector<double> res(vars->N/vars->T+1);
//     // Calculate the integral for this thread
//     for (long i = vars->thread_id; i < vars->N; i += vars->T) {
//         res.push_back(vars->w * f((vars->x) * (i + 0.5)));
//     }
//     // Add this thread's result to the global result
//     double thread_result = ksum(res);
//     pthread_mutex_lock(&(vars->mutex));
//     *(vars->global_res) += thread_result;
//     pthread_mutex_unlock(&(vars->mutex));

//     return 0;
// }

int main(int num_args, char** args) {

    auto start = std::chrono::high_resolution_clock::now();

    if (num_args != 5){
        cout << "Wrong number of arguments" << endl;
        cout << num_args << endl;
        return 0;
    }


    // Set the integral limits
    double a = stod(args[1]);
    double b = stod(args[2]);
    int num_samples = stoi(args[3]);
    int num_threads = stoi(args[4]);

    // vector<double>* result = new vector<double>(num_samples);    // Result of the integral calculation
    double x(num_samples);
    double w;

    // Initialize the mutex
    pthread_mutex_t mutex; // Mutex to protect the shared result variable
    pthread_mutex_init(&mutex, NULL);


    // Initialize the Gaussian quadrature weights and points
    // TODO: Is this the correct method? wiki looks more complicated
    w = (b - a) / num_samples;
    x = a + (b - a) / num_samples;
    // for (int i = 0; i < num_samples; i++) {
    //     x[i] = a + (b - a) * (i + 0.5) / num_samples;
    // }

    // Create the threads
    pthread_t threads[num_threads];
    parameters param_list[num_threads];
    double global_res = 0;
    for (long i = 0; i < num_threads; i++) {
        param_list[i] = {i, num_samples, num_threads, mutex, x, w, &global_res, a, b};
        pthread_create(&threads[i], NULL, integrate, (void*) &(param_list[i]));
    }

    // Wait for all threads to finish
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // double final_result = sum(*result);
    // double final_result = ksum(*result);
    // double final_result = pairwiseSum(*result);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    cout << "Threads: " << num_threads << ", "
        << "Duration: " << duration.count() << " milliseconds, "
        << "Result: " << fixed << setprecision(18) << global_res << endl;
    // cout << global_res << endl;

    pthread_mutex_destroy(&mutex);

    return 0;
}
