#include <pthread.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>
#include <immintrin.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
// #include <thread>

using namespace std;

/*
compile with:
    g++ -mavx -pthread -O3 -Wall -o integrate main.cpp
run with:
    ./integrate a b n_samples n_threads
output:
    result


integral_0^1 sin(x)/x dx | 18 digits 
= 0.946083070367183015
*/


struct parameters{
    long thread_id;
    ulong* N; // Number of samples
    int T; // Number of threads
    double* global_res;
    double a;
    double b;
    bool* running;
};


bool running;
pthread_mutex_t mutex;
pthread_mutexattr_t myMutexAttr;

// Function to be integrated
double f(double x) {
    // return 1;
    return sin(x)/x;
}


void* integrate(void* params) {
    struct parameters* vars = (struct parameters*) params;
    // uint seed = (uint) vars->thread_id + 123;
    // __m256d w = _mm256_set1_pd(vars->w);
    __m256d range = _mm256_set1_pd((vars->b - vars->a) / (double)(*vars->N));
    __m256d h = _mm256_set_pd(3.5, 2.5, 1.5, 0.5);
    __m256d a = _mm256_set1_pd(vars->a);
    __m256d sum = _mm256_setzero_pd();
    __m256d error = _mm256_setzero_pd();
    for (ulong i = 4*(vars->thread_id); i < *vars->N; i += 4*(vars->T)) {
        __m256d xi = _mm256_set1_pd(i);
        __m256d xih = _mm256_add_pd(xi, h);
        __m256d xihr = _mm256_mul_pd(xih, range);
        __m256d xihra = _mm256_add_pd(xihr, a);
        __m256d fi = _mm256_set_pd(f(xihra[3]), f(xihra[2]), f(xihra[1]), f(xihra[0]));
        __m256d weighted = _mm256_mul_pd(range, fi);
        __m256d y = _mm256_sub_pd(weighted, error);
        __m256d temp = _mm256_add_pd(sum, y);
        error = _mm256_sub_pd(_mm256_sub_pd(temp, sum), y);
        sum = temp;
    }
    __m256d sum2 = _mm256_hadd_pd(sum, sum);
    __m128d sum3 = _mm_add_pd(_mm256_extractf128_pd(sum2, 1), _mm256_castpd256_pd128(sum2));
    double thread_result = _mm_cvtsd_f64(sum3);
    pthread_mutex_lock(&(mutex));
    *(vars->global_res) += thread_result;
    pthread_mutex_unlock(&(mutex));

    return 0;
}


// void* integrate(void* params) {
//     struct parameters* vars = (struct parameters*) params;

//     double thread_result = 0.0;
//     double error = 0.0;
//     double scale = (vars->b - vars->a)/(*(vars->N));
//     for (ulong i = vars->thread_id; i < (*(vars->N)); i += vars->T) {
//         double y = vars->w * f(scale * (i + 0.5) + vars->a) - error;
//         double temp = thread_result + y;
//         error = (temp - thread_result) - y;
//         thread_result = temp;
//     }
//     // Add this thread's result to the global result
//     pthread_mutex_lock(&(mutex));
//     *(vars->global_res) += thread_result;
//     pthread_mutex_unlock(&(mutex));

//     return 0;
// }


void* timed_integrate(void* params) {
    struct parameters* vars = (struct parameters*) params;
    uint seed = (uint) vars->thread_id + 123;
    // Calculate the integral for this thread
    double thread_result = 0.0;
    double error = 0.0;
    double scale = (vars->b - vars->a)/((double)RAND_MAX);
    ulong thread_samples = 0;
    while (running) {
        thread_samples+=1;
        double y = f(scale * (double)rand_r(&seed) + vars->a) - error;
        double temp = thread_result + y;
        error = (temp - thread_result) - y;
        thread_result = temp;
    }
    pthread_mutex_lock(&(mutex));
    *(vars->N) += thread_samples;
    *(vars->global_res) += thread_result;
    pthread_mutex_unlock(&(mutex));

    return 0;
}


int main(int num_args, char** args) {

    auto start = std::chrono::high_resolution_clock::now();

    if (num_args != 5){
        cout << "Wrong number of arguments" << endl;
        cout << num_args << endl;
        return 0;
    }

    string time_string = args[3];
    bool timed = (time_string[time_string.length() - 1] == 's');
    ulong total_samples = 0;

    // Set the integral limits
    double a = stod(args[1]);
    double b = stod(args[2]);
    ulong num_samples;
    if (timed) {
        num_samples = stoi(time_string.substr(0, time_string.length() - 1));
    } else {
        num_samples = stoi(args[3]);
    }
    int num_threads = stoi(args[4]);


    // Initialize the mutex
    pthread_mutexattr_init(&myMutexAttr);
    pthread_mutexattr_setpshared(&myMutexAttr, PTHREAD_PROCESS_SHARED);
    pthread_mutex_init(&mutex, &myMutexAttr);


    // Create the threads
    pthread_t threads[num_threads];
    parameters param_list[num_threads];
    double global_res = 0.0;
    bool run = true;
    running = true;


    if (!timed){
        for (long i = 0; i < num_threads; i++) {
            param_list[i] = {i, &num_samples, num_threads, &global_res, a, b, &run};
            pthread_create(&threads[i], NULL, integrate, (void*) &(param_list[i]));
        }
    }
    else {
        //  timed
        for (int i = 0; i < num_threads; i++) {
            param_list[i] = {i, &total_samples, num_threads, &global_res, a, b, &run};
            pthread_create(&threads[i], NULL, timed_integrate, (void*) &(param_list[i]));
        }
        sleep(num_samples);
        // std::this_thread::sleep_for (std::chrono::seconds(num_samples));
        running = false;
    }
    

    // Wait for all threads to finish
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    if (timed){
        global_res *= (b-a)/total_samples;
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    // cout << "Threads: " << num_threads << ", "
    //     << "Duration: " << duration.count() << " milliseconds, "
    //     << "Result: " << fixed << setprecision(18) << global_res << endl;

    cout << num_threads << ","
        << duration.count() << ","
        << fixed << setprecision(18) << global_res << endl;

    // cout << fixed << setprecision(18) << global_res << endl;

    pthread_mutexattr_destroy(&myMutexAttr);
    pthread_mutex_destroy(&mutex);

    return 0;
}
