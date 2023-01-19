/*
compile with:
    g++ -pthread -o integrate main.cpp
run with:
    ./integrate a b n_samples n_threads
output:
    result, time taken
*/

// TODO Take command-line input
// TODO Time calculation
// TODO Optimize for speed

#include <pthread.h>
#include <iostream>
#include <cmath>

using namespace std;

const int N = 1000; // Number of samples
const int T = 10;   // Number of threads

double result;     // Result of the integral calculation
pthread_mutex_t mutex; // Mutex to protect the shared result variable

// Function to be integrated
double f(double x) {
    return sin(x)/x;
}

// Gaussian quadrature weights and points
double x[N], w[N];

// Thread function
void* integrate(void* id) {
    long thread_id = (long) id;
    double thread_result = 0;

    // Calculate the integral for this thread
    for (long i = thread_id; i < N; i += T) {
        thread_result += w[i] * f(x[i]);
    }

    // Add this thread's result to the global result
    pthread_mutex_lock(&mutex);
    result += thread_result;
    pthread_mutex_unlock(&mutex);

    return 0;
}

int main(int num_args, char** args) {

    if (num_args != 5){
        cout << "Wrong number of arguments" << endl;
        return 0;
    }

    // Set the integral limits
    double a = stod(args[1]);
    double b = stod(args[2]);
    int num_samples = stoi(args[3]);
    int num_threads = stoi(args[4]);


    // Initialize the mutex
    pthread_mutex_init(&mutex, NULL);


    // Initialize the Gaussian quadrature weights and points
    // TODO: Is this the correct method? wiki looks more complicated
    for (int i = 0; i < N; i++) {
        x[i] = a + (b - a) * (i + 0.5) / N;
        w[i] = (b - a) / N;
    }

    // Create the threads
    pthread_t threads[T];
    for (long i = 0; i < T; i++) {
        pthread_create(&threads[i], NULL, integrate, (void*) i);
    }

    // Wait for all threads to finish
    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }

    // Print the result
    cout << "Result: " << result << endl;

    // Destroy the mutex
    pthread_mutex_destroy(&mutex);

    return 0;
}
