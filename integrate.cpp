#include <pthread.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
// #include <chrono>
#include <random>

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

void *compute_integral(void *threadid)
{
    int id = (long)threadid;
    double thread_sum = 0.0;
    for (int i = id * thread_samples; i < (id + 1) * thread_samples; i++)
    {
        double x = dist(engine);
        if (x != 0)
        {
            thread_sum += sin(x) / x;
        }
    }
    pthread_mutex_lock(&mutex);
    total_sum += thread_sum;
    pthread_mutex_unlock(&mutex);
    pthread_exit(NULL);
}

int main(int argc, char *argv[])
{
    engine = std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());
    dist = std::uniform_real_distribution<double>(a, b);
    thread_samples = num_samples / num_threads;
    pthread_t threads[num_threads];
    pthread_mutex_init(&mutex, NULL);
    for (int i = 0; i < num_threads; i++)
    {
        int rc = pthread_create(&threads[i], NULL, compute_integral, (void *)i);
        if (rc)
        {
            std::cerr << "Error: Unable to create thread" << std::endl;
            exit(-1);
        }
    }
    for (int i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }
    pthread_mutex_destroy(&mutex);
    double result = (b - a) * total_sum / num_samples;
    std::cout << "Result: " << result << std::endl;
    pthread_exit(NULL);
    return 0;
}

struct parameters
{
    long thread_id;
    int N;                 // Number of samples
    int T;                 // Number of threads
    pthread_mutex_t mutex; // Mutex to protect the shared result variable
    // Gaussian quadrature weights and points
    vector<double> x;
    vector<double> w;
    vector<double> *result;
};

// Function to be integrated
double f(double x)
{
    return sin(x) / x;
}

double sum(vector<double> &values)
{
    double sum = 0;
    for (ulong i = 0; i < values.size(); i++)
    {
        sum += values[i];
    }
    return sum;
}

double ksum(const vector<double> &values)
{
    double sum = 0.0;
    double error = 0.0;
    for (double value : values)
    {
        double y = value - error;
        double temp = sum + y;
        error = (temp - sum) - y;
        sum = temp;
    }
    return sum;
}

double pairwiseSum(vector<double> &values)
{
    // Recursively add pairs of numbers together
    if (values.size() == 1)
    {
        return values[0];
    }
    if (values.size() == 2)
    {
        return values[0] + values[1];
    }
    std::vector<double> newValues;
    for (ulong i = 0; i < values.size() - 1; i += 2)
    {
        newValues.push_back(values[i] + values[i + 1]);
    }
    if (values.size() % 2 != 0)
    {
        newValues.push_back(values.back());
    }
    return pairwiseSum(newValues);
}

void *integrate(void *params)
{
    struct parameters *vars = (struct parameters *)params;
    // double thread_result = 0;

    // Calculate the integral for this thread
    for (long i = vars->thread_id; i < vars->N; i += vars->T)
    {
        // pthread_mutex_lock(&(vars->mutex));
        (*(vars->result))[i] = vars->w[i] * f(vars->x[i]);
        // cout << (*(vars->result))[i] << flush;
        // pthread_mutex_unlock(&(vars->mutex));
        // thread_result += vars->w[i] * f(vars->x[i]);
    }

    // // Add this thread's result to the global result
    // pthread_mutex_lock(&(vars->mutex));
    // *(vars->result) += thread_result;
    // pthread_mutex_unlock(&(vars->mutex));

    return 0;
}

int main(int argc, char *argv[])
{

    if (argc != 5)
    {
        cout << "Wrong number of arguments" << endl;
        cout << argc << endl;
        return 0;
    }

    auto start = std::chrono::high_resolution_clock::now();

    std::mt19937 engine;
    std::uniform_real_distribution<double> dist;

    // Set the integral limits
    double a = stod(argv[1]);
    double b = stod(argv[2]);
    int num_samples = stoi(argv[3]);
    int num_threads = stoi(argv[4]);

    vector<double> *result = new vector<double>(num_samples); // Result of the integral calculation
    vector<double> x(num_samples);
    vector<double> w(num_samples);

    // Initialize the mutex
    pthread_mutex_t mutex; // Mutex to protect the shared result variable
    pthread_mutex_init(&mutex, NULL);

    // Initialize the Gaussian quadrature weights and points
    // TODO: Is this the correct method? wiki looks more complicated
    for (int i = 0; i < num_samples; i++)
    {
        x[i] = a + (b - a) * (i + 0.5) / num_samples;
        w[i] = (b - a) / num_samples;
    }

    // Create the threads
    pthread_t threads[num_threads];
    parameters param_list[num_threads];
    for (long i = 0; i < num_threads; i++)
    {
        param_list[i] = {i, num_samples, num_threads, mutex, x, w, result};
        pthread_create(&threads[i], NULL, integrate, (void *)&(param_list[i]));
    }

    // Wait for all threads to finish
    for (int i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    // double final_result = sum(*result);
    double final_result = ksum(*result);
    // double final_result = pairwiseSum(*result);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Duration: " << duration.count() << " milliseconds" << std::endl;

    // Print the result
    // for (int i=0; i <10; i++){
    //     cout << (*result)[i] << endl;
    // }
    cout << "Result: " << fixed << setprecision(18) << final_result << endl;
    // cout << "Result: " << fixed << setprecision(18) << result << endl;

    pthread_mutex_destroy(&mutex);
    delete result;

    return 0;
}
