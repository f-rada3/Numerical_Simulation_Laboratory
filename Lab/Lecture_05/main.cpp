#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

// funzioni
double error(const vector<double> &sum_prog, const vector<double> &sum_prog2, int n) {
    if (n == 0) return 0;
    return sqrt((sum_prog2[n] - pow(sum_prog[n], 2)) / n);
}

double wave_function_1s(double x, double y, double z) {
   return 1 / sqrt(M_PI) * exp(-sqrt(x * x + y * y + z * z));
}

double wave_function_2p(double x, double y, double z) {
   return 1 / 8. * sqrt(2 / M_PI) * exp(-sqrt(x * x + y * y + z * z) / 2.) * z;
}

double acceptance_ratio(double prob_old, double prob_new) {
   return min(1., prob_new / prob_old);
}

int main() {
    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) Primes >> p1 >> p2;
    else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    // parameters
    int M = 1000000, N = 100, L = M / N;
    // step[0] = half the step for 1S with uniform distribution
    // step[1] = standard deviation for 1S with Gaussian distribution
    // step[2] = half the step for 2P with uniform distribution
    // step[3] = standard deviation for 2P with Gaussian distribution
    vector<double> step = {1.22, 0.7, 2.8, 1.6};
    
    // Posizioni iniziali fissate
    double x_1s_near = 0.1, y_1s_near = 0.1, z_1s_near = 0.1;
    double x_1s_far = 10.0, y_1s_far = 10.0, z_1s_far = 10.0;
    double x_2p_near = 0.1, y_2p_near = 0.1, z_2p_near = 0.1;
    double x_2p_far = 10.0, y_2p_far = 10.0, z_2p_far = 10.0;

    cout << "\n=================================================================================\n";
    cout << "STARTING CALCULATIONS FOR ALL DISTRIBUTIONS AND POSITIONS\n";
    cout << "=================================================================================\n\n";

    // Loop sulle due distribuzioni
    for (int dist_choice = 1; dist_choice <= 2; dist_choice++) {
        string suffix = (dist_choice == 1) ? "_uniform" : "_gauss";
        
        cout << "Processing distribution: " << (dist_choice == 1 ? "UNIFORM" : "GAUSSIAN") << endl;

        double x_1s = x_1s_near, y_1s = y_1s_near, z_1s = z_1s_near;
        double x_1s_far_curr = x_1s_far, y_1s_far_curr = y_1s_far, z_1s_far_curr = z_1s_far;
        double x_2p = x_2p_near, y_2p = y_2p_near, z_2p = z_2p_near;
        double x_2p_far_curr = x_2p_far, y_2p_far_curr = y_2p_far, z_2p_far_curr = z_2p_far;
        
        int count_accept_1s = 0, count_accept_1s_far = 0, count_accept_2p = 0, count_accept_2p_far = 0;

        // Output files
        ofstream pos_output_1s("positions_1s" + suffix + ".data");
        ofstream pos_output_2p("positions_2p" + suffix + ".data");
        ofstream pos_output_1s_far("positions_1s_far" + suffix + ".data");
        ofstream pos_output_2p_far("positions_2p_far" + suffix + ".data");
        ofstream avg_output_1s("output_average_1s" + suffix + ".data");
        ofstream avg_output_2p("output_average_2p" + suffix + ".data");
        ofstream avg_output_1s_far("output_average_1s_far" + suffix + ".data");
        ofstream avg_output_2p_far("output_average_2p_far" + suffix + ".data");
        
        vector<double> ave1_1s, ave2_1s, sum_prog_1s, sum_prog2_1s, error_prog_1s;
        vector<double> ave1_1s_far, ave2_1s_far, sum_prog_1s_far, sum_prog2_1s_far, error_prog_1s_far;
        vector<double> ave1_2p, ave2_2p, sum_prog_2p, sum_prog2_2p, error_prog_2p;
        vector<double> ave1_2p_far, ave2_2p_far, sum_prog_2p_far, sum_prog2_2p_far, error_prog_2p_far;

       
       // 1S Orbital
        for (int i = 0; i < N; i++) {
            double sum_1s = 0;
            double sum_1s_far = 0;
            
            for (int j = 0; j < L; j++) {

                double x_new = 0.0, y_new = 0.0, z_new = 0.0, x_new_far = 0.0, y_new_far = 0.0, z_new_far = 0.0;
                
                if (dist_choice == 1) {
                    x_new = x_1s + rnd.Rannyu(-step[0], step[0]);
                    y_new = y_1s + rnd.Rannyu(-step[0], step[0]);
                    z_new = z_1s + rnd.Rannyu(-step[0], step[0]);
                
                    x_new_far = x_1s_far_curr + rnd.Rannyu(-step[0], step[0]);
                    y_new_far = y_1s_far_curr + rnd.Rannyu(-step[0], step[0]);
                    z_new_far = z_1s_far_curr + rnd.Rannyu(-step[0], step[0]);
                }
                else {
                    x_new = x_1s + rnd.Gauss(0, step[1]);
                    y_new = y_1s + rnd.Gauss(0, step[1]);
                    z_new = z_1s + rnd.Gauss(0, step[1]);
                
                    x_new_far = x_1s_far_curr + rnd.Gauss(0, step[1]);
                    y_new_far = y_1s_far_curr + rnd.Gauss(0, step[1]);
                    z_new_far = z_1s_far_curr + rnd.Gauss(0, step[1]);
                }
                
                double acceptance = acceptance_ratio(pow(wave_function_1s(x_1s, y_1s, z_1s), 2), pow(wave_function_1s(x_new, y_new, z_new), 2));
                double acceptance_far = acceptance_ratio(pow(wave_function_1s(x_1s_far_curr, y_1s_far_curr, z_1s_far_curr), 2), pow(wave_function_1s(x_new_far, y_new_far, z_new_far), 2));

                if (rnd.Rannyu() <= acceptance) {
                   x_1s = x_new; y_1s = y_new; z_1s = z_new;
                   count_accept_1s++;
                }
                if (rnd.Rannyu() <= acceptance_far) {
                    x_1s_far_curr = x_new_far; y_1s_far_curr = y_new_far; z_1s_far_curr = z_new_far;
                    count_accept_1s_far++;
                }
                pos_output_1s << x_1s << "\t" << y_1s << "\t" << z_1s << endl;
                pos_output_1s_far << x_1s_far_curr << "\t" << y_1s_far_curr << "\t" << z_1s_far_curr << endl;

                sum_1s += sqrt(x_1s * x_1s + y_1s * y_1s + z_1s * z_1s);
                sum_1s_far += sqrt(x_1s_far_curr * x_1s_far_curr + y_1s_far_curr * y_1s_far_curr + z_1s_far_curr * z_1s_far_curr);
            }
            ave1_1s.push_back(sum_1s / L);
            ave2_1s.push_back(ave1_1s[i] * ave1_1s[i]);
            ave1_1s_far.push_back(sum_1s_far / L);
            ave2_1s_far.push_back(ave1_1s_far[i] * ave1_1s_far[i]);
        }

        // Data blocking for 1S
        for (int i = 0; i < N; i++) {
            sum_prog_1s.push_back(0);
            sum_prog2_1s.push_back(0);
            sum_prog_1s_far.push_back(0);
            sum_prog2_1s_far.push_back(0);
            for (int j = 0; j < i + 1; j++) {
                sum_prog_1s[i] += ave1_1s[j];
                sum_prog2_1s[i] += ave2_1s[j];
                sum_prog_1s_far[i] += ave1_1s_far[j];
                sum_prog2_1s_far[i] += ave2_1s_far[j];
            }
            sum_prog_1s[i] /= (i + 1);
            sum_prog2_1s[i] /= (i + 1);
            error_prog_1s.push_back(error(sum_prog_1s, sum_prog2_1s, i));
            sum_prog_1s_far[i] /= (i + 1);
            sum_prog2_1s_far[i] /= (i + 1);
            error_prog_1s_far.push_back(error(sum_prog_1s_far, sum_prog2_1s_far, i));
            avg_output_1s << (i + 1) * L << "\t" << sum_prog_1s[i] << "\t" << error_prog_1s[i] << "\t" << endl;
            avg_output_1s_far << (i + 1) * L << "\t" << sum_prog_1s_far[i] << "\t" << error_prog_1s_far[i] << "\t" << endl;
        }

    
       // 2P Orbital
        for (int i = 0; i < N; i++) {
            double sum_2p = 0;
            double sum_2p_far = 0;

            for (int j = 0; j < L; j++) {

                double x_new_2p = 0.0, y_new_2p = 0.0, z_new_2p = 0.0, x_new_2p_far = 0.0, y_new_2p_far = 0.0, z_new_2p_far = 0.0;

                if (dist_choice == 1) {
                    x_new_2p = x_2p + rnd.Rannyu(-step[2], step[2]);
                    y_new_2p = y_2p + rnd.Rannyu(-step[2], step[2]);
                    z_new_2p = z_2p + rnd.Rannyu(-step[2], step[2]);

                    x_new_2p_far = x_2p_far_curr + rnd.Rannyu(-step[2], step[2]);
                    y_new_2p_far = y_2p_far_curr + rnd.Rannyu(-step[2], step[2]);
                    z_new_2p_far = z_2p_far_curr + rnd.Rannyu(-step[2], step[2]);
                }
                else {
                    x_new_2p = x_2p + rnd.Gauss(0, step[3]);
                    y_new_2p = y_2p + rnd.Gauss(0, step[3]);
                    z_new_2p = z_2p + rnd.Gauss(0, step[3]);

                    x_new_2p_far = x_2p_far_curr + rnd.Gauss(0, step[3]);
                    y_new_2p_far = y_2p_far_curr + rnd.Gauss(0, step[3]);
                    z_new_2p_far = z_2p_far_curr + rnd.Gauss(0, step[3]);
                }

                double acceptance_2p = acceptance_ratio(pow(wave_function_2p(x_2p, y_2p, z_2p), 2), pow(wave_function_2p(x_new_2p, y_new_2p, z_new_2p), 2));
                double acceptance_2p_far = acceptance_ratio(pow(wave_function_2p(x_2p_far_curr, y_2p_far_curr, z_2p_far_curr), 2), pow(wave_function_2p(x_new_2p_far, y_new_2p_far, z_new_2p_far), 2));

                if (rnd.Rannyu() <= acceptance_2p) {
                   x_2p = x_new_2p; y_2p = y_new_2p; z_2p = z_new_2p;
                   count_accept_2p++;
                }
                if (rnd.Rannyu() <= acceptance_2p_far) {
                    x_2p_far_curr = x_new_2p_far; y_2p_far_curr = y_new_2p_far; z_2p_far_curr = z_new_2p_far;
                    count_accept_2p_far++;
                }

                pos_output_2p << x_2p << "\t" << y_2p << "\t" << z_2p << endl;
                pos_output_2p_far << x_2p_far_curr << "\t" << y_2p_far_curr << "\t" << z_2p_far_curr << endl;

                sum_2p += sqrt(x_2p * x_2p + y_2p * y_2p + z_2p * z_2p);
                sum_2p_far += sqrt(x_2p_far_curr * x_2p_far_curr + y_2p_far_curr * y_2p_far_curr + z_2p_far_curr * z_2p_far_curr);
            }

            ave1_2p.push_back(sum_2p / L);
            ave2_2p.push_back(ave1_2p[i] * ave1_2p[i]);
            ave1_2p_far.push_back(sum_2p_far / L);
            ave2_2p_far.push_back(ave1_2p_far[i] * ave1_2p_far[i]);
        }

        // Data blocking for 2P
        for (int i = 0; i < N; i++) {
            sum_prog_2p.push_back(0);
            sum_prog2_2p.push_back(0);
            sum_prog_2p_far.push_back(0);
            sum_prog2_2p_far.push_back(0);
            for (int j = 0; j < i + 1; j++) {
                sum_prog_2p[i] += ave1_2p[j];
                sum_prog2_2p[i] += ave2_2p[j];
                sum_prog_2p_far[i] += ave1_2p_far[j];
                sum_prog2_2p_far[i] += ave2_2p_far[j];
            }
            sum_prog_2p[i] /= (i + 1);
            sum_prog2_2p[i] /= (i + 1);
            error_prog_2p.push_back(error(sum_prog_2p, sum_prog2_2p, i));
            sum_prog_2p_far[i] /= (i + 1);
            sum_prog2_2p_far[i] /= (i + 1);
            error_prog_2p_far.push_back(error(sum_prog_2p_far, sum_prog2_2p_far, i));
            avg_output_2p << (i + 1) * L << "\t" << sum_prog_2p[i] << "\t" << error_prog_2p[i] << "\t" << endl;
            avg_output_2p_far << (i + 1) * L << "\t" << sum_prog_2p_far[i] << "\t" << error_prog_2p_far[i] << "\t" << endl;
        }

        int last_block = N - 1;
        cout << "---Distribution: " << (dist_choice == 1 ? "UNIFORM" : "GAUSSIAN") << "---\n";
        
        cout << "  Orbital 1s (starting from (0.1, 0.1, 0.1)):\n";
        cout << "   <r> = " << sum_prog_1s[last_block] << " ± " << error_prog_1s[last_block] << " a₀" << endl;
        cout << "  Acceptance rate 1s: " << (double)count_accept_1s / M << endl;

        cout << "\n  Orbital 1s (starting from (10, 10, 10)):\n";
        cout << "   <r> = " << sum_prog_1s_far[last_block] << " ± " << error_prog_1s_far[last_block] << " a₀" << endl;
        cout << "  Acceptance rate 1s: " << (double)count_accept_1s_far / M << endl;

        cout << "\n  Orbital 2p (starting from (0.1, 0.1, 0.1)):\n";
        cout << "   <r> = " << sum_prog_2p[last_block] << " ± " << error_prog_2p[last_block] << " a₀" << endl;
        cout << "  Acceptance rate 2p: " << (double)count_accept_2p / M << endl;

        cout << "\n  Orbital 2p (starting from (10, 10, 10)):\n";
        cout << "   <r> = " << sum_prog_2p_far[last_block] << " ± " << error_prog_2p_far[last_block] << " a₀" << endl;
        cout << "  Acceptance rate 2p: " << (double)count_accept_2p_far / M << "\n" << endl;

        // Chiusura file
        pos_output_1s.close();
        pos_output_2p.close();
        pos_output_1s_far.close();
        pos_output_2p_far.close();
        avg_output_1s.close();
        avg_output_2p.close();
        avg_output_1s_far.close();
        avg_output_2p_far.close();
    }

    cout << "\n";
    cout << "All calculations completed successfully.\n";
    cout << "Generated output files:\n";
    cout << "  - positions_1s_uniform.data, positions_1s_gauss.data\n";
    cout << "  - positions_1s_far_uniform.data, positions_1s_far_gauss.data\n";
    cout << "  - positions_2p_uniform.data, positions_2p_gauss.data\n";
    cout << "  - positions_2p_far_uniform.data, positions_2p_far_gauss.data\n";
    cout << "  - output_average_1s_uniform.data, output_average_1s_gauss.data\n";
    cout << "  - output_average_1s_far_uniform.data, output_average_1s_far_gauss.data\n";
    cout << "  - output_average_2p_uniform.data, output_average_2p_gauss.data\n";
    cout << "  - output_average_2p_far_uniform.data, output_average_2p_far_gauss.data\n";
    cout << "\n";

    rnd.SaveSeed();
    return 0;
}