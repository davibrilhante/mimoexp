#include <iostream>                                                             
#include <ctime>                                                                
#include <string>                                                               
#include <random>                                                               
#include <cmath>                                                                
#include <complex>                                                              
#include <numeric>

using namespace std;

complex<double> complexMean(vector<complex<double>> input)
{
    //calculates complex mean
    //https://en.wikipedia.org/wiki/Complex_random_variable#Expectation
    complex<double> accum;
    for (int i=0; i<input.size(); i++)
    {
        accum += input[i];
    }

    double realmean = real(accum)/input.size();
    double imagmean = imag(accum)/input.size();

    complex<double> m(realmean, imagmean);
    return m;
}

double variance(vector<complex<double>> input)
{
    //Calculates complex variance
    //https://en.wikipedia.org/wiki/Complex_random_variable#Variance_and_pseudo-variance
    complex<double> mean = complexMean(input);
    double var=0;

    for (int i=0; i<input.size(); i++)
    {
        complex<double> factor = input[i] - mean;
        var += abs(factor)*abs(factor);
    }
    return var/input.size();
}

class Channel
{
    public:
        virtual void realization() = 0;

        complex<double>** sample;


    protected:
        int tx_channel_len;
        int rx_channel_len;
};

class MimoRayleighChannel : public Channel
{
    public:
        void realization()
        {
            sample = (complex<double>**) malloc(rx_channel_len*sizeof(complex<double>));  

            //normal distribution constructor
            std::default_random_engine generator;                                       
            std::normal_distribution<double> distribution(0.0, 1);

            //Channel sampling       
            for (int i=0; i<rx_channel_len; i++)
            {                                        
                sample[i] = (complex<double>*) malloc(tx_channel_len*sizeof(complex<double>));
                for (int j=0;j<tx_channel_len;j++)
                {
                    double real = distribution(generator);                       
                    double imag = distribution(generator);                       

                    complex<double> noise_sample(real, imag);                               

                    sample[i][j]=noise_sample;
                }
            }

            //Normalization
            vector <complex<double>> column;
            for (int j=0;j<tx_channel_len;j++)
            {
                for (int i=0;i<rx_channel_len;i++)
                {
                    column.push_back(sample[i][j]);
                }

                double var = variance(column);

                for (int i=0;i<rx_channel_len;i++)
                {
                    sample[i][j] = sqrt(1/var)*sample[i][j];
                }
            }
        }                                                                           
        

        void setDimension(int n_r, int n_t)
        {
            tx_channel_len = n_t;
            rx_channel_len = n_r;
        }
};

class dummyChannel : public Channel
{
    public:

        void realization()
        {
            sample = (complex<double>**) malloc(rx_channel_len*sizeof(complex<double>));

            for (int i=0; i<tx_channel_len; i++)
            {
                sample[i] = (complex<double>*) malloc(tx_channel_len*sizeof(complex<double>));
                for (int j=0; j<rx_channel_len; j++)
                {
                    if (i == j)
                    {
                        double real = gain;
                        double imag = gain;
                        complex<double> channel_component(real, imag);
                        sample[i][j] = channel_component;
                    }
                    else{
                        sample[i][j] = 0.0;
                    }
                }
            }
        }

        void setDimension(int n_r, int n_t)
            {
                tx_channel_len = n_t;
                rx_channel_len = n_r;
            }

        void setGain(double G)
        {
            gain = G;
        }

    protected:
        double gain;
};

int main(){                                                                     
    float gain = 0.5;                                                           
    int n_tx = 2;
    int n_rx = 2;

    MimoRayleighChannel channel;                                                            
    channel.setDimension(n_tx, n_rx);
    //channel.setGain(gain);

    channel.realization();                                                        

    for (int i=0; i<n_tx; i++){                                                
        for (int j=0; j<n_rx; j++){                                            
            cout << channel.sample[i][j] <<" ";                                  
        }                                                                       
        cout << endl;                                                           
    }                                                                           
    return 0;                                                                   
}
