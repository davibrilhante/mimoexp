#include <iostream>                                                             
#include <ctime>                                                                
#include <string>                                                               
#include <random>                                                               
#include <cmath>                                                                
#include <complex>                                                              

using namespace std;

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
            1;
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
            sample = (complex<double>**) malloc(tx_channel_len*sizeof(complex<double>));

            for (int i=0; i<tx_channel_len; i++)
            {
                sample[i] = (complex<double>*) malloc(rx_channel_len*sizeof(complex<double>));
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

    dummyChannel channel;                                                            
    channel.setDimension(n_tx, n_rx);
    channel.setGain(gain);

    channel.realization();                                                        

    for (int i=0; i<n_tx; i++){                                                
        for (int j=0; j<n_rx; j++){                                            
            cout << channel.sample[i][j] <<" ";                                  
        }                                                                       
        cout << endl;                                                           
    }                                                                           
    return 0;                                                                   
}
