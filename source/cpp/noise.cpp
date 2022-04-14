#include <iostream>
#include <ctime>
#include <string>
#include <random>
#include <cmath>
#include <complex>

using namespace std;


class Noise
{
    public:
        virtual void noiseSample() = 0;

        void setSigma(double s)
        {
            sigma = s;
        }

        void setLength(double n)
        {
            length = n;
        }

        complex<double>** samples;


    protected:
        double sigma;
        int length;
};



//Generates a zero-mean standard deviation sigma random random gaussian noise
//It is non-zero only in the main diagonal. Assumes IID noise components
class awgnNoise : public Noise
{
    public:
        //complex<double>** noiseSample(int n_samples)
        void noiseSample() 
        {
            //complex<double>** noise=0;                                                  
            samples = (complex<double>**) malloc(length*sizeof(complex<double>));  
            double sample_var = sigma/2;                                                

            std::default_random_engine generator;                                       
            std::normal_distribution<double> distribution(0.0, sigma);

            for (int i=0; i<length; i++){                                        
                samples[i] = (complex<double>*) malloc(length*sizeof(complex<double>));

                double real = sample_var*distribution(generator);                       
                double imag = sample_var*distribution(generator);                       
                complex<double> noise_sample(real, imag);                               

                samples[i][i] = noise_sample;

                for (int j=0;j<length;j++){                                      
                    if (i!=j){                                                          
                        samples[i][j]=0.0;                                                
                    }                                                                   
                }                                                                       
            }                                                                           
            //return noise;
        }

};







int main(){
	int n_ant = 4;
	float sigma = 1.0;

    awgnNoise noise;
    noise.setSigma(sigma);
    noise.setLength(n_ant);

    noise.noiseSample();

	for (int i=0; i<n_ant; i++){
		for (int j=0; j<n_ant; j++){
			cout << noise.samples[i][j] <<" ";
		}
		cout << endl;
	}
	return 0;
}
