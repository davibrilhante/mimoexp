#include "modulators.h"

#include <stdlib.h>
#include <iostream>                                                             
#include <ctime>                                                                
#include <string>                                                               
#include <random>                                                               
#include <cmath>                                                                
#include <complex>                                                              


#define PI 3.14159265

using namespace std;
using std::string;
//string::size_type sz;

vector<string> chunkBits(int bitspersymbol, string &input, char pad)
{
    int n_chunks = input.length()/bitspersymbol;
    vector<string> output;
    string res = "";

    for (int i=0; i<n_chunks; i++)
    {
        for (int j=0; j<bitspersymbol; j++)
        {
            res += input[i*bitspersymbol+j];
            if (res.length() == bitspersymbol)
            {
                output.push_back(res);
                res = "";
            }
        }
    }

    if (res != "")
    {
        while (res.length() < bitspersymbol)
        {
            res += pad;
        }
        output.push_back(res);
    }
    
    return output;
}

Modulator::Modulator(){};

void Modulator::setOrder(int O)
{
    mod_order = O;
    bitspersymbol = sqrt(mod_order);

}

ModulatorBPSK::ModulatorBPSK(){
    mod_order = 2;
    bitspersymbol = 1;
};

//class ModulatorBPSK : public Modulator
complex<double>* ModulatorBPSK::modulate(string &input)
{
    const char* c_bit = input.c_str();

    int data_len = input.length();
    complex<double>* output = (complex<double>*) malloc(data_len*sizeof(complex<double>));

    for (int i=0; i<data_len; i++)
    {
        int i_bit = c_bit[i] - '0'; 
        output[i] = cos(PI*(1 - i_bit));
    }

    return output;
};

//class ModulatorQPSK : public Modulator
ModulatorQPSK::ModulatorQPSK()
{
    mod_order = 4;
    bitspersymbol = 2;
}

complex<double>* ModulatorQPSK::modulate(string &input)
{
    vector<string> saneInput = chunkBits(2, input, '0');

    int output_len = saneInput.size();

    complex<double>* output = (complex<double>*) malloc(output_len*sizeof(complex<double>));

    for (int i=0; i<output_len; i++)
    {
        int i_bit = stoi(saneInput[i], nullptr, 2); 
        double I_sample = cos((PI/4)*(2*i_bit + 1));
        double Q_sample = sin((PI/4)*(2*i_bit + 1));
        complex<double> sample(I_sample, Q_sample);
        output[i] = sample;
    }

    return output;
}



//class ModulatorMPSK : public Modulator
ModulatorMPSK::ModulatorMPSK(){};

complex<double>* ModulatorMPSK::modulate(string &input)
{
    vector<string> saneInput = chunkBits(bitspersymbol, input, '0');

    int output_len = saneInput.size();

    complex<double>* output = (complex<double>*) malloc(output_len*sizeof(complex<double>));

    for (int i=0; i<output_len; i++)
    {
        int i_bit = stoi(saneInput[i], nullptr, 2); 
        double I_sample = cos((2*PI/mod_order)*(i_bit));
        double Q_sample = sin((2*PI/mod_order)*(i_bit));
        complex<double> sample(I_sample, Q_sample);
        output[i] = sample;
    }

    return output;
}


/*
int main()
{
    string data = "0001101100011011";
    int mod_order = 4;

    ModulatorQPSK modulator;
    modulator.setOrder(mod_order);

    complex<double> * output = modulator.modulate(data);

    int tx_symbols = data.length()/sqrt(mod_order);
    
    for (int i=0; i < tx_symbols; i++)
    {
        cout<<output[i]<<endl;
    }

    return 0;
}*/
