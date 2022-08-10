#include "modulators.h"
#include "demodulators.h"

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

Demodulator::Demodulator(){};

void Demodulator::setOrder(int O)
{
    mod_order = O;
    bitspersymbol = log2(mod_order);

}

DemodulatorBPSK::DemodulatorBPSK()
{
    mod_order = 2;
    bitspersymbol = 1;
};


string DemodulatorBPSK::demodulate(complex<double>* input, int data_len=1)
{
    string output;

    for (int i=0; i<data_len; i++)
    {
        int data = max(0.0, real(input[i]));
        output.append(to_string(data));
    }

    return output;
}


DemodulatorQPSK::DemodulatorQPSK()
{
    mod_order = 4;
    bitspersymbol = 2;
};


string DemodulatorQPSK::demodulate(complex<double>* input, int data_len=1)
{
    string output;

    for (int i=0; i<data_len; i++)
    {
        //STUB!!!
        int data = max(0.0, real(input[i]));
        output.append(to_string(data));
    }

    return output;
}

/*
int main()
{
    string data = "0001101100011011";

    ModulatorBPSK modulator;
    DemodulatorBPSK demodulator;

    complex<double> * signal = modulator.modulate(data);

    string output = demodulator.demodulate(signal, data.length());

    cout<<data<<endl;
    cout<<output<<endl;
    
    return 0;
}*/
