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

vector<string> serialToParallel(string &input, int bitspersymbol, int n_outputs)
{
    int nbatches = input.length()/n_outputs;
    vector<string> output;
    const char pad = '0';
    string temp = "";

    for (int i=0; i<nbatches; i++)
    {
        for (int j=0; j<n_outputs; j++)
        {
            for (int k=0; k<bitspersymbol; k++)
            {
                //temp.append(input[i*n_outputs*bitspersymbol+j*bitspersymbol+k]);
                temp += input[i*n_outputs*bitspersymbol+j*bitspersymbol+k];
            }
        }
        if (temp.length() == n_outputs)
        {
            output.push_back(temp);
            temp = "";
        }
    }

    if (temp != "")
    {
        while (temp.length() < n_outputs*bitspersymbol)
        {
            temp += pad;
        }
        output.push_back(temp);
    }
    
    return output;
}

/*
string ParallelToSerial(vector<string> input, int n_inputs)
{
    int a = 1;
}*/

class Precoder
{
    public:
        virtual void precode()=0;
};

class Detector
{
    public:
        virtual void detect()=0;
};

class ZeroForcingDetector : public Detector
{
    public:
        void detect()
        {
            int a = 1;
        }

};



class MimoEncoder
{
    public:
        virtual complex<double>** encode(string &input)=0;
        
        virtual void setModulator(Modulator* m)=0;

        void setOutputs(int n)
        {
            n_outputs = n;
        }

    protected:
        Modulator * modulator;
        int n_outputs;

};


class VBlastEncoder : public MimoEncoder
{
    public:
        complex<double>** encode(string &input)
        {
            int bitspersymbol = modulator->getBitspersymbol();
            vector<string> parallelized = serialToParallel(input, bitspersymbol, n_outputs);
            complex<double>** output = (complex<double>**) malloc(parallelized.size()*n_outputs*sizeof(complex<double>));

            for (int i=0; i<parallelized.size(); i++)
            {
                /*
                output[i] = (complex<double>*) malloc(n_outputs*sizeof(complex<double>));

                for (int j=0; j<n_outputs; j++)
                {
                    string chunk = parallelized[i].substr(j*bitspersymbol, (j+1)*bitspersymbol-1);
                    output[i][j] = modulator->modulate(chunk);
                }*/
                output[i] = modulator->modulate(parallelized[i]);
            }
            return output;
        }

        void setModulator(Modulator* m)
        {
            modulator = m;
        }

        int getnOutputs()
        {
            return n_outputs;
        }
};

class MimoDecoder
{
    public:
        virtual string decode(complex<double>*input_symbol, int n_tx)=0;
        virtual void setDemodulator(Demodulator* d)=0;
        virtual void setDetector(Detector* d)=0;

        void setInputs(int n)
        {
            n_inputs = n;
        }

    protected:
        Demodulator *demodulator;
        Detector *detector;
        vector<complex<double>> buffer;
        int n_inputs;
};

class VBlastDecoder : public MimoDecoder
{
        public:                                                                     
            string decode(complex<double>*input_symbol, int n_tx)
            {
                if (buffer.size() < n_inputs)
                {
                    for (int i=0; i<n_tx; i++)
                    {
                        buffer.push_back(input_symbol[i]);
                    }

                    return ""; 
                }
                else
                {
                    string output;
                    for (int i=0; i<n_inputs; i++)
                    {
                        output += demodulator->demodulate(buffer[i]);
                    }
                    buffer.empty();
                    return output;
                }
                
            }

            void setDemodulator(Demodulator* d)
            {
                demodulator = d;
            }

            void setDetector(Detector* d)
            {
                detector = d;
            }

            int getnInputs()
            {
                return n_inputs;
            }

};
int main()
{
    string data = "0001101100011011"; 

    Modulator * modulator = new ModulatorBPSK();
    Demodulator * demodulator = new DemodulatorBPSK();

    VBlastEncoder encoder;
    encoder.setModulator(modulator);

    VBlastDecoder decoder;
    decoder.setDemodulator(demodulator);
    
    complex<double>** data_encoded = encoder.encode(data);

    int n_chunks = data.length()/(decoder.getnInputs()*modulator->getBitspersymbol());

    for (int i=0; i<n_chunks; i++)
    {
        string data_decoded = decoder.decode(data_encoded[i], encoder.getnOutputs());
    }


    delete modulator;

    return 0;
}
