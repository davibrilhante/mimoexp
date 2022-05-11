#include <complex>
#include <vector>

#ifndef DEMODULATORS_H
#define DEMODULATORS_H

class Demodulator
{
    public:
	Demodulator();
        virtual std::string demodulate(std::complex<double>* input, int data_len) = 0; 

        void setOrder(int O);

    protected:
        int mod_order;
        int bitspersymbol;
};

class DemodulatorBPSK : public Demodulator
{
    public:
	DemodulatorBPSK();
	std::string demodulate(std::complex<double>* input, int data_len);
};

class DemodulatorQPSK : public Demodulator
{
    public:
	DemodulatorQPSK();
	std::string demodulate(std::complex<double>* input, int data_len);
};

#endif
