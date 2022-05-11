#include <complex>
#include <vector>

#ifndef MODULATORS_H
#define MODULATORS_H


std::vector<std::string> chunkBits(int bitspersymbol, std::string &input, char pad);

class Modulator
{

    public:
	Modulator();
        virtual std::complex<double>* modulate(std::string &input) = 0; 

        void setOrder(int O);

    protected:
        int mod_order;
        int bitspersymbol;
};

class ModulatorBPSK : public Modulator
{
    public:                                                                     
	ModulatorBPSK();
	std::complex<double>* modulate(std::string &input);
};

//class ModulatorQPSK : public Modulator{};
class ModulatorQPSK : public Modulator
{
    public:                                                                     
	ModulatorQPSK();
	std::complex<double>* modulate(std::string &input);
};

//class ModulatorMPSK : public Modulator{};
class ModulatorMPSK : public Modulator
{
    public:                                                                     
	ModulatorMPSK();
	std::complex<double>* modulate(std::string &input);
};

#endif
