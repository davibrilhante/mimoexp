
#include <complex>
#include <vector>

#ifndef DEMODULATORS_H
#define DEMODULATORS_H

std::vector<std::string> serialToParallel(std::string &input);

std::string ParallelToSerial(std::vector<std::string> input);

class Precoder
{
    public:
	Precoder();
        virtual void precode()=0;
};

class Detector
{
    public:
	Detector();
        virtual void detect()=0;
};


class MimoEncoder
{
    public:
	MimoEncoder();
        virtual std::complex<double>** encode(std::string &input)=0;
};

class MimoDecoder
{
    public:
	MimoDecoder();
        virtual std::tring* function(std::complex<double>* input symbol)=0;
};

class VBlastEncoder : public MimoEncoder
{
    public:
	VBlastEncoder();
	std::complex<double>** encode(std::string &input);
};

class VBlastDecoder : public MimoDecoder
{
        public:
	VBlastDecoder();
	std::string* decode(std::complex<double>* input_symbol);
};

class ZeroForcingDetector : public Detector
{
    public:
	ZeroForcingDetector();
        void detect();
        

};

#endif
