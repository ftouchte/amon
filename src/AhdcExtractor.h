#ifndef AHDC_EXTRACTOR_H
#define AHDC_EXTRACTOR_H

#include <vector>
#include <map>
#include <string>

/**
 * @class AhdcExtractor
 *
 * @author Felix Touchte Codjo
 */
class AhdcExtractor {
	public :
		float samplingTime; ///< time between two ADC bins
		int sparseSample = 0; ///< used to defined binOffset
		short adcOffset = 0; ///< pedestal or noise level
		long timeStamp = 0; ///< timeStamp timing informations (used to make fine corrections)
		float fineTimeStampResolution = 0; ///< precision of dream clock (usually 8)
		static const short ADC_LIMIT = 4095; ///< Maximum value of ADC : 2^12-1
		float amplitudeFractionCFA; ///< amplitude fraction between 0 and 1
		int binDelayCFD; ///< CFD delay parameter
		float fractionCFD; ///< CFD fraction parameter between 0 and 1

	public :
		int binMax; ///< Bin of the max ADC over the pulse
                int binOffset; ///< Offset due to sparse sample
                float adcMax; ///< Max value of ADC over the pulse (fitted)
                float timeMax; ///< Time of the max ADC over the pulse (fitted)
                float integral; ///< Sum of ADCs over the pulse (not fitted)

		std::vector<short> samplesCorr; ///< Waveform after offset (pedestal) correction
                int binNumber; ///< Number of bins in one waveform

                float leadingEdgeTime; ///< moment when the signal reaches a Constant Fraction of its Amplitude uphill (fitted)
                float trailingEdgeTime; ///< moment when the signal reaches a Constant Fraction of its Amplitude downhill (fitted)
                float timeOverThreshold; ///< is equal to (trailingEdgeTime - leadingEdgeTime)
                float constantFractionTime; ///< time extracted using the Constant Fraction Discriminator (CFD) algorithm (fitted)

		/** @brief Default constructor */
		AhdcExtractor() = default;

		/** @brief Constructor */
		AhdcExtractor(float _samplingTime,float _amplitudeFractionCFA, int _binDelayCFD, float _fractionCFD) :
			samplingTime(_samplingTime), amplitudeFractionCFA(_amplitudeFractionCFA), binDelayCFD(_binDelayCFD), fractionCFD(_fractionCFD) {}

		/** @brief Destructor */
		~AhdcExtractor(){;}

		/**
		 * This method extracts relevant informations from the digitized signal
		 * (the samples) and store them in a Pulse
		 *
		 * @param samples ADC samples
		 */
		std::map<std::string,double> extract(const std::vector<short> samples);

	private :
		/**
		* This method subtracts the pedestal (noise) from samples and stores it in : `samplesCorr`
		* It also computes a first value for : `adcMax`, `binMax`, `timeMax` and `integral`
		* This code is inspired by the one of coatjava/.../MVTFitter.java
		*/
		void waveformCorrection();

		/**
		 * This method gives a more precise value of the max of the waveform by computing the average of five points around the binMax
		 * It is an alternative to fitParabolic()
		 * The suitability of one of these fits can be the subject of a study
		 * Remark : This method updates adcMax but doesn't change timeMax
		 */
		void fitAverage();

		/**
		 * @brief Alternative to `fitAverage`
		 */
		void fitParabolic();

		/**
		 * From MVTFitter.java
		 * Make fine timestamp correction (using dream (=electronic chip) clock)
		 *
		 * Parameter dependency :
		 * - `timeStamp`
		 * - `fineTimeStampResolution`
		 */
		void fineTimeStampCorrection();

		/**
		 * This method determines the moment when the signal reaches a Constant Fraction of its Amplitude (i.e amplitudeFraction*adcMax)
		 * It fills the attributs : `leadingEdgeTime`, `trailingEdgeTime`, `timeOverThreshold`
		 *
		 * Parameter dependency :
		 * - `samplingTime` time between 2 ADC bins
		 * - `amplitudeFraction` amplitude fraction between 0 and 1
		 */
		void computeTimeAtConstantFractionAmplitude();

		/**
		 * This methods extracts a time using the Constant Fraction Discriminator (CFD) algorithm
		 * It fills the attribut : `constantFractionTime`
		 *
		 * Parameter dependency :
		 * - `samplingTime` time between 2 ADC bins
		 * - `fractionCFD` CFD fraction parameter between 0 and 1
		 * - `binDelayCFD` CFD delay parameter
		 */
		void computeTimeUsingConstantFractionDiscriminator();
	public:
		std::vector<float> samplesCFD; ///< samples corresponding to the CFD signal
		
		/** @brief Show the extracted data in a figure
		 * 
		 * @return generate pdf file
		 */
		void Show(const char * filename);

		/** @brief Show the CFD signal 
		 *
		 * @return generate a pdf file
		 */
		void ShowCFD(const char * filename);
};

#endif
