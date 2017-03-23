/*
* test_mcarray.cpp
*
* Copyright 2016 (c) Jordi Adell
* Created on: 2015
*       Author: Jordi Adell - adellj@gmail.com
*
* This file is part of MCARRAY
*
* MCARRAY is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* MCARRAY is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with MCARRAY.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <gtest/gtest.h>

#include <mcarray/ArrayModules.h>
#include <mcarray/SoundLocalisationCallback.h>
#include <mcarray/SoundLocalisationCallback.h>
#include <mcarray/MultibandBinarualLocalisation.h>
#include <mcarray/SourceSeparationAndLocalisation.h>
#include <mcarray/BinauralMaskingImpl.h>
#include <mcarray/SourceLocalisation.h>
#include <mcarray/FastBinauralMasking.h>
#include <mcarray/microhponeArrayHelpers.h>
#include <mcarray/ArrayDescription.h>

#include <dspone/algorithm/fft.h>
#include <dspone/filter/BandPassFIRFilter.h>
#include <dspone/algorithm/signalPower.h>

#include <wipp/wipputils.h>
#include <wipp/wippsignal.h>
#include <wipp/wippstats.h>

#include <fstream>

namespace mca {
namespace test {



inline std::string getAudioPath()
{
  std::string audio = "./audiofiles/";
  DEBUG_STREAM(audio);
  return audio;
}

inline std::string getTmpPath()
{
  std::string tmp = "./tmp/";
  return tmp;
}

//template <class T> std::string  saveBufferToFile(T* buffer, int length, std::string file)
//{
//  std::ofstream ofs;
//  util::file::open(getTmpPath()+"/"+file,ofs,true,false);
//  for (int i=0; i<length; ++i)
//    ofs << buffer[i] << std::endl;
//  return util::context::etcToAbs(getTmpPath()+"/"+file);
//}


// Helpers declaration
void testTemporalMaskingCore(dsp::ShortTimeProcess *masking, int sampleRate);
void testSpatialMaskingCore(dsp::ShortTimeProcess *masking);
void testLocalisationAndSeparationCore(std::string file, dsp::ShortTimeProcess *ls, std::vector<unsigned int> channels,
				       int sampleRate, bool useAntiAliasingFilter);
void assertLocalisationMessage(double doa, double prob, double power);
void testLocalisationCore(std::string file, dsp::ShortTimeAnalysis *loc, std::vector<unsigned int> channels,
			  int sampleRate, bool useAntiAliasingFilter);



// localisation callbacks to evaluate localisation results

class TestLocalisationCallback : public LocalisationCallback
{
    const double _mindoa;
    const double _maxdoa;
  public:
    TestLocalisationCallback(double mindoa, double maxoda) : _mindoa(mindoa), _maxdoa(maxoda) {}
    virtual void setDOA(SignalPtr doa, SignalPtr prob, double power, int numOfSources)
    {
      DEBUG_STREAM("DOA: " << doa[0] << " [" << _mindoa << ", " << _maxdoa << "]" << ", P:" << power << ", p: " << prob[0]);
      //            CPPUNIT_ASSERT(doa[0] >= _mindoa);
      EXPECT_GE(doa[0], _mindoa);
      //            CPPUNIT_ASSERT(doa[0] <= _maxdoa);
      EXPECT_LE(doa[0], _maxdoa);
    }
};

class TestMeanLocalisationCallback : public LocalisationCallback
{
    const double _mean;
    const double _stddev;
    int _count;
    double _sum;
    double _sum2;
  public:
    TestMeanLocalisationCallback(double mean, double stddev) :
      _mean(mean), _stddev(stddev), _count(0), _sum(0), _sum2(0)
    {
    }

    ~TestMeanLocalisationCallback()
    {
      double mean = 0;
      double var = 0;
      //      double stddev = 0;


      if (_count != 0)
      {
	mean = _sum/_count;
	var = _sum2/_count - mean*mean;
	//	stddev = sqrt(var);
      }

      DEBUG_STREAM("*** Results: mean " << mean << " stddev " << stddev);

      //            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("mean out of range", _mean, mean, 10);
      //            CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("standard deviation out of range", _stddev, stddev, 8);
    }

    virtual void setDOA(SignalPtr doa, SignalPtr prob, double power, int numOfSources)
    {
      //            DEBUG_STREAM("DOA: " << doa[0] << " [" << _mean << ", " << _stddev << "]" << ", P:" << power << ", p: " << prob[0]);
      _sum += doa[0];
      _sum2 += doa[0]*doa[0];
      ++_count;
    }

};

class TestInteractiveLocalisationCallback : public LocalisationCallback
{

  private:
    std::string _filename;
    std::string _extension;

  public:
    TestInteractiveLocalisationCallback(std::string filename)
    {
      _filename = filename.substr(0, filename.find_last_of("."));
      _extension = filename.substr(filename.find_last_of("."));
    }

    virtual void setDOA(SignalPtr doa, SignalPtr prob, double power, int numOfSources)
    {
      std::ofstream fs;
      std::ofstream fs_source;
      fs.open(_filename+_extension, std::ofstream::app);

      for (unsigned int s = 0; s < numOfSources; ++s)
      {
	DEBUG_STREAM("Source " << s << " --> DOA: " << doa[s] << ", P: " << power << ", p: " << prob[s]);
	std::string filename = _filename + "_source_" + std::to_string(s) + _extension;
	fs_source.open(filename, std::ofstream::app);
	fs_source << doa[s] << "\t" << power << "\t" << prob[s] << std::endl;
	fs_source.close();
      }

      for (unsigned int s = 0; s < numOfSources; ++s)
	fs << doa[s] << "\t";
      fs << std::endl;
      fs.close();
    }

};

// This callback only takes into account the first source
class TestInteractiveLocalisationStatisticsCallback : public LocalisationCallback
{
    const int _expectedDOA;
    double _mean;
    double _var;
    double _mse;
    int _counter;
    int _accepted;
    float _confidence;
    int _confidenceInterval;
    std::string _filename;

  public:
    TestInteractiveLocalisationStatisticsCallback(int expectedDOA, int confidenceInterval, std::string filename) :
      _expectedDOA(expectedDOA), _mean(0), _var(0), _mse(0), _counter(0), _accepted(0), _confidence(0),
      _confidenceInterval(confidenceInterval), _filename(filename) {}


    virtual void setDOA(SignalPtr doa, SignalPtr prob, double power, int numOfSources)
    {
      float prev_mean = _mean;

      _mean = (_mean*_counter + doa[0]) / (_counter + 1);
      _var = ((_var + prev_mean*prev_mean)*_counter + doa[0]*doa[0]) / (_counter+1) - _mean*_mean;
      if(_var < 0) // Var could be negative due to numerical precision
	_var = 0;
      _mse = (_mse*_counter + (doa[0]-_expectedDOA)*(doa[0]-_expectedDOA)) / (_counter + 1);

      if (abs(doa[0]-_expectedDOA) <= _confidenceInterval)
	_accepted++;

      _confidence = 100*_accepted/(_counter+1);

      _counter++;

      DEBUG_STREAM("DOA: " << doa[0] << ", P:" << power << ", p: " << prob[0]);
      DEBUG_STREAM("Expected DOA: " << _expectedDOA << ", Mean DOA: " << _mean << ", Var: " << _var << ", MSE: " << _mse << ", Confidence (+-"  << _confidenceInterval << "º): " << _confidence << "%");

      std::ofstream fs;
      fs.open(_filename, std::ofstream::app);
      fs << _expectedDOA << "\t" << _mean << "\t" << _var << "\t" << _mse << "\t" << _confidence << std::endl;
      fs.close();
    }
};


//actual test functions.

TEST(MicrophoneArrayTest, testTemporalMasking)
{
  int sampleRate = 16000;

  FastBinauralMasking fastmasking(sampleRate, 0.086, 500, 5000, BinauralMasking::FULL);
  {
    testTemporalMaskingCore(&fastmasking, sampleRate);
  }

  BinauralMaskingImpl masking(sampleRate, 0.086, 500, 5000, BinauralMaskingImpl::FULL);
  {
    testTemporalMaskingCore(&masking, sampleRate);
  }
}

TEST(MicrophoneArrayTest, testSpatialMasking)
{
  FastBinauralMasking fastmasking(16000, 0.086, 500, 5000, BinauralMasking::FULL);
  {
    testSpatialMaskingCore(&fastmasking);
  }

  BinauralMaskingImpl masking(16000, 0.086, 500, 5000, BinauralMaskingImpl::FULL);
  {
    testSpatialMaskingCore(&masking);
  }
  //::testing::UnitTest::elapsed_time()
}



TEST(MicrophoneArrayTest, testBinauralLocalisation)
{

  std::unique_ptr<LocalisationCallback> callback;

  // *********** Freq GCC Binaural Localisation Test *********** //
  {
    DEBUG_STREAM("Freq GCC binaural localisation");

    std::unique_ptr<FreqGCCBinauralLocalisation> fbl;

    std::string filename;
    std::string path = getAudioPath() + "/doa/recs_28_03_2013_44100Hz/angle_";

    // *** Simulation Parameters *** //
    int initialDOA = -60;
    int finalDOA = 60;
    int step = 10;
    int tolerance = 35;
    int sampleRate = 44100;
    double microDistance = 0.089;
    bool useAntiAliasingFilter = false;
    bool usePowerFloor = true; // tested with real audio

    ArrayDescription adesc;
    adesc.pushPosition(0, 0, 0);
    adesc.pushPosition(microDistance, 0, 0);

    for (int doa=initialDOA; doa<=finalDOA; doa+=step)
    {
      filename = path + std::to_string(doa) + ".raw";

      fbl.reset(new FreqGCCBinauralLocalisation(sampleRate, adesc, usePowerFloor));
      callback.reset(new TestMeanLocalisationCallback(doa, 0));

      fbl->setCallback(callback.get());
      testLocalisationCore(filename, fbl.get(), {0,1}, sampleRate, useAntiAliasingFilter);
      callback.reset();
    }
  }

  // *********** Temporal GCC Binaural Localisation Test *********** //
  {
    DEBUG_STREAM("Temporal GCC binaural localisation");

    std::unique_ptr<TemporalGCCBinauralLocalisation> tbl;

    std::vector<std::string> filenames;
    std::vector<std::vector<double> > acceptedRange;
    filenames.push_back(getAudioPath() + "/right90.raw");
    acceptedRange.push_back({-90, -30});
    filenames.push_back(getAudioPath() + "/right45.raw");
    acceptedRange.push_back({-90, 0});
    filenames.push_back(getAudioPath() + "/front.raw");
    acceptedRange.push_back({-20, 20});
    filenames.push_back(getAudioPath() + "/left45.raw");
    acceptedRange.push_back({0, 90});
    filenames.push_back(getAudioPath() + "/left90.raw");
    acceptedRange.push_back({30, 90});

    int sampleRate = 44100;
    float microDistance = 0.086;
    bool useAntiAliasingFilter = true;

    ArrayDescription adesc;
    adesc.pushPosition(0, 0, 0);
    adesc.pushPosition(microDistance, 0, 0);

    for (unsigned int i = 0; i < filenames.size(); ++i)
    {
      tbl.reset(new TemporalGCCBinauralLocalisation(sampleRate, adesc));
      callback.reset(new TestLocalisationCallback(acceptedRange[i][0], acceptedRange[i][1]));
      tbl->setCallback(callback.get());
      testLocalisationCore(filenames[i], tbl.get(), {0,1}, sampleRate, useAntiAliasingFilter);
    }
  }


}

TEST(MicrophoneArrayTest, testMultibandBinauralLocalisation)
{

  // *** Simulation Parameters *** //
  int initialDOA = -90;
  int finalDOA = 90;
  int step = 10;
  int tolerance = 15;
  int sampleRate = 48000;
  float microDistance = 0.086;
  int numBins = 25;
  bool usePowerFloor = false; // tested with sines
  bool useAntiAliasingFilter = false;

  std::string filename;
  std::string path = getAudioPath() + "/doa/sine_f_1000_fs_48000_microdistance_0.086_duration_1/sine_angle_";

  std::unique_ptr<TestLocalisationCallback> callback;
  std::unique_ptr<MultibandBinarualLocalisation> mbl;

  ArrayDescription adesc;
  adesc.pushPosition(0, 0, 0);
  adesc.pushPosition(microDistance, 0, 0);

  // *********** Multiband Binaural Localisation Test *********** //

  DEBUG_STREAM("Multi-band binaural localisation");

  for (int doa=initialDOA; doa<=finalDOA; doa+=step)
  {
    filename = path + std::to_string(doa) + ".raw";

    callback.reset(new TestLocalisationCallback(doa-tolerance, doa+tolerance));
    mbl.reset(new MultibandBinarualLocalisation(sampleRate, adesc, numBins, usePowerFloor));

    mbl->setCallback(callback.get());
    testLocalisationCore(filename, mbl.get(), {0,1}, sampleRate, useAntiAliasingFilter);
  }
}

TEST(MicrophoneArrayTest, testBeamformingSoundLocalisation)
{
  // *** Simulation Parameters *** //
  int initialDOA = -80;
  int finalDOA = 80;
  int step = 10;
  int tolerance = 7;
  int sampleRate = 48000;
  int numOfSources = 1;
  bool usePowerFloor = false; // tested with sines
  bool useAntiAliasingFilter = false;

  ArrayDescription microPositions =
      ArrayDescription::make_linear_array_description({0, 0.07, 0.175, 0.21}); //Reem C

  std::vector<unsigned int> channels;
  for (unsigned int c = 0; c < microPositions.size(); ++c)
    channels.push_back(c);

  std::string filename;
  std::string path = getAudioPath() + "/doa/sine_4ch_ReemC/sines_f_1000_fs_48000/sine_DOA_";

  std::unique_ptr<TestLocalisationCallback> callback;
  std::unique_ptr<SourceLocalisation> ssl;

  // *********** Beamforming Sound Localisation Test *********** //

  DEBUG_STREAM("Beamforming Sound Localisation");

  for (int doa=initialDOA; doa<=finalDOA; doa+=step)
  {
    filename = path + std::to_string(doa) + ".raw";

    callback.reset(new TestLocalisationCallback(doa-tolerance, doa+tolerance));
    ssl.reset(new SourceLocalisation(sampleRate, microPositions, numOfSources, usePowerFloor));

    ssl->setCallback(callback.get());
    testLocalisationCore(filename, ssl.get(), channels, sampleRate, useAntiAliasingFilter);
  }
}

TEST(MicrophoneArrayTest, testInteractiveLocalisation)
{

  return;

  // *** Select file path *** //
  std::string filename = "";

  // *** Simulation Parameters *** //
  int tolerance = 10;
  int doa = 0;
  int sampleRate = 48000;
  int numBins = 15; // Used in algorithm 3 (multiband)
  bool usePowerFloor = true;
  bool useAntiAliasingFilter = false;
  std::string outputFilename = "/home/student/output.txt";

  int algorithm, cb;
  std::string s;
  std::cout << "Select the localisation algorithm to be used:\n";
  std::cout << "\t(1) Temporal GCC Binaural\n";
  std::cout << "\t(2) Freq GCC Binaural\n";
  std::cout << "\t(3) Multiband GCC Binaural\n";
  std::cout << "\t(4) Steering Beamforming\n";
  std::cin >> s;
  algorithm = std::stoi(s);
  std::cout << "\nSelect the callback to be used:\n";
  std::cout << "\t(1) Interactive callback (print DOA on screen and store in output file)\n";
  std::cout << "\t(2) Statistics callback (compute statistics and store them in output file)\n";
  std::cin >> s;
  cb = std::stoi(s);
  std::cout << "\nUse power floor? (y/n)\n";
  std::cin >> s;
  if (s == "n")
    usePowerFloor = false;
  else
    usePowerFloor = true;
  std::cout << "\nUse antialiasing filter? (y/n)\n";
  std::cin >> s;
  if (s == "n")
    useAntiAliasingFilter = false;
  else
    useAntiAliasingFilter = true;


  // Used in Binaural algorithms (1...3)
  float microDistance = 0.089;
  ArrayDescription twoMicsDescription =
      ArrayDescription::make_linear_array_description({0, microDistance});

  // Used in beamforming algorithm (a number of microphones)
  ArrayDescription microPositions =
      ArrayDescription::make_linear_array_description({0, 0.07, 0.175, 0.21});
  unsigned int nchannels = microPositions.size();
  std::vector<unsigned int> channels;
  for (unsigned int c = 0; c < nchannels; ++c)
    channels.push_back(c);
  int numOfSources = 2;

  std::unique_ptr<dsp::ShortTimeAnalysis> bl;
  std::unique_ptr<LocalisationCallback> callback;

  if (cb == 1) // Standard callback (print the DOA + assert)
    callback.reset(new TestInteractiveLocalisationCallback(outputFilename));
  else if (cb == 2) // Error callback (compute statistics and store them in an output file.
    callback.reset(new TestInteractiveLocalisationStatisticsCallback(doa, tolerance, outputFilename));

  if (algorithm == 1) // Temporal GCC Binaural
  {
    bl.reset(new TemporalGCCBinauralLocalisation(sampleRate, twoMicsDescription));
    dynamic_cast<TemporalGCCBinauralLocalisation*>(bl.get())->setCallback(callback.get());
  }
  else if (algorithm == 2) // Freq GCC Binaural
  {
    std::cout << "\nNum of bins: ";
    bl.reset(new FreqGCCBinauralLocalisation(sampleRate, twoMicsDescription, usePowerFloor));
    dynamic_cast<FreqGCCBinauralLocalisation*>(bl.get())->setCallback(callback.get());
  }
  else if (algorithm == 3) // Multiband Binaural
  {
    bl.reset(new MultibandBinarualLocalisation(sampleRate, twoMicsDescription, numBins, usePowerFloor));
    dynamic_cast<MultibandBinarualLocalisation*>(bl.get())->setCallback(callback.get());
  }
  else if (algorithm == 4) // Steering beamforming
  {
    bl.reset(new SourceLocalisation(sampleRate, microPositions, numOfSources, usePowerFloor));
    dynamic_cast<SourceLocalisation*>(bl.get())->setCallback(callback.get());
  }

  testLocalisationCore(filename, bl.get(), channels, sampleRate, useAntiAliasingFilter);
}


TEST(MicrophoneArrayTest, testArrayDescription)
{
  ArrayDescription description;
  ArrayDescription::ElementId l, cl, cr, r;

  EXPECT_EQ(description.size(), 0);
  l  = description.pushPosition(0       /*0cm */,  2, 3, "left");
  EXPECT_EQ(description.size(), 1);
  cl = description.pushPosition(0.035*2 /*7cm*/,   2, 3, "center-left");
  EXPECT_EQ(description.size(), 2);
  cr = description.pushPosition(0.035*5 /*17.5cm*/, 2, 3, "center-right");
  EXPECT_EQ(description.size(), 3);
  r  = description.pushPosition(0.035*6 /*21cm*/,   2, 3, "right");
  EXPECT_EQ(description.size(), 4);

  EXPECT_EQ (description.getName(l), "left");
  EXPECT_EQ (description.getName(cl), "center-left");
  EXPECT_EQ (description.getName(cr), "center-right");
  EXPECT_EQ (description.getName(r),  "right");

  EXPECT_DOUBLE_EQ(description.getX(l),  0.000);
  EXPECT_DOUBLE_EQ(description.getX(cl), 0.070);
  EXPECT_DOUBLE_EQ(description.getX(cr), 0.175);
  EXPECT_DOUBLE_EQ(description.getX(r),  0.210);

  EXPECT_DOUBLE_EQ(description.getY(l),  2);
  EXPECT_DOUBLE_EQ(description.getY(cl), 2);
  EXPECT_DOUBLE_EQ(description.getY(cr), 2);
  EXPECT_DOUBLE_EQ(description.getY(r),  2);

  EXPECT_DOUBLE_EQ(description.getZ(l),  3);
  EXPECT_DOUBLE_EQ(description.getZ(cl), 3);
  EXPECT_DOUBLE_EQ(description.getZ(cr), 3);
  EXPECT_DOUBLE_EQ(description.getZ(r),  3);

  EXPECT_DOUBLE_EQ(description.distance(0,0),  0.000);
  EXPECT_DOUBLE_EQ(description.distance(0,1),  0.070);
  EXPECT_DOUBLE_EQ(description.distance(0,2),  0.175);
  EXPECT_DOUBLE_EQ(description.distance(0,3),  0.210);

  EXPECT_DOUBLE_EQ(description.distance(1,0),  0.070);
  EXPECT_DOUBLE_EQ(description.distance(1,1),  0.000);
  EXPECT_DOUBLE_EQ(description.distance(1,2),  0.105);
  EXPECT_DOUBLE_EQ(description.distance(1,3),  0.140);

  EXPECT_DOUBLE_EQ(description.distance(2,0),  0.175);
  EXPECT_DOUBLE_EQ(description.distance(2,1),  0.105);
  EXPECT_DOUBLE_EQ(description.distance(2,2),  0.000);
  EXPECT_DOUBLE_EQ(description.distance(2,3),  0.035);

  EXPECT_DOUBLE_EQ(description.distance(3,0),  0.210);
  EXPECT_DOUBLE_EQ(description.distance(3,1),  0.140);
  EXPECT_DOUBLE_EQ(description.distance(3,2),  0.035);
  EXPECT_DOUBLE_EQ(description.distance(3,3),  0.000);

  EXPECT_DOUBLE_EQ(description.maxDistance(),  0.210);

  EXPECT_DOUBLE_EQ(description.distance("left","center-left"),         0.070);
  EXPECT_DOUBLE_EQ(description.distance("left","center-right"),        0.175);
  EXPECT_DOUBLE_EQ(description.distance("center-right","center-left"), 0.105);
  EXPECT_DOUBLE_EQ(description.distance("right","center-left"),        0.140);

}



// Helpers implementation
void testLocalisationCore(std::string file, dsp::ShortTimeAnalysis *loc, std::vector<unsigned int> channels,
			  int sampleRate, bool useAntiAliasingFilter)
{

  //FileAudioRecorder recorder(file, sampleRate, channels.size(), false, false);
  dsp::BandPassFIRFilter filter(1024, 100.0/sampleRate, 1900.0/sampleRate);

  DEBUG_STREAM("====================================================");
  DEBUG_STREAM("Trying to localise sound in " << file);

  int samplesProvided = 0;
  int buffersize = 3*loc->getFrameSize();

  SignalVector input;
  SignalVector filtered;
  for (unsigned int c = 0; c < channels.size(); ++c)
  {
    input.push_back(SignalPtr(new BaseType[buffersize]));
    filtered.push_back(SignalPtr(new BaseType[buffersize]));
  }

  //  recorder.startAudioRecord();

  do{

    //    recorder.getAudioChannels(input, buffersize, &samplesProvided, channels);
    samplesProvided = -1;
    if (samplesProvided > 0)
    {
      if (useAntiAliasingFilter)
      {
	for (unsigned int c = 0; c < channels.size(); ++c)
	  filter.filterBuffer(input[c].get(), filtered[c].get(), samplesProvided);
	loc->process(filtered, samplesProvided);
      }
      else
      {
	loc->process(input, samplesProvided);
      }
      DEBUG_STREAM(samplesProvided << " samples processed");
    }
  }while(samplesProvided > -1);

  //  recorder.stopAudioRecord();
}

void testBeamformingSeparation()
{

  std::vector<std::string> filenames;
  std::vector<double> freqSignal1;
  std::vector<double> freqSignal2;
  std::vector<double> doaSignal1;
  std::vector<double> doaSignal2;

  filenames.push_back(getAudioPath() + "/doa/sine_4ch_ReemC/beamformingSeparationTest/sine_f_800_fs_48000_SIR_0_freq_interf_2000_doa_interf_-45/sine_DOA_45.raw");
  freqSignal1.push_back(800);
  freqSignal2.push_back(2000);
  doaSignal1.push_back(45*M_PI/180);
  doaSignal2.push_back(-45*M_PI/180);

  filenames.push_back(getAudioPath() + "/doa/sine_4ch_ReemC/beamformingSeparationTest/sine_f_1000_fs_48000_SIR_0_freq_interf_4000_doa_interf_-20/sine_DOA_20.raw");
  freqSignal1.push_back(1000);
  freqSignal2.push_back(4000);
  doaSignal1.push_back(20*M_PI/180);
  doaSignal2.push_back(-20*M_PI/180);

  filenames.push_back(getAudioPath() + "/doa/sine_4ch_ReemC/beamformingSeparationTest/sine_f_1000_fs_48000_SIR_0_freq_interf_1500_doa_interf_10/sine_DOA_80.raw");
  freqSignal1.push_back(1000);
  freqSignal2.push_back(1500);
  doaSignal1.push_back(80*M_PI/180);
  doaSignal2.push_back(10*M_PI/180);

  // *** Simulation Parameters *** //
  int sampleRate = 48000;
  int fftOrder = 11;
  int fftLength = pow(2, fftOrder);
  int fftCCSLength = fftLength + 2;
  int samplesProvided = 0;

  // Array Description
  ArrayDescription microPositions =
      ArrayDescription::make_linear_array_description({0, 0.07, 0.175, 0.21}); //Reem C
  std::vector<unsigned int> channels = {0, 1, 2, 3};
  unsigned int nchannels = channels.size();

  Beamformer beamformer(sampleRate, microPositions, fftCCSLength, nchannels);
  dsp::FFT fft(fftOrder);


  SignalVector16s inputSamples; // time-domain input samples
  SignalPtr auxInput; //used to convert from BaseType16s to BaseType
  SignalVector analysisFrames; // freq-domain input
  SignalPtr outputFrames; // freq-domain output
  SignalPtr outputSamples; // time-domain output samples (BaseType64)
  SignalPtr32 outputSamples32; // time-domain output samples (BaseType32)
  SignalPtr magnitude;
  BaseType32 zeroCrossRate; // Zero Crossing Rate

  for (unsigned int c = 0; c < nchannels; ++c)
  {
    inputSamples.push_back(SignalPtr16s(new BaseType16s[fftLength]));
    analysisFrames.push_back(SignalPtr(new BaseType[fftCCSLength]));
  }
  outputFrames.reset(new BaseType[fftCCSLength]);
  outputSamples.reset(new BaseType[fftCCSLength]);

  outputSamples32.reset(new BaseType32[fftCCSLength]);
  auxInput.reset(new BaseType[fftLength]);
  magnitude.reset(new BaseType[fftCCSLength/2]);


  for (unsigned int f = 0; f < filenames.size(); ++f)
  {

    DEBUG_STREAM("====================================================");
    DEBUG_STREAM("Trying to localise sound in " << filenames[f]);

    DEBUG_STREAM("Signal 1: DOA = " << toDegrees(doaSignal1[f])  << ". Freq = " << freqSignal1[f]);
    DEBUG_STREAM("Signal 2: DOA = " << toDegrees(doaSignal2[f])  << ". Freq = " << freqSignal2[f]);

    //recorder(filenames[f], sampleRate, nchannels, false, false);
    //recorder.startAudioRecord();

    do{

      //      recorder.getAudioChannels(inputSamples, fftLength, &samplesProvided, channels);

      if (samplesProvided == fftLength)
      {
	double inputPowerSignal1 = 0;
	double inputPowerSignal2 = 0;
	double power;

	// Compute FFT (analysisFrames) and mean power of input signals
	for (unsigned int ch = 0; ch < nchannels; ++ch)
	{
	  wipp::copyBuffer(inputSamples[ch].get(), auxInput.get(), fftLength); // conversion from SignalTypePtr to BaseType
	  fft.fwdTransform(auxInput.get(), analysisFrames[ch].get());
	  wipp::magnitude(reinterpret_cast<wipp::wipp_complex_t*>(analysisFrames[ch].get()), magnitude.get(), fftCCSLength/2);
	  // Find the delta around the position where it should be placed.

	  wipp::max(&(magnitude[(int)(freqSignal1[f]/sampleRate*fftLength)-4]), 8, &power);
	  inputPowerSignal1 += power;
	  wipp::max(&(magnitude[(int)(freqSignal2[f]/sampleRate*fftLength)-4]), 8, &power);
	  inputPowerSignal2 += power;
	}
	inputPowerSignal1 /= nchannels;
	inputPowerSignal2 /= nchannels;

	double outputPowerSignal1;
	double outputPowerSignal2;
	double attenuation;
	double computedFreq;

	// Pointing to Signal 1
	beamformer.processFrame(analysisFrames, outputFrames, doaSignal1[f]);
	wipp::magnitude(reinterpret_cast<wipp::wipp_complex_t*>(outputFrames.get()), magnitude.get(), fftCCSLength/2);
	wipp::max(&(magnitude[(int)(freqSignal1[f]/sampleRate*fftLength)-4]), 8, &outputPowerSignal1);
	wipp::max(&(magnitude[(int)(freqSignal2[f]/sampleRate*fftLength)-4]), 8, &outputPowerSignal2);
	attenuation = 20*log10((outputPowerSignal1/inputPowerSignal1) / (outputPowerSignal2/inputPowerSignal2));

	fft.invTransfrom(outputSamples.get(), outputFrames.get());
	wipp::copyBuffer(outputSamples.get(), outputSamples32.get(), fftLength);
	//ippsZeroCrossing_32f(outputSamples32.get(), fftLength, &zeroCrossRate, ippZCR);
	computedFreq = zeroCrossRate / fftLength * sampleRate / 2;

	DEBUG_STREAM("Pointing to Signal 1... attenuation = " << attenuation  << " dB, freq = " << computedFreq << " Hz (expected = " << freqSignal1[f] << " Hz)");
	if(attenuation >= 5.5)
	{
	  ERROR_STREAM("Attenuation is less than 5.5dB.");
	}
	EXPECT_GE(attenuation, 5.5);

	//                CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Signal 1 frequency out of expected interval", freqSignal1[f], computedFreq, 45);

	if ( fabs(freqSignal1[f] - computedFreq) < 45)
	{
	  ERROR_STREAM("Signal 1 frequency out of expected interval");
	}
	EXPECT_DOUBLE_EQ(freqSignal1[f], computedFreq);


	// Pointing to Signal 2
	beamformer.processFrame(analysisFrames, outputFrames, doaSignal2[f]);
	wipp::magnitude(reinterpret_cast<wipp::wipp_complex_t*>(outputFrames.get()), magnitude.get(), fftCCSLength/2);
	wipp::max(&(magnitude[(int)(freqSignal1[f]/sampleRate*fftLength)-4]), 8, &outputPowerSignal1);
	wipp::max(&(magnitude[(int)(freqSignal2[f]/sampleRate*fftLength)-4]), 8, &outputPowerSignal2);
	attenuation = 20*log10((outputPowerSignal2/inputPowerSignal2) / (outputPowerSignal1/inputPowerSignal1));

	fft.invTransfrom(outputSamples.get(), outputFrames.get());
	wipp::copyBuffer(outputSamples.get(), outputSamples32.get(), fftLength);
	//ippsZeroCrossing_32f(outputSamples32.get(), fftLength, &zeroCrossRate, ippZCR);
	computedFreq = zeroCrossRate / fftLength * sampleRate / 2;

	DEBUG_STREAM("Pointing to Signal 2... attenuation = " << attenuation  << " dB, freq = " << computedFreq << " Hz (expected = " << freqSignal2[f] << " Hz)");
	if(attenuation >= 5.5)
	{
	  ERROR_STREAM("Attenuation is less than 5.5dB.");
	}
	EXPECT_GE(attenuation, 5.5);
	//                CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Signal 2 frequency out of expected interval", freqSignal2[f], computedFreq, 45);
	if ( fabs(freqSignal2[f] - computedFreq) < 45)
	{
	  ERROR_STREAM("Signal 1 frequency out of expected interval");
	}
	EXPECT_DOUBLE_EQ(freqSignal2[f], computedFreq);


	DEBUG_STREAM(samplesProvided << " samples processed");
      }
    }while(samplesProvided == fftLength);

    //recorder.stopAudioRecord();
  }
}

void testInteractiveLocalisationAndSeparation()
{

  return;

  // *** Select file path *** //
  std::string filename = "/home/student/Audios/ReemC/recs_two_sources_44100/angles_0_and_80.wav";

  // *** Simulation Parameters *** //
  bool usePowerFloor = true;
  bool useAntiAliasingFilter = false;
  int sampleRate = 44100;
  int numOfSources = 2;
  std::string outputFilename = "/home/student/output.txt";

  // Array Description
  ArrayDescription microPositions =
      ArrayDescription::make_linear_array_description({0.29, 0.348, 0.45, 0.485}); // Experiments
  std::vector<unsigned int> channels = {0, 1, 2, 3};

  std::unique_ptr<SourceSeparationAndLocalisation> ls;
  ls.reset(new SourceSeparationAndLocalisation(sampleRate, microPositions, numOfSources, usePowerFloor));

  std::unique_ptr<TestInteractiveLocalisationCallback> callback;
  callback.reset(new TestInteractiveLocalisationCallback(outputFilename));

  ls->setCallback(callback.get());
  testLocalisationAndSeparationCore(filename, ls.get(), channels, sampleRate, useAntiAliasingFilter);
}

void testLocalisationAndSeparationCore(std::string file, dsp::ShortTimeProcess *ls, std::vector<unsigned int> channels,
				       int sampleRate, bool useAntiAliasingFilter)
{
  // FileAudioRecorder recorder(file, sampleRate, channels.size(), false, false);
  dsp::BandPassFIRFilter filter(1024, 100.0/sampleRate, 1900.0/sampleRate);

  DEBUG_STREAM("====================================================");
  DEBUG_STREAM("Trying to localise sound in " << file);

  int inputSamplesProvided = 0;
  int outputSamplesProvided = 0;
  const int buffersize = 3*ls->getFrameSize();
  const int outbuffersize = buffersize + ls->getMaxLatency();

  std::ofstream ofs;
  ofs.open(file.substr(0, file.find_last_of(".")) + "_output.raw", std::ofstream::app);

  SignalVector16s input, filtered, output;
  for (unsigned int c = 0; c < channels.size(); ++c)
  {
    input.push_back(SignalPtr16s(new BaseType16s[buffersize]));
    filtered.push_back(SignalPtr16s(new BaseType16s[buffersize]));
    output.push_back(SignalPtr16s(new BaseType16s[outbuffersize]));
  }

  //recorder.startAudioRecord();

  do{

    //recorder.getAudioChannels(input, buffersize, &inputSamplesProvided, channels);

    if (inputSamplesProvided > 0)
    {
      if (useAntiAliasingFilter)
      {
	for (unsigned int c = 0; c < channels.size(); ++c)
	  filter.filterBuffer(input[c].get(), filtered[c].get(), inputSamplesProvided);
	outputSamplesProvided = ls->process(filtered, inputSamplesProvided, output, outbuffersize);
      }
      else
      {
	outputSamplesProvided = ls->process(input, inputSamplesProvided, output, outbuffersize);
      }

      for (int pos = 0; pos < outputSamplesProvided; ++pos)
      {
	for (unsigned int channel = 0; channel < channels.size(); ++channel)
	{
	  //	  BaseType16s *ptr = output[channel].get();
	  //	  ofs.write(reinterpret_cast<char*>(&ptr[pos]), bytesPerFrame());
	}
      }

      DEBUG_STREAM(inputSamplesProvided << " samples processed");
    }
  }while(inputSamplesProvided > -1);

  //  recorder.stopAudioRecord();
}

void testSpatialMaskingCore(dsp::ShortTimeProcess *masking)
{

  int buffersize = 5*1024;
  int outbuffersize = buffersize + masking->getMaxLatency();
  BaseType16s interestLeft[buffersize];
  BaseType16s interestRight[buffersize];
  BaseType16s interferenceLeft[buffersize];
  BaseType16s interferenceRight[buffersize];
  BaseType16s signalLeft[buffersize];
  BaseType16s signalRight[buffersize];
  BaseType16s outputLeft[outbuffersize];
  BaseType16s outputRight[outbuffersize];
  BaseType16s filteredLeft[outbuffersize];
  BaseType16s filteredRight[outbuffersize];

  int magn = 5000;
  int delay = 6;
  float sFreq = 0.1, iFreq=0.3;
  float phase = 0;

  dsp::BandPassFIRFilter signalFilter(256, 0.05, 0.15);
  dsp::BandPassFIRFilter interfFilter(256, 0.25, 0.35);

//  ippsTone_Direct_16s(interestLeft,  buffersize, magn, sFreq, &phase, ippAlgHintNone);
  wipp::tone(interestLeft, buffersize, magn, sFreq, phase);
  wipp::copyBuffer(interestLeft, interestRight, buffersize);

  phase=0;
  wipp::setZeros(interferenceRight, buffersize);
  wipp::tone(interferenceLeft, buffersize, magn, iFreq, phase);
//  ippsTone_Direct_16s(interferenceLeft,  buffersize, magn, iFreq, &phase, ippAlgHintNone);
  wipp::copyBuffer(&interferenceLeft[delay], interferenceRight, buffersize-delay);

  wipp::copyBuffer(interestLeft,  signalLeft,  buffersize);
  wipp::copyBuffer(interestRight, signalRight, buffersize);
  wipp::add(interferenceLeft,  signalLeft,  buffersize);
  wipp::add(interferenceRight, signalRight, buffersize);

  std::vector<BaseType16s*> signal, output;
  signal.push_back(signalLeft);
  signal.push_back(signalRight);
  output.push_back(outputLeft);
  output.push_back(outputRight);

  int processedSamples = masking->process(signal, buffersize, output, outbuffersize);

  signalFilter.filterBuffer(outputLeft, filteredLeft, processedSamples);

  float power = dsp::SignalPower::logPower(filteredLeft, processedSamples);

  EXPECT_LE(fabs(70-power), 10);

  //          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Important signal has been removed in binaural masking", 70, power, 10);
  signalFilter.filterBuffer(outputLeft, filteredRight, processedSamples);
  power = dsp::SignalPower::logPower(filteredRight, processedSamples);
  EXPECT_LE(fabs(70-power), 10);
  //          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Important signal has been removed in binaural masking", 70, power, 10);

  power = dsp::SignalPower::logPower(filteredLeft, processedSamples);
  EXPECT_LE(fabs(70-power), 10);
  //          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Interferring signal has not been properly removed in binaural masking", 70, power, 10);
  interfFilter.filterBuffer(outputLeft, filteredRight, processedSamples);
  power = dsp::SignalPower::logPower(filteredRight, processedSamples);
  EXPECT_LE(fabs(70-power), 10);
  //          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Interferring signal has not been properly removed in binaural masking", 70, power, 10);
}

void testTemporalMaskingCore(dsp::ShortTimeProcess *masking, int sampleRate)
{
  const int magn = 5000;
  const int delay = 0.1*sampleRate;
  const int tonestep = 0.1*sampleRate;
  const float freqstep = 0.01;
  float sFreq = 0.01;
  const int nreverb = 5;
  float phase = 0;
  int interestStart = 0;
  float interestFreq = 0.2;


  dsp::BandPassFIRFilter signalFilter(256, 0.19, 0.21);
  dsp::BandPassFIRFilter interfFilter(256, 0.15, 0.20);

  int buffersize = 50*1024;
  int outbuffersize = buffersize + masking->getMaxLatency();

  BaseType16s tone_buffer[buffersize];

  BaseType16s interferenceLeft[buffersize];
  BaseType16s interferenceRight[buffersize];
  BaseType16s signalLeft[buffersize];
  BaseType16s signalRight[buffersize];
  BaseType16s outputLeft[outbuffersize];
  BaseType16s outputRight[outbuffersize];
  BaseType16s filteredLeft[outbuffersize];
  BaseType16s filteredRight[outbuffersize];

  wipp::setZeros(tone_buffer, buffersize);

  int i = 0;
  for (;
       i < buffersize-tonestep && sFreq<0.5;
       i += tonestep, sFreq += freqstep)
  {
    if (interestFreq - freqstep/2 < sFreq && sFreq < interestFreq + freqstep/2)
      interestStart = i;
    //    ippsTone_Direct_16s(&tone_buffer[i],  tonestep, magn, sFreq, &phase, ippAlgHintNone);
    wipp::tone(&tone_buffer[i], tonestep, magn, sFreq, phase);
  }

  DEBUG_STREAM("Filled up to " << i << " samples of " << buffersize << " and up to " << sFreq << "1/s");
  DEBUG_STREAM("Startindex: " << interestStart);

  std::vector<BaseType16s*> signal, output;
  signal.push_back(signalLeft);
  signal.push_back(signalRight);
  output.push_back(outputLeft);
  output.push_back(outputRight);

  wipp::setZeros(signalLeft, buffersize);
  wipp::setZeros(signalRight, buffersize);

  for (i = 0; i <= nreverb; ++i)
  {
    wipp::divC(2, tone_buffer, buffersize);
    wipp::add(tone_buffer, &signalLeft[delay*i], buffersize-delay*i);
  }
  wipp::copyBuffer(signalLeft, signalRight, buffersize);

  //          saveBufferToFile(signal[0], buffersize, "sig.txt");
  int outLength = masking->process(signal, buffersize, output, outbuffersize);
  //          saveBufferToFile(output[0], outLength, "out.txt");


  signalFilter.filterBuffer(signal[0], filteredLeft, outLength);
  signalFilter.filterBuffer(signal[1], filteredRight, outLength);
  DEBUG_STREAM("Interest Band L: " << dsp::SignalPower::logPower(&filteredLeft[interestStart], tonestep) << " R: " << dsp::SignalPower::logPower(&filteredRight[interestStart], tonestep));

  interfFilter.filterBuffer(signal[0], interferenceLeft, outLength);
  interfFilter.filterBuffer(signal[1], interferenceRight, outLength);
  DEBUG_STREAM("Interference Band L: " << dsp::SignalPower::logPower(&interferenceLeft[interestStart], tonestep) << " R: " << dsp::SignalPower::logPower(&interferenceRight[interestStart], tonestep));



  double powerDB = dsp::SignalPower::logPower(&filteredLeft[interestStart], tonestep) -
		   dsp::SignalPower::logPower(&interferenceLeft[interestStart], tonestep);
  EXPECT_LT(fabs(powerDB - 2), 0.5);
  //          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Something happened to the way we contruct the transmitted signal", 2, powerDB, 0.5);



  powerDB = dsp::SignalPower::logPower(&filteredRight[interestStart], tonestep) - dsp::SignalPower::logPower(&interferenceRight[interestStart], tonestep);
  EXPECT_LT(fabs(powerDB-2), 0.5);
  //          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Something happened to the way we contruct the transmitted signal", 2,powerDB,0.5);


  signalFilter.filterBuffer(output[0], filteredLeft, outLength);
  signalFilter.filterBuffer(output[1], filteredRight, outLength);
  DEBUG_STREAM("Interest Band L: " << dsp::SignalPower::logPower(&filteredLeft[interestStart], tonestep) << " R: " << dsp::SignalPower::logPower(&filteredRight[interestStart], tonestep));

  interfFilter.filterBuffer(output[0], interferenceLeft, outLength);
  interfFilter.filterBuffer(output[1], interferenceRight, outLength);
  DEBUG_STREAM("Interference Band L: " << dsp::SignalPower::logPower(&interferenceLeft[interestStart], tonestep) << " R: " << dsp::SignalPower::logPower(&interferenceRight[interestStart], tonestep));

  //          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Reverb has not been removed as intended", 5,
  //                                               dsp::SignalPower::logPower(&filteredRight[interestStart], tonestep) - dsp::SignalPower::logPower(&interferenceRight[interestStart], tonestep), 1);
  //          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Reverb has not been removed as intended", 5,
  //                                              dsp::SignalPower::logPower(&filteredLeft[interestStart], tonestep) - dsp::SignalPower::logPower(&interferenceLeft[interestStart], tonestep), 1);
  powerDB = dsp::SignalPower::logPower(&filteredRight[interestStart], tonestep) - dsp::SignalPower::logPower(&interferenceRight[interestStart], tonestep);
  EXPECT_LT(fabs(powerDB-5), 1);
  powerDB = dsp::SignalPower::logPower(&filteredLeft[interestStart],   tonestep) - dsp::SignalPower::logPower(&interferenceLeft[interestStart],  tonestep);
  EXPECT_LT(fabs(powerDB-5), 1);
}

void assertLocalisationMessage(double doa, double prob, double power)
{
  DEBUG_STREAM("DOA: " << doa << " PROB: " << prob << " POWER: " << power);
  EXPECT_GT(doa, -90);
  EXPECT_LT(doa,  90);
  EXPECT_GT(prob,  0);
  EXPECT_GT(power, 0);
}

}
}



int
main(int argc, char *argv[])
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

