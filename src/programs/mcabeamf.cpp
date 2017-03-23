#include <mcarray/mcalogger.h>
#include <mcarray/ArrayDescription.h>
#include <mcarray/SourceSeparationAndLocalisation.h>
#include <mcarray/ArrayModules.h>

#include <dspone/rt/ShortTimeFourierTransform.h>

#ifdef DSPONE_GUI
#include <dspone/plot/dspgui.h>
#endif

#include <sndfile.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <functional>

using namespace std::placeholders;

void usage(int argc, char *argv[])
{
  std::cout << "Use: " << argv[0] << " [options]" << std::endl;
  std::cout << "  options:" << std::endl;
  std::cout << "      -i file    Input file" << std::endl;
  std::cout << "      -o file    Output file" << std::endl;
  std::cout << "      -a file    Array description [not implemented]" << std::endl;
  std::cout << "      -d file    Print DOA in file [use - for stdout]" << std::endl;
  std::cout << "      -h         This help message" << std::endl;
  std::cout << std::endl;
}


class dumystft : public dsp::STFT
{
  public:
    dumystft(int channels) : STFT(channels, dsp::STFT::_defaultFFTOrder)
    {

    }

    virtual ~dumystft() {}

  private:
    virtual void processParametrisation(std::vector<double *> &analysisFrames, int analysisLength,
					std::vector<double *> &dataChannels, int dataLength)
    {

    }
};


class McBeamCallback : public mca::LocalisationCallback
{
  public:
    McBeamCallback(std::ostream &os) : os_(os) {}
    virtual ~McBeamCallback() {}

    virtual void setDOA(mca::SignalPtr doas,
			mca::SignalPtr probs,
			double power,
			int numOfSources)
    {
      for (int i = 0; i < numOfSources; ++i)
      {
	os_ << "[DOA: " << doas[i] << ", p=" << probs[i] << ", P=" << power << "] ";
      }
      os_ << std::endl;
    }

  private:
    std::ostream &os_;
};



void process_file(dsp::ShortTimeProcess &processor,
		  SNDFILE *snd_file_in, SNDFILE *snd_file_out,
		  const SF_INFO &sf_info)
{

  int nchannels = sf_info.channels;
  int nframes = 1 << 10;
  int buffer_length = nframes * nchannels;
  int out_buffer_length = (nframes + processor.getMaxLatency()) * nchannels;
  double buffer[buffer_length];
  double output_buffer[out_buffer_length];

  std::vector<double*> raw_audio_channels, raw_output;
  mca::SignalVector audio_channels, output;

  for (int c = 0; c < nchannels; ++c)
  {
    audio_channels.push_back(mca::SignalPtr(new double[nframes]));
    raw_audio_channels.push_back(audio_channels.back().get());
    output.push_back(mca::SignalPtr(new double[nframes + processor.getMaxLatency()]));
    raw_output.push_back(output.back().get());
  }

  int read_frames = 0;
  while( (read_frames = sf_readf_double(snd_file_in, buffer, nframes)) > 0)
  {

    for (int f = 0; f < read_frames; ++f)
    {
      for (int c = 0; c < nchannels; ++c)
      {
	audio_channels[c][f] = buffer[f*nchannels + c];
      }
    }

    int procesed_frames = processor.process(raw_audio_channels, read_frames, raw_output, out_buffer_length);

    for (int f = 0; f < procesed_frames; ++f)
    {
      output_buffer[f] = output[0][f];
    }

    sf_writef_double(snd_file_out, output_buffer, procesed_frames);
  }

}



int main(int argc, char *argv[])
{

  int c;
  const char shortopts[] = "i:o:a:d:gh";
  extern char *optarg;
  std::string input_file, output_file, array_description_file, doa_file;
  bool use_gui = false;

  while ( (c = getopt(argc, argv, shortopts)) != -1)
  {
    switch (c) {
      case 'i':
	input_file = optarg;
      break;
      case 'o':
	output_file = optarg;
      break;
      case 'a':
	array_description_file = optarg;
      break;
      case 'd':
	doa_file = optarg;
      break;
      case 'g':
	use_gui = true;
      break;
      default:
	std::cerr << "Unknon option '" << reinterpret_cast<const char*>(&c) << "'" << std::endl;
      case 'h':
	usage(argc, argv);
	exit(1);
      break;
    }
  }

  if (input_file.empty() || output_file.empty())
  {
    usage(argc, argv);
    exit(1);
  }

  std::ostream *os = NULL;
  std::unique_ptr<std::ofstream> fos;
  if (!doa_file.empty())
  {
    if(doa_file == "-")
      os = &(std::cout);
    else
    {
      fos.reset(new std::ofstream(doa_file));
      os = fos.get();
    }
  }

  //  std::vector<double> mics_positions = {0, 1, 3.5, 4.5};
  std::vector<double> mics_positions = {-2.25, -1.25, 1.25, 2.25};

  mca::ArrayDescription array_description = mca::ArrayDescription::make_linear_array_description(mics_positions);

  SF_INFO sf_info_in, sf_info_out;
  SNDFILE *snd_file_in  = sf_open(input_file.c_str(), SFM_READ, &sf_info_in);
  sf_info_out = sf_info_in;
  sf_info_out.channels = 1;
  SNDFILE *snd_file_out = sf_open(output_file.c_str(), SFM_WRITE, &sf_info_out);

  int sampleRate = sf_info_in.samplerate;
  McBeamCallback callback(*os);
  mca::SourceSeparationAndLocalisation sss(sampleRate, array_description, 1, false);


#ifdef DSPONE_GUI
  dsp::DspGui::thread_processing_run_t f_process = std::bind(process_file, _1, snd_file_in, snd_file_out, sf_info_in);
  dsp::DspGui gui(sss, f_process);
#else
  auto f_process = std::bind(process_file, _1, snd_file_in, snd_file_out, sf_info_in);
#endif

  if (os) sss.setCallback(callback);

  if (snd_file_in != NULL && snd_file_out != NULL)
  {

    dsp::ShortTimeProcess *p = reinterpret_cast<dsp::ShortTimeProcess*>(&sss);
#ifdef DSPONE_GUI
    gui.start(0);
#else
      (f_process)(sss);
#endif

    sf_close(snd_file_out);
    sf_close(snd_file_in);

  }
  else
  {
    const char *error_message = sf_strerror(snd_file_in);
    std::cout << "error" << std::endl;
    ERROR_STREAM("While reading file " << input_file << " (" << error_message << ")");
  }

}
