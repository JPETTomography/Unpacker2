#include "Unpacker2D.h"
#include "EventIII.h"
#include "TDCChannel.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

UInt_t ReverseHexDJ(UInt_t n)
{
  UInt_t a, b, c, d, e;
  a = n & 0x000000ff;
  b = n & 0x0000ff00;
  c = n & 0x00ff0000;
  d = n & 0xff000000;

  a <<= 8;
  b >>= 8;
  c <<= 8;
  d >>= 8;

  e = a | b | c | d;

  return e;
}

UInt_t ReverseHexTDC(UInt_t n)
{
  UInt_t a, b, c, d, e;
  a = n & 0x000000ff;
  b = n & 0x0000ff00;
  c = n & 0x00ff0000;
  d = n & 0xff000000;

  a <<= 24;
  b <<= 8;
  c >>= 8;
  d >>= 24;

  e = a | b | c | d;

  return e;
}

void Unpacker2D::BuildEvent(EventIII* e, std::map<UInt_t, std::vector<UInt_t>>* m, std::map<UInt_t, UInt_t>* refTimes)
{
  UInt_t data;
  UInt_t fine;
  UInt_t coarse;
  UInt_t rising;
  double refTime = 0;
  double time = 0;
  double time_t = 0;

  map<UInt_t, double> refTimes_local;

  // copy reference channels into the tree
  map<UInt_t, UInt_t>::iterator ref_it;
  for (ref_it = refTimes->begin(); ref_it != refTimes->end(); ref_it++)
  {
    TDCChannel* tc = e->AddTDCChannel(kReferenceChannel + ref_it->first);
    data = ref_it->second;
    fine = data & 0x7f;
    coarse = (data >> 7) & 0xffff;

    if (useTDCcorrection)
    {
      double correction = 0.0;
      auto search = fTDCCorrections.find(ref_it->first);
      if (search != fTDCCorrections.end())
      {
        correction = fTDCCorrections.at(ref_it->first)->GetBinContent(fine + 1);
      }
      time = (coarse * kCoarseConstant) - correction;
    }
    else
    {
      time = (coarse * kCoarseConstant) - ((fine - 4) * kFineConstant);
    }

    tc->AddLead(time);
    refTimes_local[ref_it->first] = time;
  }

  map<UInt_t, vector<UInt_t>>::iterator m_it;
  for (m_it = m->begin(); m_it != m->end(); m_it++)
  {
    TDCChannel* tc = e->AddTDCChannel(m_it->first);

    for (int i = 0; i < m_it->second.size(); i++)
    {
      refTime = refTimes_local.find((int)((m_it->first - kModularOffset) / kNumberOfChannels) * kNumberOfChannels + kModularOffset)->second;
      data = m_it->second[i];
      fine = data & 0x7f;
      coarse = (data >> 7) & 0xffff;
      rising = (data >> 31);

      if (useTDCcorrection)
      {
        double correction = 0.0;
        auto search = fTDCCorrections.find(m_it->first);
        if (search != fTDCCorrections.end())
        {
          correction = fTDCCorrections.at(m_it->first)->GetBinContent(fine + 1);
        }
        time = (coarse * kCoarseConstant) - correction;
      }
      else
      {
        time = coarse * kCoarseConstant - (fine - 4) * kFineConstant;
      }

      if (rising == 0)
      {
        if (time - refTime < 0)
        {
          time_t = time + ((0x10000 * kCoarseConstant) - refTime);
        }
        else
        {
          time_t = time - refTime;
        }
        tc->AddLead(time_t);
      }
      else
      {
        if (time - refTime < 0)
        {
          time_t = time + ((0x10000 * kCoarseConstant) - refTime);
        }
        else
        {
          time_t = time - refTime;
        }
        tc->AddTrail(time_t);
      }
    }
  }
}

void Unpacker2D::ParseConfigFile()
{
  // parsing xml config file
  boost::property_tree::ptree tree;

  try
  {
    boost::property_tree::read_xml(fConfigFile, tree);
  }
  catch (boost::property_tree::xml_parser_error e)
  {
    cerr << "ERROR: Failed to read config file" << endl;
    exit(0);
  }

  // get the config options from the config file
  try
  {
    if (tree.get<string>("READOUT.DEBUG") == "ON")
    {
      debugMode = true;
    }
  }
  catch (exception e)
  {
    cerr << "ERROR: Incorrect config file structure" << endl;
    exit(0);
  }

  if (debugMode)
    cerr << "DEBUG mode on" << endl;

  // get the first data source entry in the config file
  boost::property_tree::ptree readoutTree = tree.get_child("READOUT");
  string type;
  string address_s;
  UInt_t address;
  string hubAddress;
  string correctionFile;
  int channels = 0;
  int offset = 0;
  string measurementType("");
  highest_channel_number = kModularOffset + 2 * kNumberOfModules * kNumberOfChannels;

  // iterate through entries and create appropriate unpackers
  for (const auto& readoutEntry : readoutTree)
  {
    // read out values from xml entry
    if ((readoutEntry.first) == "DATA_SOURCE")
    {
      type = (readoutEntry.second).get<string>("TYPE");
      address_s = (readoutEntry.second).get<string>("TRBNET_ADDRESS");
      hubAddress = (readoutEntry.second).get<string>("HUB_ADDRESS");
      correctionFile = (readoutEntry.second).get<string>("CORRECTION_FILE");

      if (correctionFile.compare("raw") != 0)
      {
        cerr << "WARNING: The TDC correction file path was set in the XML config file of the Unpacker!" << endl;
        cerr << "This file path should be defined in the user parameters JSON file instead." << endl;
        cerr << "The setting from the XML file fill be ignored!" << endl;
      }

      if (type == "TRB3_S" || type == "DJPET_ENDP")
      {
        // create additional unpackers for internal modules
        boost::property_tree::ptree modulesTree = (readoutEntry.second).get_child("MODULES");
        for (const auto& module : modulesTree)
        {
          type = (module.second).get<string>("TYPE");
          address_s = (module.second).get<string>("TRBNET_ADDRESS");
          address = std::stoul(address_s, 0, 16);
          channels = (module.second).get<int>("NUMBER_OF_CHANNELS");
          offset = (module.second).get<int>("CHANNEL_OFFSET");
          measurementType = (module.second).get<string>("MEASUREMENT_TYPE");
          tdc_offsets[address] = offset;
        }
      }
      else
      {
        cerr << "Incorrect configuration in the xml file!" << endl;
        cerr << "The DATA_SOURCE entry is missing!" << endl;
      }
    }
  }
}

void Unpacker2D::UnpackSingleStep(string inputFile, string inputPath, string outputPath, string configFile, int numberOfEvents, int refChannelOffset,
                                  std::string tdcCalibFile)
{
  fInputFile = inputFile;
  fInputFilePath = inputPath;
  fOutputFilePath = outputPath;
  fConfigFile = configFile;
  fTDCCalibFile = tdcCalibFile;
  fEventsToAnalyze = numberOfEvents;
  fRefChannelOffset = refChannelOffset;
  ParseConfigFile();

  if (!fTDCCalibFile.empty())
  {
    useTDCcorrection = loadTDCcalibFile(fTDCCalibFile.c_str());
  }

  DistributeEventsSingleStep();
}

bool Unpacker2D::loadTDCcalibFile(const char* calibFileName)
{
  TFile* file = new TFile(calibFileName, "READ");
  TDirectory* dir = gDirectory->GetDirectory("Rint:/");

  if (file->IsOpen())
  {
    for (int i = kModularOffset; i < highest_channel_number; i++)
    {
      TH1F* tmp = dynamic_cast<TH1F*>(file->Get(Form("correction%d", i)));
      if (tmp)
      {
        fTDCCorrections[i] = dynamic_cast<TH1F*>(tmp->Clone(tmp->GetName()));
        fTDCCorrections[i]->SetDirectory(dir);
      }
    }

    if (debugMode)
    {
      cerr << "Loaded TDC nonlinearity corrections." << endl;
    }
  }
  else
  {
    if (debugMode)
    {
      cerr << "The TDC calibration file " << calibFileName << " could not be properly opened." << endl;
      cerr << "TDC nonlinearity correction will not be used!" << endl;
    }
    file->Close();
    delete file;
    return false;
  }
  file->Close();
  delete file;
  return true;
}

void Unpacker2D::DistributeEventsSingleStep()
{
  string fileName = fInputFilePath + fInputFile;
  ifstream* file = new ifstream(fileName.c_str());

  if (file->is_open())
  {
    EventIII* eventIII = new EventIII();

    // open a new file
    string newFileName = fOutputFilePath + fInputFile + ".root";

    TFile* newFile = new TFile(newFileName.c_str(), "RECREATE");
    TTree* newTree = new TTree("T", "Tree");
    newTree->Branch("eventIII", "EventIII", &eventIII, 64000, 99);

    UInt_t data4;
    UInt_t nEvents = 0;
    UInt_t eventSize = 0;
    UInt_t queueSize = 0;
    UInt_t subDecoding = 0;
    UInt_t trgId = 0;
    UInt_t concId = 0;
    UInt_t ftabId = 0;
    UInt_t ftabSize = 0;
    UInt_t ftabTrgn = 0;
    UInt_t ftabDbg = 0;
    UInt_t dataCtr = 0;
    UInt_t channel = 0;
    UInt_t currentOffset = 0;
    UInt_t ftabWords = 0;
    UInt_t prevTrgId = 0;
    UInt_t prevSTrgId = 0;

    map<UInt_t, UInt_t>::iterator offsets_it;
    map<UInt_t, vector<UInt_t>> tdc_channels;
    map<UInt_t, UInt_t> refTimes;
    map<UInt_t, UInt_t> refTimes_previous;

    bool trgSequence = false;
    bool trgSequenceHold = false;
    bool missing_ref = false;
    bool badIdSkip = false;

    int built = 0;
    int skip = 0;
    int missingRefCtr = 0;
    int badIdCtr = 0;

    short endpointsCtr = 0;
    short refTimesCtr = 0;

    // skip the first entry (EB file headers)
    file->ignore(32);

    while (!file->eof())
    {
      if (nEvents % 10000 == 0)
      {
        printf("%d\n", nEvents);
      }
      eventIII->Clear();

      // EVENT HEADERS FROM EB
      // size
      file->read((char*)&data4, 4);
      eventSize = data4 / 4; // in 32b words

      nEvents++;

      // skip bad entries
      if (eventSize < 0x10)
      {
        file->ignore(28);
        newTree->Fill();
        continue;
      }

      // skip event headers
      file->ignore(7 * 4);

      endpointsCtr = 0;
      refTimesCtr = 0;

      refTimes.clear();
      
      //*****************
      //  LOOP OVER CONCENTRATORS IN THE QUEUE
      //*****************
      while (!file->eof())
      {
        // QUEUE FROM CONCENTRATOR SIZE
        file->read((char*)&data4, 4);
        data4 = ReverseHex(data4);
        queueSize = data4 / 4; // queue size in 32 bit words
        eventSize -= queueSize;
        file->ignore(4); // queue decoding

        // CONCENTRATOR ID
        file->read((char*)&data4, 4);
        data4 = ReverseHex(data4);
        concId = data4;
        queueSize--;

        // TRIGGER NUMBER FROM CONCENTRATOR
        file->read((char*)&data4, 4);
        data4 = ReverseHex(data4);
        trgId = data4;
        queueSize--;

        //*****************
        // LOOP OVER FTABS FROM CONCENTRATOR
        //*****************
        while (!file->eof())
        {

          // FTAB HEADERS ID AND SIZE
          file->read((char*)&data4, 4);
          data4 = ReverseHexDJ(data4);
          ftabSize = data4 >> 16;
          ftabId = data4 & 0xffff;
          queueSize--;
          endpointsCtr += 1;
          ftabWords = ftabSize - 2;

          // FTAB TRG NR
          file->read((char*)&data4, 4);
          data4 = ReverseHexDJ(data4);
          ftabTrgn = data4 & 0xffff;
          queueSize--;
          offsets_it = tdc_offsets.find(ftabId);
          badIdSkip = false;

          // MISSING FTAB ID
          if (offsets_it == tdc_offsets.end())
          {
            badIdSkip = true;
            badIdCtr++;
          }
          currentOffset = offsets_it->second;

          if (nEvents > 0 && ftabId != 0xffff)
          {
            refTimes_previous[currentOffset] = refTimes[currentOffset];
          }

          bool gotRef = false;

          //*****************
          // LOOP OVER DATA FROM FTAB
          //*****************
          while (!file->eof())
          {
            file->read((char*)&data4, 4);
            data4 = ReverseHexTDC(data4);
            queueSize--;
            ftabWords--;

            if (ftabWords == 1)
            {
              file->read((char*)&data4, 4);
              file->read((char*)&data4, 4);

              if (data4 == 0x05050505)
              {
                // padding found
                queueSize -= 2;
                break;
              }
              else
              {
                // regular ending
                file->seekg(file->tellg() - 4);
                queueSize--;
                break;
              }
            }

            // FTAB TDC HIT
            if ((data4 >> 24) != 0xfc && ftabId != 0xffff)
            {
              channel = (data4 >> 24) & 0x7f;
              if (channel == kReferenceChannel)
              {
                // reference channel
                refTimes[currentOffset] = data4;
                gotRef = true;
                refTimesCtr++;
              }
              else
              {
                // regular channel
                tdc_channels[channel + currentOffset].push_back(data4);
              }
            }
          }

          if (!gotRef)
          {
            missingRefCtr++;
          }

          // conditions to exit the current concentrator
          if (queueSize <= 4)
          {
            break;
          }
        } // loop over ftabs in queue

        if (eventSize <= 8)
        {
          break;
        }

      } // end of loop over concentrator in queue

      if (!badIdSkip && endpointsCtr == refTimesCtr)
      {
        // copy all times into the tree
        BuildEvent(eventIII, &tdc_channels, &refTimes);
        newTree->Fill();
      }

      tdc_channels.clear();

      if (nEvents == fEventsToAnalyze)
      {
        printf("Max timeslots reached\n");
        break;
      }

    } // end of loop over events

    newFile->Write();
    delete newTree;
    file->close();
  }
  else
  {
    cerr << "ERROR: failed to open data file" << endl;
  }
}
