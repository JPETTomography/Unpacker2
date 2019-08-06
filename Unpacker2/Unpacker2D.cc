#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "EventIII.h"
#include "TDCChannel.h"
#include "Unpacker2D.h"

using namespace std;

UInt_t ReverseHexDJ(UInt_t n) {
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

UInt_t ReverseHexTDC(UInt_t n) {
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

void BuildEvent(EventIII *e, map<UInt_t, vector<UInt_t>> *m) {

  UInt_t data;
  UInt_t fine;
  UInt_t coarse;
  UInt_t rising;
  double time = 0;

  map<UInt_t, vector<UInt_t>>::iterator m_it;
  for (m_it = m->begin(); m_it != m->end(); m_it++) {
    TDCChannel *tc = e->AddTDCChannel(m_it->first);
    for (int i = 0; i < m_it->second.size(); i++) {

      data = m_it->second[i];

      fine = data & 0xff;
      coarse = (data >> 8) & 0xffff;
      rising = (data >> 31);

      time = coarse * 3.3 + fine * 0.025;

      if (rising == 1) {
        tc->AddLead(time);
      } else {
        tc->AddTrail(time);
      }
    }
  }

  m->clear();
}

// parsing xml config file
void Unpacker2D::ParseConfigFile() {
  boost::property_tree::ptree tree;

  try {
    boost::property_tree::read_xml(fConfigFile, tree);
  } catch (boost::property_tree::xml_parser_error e) {
    cerr << "ERROR: Failed to read config file" << endl;
    exit(0);
  }

  // get the config options from the config file
  try {
    if (tree.get<string>("READOUT.DEBUG") == "ON") {
      debugMode = true;
    }
  } catch (exception e) {
    cerr << "ERROR: Incorrect config file structure" << endl;
    exit(0);
  }

  if (debugMode == true)
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
  // int resolution = 0;
  // int referenceChannel = 0;
  string measurementType("");
  highest_channel_number = 0;

  // iterate through entries and create appropriate unpackers
  for (const auto &readoutEntry : readoutTree) {
    // read out values from xml entry
    if ((readoutEntry.first) == "DATA_SOURCE") {
      type = (readoutEntry.second).get<string>("TYPE");
      address_s = (readoutEntry.second).get<string>("TRBNET_ADDRESS");
      hubAddress = (readoutEntry.second).get<string>("HUB_ADDRESS");
      // referenceChannel = (readoutEntry.second).get<int>("REFERENCE_CHANNEL");
      correctionFile = (readoutEntry.second).get<string>("CORRECTION_FILE");

      if (correctionFile.compare("raw") != 0) {
        cerr << "WARNING: The TDC correction file path was set in the XML "
                "config file of the Unpacker!"
             << endl;
        cerr << "This file path should be defined in the user parameters JSON "
                "file instead."
             << endl;
        cerr << "The setting from the XML file fill be ignored!" << endl;
      }

      if (type == "TRB3_S" || type == "DJPET_ENDP") {

        // create additional unpackers for internal modules
        boost::property_tree::ptree modulesTree =
            (readoutEntry.second).get_child("MODULES");
        for (const auto &module : modulesTree) {
          type = (module.second).get<string>("TYPE");
          address_s = (module.second).get<string>("TRBNET_ADDRESS");
          address = std::stoul(address_s, 0, 16);
          channels = (module.second).get<int>("NUMBER_OF_CHANNELS");
          offset = (module.second).get<int>("CHANNEL_OFFSET");
          // resolution = (module.second).get<int>("RESOLUTION");
          measurementType = (module.second).get<string>("MEASUREMENT_TYPE");

          tdc_offsets[address] = offset;
          if (offset + channels > highest_channel_number) {
            highest_channel_number = offset + channels;
          }
        }
      } else {
        cerr << "Incorrect configuration in the xml file!" << endl;
        cerr << "The DATA_SOURCE entry is missing!" << endl;
      }
    }
  }
}

void Unpacker2D::UnpackSingleStep(
  string inputFile, string inputPath, string outputPath,
  string configFile, int numberOfEvents, int refChannelOffset,
  string totCalibFile, std::string tdcCalibFile)
{
  fInputFile = inputFile;
  fInputFilePath = inputPath;
  fOutputFilePath = outputPath;
  fConfigFile = configFile;
  fTOTCalibFile = totCalibFile;
  fTDCCalibFile = tdcCalibFile;

  fEventsToAnalyze = numberOfEvents;
  fRefChannelOffset = refChannelOffset;
  // tdc_offsets[0xa110] = 0;

  ParseConfigFile();

  DistributeEventsSingleStep();
}

void Unpacker2D::DistributeEventsSingleStep() {
  string fileName = fInputFilePath+fInputFile;
  ifstream *file = new ifstream(fileName.c_str());

  if (file->is_open()) {

    EventIII *eventIII = new EventIII();

    string newFileName = fOutputFilePath + fInputFile + ".root";
    TFile *newFile = new TFile(newFileName.c_str(), "RECREATE");
    TTree *newTree = new TTree("T", "Tree");

    newTree->Branch("eventIII", "EventIII", &eventIII, 64000, 99);

    UInt_t data4;
    Int_t nBytes;
    UInt_t nEvents = 0;

    UInt_t queueSize = 0;
    UInt_t queueDecoding = 0;
    UInt_t subSize = 0;
    UInt_t subDecoding = 0;
    UInt_t ftabId = 0;
    UInt_t ftabSize = 0;
    UInt_t ftabTrgn = 0;
    UInt_t ftabDbg = 0;
    UInt_t dataCtr = 0;
    UInt_t channel = 0;
    UInt_t currentOffset = 0;

    map<UInt_t, UInt_t>::iterator offsets_it;
    map<UInt_t, vector<UInt_t>> tdc_channels;

    // skip the first entry
    file->ignore(32);

    while (!file->eof()) {
      nBytes = 0;

      if (nEvents % 10000 == 0) {
        printf("%d\n", nEvents);
      }

      eventIII->Clear();

      // queue headers
      // size
      file->read((char *)&data4, 4);
      nBytes += 4;
      queueSize = data4;

      // skip bad entries
      if (queueSize == 0x20) {
        file->ignore(28);
        nBytes += 28;
        if (debugMode)
          printf("Skipping too small queue\n");
        continue;
      }

      // decoding
      file->read((char *)&data4, 4);
      nBytes += 4;
      queueDecoding = data4;

      // skip some headers
      file->ignore(8);
      nBytes += 8;

      // subevent
      // sub size
      file->read((char *)&data4, 4);
      nBytes += 4;
      subSize = data4 & 0xffff;

      // sub decoding
      file->read((char *)&data4, 4);
      nBytes += 4;
      subDecoding = data4;

      // skip some headers
      file->ignore(24);
      nBytes += 24;

      // ftab header
      // ftab size and id
      file->read((char *)&data4, 4);
      nBytes += 4;
      data4 = ReverseHexDJ(data4);
      ftabSize = data4 >> 16;
      ftabId = data4 & 0xffff;

      // ftab trigger number and debug
      file->read((char *)&data4, 4);
      nBytes += 4;
      data4 = ReverseHexDJ(data4);
      ftabDbg = data4 >> 16;
      ftabTrgn = data4 & 0xffff;

      offsets_it = tdc_offsets.find(ftabId);
      if (offsets_it == tdc_offsets.end()) {
        printf("Wrong ftab ID\n");
        break;
      }
      currentOffset = offsets_it->second;

      if (nEvents == fEventsToAnalyze) {
        printf("Max timeslots reached\n");
        break;
      }

      // ftab data
      while (!file->eof()) {
        file->read((char *)&data4, 4);
        nBytes += 4;

        data4 = ReverseHexTDC(data4);

        if ((data4 >> 24) != 0xfc) {
          channel = (data4 >> 24) & 0x7f;

          tdc_channels[channel + currentOffset].push_back(data4);
        }

        if (ftabSize % 2 == 0) {
          if (dataCtr == ftabSize - 5) {
            dataCtr = 0;
            // skip trailing
            file->ignore(8);
            nBytes += 8;

            BuildEvent(eventIII, &tdc_channels);

            newTree->Fill();
            break;
          }
        } else {
          if (dataCtr == ftabSize - 5) {
            dataCtr = 0;
            // skip trailing
            file->ignore(12);
            nBytes += 12;

            BuildEvent(eventIII, &tdc_channels);

            newTree->Fill();
            break;
          }
        }

        dataCtr++;
      }

      nEvents++;
    }

    newFile->Write();
    delete newTree;
    file->close();
  } else {
    cerr << "ERROR: failed to open data file" << endl;
  }
}
