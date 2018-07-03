/**
 * This file contains all the necessary enumerations and structs used throughout the simulation.
 * It also defines 3 very important classed; the Block, Chunk and Blockchain.
 */


#ifndef BITCOIN_H
#define BITCOIN_H

#include <vector>
#include <map>
#include "ns3/address.h"
#include <algorithm>

namespace ns3 {

/**
 * The bitcoin message types that have been implemented.
 */
enum Messages
{
  INV,              //0
  GET_DATA,         //1
  TX,
  FILTER,
  MODE,
  BLOCK
};

enum ProtocolType
{
  STANDARD_PROTOCOL,           //DEFAULT
  FILTERS_ON_LINKS
};

enum ModeType
{
  REGULAR,           //DEFAULT
  BLOCKS_ONLY,
  SPY
};

typedef struct {
  int nodeId;
  int txHash;
  int txTime;
} txRecvTime;

/**
 * The struct used for collecting node statistics.
 */
typedef struct {
  int      nodeId;
  long     invReceivedBytes;
  long     invSentBytes;
  long     invReceivedMessages;
  long     invSentMessages;
  long     getDataReceivedBytes;
  long     getDataSentBytes;
  long     getDataReceivedMessages;
  long     getDataSentMessages;
  long txCreated;
  int      connections;

  int      blocksRelayed;
  double firstSpySuccess;

  int txReceived;
  int systemId;

  std::vector<txRecvTime> txReceivedTimes;
  std::vector<uint64_t> invSentTimes;

  int ignoredFilters;
} nodeStatistics;



typedef struct {
  double downloadSpeed;
  double uploadSpeed;
} nodeInternetSpeeds;

}// Namespace ns3

#endif /* BITCOIN_H */
