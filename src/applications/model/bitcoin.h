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
  FILTER_REQUEST,
  MODE,
  BLOCK,
  UPDATE_FILTER_BEGIN,
  UPDATE_FILTER_END,
};

enum ProtocolType
{
  STANDARD_PROTOCOL,           //DEFAULT
  FILTERS_ON_INCOMING_LINKS,
  OUTGOING_FILTERS,
  PREFERRED_DESTINATIONS,
  DANDELION_LIKE
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
  long     txCreated;
  int      connections;

  int      blocksRelayed;
  double firstSpySuccess;

  int txReceived;
  int systemId;

  std::vector<txRecvTime> txReceivedTimes;
  int ignoredFilters;
} nodeStatistics;

typedef struct {
  long numUsefulInvReceived;
  long numUselessInvReceived;
  long numGetDataReceived;
  long numGetDataSent;
  double connectionLength;
  double usefulInvRate;
} peerStatistics;

typedef struct {
  double downloadSpeed;
  double uploadSpeed;
} nodeInternetSpeeds;

#define FILTER_BASE_NUMBERING 1000

}// Namespace ns3

#endif /* BITCOIN_H */
