/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <fstream>
#include <time.h>
#include <sys/time.h>
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/point-to-point-layout-module.h"
#include "ns3/mpi-interface.h"
#include <climits>

#define MPI_TEST

#ifdef NS3_MPI
#include <mpi.h>
#endif

using namespace ns3;

double get_wall_time();
int GetNodeIdByIpv4 (Ipv4InterfaceContainer container, Ipv4Address addr);
void PrintStatsForEachNode (nodeStatistics *stats, int totalNodes, int publicIPNodes, int blocksOnlyPrivateIpNodes, int txToCreate);
void PrintBitcoinRegionStats (uint32_t *bitcoinNodesRegions, uint32_t totalNodes);
void CollectTxData(nodeStatistics *stats, int totalNoNodes, int txToCreate,
   int systemId, int systemCount, int nodesInSystemId0, BitcoinTopologyHelper bitcoinTopologyHelper);
int PoissonDistribution(int value);
std::vector<int> generateTxCreateList(int n, int nodes);


NS_LOG_COMPONENT_DEFINE ("MyMpiTest");

int
main (int argc, char *argv[])
{
  std::cout << "Start \n";

  bool nullmsg = false;
  bool testScalability = false;
  int invTimeoutMins = -1;
  double tStart = get_wall_time(), tStartSimulation, tFinish;
  const int secsPerMin = 60;
  const uint16_t bitcoinPort = 8333;
  int start = 0;
//
  int totalNoNodes = 16;
  int minConnectionsPerNode = -1;
  int maxConnectionsPerNode = -1;

  uint32_t protocol;


  uint64_t txToCreate = 1024;
  int publicIPNodes, blocksOnlyPrivateIpNodes;

  double stop;


//
  Ipv4InterfaceContainer                               ipv4InterfaceContainer;
  std::map<uint32_t, std::vector<Ipv4Address>>         nodesConnections;
  std::map<uint32_t, std::map<Ipv4Address, double>>    peersDownloadSpeeds;
  std::map<uint32_t, std::map<Ipv4Address, double>>    peersUploadSpeeds;
  std::map<uint32_t, nodeInternetSpeeds>               nodesInternetSpeeds;
  int                                                  nodesInSystemId0 = 0;

  int netGroups = 0;

  Time::SetResolution (Time::NS);

  int r = 1;

  CommandLine cmd;
  cmd.AddValue ("nullmsg", "Enable the use of null-message synchronization", nullmsg);
  cmd.AddValue ("nodes", "The total number of nodes in the network", totalNoNodes);
  cmd.AddValue ("minConnections", "The minConnectionsPerNode of the grid", minConnectionsPerNode);
  cmd.AddValue ("maxConnections", "The maxConnectionsPerNode of the grid", maxConnectionsPerNode);
  cmd.AddValue ("invTimeoutMins", "The inv block timeout", invTimeoutMins);
  cmd.AddValue ("test", "Test the scalability of the simulation", testScalability);

  cmd.AddValue ("txToCreate", "The number of transactions each the network should generate", txToCreate);

  cmd.AddValue ("publicIPNodes", "How many nodes has public IP", publicIPNodes);

  cmd.AddValue ("blocksOnlyPrivateIPNodes", "How many nodes with private IP run blocksOnly", blocksOnlyPrivateIpNodes);

  cmd.AddValue ("protocol", "Used protocol: 0 — Default, 1 — Filters on links", protocol);
  cmd.AddValue ("netGroups", "How many groups each node has", netGroups);
  cmd.AddValue ("r", "incoming_inv_interval/outgoing_inv_interval", r);


  cmd.Parse(argc, argv);

  assert(netGroups > 0);

  // TODO Configure
  uint averageBlockGenInterval = 10 * 60;
  uint targetNumberOfBlocks = 5000;

  stop = targetNumberOfBlocks * averageBlockGenInterval / 60; // minutes
  nodeStatistics *stats = new nodeStatistics[totalNoNodes];

  #ifdef MPI_TEST
    // Distributed simulation setup; by default use granted time window algorithm.
    if(nullmsg)
      {
        GlobalValue::Bind ("SimulatorImplementationType",
                           StringValue ("ns3::NullMessageSimulatorImpl"));
      }
    else
      {
        GlobalValue::Bind ("SimulatorImplementationType",
                           StringValue ("ns3::DistributedSimulatorImpl"));
      }

    // Enable parallel simulator with the command line arguments
    MpiInterface::Enable (&argc, &argv);
    uint32_t systemId = MpiInterface::GetSystemId ();
    uint32_t systemCount = MpiInterface::GetSize ();
  #else
    uint32_t systemId = 0;
    uint32_t systemCount = 1;
  #endif


  LogComponentEnable("BitcoinNode", LOG_LEVEL_INFO);

  BitcoinTopologyHelper bitcoinTopologyHelper (systemCount, totalNoNodes, publicIPNodes, minConnectionsPerNode,
                                               maxConnectionsPerNode, systemId);

  // Install stack on Grid
  InternetStackHelper stack;
  bitcoinTopologyHelper.InstallStack (stack);


  // Assign Addresses to Grid
  bitcoinTopologyHelper.AssignIpv4Addresses (Ipv4AddressHelperCustom ("1.0.0.0", "255.255.255.0", false));
  ipv4InterfaceContainer = bitcoinTopologyHelper.GetIpv4InterfaceContainer();
  nodesConnections = bitcoinTopologyHelper.GetNodesConnectionsIps();
  peersDownloadSpeeds = bitcoinTopologyHelper.GetPeersDownloadSpeeds();
  peersUploadSpeeds = bitcoinTopologyHelper.GetPeersUploadSpeeds();
  nodesInternetSpeeds = bitcoinTopologyHelper.GetNodesInternetSpeeds();


  std::cout << "Total nodes: " << totalNoNodes << "\n";
  //Install simple nodes
  BitcoinNodeHelper bitcoinNodeHelper ("ns3::TcpSocketFactory", InetSocketAddress (Ipv4Address::GetAny (), bitcoinPort),
                                        nodesConnections[0], peersDownloadSpeeds[0],  peersUploadSpeeds[0], nodesInternetSpeeds[0], stats, r);
  ApplicationContainer bitcoinNodes;


  int startedblocksOnlyPrivateIpNodes;

  auto txCreateList = generateTxCreateList(txToCreate, totalNoNodes);

  std::map<int, int> nodeSystemIds;

  for(auto &node : nodesConnections)
  {
    Ptr<Node> targetNode = bitcoinTopologyHelper.GetNode (node.first);

  	if (systemId == targetNode->GetSystemId())
  	{
      bitcoinNodeHelper.SetPeersAddresses (node.second);
      bitcoinNodeHelper.SetPeersDownloadSpeeds (peersDownloadSpeeds[node.first]);
      bitcoinNodeHelper.SetPeersUploadSpeeds (peersUploadSpeeds[node.first]);
      bitcoinNodeHelper.SetNodeInternetSpeeds (nodesInternetSpeeds[node.first]);

      // bitcoinNodeHelper.SetProperties(1, ProtocolType(protocol), REGULAR, netGroups, systemId);
      bitcoinNodeHelper.SetProperties(txCreateList[targetNode->GetId()], ProtocolType(protocol), REGULAR, netGroups, systemId);
  	  bitcoinNodeHelper.SetNodeStats (&stats[node.first]);
      bitcoinNodes.Add(bitcoinNodeHelper.Install (targetNode));

      if (systemId == 0)
        nodesInSystemId0++;
  	}
  }

  bitcoinNodes.Start (Seconds (start));
  bitcoinNodes.Stop (Minutes (stop));

  tStartSimulation = get_wall_time();


  if (systemId == 0) {
    std::cout << "start: " << start << "\n";
    std::cout << "stop: " << stop << "\n";
    std::cout << "The applications have been setup.\n";
    std::cout << "Setup time = " << tStartSimulation - tStart << "s\n";
    std::cout << "Total nodes: " << totalNoNodes << "\n";
  }

  Simulator::Stop (Minutes (stop + 0.1));
  Simulator::Run ();
  Simulator::Destroy ();

  #ifdef MPI_TEST

    int            blocklen[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                   1, 1, 1};
    MPI_Aint       disp[15];
    MPI_Datatype   dtypes[15] = {MPI_INT, MPI_LONG, MPI_LONG, MPI_LONG, MPI_LONG, MPI_LONG, MPI_LONG, MPI_LONG, MPI_LONG, MPI_LONG, MPI_INT, MPI_INT,
                                 MPI_DOUBLE, MPI_INT, MPI_INT};
    MPI_Datatype   mpi_nodeStatisticsType;

    disp[0] = offsetof(nodeStatistics, nodeId);
    disp[1] = offsetof(nodeStatistics, invReceivedBytes);
    disp[2] = offsetof(nodeStatistics, invSentBytes);
    disp[3] = offsetof(nodeStatistics, invReceivedMessages);
    disp[4] = offsetof(nodeStatistics, invSentMessages);
    disp[5] = offsetof(nodeStatistics, getDataReceivedBytes);
    disp[6] = offsetof(nodeStatistics, getDataSentBytes);
    disp[7] = offsetof(nodeStatistics, getDataReceivedMessages);
    disp[8] = offsetof(nodeStatistics, getDataSentMessages);
    disp[9] = offsetof(nodeStatistics, txCreated);
    disp[10] = offsetof(nodeStatistics, connections);
    disp[11] = offsetof(nodeStatistics, blocksRelayed);
    disp[12] = offsetof(nodeStatistics, firstSpySuccess);
    disp[13] = offsetof(nodeStatistics, txReceived);
    disp[14] = offsetof(nodeStatistics, systemId);

    MPI_Type_create_struct (15, blocklen, disp, dtypes, &mpi_nodeStatisticsType);
    MPI_Type_commit (&mpi_nodeStatisticsType);

    if (systemId != 0 && systemCount > 1)
    {
      for(int i = 0; i < totalNoNodes; i++)
      {
        Ptr<Node> targetNode = bitcoinTopologyHelper.GetNode (i);

    	  if (systemId == targetNode->GetSystemId())
    	  {
            MPI_Send(&stats[i], 1, mpi_nodeStatisticsType, 0, 8888, MPI_COMM_WORLD);
    	  }
      }
    }
    else if (systemId == 0 && systemCount > 1)
    {
      int count = nodesInSystemId0;

    	while (count < totalNoNodes)
    	{
    	  MPI_Status status;
        nodeStatistics recv;

    	  MPI_Recv(&recv, 1, mpi_nodeStatisticsType, MPI_ANY_SOURCE, 8888, MPI_COMM_WORLD, &status);
        stats[recv.nodeId].nodeId = recv.nodeId;
        stats[recv.nodeId].connections = recv.connections;
        stats[recv.nodeId].txCreated = recv.txCreated;
        stats[recv.nodeId].invSentMessages = recv.invSentMessages;
        stats[recv.nodeId].getDataReceivedMessages = recv.getDataReceivedMessages;
        stats[recv.nodeId].firstSpySuccess = recv.firstSpySuccess;
        stats[recv.nodeId].txReceived = recv.txReceived;
        stats[recv.nodeId].systemId = recv.systemId;
  	    count++;
      }
    }

    CollectTxData(stats, totalNoNodes, txToCreate, systemId, systemCount, nodesInSystemId0, bitcoinTopologyHelper);

  #endif


  if (systemId == 0)
  {
    tFinish=get_wall_time();

    PrintStatsForEachNode(stats, totalNoNodes, publicIPNodes, blocksOnlyPrivateIpNodes, txToCreate);


    std::cout << "\nThe simulation ran for " << tFinish - tStart << "s simulating "
              << stop << "mins. Performed " << stop * secsPerMin / (tFinish - tStart)
              << " faster than realtime.\n" << "Setup time = " << tStartSimulation - tStart << "s\n"
              <<"It consisted of " << totalNoNodes << " nodes ( with minConnectionsPerNode = "
              << minConnectionsPerNode << " and maxConnectionsPerNode = " << maxConnectionsPerNode
              << "\n" << "Protocol Type: " << protocol << "\n";

  }

  #ifdef MPI_TEST

    // Exit the MPI execution environment
    MpiInterface::Disable ();


  #else
     NS_FATAL_ERROR ("Can't use distributed simulator without MPI compiled in");
   #endif

  delete[] stats;

  return 0;
//
// #else
//   NS_FATAL_ERROR ("Can't use distributed simulator without MPI compiled in");
// #endif
}

double get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int GetNodeIdByIpv4 (Ipv4InterfaceContainer container, Ipv4Address addr)
{
  for (auto it = container.Begin(); it != container.End(); it++)
  {
	int32_t interface = it->first->GetInterfaceForAddress (addr);
	if ( interface != -1)
      return it->first->GetNetDevice (interface)-> GetNode()->GetId();
  }
  return -1; //if not found
}

void PrintStatsForEachNode (nodeStatistics *stats, int totalNodes, int publicIPNodes, int blocksOnlyPrivateIpNodes, int txToCreate)
{
  double totalUsefulInvSentRatePublicIPNode = 0;
  double totalUsefulInvSentRatePrivateIPNode = 0;
  double totalUsefulInvReceivedRate = 0;
  double totaluselessInvSentMegabytesPublicIPNode = 0;

  std::map<int, std::vector<double>> allTxRelayTimes;

  for (int it = 0; it < totalNodes; it++ )
  {
    std::cout << "\nNode " << stats[it].nodeId << " statistics:\n";
    std::cout << "Connections = " << stats[it].connections << "\n";
    std::cout << "Transactions created = " << stats[it].txCreated << "\n";
    std::cout << "Inv sent = " << stats[it].invSentMessages << "\n";
    // std::cout << "Inv received = " << stats[it].invReceivedMessages << "\n";
    // std::cout << "GetData sent = " << stats[it].getDataSentMessages << "\n";
    std::cout << "GetData received = " << stats[it].getDataReceivedMessages << "\n";
    std::cout << "Tx received = " << stats[it].txReceived << "\n";

    double usefulInvSentRate = 0, usefulInvReceivedRate = 0, invSentMegabytes = 0;
    if (stats[it].invSentMessages != 0)
      usefulInvSentRate = double(stats[it].getDataReceivedMessages) / stats[it].invSentMessages;
    usefulInvReceivedRate = double(stats[it].getDataSentMessages) / stats[it].invReceivedMessages;
    invSentMegabytes = double(stats[it].invSentBytes) / 1024 / 1024;

    std::cout << "Inv sent megabytes = " << invSentMegabytes << "\n";
    std::cout << "Useless inv sent megabytes = " << (1.0-usefulInvSentRate) * invSentMegabytes << "\n";

    std::cout << "Useful inv sent rate = " << usefulInvSentRate << "\n";
    std::cout << "Useful inv received rate = " << usefulInvReceivedRate << "\n";

    if (it < publicIPNodes) {
      totalUsefulInvSentRatePublicIPNode += usefulInvSentRate;
      totaluselessInvSentMegabytesPublicIPNode += (1.0-usefulInvSentRate) * invSentMegabytes;
    }

    if (it >= publicIPNodes)
      totalUsefulInvSentRatePrivateIPNode += usefulInvSentRate;


    totalUsefulInvReceivedRate += usefulInvReceivedRate;

    for (int txCount = 0; txCount < stats[it].txReceived; txCount++)
    {
      txRecvTime txTime = stats[it].txReceivedTimes[txCount];
      allTxRelayTimes[txTime.txHash].push_back(txTime.txTime);
    }

  }

  std::vector<double> fiftyPercentRelayTimes;
  std::vector<double> seventyFivePercentRelayTimes;
  std::vector<double> ninetyNinePercentRelayTimes;
  std::vector<double> fullRelayTimes;

  for (std::map<int, std::vector<double>>::iterator txTimes=allTxRelayTimes.begin();
    txTimes!=allTxRelayTimes.end(); ++txTimes)
  {
    std::vector<double> relayTimes = txTimes->second;
    std::sort(relayTimes.begin(), relayTimes.end());

    int wasRelayedTimes = relayTimes.size();

    if (wasRelayedTimes < (totalNodes - blocksOnlyPrivateIpNodes) * 0.5) {
      continue;
    }

    fiftyPercentRelayTimes.push_back(relayTimes.at(wasRelayedTimes / 2) - relayTimes.front());

    if (wasRelayedTimes < (totalNodes - blocksOnlyPrivateIpNodes) * 0.75) {
      continue;
    }
    seventyFivePercentRelayTimes.push_back(relayTimes.at(wasRelayedTimes * 3 / 4) - relayTimes.front());

    if (wasRelayedTimes < (totalNodes - blocksOnlyPrivateIpNodes) * 0.99) {
      continue;
    }

    ninetyNinePercentRelayTimes.push_back(relayTimes.at(wasRelayedTimes * 99 / 100) - relayTimes.front());

    if (wasRelayedTimes < (totalNodes - blocksOnlyPrivateIpNodes)) {
      continue;
    }

    fullRelayTimes.push_back(relayTimes.back() - relayTimes.front());
  }

  std::cout << "Average 50% relay time: " << accumulate(fiftyPercentRelayTimes.begin(), fiftyPercentRelayTimes.end(), 0.0) / fiftyPercentRelayTimes.size() << "\n";
  std::cout << "Average 75% relay time: " << accumulate(seventyFivePercentRelayTimes.begin(), seventyFivePercentRelayTimes.end(), 0.0) / seventyFivePercentRelayTimes.size() << "\n";
  std::cout << "Average 99% relay time: " << accumulate(ninetyNinePercentRelayTimes.begin(), ninetyNinePercentRelayTimes.end(), 0.0) / ninetyNinePercentRelayTimes.size() << "\n";
  std::cout << "Average 100% relay time: " << accumulate(fullRelayTimes.begin(), fullRelayTimes.end(), 0.0) / fullRelayTimes.size() << "\n";
  std::cout << "Generated transactions: " << allTxRelayTimes.size() << "\n";
  std::cout << "Fully relayed transactions: " << fullRelayTimes.size() << "\n";


  // std::cout << "Average useful inv sent rate (private IP nodes) = " << totalUsefulInvSentRatePrivateIPNode / (totalNodes - publicIPNodes) << "\n";
  // std::cout << "Average useful inv received rate (all) = " << totalUsefulInvReceivedRate / totalNodes << "\n";

  // if (publicIPNodes > 0) {
  //   std::cout << "Average useful inv sent rate (public IP nodes) =" << totalUsefulInvSentRatePublicIPNode / publicIPNodes << "\n";
  //   std::cout << "Average useless inv megabytes sent (public IP) = " << totaluselessInvSentMegabytesPublicIPNode / publicIPNodes << "\n";
  // }
  //

}

int PoissonDistribution(int value) {
    const uint64_t range_from  = 0;
    const uint64_t range_to    = 1ULL << 48;
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<uint64_t>  distr(range_from, range_to);
    auto bigRand = distr(generator);
    return (int)(log1p(bigRand * -0.0000000000000035527136788 /* -1/2^48 */) * value * -1 + 0.5);
}

std::vector<int> generateTxCreateList(int n, int nodes) {
  std::vector<int> result;
  int averageTxPerNode = n / nodes;
  int alreadyAssigned = 0;
  for (int i = 0; i < nodes - 1; i++) {
    int txToCreate = PoissonDistribution(averageTxPerNode);
    result.push_back(txToCreate);
    alreadyAssigned += txToCreate;
    if (alreadyAssigned > n) {
      result[i] -= (n - alreadyAssigned);
      alreadyAssigned -= (n - alreadyAssigned);
      for (int j = i; j < nodes; j++)
        result.push_back(0);
      break;
    }
  }

  result.push_back(std::max(0, n - alreadyAssigned));
  return result;
}


void CollectTxData(nodeStatistics *stats, int totalNoNodes, int txToCreate,
  int systemId, int systemCount, int nodesInSystemId0, BitcoinTopologyHelper bitcoinTopologyHelper)
{
#ifdef MPI_TEST
  int            blocklen[3] = {1, 1, 1};
  MPI_Aint       disp[3];
  MPI_Datatype   dtypes[3] = {MPI_INT, MPI_INT, MPI_INT};
  MPI_Datatype   mpi_txRecvTime;

  disp[0] = offsetof(txRecvTime, nodeId);
  disp[1] = offsetof(txRecvTime, txHash);
  disp[2] = offsetof(txRecvTime, txTime);

  MPI_Type_create_struct (3, blocklen, disp, dtypes, &mpi_txRecvTime);
  MPI_Type_commit (&mpi_txRecvTime);

  if (systemId != 0 && systemCount > 1)
  {
    for(int i = nodesInSystemId0; i < totalNoNodes; i++)
    {
      Ptr<Node> targetNode = bitcoinTopologyHelper.GetNode(i);
      if (systemId == targetNode->GetSystemId())
      {
          for (int j = 0; j < stats[i].txReceived; j++) {
            MPI_Send(&stats[i].txReceivedTimes[j], 1, mpi_txRecvTime, 0, 9999, MPI_COMM_WORLD);
            auto recv = stats[i].txReceivedTimes[j];
          }
      }
    }
  }
  else if (systemId == 0 && systemCount > 1)
  {
    int count = 0;
    MPI_Status status;
    txRecvTime recv;

    while (count < totalNoNodes)
    {
      if (stats[count].systemId == 0) {
        count++;
        continue;
      }
      for (int j = 0; j < stats[count].txReceived; j++)
       {
          MPI_Recv(&recv, 1, mpi_txRecvTime, MPI_ANY_SOURCE, 9999, MPI_COMM_WORLD, &status);
          stats[recv.nodeId].txReceivedTimes.push_back(recv);
      }
      count++;
    }
  }
#endif
}
