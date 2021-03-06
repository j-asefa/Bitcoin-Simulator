/**
 * This file contains the definitions of the functions declared in bitcoin-node.h
 */

#include "ns3/address.h"
#include "ns3/address-utils.h"
#include "ns3/log.h"
#include "ns3/inet-socket-address.h"
#include "ns3/inet6-socket-address.h"
#include "ns3/node.h"
#include "ns3/socket.h"
#include "ns3/udp-socket.h"
#include "ns3/simulator.h"
#include "ns3/socket-factory.h"
#include "ns3/packet.h"
#include "ns3/trace-source-accessor.h"
#include "ns3/udp-socket-factory.h"
#include "ns3/tcp-socket-factory.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "bitcoin-node.h"
#include "../helper/bitcoin-node-helper.h"
#include <random>
#include <math.h>

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("BitcoinNode");

NS_OBJECT_ENSURE_REGISTERED (BitcoinNode);

int invIntervalSeconds = 10;

TypeId
BitcoinNode::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::BitcoinNode")
    .SetParent<Application> ()
    .SetGroupName("Applications")
    .AddConstructor<BitcoinNode> ()
    .AddAttribute ("Local",
                   "The Address on which to Bind the rx socket.",
                   AddressValue (),
                   MakeAddressAccessor (&BitcoinNode::m_local),
                   MakeAddressChecker ())
    .AddAttribute ("Protocol",
                   "The type id of the protocol to use for the rx socket.",
                   TypeIdValue (UdpSocketFactory::GetTypeId ()),
                   MakeTypeIdAccessor (&BitcoinNode::m_tid),
                   MakeTypeIdChecker ())
    .AddAttribute ("InvTimeoutMinutes",
				   "The timeout of inv messages in minutes",
                   TimeValue (Minutes (20)),
                   MakeTimeAccessor (&BitcoinNode::m_invTimeoutMinutes),
                   MakeTimeChecker())
    .AddTraceSource ("Rx",
                     "A packet has been received",
                     MakeTraceSourceAccessor (&BitcoinNode::m_rxTrace),
                     "ns3::Packet::AddressTracedCallback")
  ;
  return tid;
}

BitcoinNode::BitcoinNode (void) : m_bitcoinPort (8333), m_secondsPerMin(60), m_countBytes (4), m_bitcoinMessageHeader (90),
                                  m_inventorySizeBytes (36), m_getHeadersSizeBytes (72), m_headersSizeBytes (81),
                                  m_averageTransactionSize (522.4), m_txToCreate(0), m_mode(REGULAR),
                                  m_netGroups(2)
{
  NS_LOG_FUNCTION (this);
  m_socket = 0;
  heardTotal = 0;
  firstTimeHops = std::vector<int>(1024);
  m_numberOfPeers = m_peersAddresses.size();
  m_numInvsSent = 0;
  txCreator = false;
}

BitcoinNode::~BitcoinNode(void)
{
  NS_LOG_FUNCTION (this);
}

Ptr<Socket>
BitcoinNode::GetListeningSocket (void) const
{
  NS_LOG_FUNCTION (this);
  return m_socket;
}


std::vector<Ipv4Address>
BitcoinNode::GetPeersAddresses (void) const
{
  NS_LOG_FUNCTION (this);
  return m_peersAddresses;
}


void
BitcoinNode::SetPeersAddresses (const std::vector<Ipv4Address> &peers)
{
  NS_LOG_FUNCTION (this);
  m_peersAddresses = peers;
  m_numberOfPeers = m_peersAddresses.size();
}


void
BitcoinNode::SetPeersDownloadSpeeds (const std::map<Ipv4Address, double> &peersDownloadSpeeds)
{
  NS_LOG_FUNCTION (this);
  m_peersDownloadSpeeds = peersDownloadSpeeds;
}


void
BitcoinNode::SetPeersUploadSpeeds (const std::map<Ipv4Address, double> &peersUploadSpeeds)
{
  NS_LOG_FUNCTION (this);
  m_peersUploadSpeeds = peersUploadSpeeds;
}

void
BitcoinNode::SetNodeInternetSpeeds (const nodeInternetSpeeds &internetSpeeds)
{
  NS_LOG_FUNCTION (this);

  m_downloadSpeed = internetSpeeds.downloadSpeed * 1000000 / 8 ;
  m_uploadSpeed = internetSpeeds.uploadSpeed * 1000000 / 8 ;
}


void
BitcoinNode::SetNodeStats (nodeStatistics *nodeStats)
{
  NS_LOG_FUNCTION (this);
  m_nodeStats = nodeStats;
}

void
BitcoinNode::SetProperties (uint64_t txToCreate, enum ProtocolType protocol, enum ModeType mode, double overlap, int netGroups, int r, int systemId, std::vector<Ipv4Address> outPeers)
{
  NS_LOG_FUNCTION (this);
  m_txToCreate = txToCreate;
  if (txToCreate > 0)
    txCreator = true;
  m_protocol = protocol;
  m_mode = mode;
  m_netGroups = netGroups;
  m_r = r;
  m_systemId = systemId;
  m_outPeers = outPeers;
  m_overlap = overlap;
  if (m_protocol == RECONCILIATION) 
  {
      for (auto peer: m_peersAddresses) 
      {
          std::vector<std::string> s;
          m_peerReconciliationSets.insert(std::pair<Ipv4Address,std::vector<std::string>>(peer, s));
          m_reconcilePeers.push_back(peer);
      }
  }
  if (m_protocol == PREFERRED_DESTINATIONS) {
      peerStatistics peerstats;
      peerstats.numUsefulInvReceived = 0;
      peerstats.numUselessInvReceived = 0;
      peerstats.numGetDataReceived = 0;
      peerstats.numGetDataSent = 0;
      peerstats.connectionLength = 0;
      peerstats.usefulInvRate = 0;
      for(auto peer : m_peersAddresses) {
        m_peerStatistics.insert(std::pair<Ipv4Address, peerStatistics>(peer, peerstats));
      }

      // set up preferred peers
      std::cout << "peer address size: " << m_peersAddresses.size() << std::endl;
      for (int i = 0; i < m_peersAddresses.size(); i++) {
          m_preferredPeers.push_back(m_peersAddresses.at(i));
      }
      std::cout << "preferred peers size initially: " << m_preferredPeers.size() << std::endl;
    }
}

void
BitcoinNode::DoDispose (void)
{
  NS_LOG_FUNCTION (this);
  m_socket = 0;

  // chain up
  Application::DoDispose ();
}


// Application Methods
void
BitcoinNode::StartApplication ()    // Called at time specified by Start
{
  NS_LOG_FUNCTION (this);
  // Create the socket if not already

  srand(time(NULL) + GetNode()->GetId());
  NS_LOG_INFO ("Node " << GetNode()->GetId() << ": download speed = " << m_downloadSpeed << " B/s");
  NS_LOG_INFO ("Node " << GetNode()->GetId() << ": upload speed = " << m_uploadSpeed << " B/s");
  NS_LOG_INFO ("Node " << GetNode()->GetId() << ": m_numberOfPeers = " << m_numberOfPeers);
  NS_LOG_INFO ("Node " << GetNode()->GetId() << ": m_invTimeoutMinutes = " << m_invTimeoutMinutes.GetMinutes() << "mins");

  NS_LOG_INFO ("Node " << GetNode()->GetId() << ": My peers are");

  for (auto it = m_peersAddresses.begin(); it != m_peersAddresses.end(); it++)
    NS_LOG_INFO("\t" << *it);

  double currentMax = 0;

  if (!m_socket)
  {
    m_socket = Socket::CreateSocket (GetNode (), m_tid);
    m_socket->Bind (m_local);
    m_socket->Listen ();
    m_socket->ShutdownSend ();
    if (addressUtils::IsMulticast (m_local))
    {
      Ptr<UdpSocket> udpSocket = DynamicCast<UdpSocket> (m_socket);
      if (udpSocket)
      {
        // equivalent to setsockopt (MCAST_JOIN_GROUP)
        udpSocket->MulticastJoinGroup (0, m_local);
      }
      else
      {
        NS_FATAL_ERROR ("Error: joining multicast on a non-UDP socket");
      }
    }
  }

  m_socket->SetRecvCallback (MakeCallback (&BitcoinNode::HandleRead, this));
  m_socket->SetAcceptCallback (
    MakeNullCallback<bool, Ptr<Socket>, const Address &> (),
    MakeCallback (&BitcoinNode::HandleAccept, this));
  m_socket->SetCloseCallbacks (
    MakeCallback (&BitcoinNode::HandlePeerClose, this),
    MakeCallback (&BitcoinNode::HandlePeerError, this));

  NS_LOG_DEBUG ("Node " << GetNode()->GetId() << ": Before creating sockets");
  for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
  {
    m_peersSockets[*i] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
    m_peersSockets[*i]->Connect (InetSocketAddress (*i, m_bitcoinPort));
  }
  NS_LOG_DEBUG ("Node " << GetNode()->GetId() << ": After creating sockets");

  m_nodeStats->nodeId = GetNode()->GetId();
  m_nodeStats->systemId = GetNode()->GetId();

  m_nodeStats->invReceivedMessages = 0;
  m_nodeStats->invSentMessages = 0;
  m_nodeStats->invReceivedBytes = 0;
  m_nodeStats->invSentBytes = 0;
  m_nodeStats->getDataReceivedMessages = 0;
  m_nodeStats->getDataSentMessages = 0;
  m_nodeStats->getDataReceivedBytes = 0;
  m_nodeStats->getDataSentBytes = 0;
  m_nodeStats->txCreated = 0;
  m_nodeStats->blocksRelayed = 0;
  m_nodeStats->connections = m_peersAddresses.size();

  m_nodeStats->firstSpySuccess = 0;
  m_nodeStats->txReceived = 0;

  m_nodeStats->systemId = m_systemId;

  if (m_protocol == FILTERS_ON_INCOMING_LINKS) 
  {
    RequestIncomingFilters();
  } 
  else if (m_protocol == PREFERRED_DESTINATIONS) 
  {
    
    // In this case we don't set up any filters, but on INV messages we send to preferred peers.
    // maybe we set up preferred peers
  } 
  else if (m_protocol == OUTGOING_FILTERS) 
  {
    ConstructOutgoingFilters();
  } 
  else if (m_protocol == DANDELION_LIKE) 
  {
    ConstructDandelionLinks();
  }
  AnnounceMode();
}

void
BitcoinNode::StopApplication ()     // Called at time specified by Stop
{
  NS_LOG_FUNCTION (this);

  for (std::vector<Ipv4Address>::iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) //close the outgoing sockets
  {
    m_peersSockets[*i]->Close ();
  }


  if (m_socket)
  {
    m_socket->Close ();
    m_socket->SetRecvCallback (MakeNullCallback<void, Ptr<Socket> > ());
  }

  NS_LOG_WARN ("\n\nBITCOIN NODE " << GetNode ()->GetId () << ":");
}

/* 
 * For each incoming peer, choose an outgoing peer != incoming to relay transactions
 * from incoming -> outgoing peers.
 */
void
BitcoinNode::ConstructDandelionLinks (void) 
{
    for (int i = 0; i < m_peersAddresses.size(); i++)
    {
        std::random_device rd; 
		std::mt19937 eng(rd()); 
		std::uniform_int_distribution<> distr(0, m_peersAddresses.size() - 1); 
		int outgoingPeer = distr(eng);

        while (m_DandelionLinks.find(m_peersAddresses.at(outgoingPeer)) != m_DandelionLinks.end()) 
        {
            outgoingPeer = distr(eng);
        }
        
        m_DandelionLinks.insert(std::pair<Ipv4Address, Ipv4Address>(m_peersAddresses.at(i), m_peersAddresses.at(outgoingPeer)));
    }
}

void
BitcoinNode::ReconcileWithNewPeer(void) {
    const uint8_t delimiter[] = "#";
    Ipv4Address peer = m_reconcilePeers.front();
    m_reconcilePeers.pop_front();
    size_t set_size = m_peerReconciliationSets[peer].size();

    rapidjson::Document reconcileData;
    rapidjson::Value value;
    value = RECONCILE_TX_REQUEST;
    reconcileData.SetObject();
    reconcileData.AddMember("message", value, reconcileData.GetAllocator());

    rapidjson::Value setSize;

    setSize.SetInt(set_size);

    reconcileData.AddMember("setSize", setSize, reconcileData.GetAllocator());

    rapidjson::StringBuffer reconcileInfo;
    rapidjson::Writer<rapidjson::StringBuffer> reconcileWriter(reconcileInfo);
    reconcileData.Accept(reconcileWriter);

    m_peersSockets[peer]->Send(reinterpret_cast<const uint8_t*>(reconcileInfo.GetString()), reconcileInfo.GetSize(), 0);
    m_peersSockets[peer]->Send(delimiter, 1, 0);
    m_reconcilePeers.push_back(peer);
}

void
BitcoinNode::UpdatePreferredPeersList(void) 
{
    double maxUsefulInvRate = 0;
    if (m_peerStatistics.empty() || m_preferredPeers.empty())
        NS_FATAL_ERROR ("Error: m_peerStatistics is empty"); // die for now. Should be set up in bitcoin-node-helper.cc
    else { 
        maxUsefulInvRate = m_peerStatistics.begin()->second.usefulInvRate; 
    }

    for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i) 
    {
        if (m_peerStatistics[*i].usefulInvRate > maxUsefulInvRate) {
           m_preferredPeers.push_back(*i); 
        }
    }

    int update = std::min((int) m_preferredPeers.size() / 2, 4); 
    std::cout << "update: " << update << std::endl;
    std::cout << "m_preferredPeers size before: " << m_preferredPeers.size() << std::endl;
    if (update != 0) {
        std::sort(m_preferredPeers.begin(), m_preferredPeers.end());
        std::vector<Ipv4Address>::iterator it = m_preferredPeers.begin();
        std::advance(it, update); 
        m_preferredPeers.erase(it,m_preferredPeers.end());
    }
    std::cout << "m_preferredPeers size after: " << m_preferredPeers.size() << std::endl;
}

Ipv4Address
BitcoinNode::ChooseFromPreferredPeers(void) 
{
    if (m_peerStatistics.empty())
        NS_FATAL_ERROR ("Error: m_peerStatistics is empty"); // die for now. Should be set up in bitcoin-node-helper.cc
    std::random_device rd; 
    std::mt19937 eng(rd()); 
    std::uniform_int_distribution<> distr(0, m_preferredPeers.size() - 1); 
    int outgoingPeer = distr(eng);
    return m_preferredPeers.at(outgoingPeer);
}

void
BitcoinNode::ConstructOutgoingFilters (void)
{
  uint32_t filterLength = (FILTER_BASE_NUMBERING % m_numberOfPeers == 0) ? FILTER_BASE_NUMBERING / m_numberOfPeers : FILTER_BASE_NUMBERING / m_numberOfPeers + 1;
  uint32_t from = 0;
  uint32_t to = m_numberOfPeers > 1 ? ceil((filterLength - 1) * (1 + m_overlap)) : filterLength - 1; // no need for overlap when we only have 1 peer
  for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
  {
      filterBegin[*i] = from;
      filterEnd[*i] = to;
      from += filterLength;
      to += filterLength;
      if (to > FILTER_BASE_NUMBERING && m_overlap == 0)
          to = FILTER_BASE_NUMBERING - 1;
      to = to % FILTER_BASE_NUMBERING;
  }
}

void
BitcoinNode::RequestIncomingFilters (void)
{
  const uint8_t delimiter[] = "#";
  uint32_t filterLength = (FILTER_BASE_NUMBERING % m_numberOfPeers == 0) ? FILTER_BASE_NUMBERING / m_numberOfPeers : FILTER_BASE_NUMBERING / m_numberOfPeers + 1;
  uint32_t from = 0;
  uint32_t to = m_numberOfPeers > 1 ? ceil((filterLength - 1) * (1 + m_overlap)) : filterLength - 1; // no need for overlap when we only have 1 peer
  for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
  {
    rapidjson::Document filterData;
    rapidjson::Value value;
    value = FILTER_REQUEST;
    filterData.SetObject();
    filterData.AddMember("message", value, filterData.GetAllocator());

    rapidjson::Value filterBegin;
    rapidjson::Value filterEnd;

    filterBegin.SetInt(from);
    filterEnd.SetInt(to);

    filterData.AddMember("filterBegin", filterBegin, filterData.GetAllocator());
    filterData.AddMember("filterEnd", filterEnd, filterData.GetAllocator());

    rapidjson::StringBuffer filterInfo;
    rapidjson::Writer<rapidjson::StringBuffer> filterWriter(filterInfo);
    filterData.Accept(filterWriter);

    m_peersSockets[*i]->Send(reinterpret_cast<const uint8_t*>(filterInfo.GetString()), filterInfo.GetSize(), 0);
    m_peersSockets[*i]->Send(delimiter, 1, 0);
    
    from += filterLength;
    to += filterLength;
    if (to > FILTER_BASE_NUMBERING && m_overlap == 0)
        to = FILTER_BASE_NUMBERING - 1;
    to = to % FILTER_BASE_NUMBERING;
  }
}

/*
 * Go through all the filters each peer announced to us and validate that it covers the full transaction range.
 * If not, coordinate with the necessary peers to cover the tx space.
 */
void
BitcoinNode::ValidateNodeFilters(void) {
  uint32_t prevFilterEnd = 0;
  uint32_t maxFilterEnd = 0;
  std::map<uint32_t, uint32_t> filterRange;
  std::map<uint32_t, Ipv4Address> invertedFilterBegin;
  std::map<uint32_t, Ipv4Address> invertedFilterEnd;

  for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); i++) 
  {
    filterRange.insert( std::pair<uint32_t, uint32_t>(filterBegin[*i], filterEnd[*i]));
    invertedFilterBegin.insert( std::pair<uint32_t, Ipv4Address>(filterBegin[*i], *i));
    invertedFilterEnd.insert( std::pair<uint32_t, Ipv4Address>(filterEnd[*i], *i));
  }

  /*
   * Have to loop again, because the map keys are sorted now. 
   * Also, we keep track of the highest "filterEnd" value out of all peers 
   * to make sure at least one peer has assigned us to cover the top end of the filter space.
   */
  for (std::map<uint32_t, uint32_t>::const_iterator i = filterRange.begin(); i != filterRange.end(); i++) {
      if (prevFilterEnd < i->first) {
        /* 
         * In this case, there is a gap in the filters. So extend the upper bound of this peer's filter so it covers the gap.
         * Note this does not cover the case where no peer has assigned us a bucket containing FILTER_BASE_NUMBERING - 1,
         * That is covered below by updating maxFilterEnd
         */
        filterBegin[invertedFilterBegin[i->first]] = prevFilterEnd;
      } 
      
      prevFilterEnd = i->second;
      if (prevFilterEnd > maxFilterEnd)
          maxFilterEnd = prevFilterEnd;
  }

  /*
   * Update our largest filterEnd value to make sure it covers the top end of the Tx space
   */
  if (maxFilterEnd < FILTER_BASE_NUMBERING - 1)
      filterEnd[invertedFilterEnd[maxFilterEnd]] = FILTER_BASE_NUMBERING - 1;
}


/* 
 * Construct and send UPDATE_FILTER_BEGIN message
 */
void
BitcoinNode::UpdateFilterBegin(Ipv4Address& peer, uint32_t newVal) 
{
    const uint8_t delimiter[] = "#";
    rapidjson::Document filterData;

    rapidjson::Value message;
    message = UPDATE_FILTER_BEGIN;
    filterData.SetObject();

    filterData.AddMember("message", message, filterData.GetAllocator());

    rapidjson::Value newFilterBegin;
    newFilterBegin.SetInt(newVal);

    filterData.AddMember("newFilterBegin", newFilterBegin, filterData.GetAllocator());
    
    rapidjson::StringBuffer filterInfo;
    rapidjson::Writer<rapidjson::StringBuffer> filterWriter(filterInfo);
    filterData.Accept(filterWriter);

    m_peersSockets[peer]->Send (reinterpret_cast<const uint8_t*>(filterInfo.GetString()), filterInfo.GetSize(), 0);
    m_peersSockets[peer]->Send(delimiter, 1, 0);
}

/*
 * Construct and send UPDATE_FILTER_END message
 */
void
BitcoinNode::UpdateFilterEnd(Ipv4Address& peer, uint32_t newVal) 
{
    const uint8_t delimiter[] = "#";
    rapidjson::Document filterData;

    rapidjson::Value message;
    message = UPDATE_FILTER_END;
    filterData.SetObject();

    filterData.AddMember("message", message, filterData.GetAllocator());

    rapidjson::Value newFilterEnd;
    newFilterEnd.SetInt(newVal);

    filterData.AddMember("newFilterEnd", newFilterEnd, filterData.GetAllocator());
    
    rapidjson::StringBuffer filterInfo;
    rapidjson::Writer<rapidjson::StringBuffer> filterWriter(filterInfo);
    filterData.Accept(filterWriter);

    m_peersSockets[peer]->Send(reinterpret_cast<const uint8_t*>(filterInfo.GetString()), filterInfo.GetSize(), 0);
    m_peersSockets[peer]->Send(delimiter, 1, 0);
}

void
BitcoinNode::AnnounceMode (void)
{
  int count = 0;
  const uint8_t delimiter[] = "#";

  for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
  {
    rapidjson::Document modeData;

    rapidjson::Value value;
    value = MODE;
    modeData.SetObject();

    modeData.AddMember("message", value, modeData.GetAllocator());

    rapidjson::Value modeValue;
    modeValue.SetInt(m_mode);


    modeData.AddMember("mode", modeValue, modeData.GetAllocator());


    rapidjson::StringBuffer modeInfo;
    rapidjson::Writer<rapidjson::StringBuffer> modeWriter(modeInfo);
    modeData.Accept(modeWriter);

    m_peersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(modeInfo.GetString()), modeInfo.GetSize(), 0);
    m_peersSockets[*i]->Send(delimiter, 1, 0);
  }

  Simulator::Schedule (Seconds(100), &BitcoinNode::ScheduleNextTransactionEvent, this);

}

int MurmurHash3Mixer(int key) // TODO a better hash function
{
  key ^= (key >> 13);
  key *= 0xff51afd7ed558ccd;
  key ^= (key >> 13);
  key *= 0xc4ceb9fe1a85ec53;
  key ^= (key >> 13);
  return key;
}


std::map<int, int64_t> lastGroupInv;

int PoissonNextSend(int averageIntervalSeconds) {
    const uint64_t range_from  = 0;
    const uint64_t range_to    = 1ULL << 48;
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<uint64_t>  distr(range_from, range_to);
    auto bigRand = distr(generator);
    return (int)(log1p(bigRand * -0.0000000000000035527136788 /* -1/2^48 */) * averageIntervalSeconds * -1 + 0.5);
}

int PoissonNextSendTo(int averageIntervalSeconds, int netGroup) {
  auto now = Simulator::Now().GetSeconds();
  if (lastGroupInv[netGroup] < now) {
    auto newDelay = PoissonNextSend(averageIntervalSeconds);
    lastGroupInv[netGroup] = now + newDelay;
    return newDelay;
  }
  return lastGroupInv[netGroup] - now;
}


void
BitcoinNode::ScheduleNextTransactionEvent (void)
{
  NS_LOG_FUNCTION (this);

  if (m_txToCreate == 0)
    return;

  auto delay = 1;
  EventId m_nextTransactionEvent = Simulator::Schedule (Seconds(delay), &BitcoinNode::EmitTransaction, this);
  Simulator::Schedule (Seconds(60), &BitcoinNode::ReconcileWithNewPeer, this);
}


void
BitcoinNode::ScheduleNextBlockEvent (void)
{
  int fixedBlockTimeGeneration = 60*100;
  Simulator::Schedule (Seconds(fixedBlockTimeGeneration), &BitcoinNode::EmitBlock, this);
}

void
BitcoinNode::EmitBlock (void)
{
  rapidjson::Document blockData;

  rapidjson::Value value;
  value = BLOCK;
  blockData.SetObject();

  blockData.AddMember("message", value, blockData.GetAllocator());

  rapidjson::Value blockValue;
  blockValue.SetInt(1);


  blockData.AddMember("block", blockValue, blockData.GetAllocator());


  rapidjson::StringBuffer blockInfo;
  rapidjson::Writer<rapidjson::StringBuffer> blockWriter(blockInfo);
  blockData.Accept(blockWriter);

  const uint8_t delimiter[] = "#";

  for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
  {
    if (peersMode[*i] != REGULAR) {
      continue;
    }

    // Block is 1MB
    int blockSize = 1024 * 1024;
    m_peersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(blockInfo.GetString()), blockSize, 0);
    m_peersSockets[*i]->Send(delimiter, 1, 0);
  }

  m_nodeStats->blocksRelayed+=1;

  ScheduleNextBlockEvent();
}


void
BitcoinNode::EmitTransaction (void)
{
  NS_LOG_FUNCTION (this);

  int nodeId = GetNode()->GetId();
  double currentTime = Simulator::Now().GetSeconds();
  std::ostringstream stringStream;

  m_nodeStats->txCreated++;

  stringStream << m_nodeStats->txCreated << "/" << nodeId;
  std::string transactionId = stringStream.str();
  AdvertiseTransactionInvWrapper(InetSocketAddress::ConvertFrom(m_local).GetIpv4(), transactionId, 0);
  SaveTxData(transactionId);

  if (m_nodeStats->txCreated >= m_txToCreate)
    return;

  ScheduleNextTransactionEvent();
}


void
BitcoinNode::HandleRead (Ptr<Socket> socket)
{
  NS_LOG_FUNCTION (this << socket);
  Ptr<Packet> packet;
  Address from;

  while ((packet = socket->RecvFrom (from)))
  {
      if (packet->GetSize () == 0)
      { //EOF
         break;
      }

      if (InetSocketAddress::IsMatchingType (from))
      {
        /**
         * We may receive more than one packets simultaneously on the socket,
         * so we have to parse each one of them.
         */
        std::string delimiter = "#";
        std::string parsedPacket;
        size_t pos = 0;
        char *packetInfo = new char[packet->GetSize () + 1];
        std::ostringstream totalStream;

        packet->CopyData (reinterpret_cast<uint8_t*>(packetInfo), packet->GetSize ());
        packetInfo[packet->GetSize ()] = '\0'; // ensure that it is null terminated to avoid bugs

        /**
         * Add the buffered data to complete the packet
         */
        totalStream << m_bufferedData[from] << packetInfo;
        std::string totalReceivedData(totalStream.str());
        NS_LOG_INFO("Node " << GetNode ()->GetId () << " Total Received Data: " << totalReceivedData);

        while ((pos = totalReceivedData.find(delimiter)) != std::string::npos)
        {
          parsedPacket = totalReceivedData.substr(0, pos);
          NS_LOG_INFO("Node " << GetNode ()->GetId () << " Parsed Packet: " << parsedPacket);

          rapidjson::Document d;
          d.Parse(parsedPacket.c_str());

          if(!d.IsObject())
          {
            NS_LOG_WARN("The parsed packet is corrupted");
            totalReceivedData.erase(0, pos + delimiter.length());
            continue;
          }

          rapidjson::StringBuffer buffer;
          rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
          d.Accept(writer);

          NS_LOG_INFO ("At time "  << Simulator::Now ().GetSeconds ()
                        << "s bitcoin node " << GetNode ()->GetId () << " received "
                        <<  packet->GetSize () << " bytes from "
                        << InetSocketAddress::ConvertFrom(from).GetIpv4 ()
                        << " port " << InetSocketAddress::ConvertFrom (from).GetPort ()
                        << " with info = " << buffer.GetString());


          switch (d["message"].GetInt())
          {
            case FILTER_REQUEST:
            {
              uint32_t a_filterBegin = d["filterBegin"].GetInt();
              uint32_t a_filterEnd   = d["filterEnd"].GetInt();
              filterBegin[InetSocketAddress::ConvertFrom(from).GetIpv4()] = a_filterBegin;
              filterEnd[InetSocketAddress::ConvertFrom(from).GetIpv4()] = a_filterEnd;
              break;
            }
            case MODE:
            {
              ModeType mode = ModeType(d["mode"].GetInt());
              peersMode[InetSocketAddress::ConvertFrom(from).GetIpv4()] = mode;
              break;
            }
            case UPDATE_FILTER_BEGIN: 
            {
                filterBegin[InetSocketAddress::ConvertFrom(from).GetIpv4()] = d["newFilterBegin"].GetInt();
                break;
            }
            case UPDATE_FILTER_END:
            {
                filterEnd[InetSocketAddress::ConvertFrom(from).GetIpv4()] = d["newFilterEnd"].GetInt();
                break;
            }
            case RECONCILE_TX_REQUEST:
            {
                size_t set =  d["setSize"].GetInt();
                rapidjson::Document reconcileData;
                reconcileData.SetObject();

                rapidjson::Value msg;
                rapidjson::Value txArray(rapidjson::kArrayType);

                rapidjson::Document::AllocatorType& allocator = reconcileData.GetAllocator();

                msg = RECONCILE_TX_RESPONSE;
                reconcileData.AddMember("message", msg, allocator);

                Ipv4Address peer = InetSocketAddress::ConvertFrom(from).GetIpv4();
                for (auto it = m_peerReconciliationSets[peer].begin(); it != m_peerReconciliationSets[peer].end(); ++it)
                {
                    rapidjson::Value txhash;
                    txhash.SetString(it->c_str(), it->length(), allocator);
                    txArray.PushBack(txhash,allocator);
                }

                // try to send as one batch
                reconcileData.AddMember("transactions", txArray, allocator);
                rapidjson::StringBuffer reconcileInfo;
                rapidjson::Writer<rapidjson::StringBuffer> reconcileWriter(reconcileInfo);
                reconcileData.Accept(reconcileWriter);

                const uint8_t delimiter[] = "#";
                m_peersSockets[peer]->Send(reinterpret_cast<const uint8_t*>(reconcileInfo.GetString()), reconcileInfo.GetSize(), 0);
                m_peersSockets[peer]->Send(delimiter, 1, 0);

                m_peerReconciliationSets[peer].clear();

                //move the reconcile initiator to back of "queue"
                m_reconcilePeers.remove(peer); 
                m_reconcilePeers.push_back(peer);
                break;
            }
            case RECONCILE_TX_RESPONSE:
            {
                std::cout << " inside reconcile response beginning\n";
                std::set<std::string> nodeBtransactions;
                Ipv4Address peer = InetSocketAddress::ConvertFrom(from).GetIpv4();
                for (rapidjson::Value::ConstValueIterator itr = d["transactions"].Begin(); itr != d["transactions"].End(); ++itr)
                {
                    if (std::find(knownTxHashes.begin(), knownTxHashes.end(), itr->GetString()) == knownTxHashes.end()) 
                    {
                        rapidjson::Document reconcileData;
                        reconcileData.SetObject();

                        rapidjson::Value txHash;
                        rapidjson::Value msg;
                        msg = GET_DATA;
                        reconcileData.AddMember("message", msg, reconcileData.GetAllocator());
                        txHash.SetString(itr->GetString(), itr->GetStringLength(), reconcileData.GetAllocator());
                        reconcileData.AddMember("txHash", txHash, reconcileData.GetAllocator());

                        rapidjson::StringBuffer reconcileInfo;
                        rapidjson::Writer<rapidjson::StringBuffer> reconcileWriter(reconcileInfo);
                        reconcileData.Accept(reconcileWriter);

                        const uint8_t delimiter[] = "#";
                        m_peersSockets[peer]->Send(reinterpret_cast<const uint8_t*>(reconcileInfo.GetString()), reconcileInfo.GetSize(), 0);
                        m_peersSockets[peer]->Send(delimiter, 1, 0);
                    }
                    nodeBtransactions.insert(itr->GetString());
                }

                for (auto it = knownTxHashes.begin(); it != knownTxHashes.end(); ++it) 
                {
                    if (std::find(nodeBtransactions.begin(), nodeBtransactions.end(), *it) == nodeBtransactions.end()) 
                    {
                        rapidjson::Document reconcileData;
                        reconcileData.SetObject();

                        rapidjson::Value txHash;
                        rapidjson::Value msg;
                        msg = INV;
                        reconcileData.AddMember("message", msg, reconcileData.GetAllocator());
                        txHash.SetString(it->c_str(), it->length(), reconcileData.GetAllocator());
                        reconcileData.AddMember("txHash", txHash, reconcileData.GetAllocator());


                        rapidjson::StringBuffer reconcileInfo;
                        rapidjson::Writer<rapidjson::StringBuffer> reconcileWriter(reconcileInfo);
                        reconcileData.Accept(reconcileWriter);

                        const uint8_t delimiter[] = "#";
                        m_peersSockets[peer]->Send(reinterpret_cast<const uint8_t*>(reconcileInfo.GetString()), reconcileInfo.GetSize(), 0);
                        m_peersSockets[peer]->Send(delimiter, 1, 0);
                    }
                }

                std::cout << " inside reconcile response end\n";
                break;
            }
            case INV:
            {
              m_nodeStats->invReceivedMessages += 1;
              std::vector<std::string>            requestTxs;
              for (int j=0; j<d["inv"].Size(); j++)
              {
                std::string   parsedInv = d["inv"][j].GetString();
                int   hopNumber = d["hop"].GetInt();
                peersKnowTx[parsedInv].push_back(InetSocketAddress::ConvertFrom(from).GetIpv4());
                if(std::find(knownTxHashes.begin(), knownTxHashes.end(), parsedInv) != knownTxHashes.end()) {
                    // TODO: this is an attack vector, because we progressively increase trust in the peer that shows us the most 
                    // INV's that we don't have. But INV's may be invalid, thus allowing an attacker to become our most trusted peer.
                    // need some way to validate.
                    if (m_protocol == PREFERRED_DESTINATIONS) {
                        m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].numUselessInvReceived++; 
                        m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].usefulInvRate = ((double) m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].numUsefulInvReceived) / m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].numUselessInvReceived; 
                        UpdatePreferredPeersList();
                    }
                    else
                        continue;
                } else {
                  if (m_protocol == PREFERRED_DESTINATIONS) {
                      m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].numUsefulInvReceived++; 
                      uint32_t useless = m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].numUselessInvReceived;
                      uint32_t useful = m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].numUsefulInvReceived;
                      m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].usefulInvRate = useless == 0 ? (double) useful : ((double) useful) / useless;
                      UpdatePreferredPeersList();
                  }
                  SaveTxData(parsedInv);
                  if (m_mode == SPY) {
                    heardTotal++;
                    firstTimeHops[hopNumber]++;
                    std::cout << knownTxHashes.size() << ", Fraction heard from creator: " << (float)firstTimeHops[0] / heardTotal << std::endl;
                    m_nodeStats->firstSpySuccess = (float)firstTimeHops[0] / heardTotal;
                  } else {
                    AdvertiseTransactionInvWrapper(from, parsedInv, hopNumber + 1);
                  }
                }
                requestTxs.push_back(parsedInv);
              }
              if (requestTxs.size() == 0)
                break;

              rapidjson::Value   value;
              rapidjson::Value   array(rapidjson::kArrayType);
              d.RemoveMember("inv");

              for (auto tx_it = requestTxs.begin(); tx_it < requestTxs.end(); tx_it++)
              {
                value.SetString(tx_it->c_str(), tx_it->size(), d.GetAllocator());
                array.PushBack(value, d.GetAllocator());
              }

              d.AddMember("transactions", array, d.GetAllocator());
              SendMessage(INV, GET_DATA, d, from);
              if (m_protocol == PREFERRED_DESTINATIONS)
                  m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].numGetDataSent++;
              break;
            }
            case GET_DATA:
            {
              NS_LOG_INFO ("GET_DATA");
              m_nodeStats->getDataReceivedMessages += 1;
              if (m_protocol == PREFERRED_DESTINATIONS)
                  m_peerStatistics[InetSocketAddress::ConvertFrom(from).GetIpv4()].numGetDataReceived++; 
              //SendMessage(GET_DATA, TX, d, from);
              break;
            }
            // case TX:
            // {
            //   rapidjson::StringBuffer buffer;
            //   rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
            //   d.Accept(writer);
            //
            //   for (int j=0; j<d["transactions"].Size(); j++)
            //   {
            //     std::string   parsedInv = d["transactions"][j].GetString();
            //     int   hopNumber = d["hop"].GetInt();            //
            //     // processing delay
            //     auto delay = 10;
            //     Simulator::Schedule(Seconds(delay), &BitcoinNode::AdvertiseNewTransactionInv, this, from, parsedInv, hopNumber + 1);
            //   }
            // }
            default:
              NS_LOG_INFO ("Default");
              break;
          }

          totalReceivedData.erase(0, pos + delimiter.length());
        }

        /**
        * Buffer the remaining data
        */

        m_bufferedData[from] = totalReceivedData;
        delete[] packetInfo;
      }
      else if (Inet6SocketAddress::IsMatchingType (from))
      {
        NS_LOG_INFO ("At time " << Simulator::Now ().GetSeconds ()
                     << "s bitcoin node " << GetNode ()->GetId () << " received "
                     <<  packet->GetSize () << " bytes from "
                     << Inet6SocketAddress::ConvertFrom(from).GetIpv6 ()
                     << " port " << Inet6SocketAddress::ConvertFrom (from).GetPort ());
      }
      m_rxTrace (packet, from);
  }
}

void
BitcoinNode::AdvertiseTransactionInvWrapper (Address from, const std::string transactionHash, int hopNumber)
{
    switch(m_protocol)
    {
        case STANDARD_PROTOCOL:
        {
            AdvertiseNewTransactionInvStandard(from, transactionHash, hopNumber);
            break;
        }
        case PREFERRED_DESTINATIONS:
        {
			AdvertiseNewTransactionInvPreferredPeers(from, transactionHash, hopNumber);
            break;
        }
        case OUTGOING_FILTERS:
        {
			AdvertiseNewTransactionInvFilters(from, transactionHash, hopNumber);
            break;
        }
        case FILTERS_ON_INCOMING_LINKS:
        {
			AdvertiseNewTransactionInvFilters(from, transactionHash, hopNumber);
            break;
        }
        case DANDELION_LIKE:
        {
			AdvertiseNewTransactionInvDandelion(from, transactionHash, hopNumber);
            break;
        }
        case RECONCILIATION:
        {
            AdvertiseNewTransactionInvStandard(from, transactionHash, hopNumber);
        }
    }
}

void
BitcoinNode::AdvertiseNewTransactionInvStandard(Address from, const std::string transactionHash, int hopNumber)
{
  NS_LOG_FUNCTION (this);

  for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
  {
    // if (peersMode[*i] == BLOCKS_ONLY) {
    //   continue;
    // }
    if (*i != InetSocketAddress::ConvertFrom(from).GetIpv4())
    {
      auto delay = 0;
      if (std::find(m_outPeers.begin(), m_outPeers.end(), *i) != m_outPeers.end())
        delay = PoissonNextSend(invIntervalSeconds);
      else
        delay = PoissonNextSend(invIntervalSeconds * 2);
      Simulator::Schedule (Seconds(delay), &BitcoinNode::SendInvToNode, this, *i, transactionHash, hopNumber);
    }
  }
}

void
BitcoinNode::AdvertiseNewTransactionInvDandelion(Address from, const std::string transactionHash, int hopNumber)
{
    NS_LOG_FUNCTION (this);
    auto delay = 0;
    Ipv4Address outPeer = m_DandelionLinks[InetSocketAddress::ConvertFrom(from).GetIpv4()];
    if (std::find(m_outPeers.begin(), m_outPeers.end(), outPeer) != m_outPeers.end()) 
        delay = PoissonNextSend(invIntervalSeconds);
    else 
        delay = PoissonNextSend(invIntervalSeconds * 2);
    Simulator::Schedule (Seconds(delay), &BitcoinNode::SendInvToNode, this, outPeer, transactionHash, hopNumber);
}

void
BitcoinNode::AdvertiseNewTransactionInvPreferredPeers(Address from, const std::string transactionHash, int hopNumber)
{
    NS_LOG_FUNCTION (this);
    auto delay1 = 0;
    auto delay2 = 0;

    if (m_preferredPeers.size() > 1) {
        Ipv4Address peer1 = ChooseFromPreferredPeers();
        // make sure we're not sending to our from peer
        while (peer1.Get() == InetSocketAddress::ConvertFrom(from).GetIpv4().Get())
            peer1 = ChooseFromPreferredPeers();
        Ipv4Address peer2 = ChooseFromPreferredPeers();

        // make sure peer 1 and 2 are different
        while (peer2.Get() == peer1.Get())
            peer2 = ChooseFromPreferredPeers();

        if (std::find(m_outPeers.begin(), m_outPeers.end(), peer1) != m_outPeers.end())
            delay1 = PoissonNextSend(invIntervalSeconds);
        else
            delay1 = PoissonNextSend(invIntervalSeconds * 2);

        if (std::find(m_outPeers.begin(), m_outPeers.end(), peer2) != m_outPeers.end())
            delay2 = PoissonNextSend(invIntervalSeconds);
        else
            delay2 = PoissonNextSend(invIntervalSeconds * 2);


        Simulator::Schedule (Seconds(delay1), &BitcoinNode::SendInvToNode, this, peer1, transactionHash, hopNumber);
        Simulator::Schedule (Seconds(delay2), &BitcoinNode::SendInvToNode, this, peer2, transactionHash, hopNumber);
    } 
}

void
BitcoinNode::AdvertiseNewTransactionInvFilters(Address from, const std::string transactionHash, int hopNumber)
{
    NS_LOG_FUNCTION (this);
    size_t bucket = std::hash<std::string>()(transactionHash) % FILTER_BASE_NUMBERING;
    for (int i = 0; i < m_peersAddresses.size(); i++)
    {
        if (m_peersAddresses.at(i).Get() != InetSocketAddress::ConvertFrom(from).GetIpv4().Get())
        {
            auto delay = 0;
            if (filterBegin[m_peersAddresses.at(i)] < filterEnd[m_peersAddresses.at(i)]) { 
                if (filterBegin[m_peersAddresses.at(i)] <= bucket && bucket < filterEnd[m_peersAddresses.at(i)]) {
                    if (std::find(m_outPeers.begin(), m_outPeers.end(), m_peersAddresses.at(i)) != m_outPeers.end())
                        delay = PoissonNextSend(invIntervalSeconds);
                    else
                        delay = PoissonNextSend(invIntervalSeconds * 2);
                    Simulator::Schedule (Seconds(delay), &BitcoinNode::SendInvToNode, this, m_peersAddresses.at(i), transactionHash, hopNumber);
                }
            } else {
                // special case where we have overlap > 0 and the filter has wrapped around modulo 1000
                if (bucket < filterEnd[m_peersAddresses.at(i)] || bucket >= filterBegin[m_peersAddresses.at(i)]) {
                    if (std::find(m_outPeers.begin(), m_outPeers.end(), m_peersAddresses.at(i)) != m_outPeers.end())
                        delay = PoissonNextSend(invIntervalSeconds);
                    else
                        delay = PoissonNextSend(invIntervalSeconds * 2);
                    Simulator::Schedule (Seconds(delay), &BitcoinNode::SendInvToNode, this, m_peersAddresses.at(i), transactionHash, hopNumber);
                }
            }
        }
    }
}

void
BitcoinNode::SendInvToNode(Ipv4Address receiver, const std::string transactionHash, int hopNumber) {
  if (std::find(peersKnowTx[transactionHash].begin(), peersKnowTx[transactionHash].end(), receiver) != peersKnowTx[transactionHash].end()) 
    return;

  rapidjson::Document inv;
  inv.SetObject();


  rapidjson::Value value;
  value.SetString("tx");
  inv.AddMember("type", value, inv.GetAllocator());

  value = INV;
  inv.AddMember("message", value, inv.GetAllocator());


  rapidjson::Value array(rapidjson::kArrayType);
  value.SetString(transactionHash.c_str(), transactionHash.size(), inv.GetAllocator());
  array.PushBack(value, inv.GetAllocator());
  inv.AddMember("inv", array, inv.GetAllocator());

  value = hopNumber;
  inv.AddMember("hop", value, inv.GetAllocator());


  rapidjson::StringBuffer invInfo;
  rapidjson::Writer<rapidjson::StringBuffer> invWriter(invInfo);
  inv.Accept(invWriter);
  const uint8_t delimiter[] = "#";
  m_peersSockets[receiver]->Send (reinterpret_cast<const uint8_t*>(invInfo.GetString()), invInfo.GetSize(), 0);
  m_peersSockets[receiver]->Send (delimiter, 1, 0);


  m_nodeStats->invSentMessages += 1;
  m_nodeStats->invSentBytes += m_bitcoinMessageHeader + m_countBytes + 1*m_inventorySizeBytes;


  peersKnowTx[transactionHash].push_back(receiver);

}

void
BitcoinNode::SendMessage(enum Messages receivedMessage,  enum Messages responseMessage, rapidjson::Document &d, Address &outgoingAddress)
{
  NS_LOG_FUNCTION (this);

  const uint8_t delimiter[] = "#";

  rapidjson::StringBuffer buffer;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);



  d["message"].SetInt(responseMessage);
  d.Accept(writer);

  Ipv4Address outgoingIpv4Address = InetSocketAddress::ConvertFrom(outgoingAddress).GetIpv4 ();
  std::map<Ipv4Address, Ptr<Socket>>::iterator it = m_peersSockets.find(outgoingIpv4Address);

  if (it == m_peersSockets.end()) //Create the socket if it doesn't exist
  {
    m_peersSockets[outgoingIpv4Address] = Socket::CreateSocket (GetNode (), TcpSocketFactory::GetTypeId ());
    m_peersSockets[outgoingIpv4Address]->Connect (InetSocketAddress (outgoingIpv4Address, m_bitcoinPort));
  }

  m_peersSockets[outgoingIpv4Address]->Send (reinterpret_cast<const uint8_t*>(buffer.GetString()), buffer.GetSize(), 0);
  m_peersSockets[outgoingIpv4Address]->Send (delimiter, 1, 0);


  switch (d["message"].GetInt())
  {
    case INV:
    {
      m_nodeStats->invSentBytes += m_bitcoinMessageHeader + m_countBytes + d["inv"].Size()*m_inventorySizeBytes;
      m_nodeStats->invSentMessages += 1;
      break;
    }
    case GET_DATA:
    {
      m_nodeStats->getDataSentBytes += m_bitcoinMessageHeader + m_countBytes + d["transactions"].Size()*m_inventorySizeBytes;
      m_nodeStats->getDataSentMessages += 1;
      break;
    }
  }
}

void BitcoinNode::SaveTxData(std::string txId) {
  int numericHash = std::hash<std::string>{}(txId);
  txRecvTime txTime;
  txTime.nodeId = GetNode()->GetId();
  txTime.txHash = numericHash;
  txTime.txTime = Simulator::Now().GetSeconds();
  m_nodeStats->txReceivedTimes.push_back(txTime);
  knownTxHashes.push_back(txId);
  m_nodeStats->txReceived++;
}


void
BitcoinNode::HandlePeerClose (Ptr<Socket> socket)
{
  NS_LOG_FUNCTION (this << socket);
}

void BitcoinNode::HandlePeerError (Ptr<Socket> socket)
{
  NS_LOG_FUNCTION (this << socket);
}


void
BitcoinNode::HandleAccept (Ptr<Socket> s, const Address& from)
{
  NS_LOG_FUNCTION (this << s << from);
  s->SetRecvCallback (MakeCallback (&BitcoinNode::HandleRead, this));
}


} // Namespace ns3
