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
#include <random>

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
BitcoinNode::SetProperties (uint64_t txToCreate, enum ProtocolType protocol, enum ModeType mode, int netGroups, int r)
{
  NS_LOG_FUNCTION (this);
  m_txToCreate = txToCreate;
  if (txToCreate > 0)
    txCreator = true;
  m_protocol = protocol;
  m_mode = mode;
  m_netGroups = netGroups;
  m_r = r;
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

  if (m_protocol == FILTERS_ON_LINKS) {
    AnnounceFilters();
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

void
BitcoinNode::AnnounceFilters (void)
{
  const uint8_t delimiter[] = "#";
  rapidjson::Document filterData;
  rapidjson::Value value;
  value = FILTER;
  filterData.SetObject();
  filterData.AddMember("message", value, filterData.GetAllocator());
  rapidjson::Value filterValue;

  int count = 0;
  for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
  {
    // filterValue.SetInt(count++ % 8);

    // Modified filters to choose odd-even
    // filterData.AddMember("nodeId", GetNode()->GetId(), filterData.GetAllocator());
    filterValue.SetInt(count++ % 2);



    filterData.AddMember("filter", filterValue, filterData.GetAllocator());
    rapidjson::StringBuffer filterInfo;
    rapidjson::Writer<rapidjson::StringBuffer> filterWriter(filterInfo);
    filterData.Accept(filterWriter);

    m_peersSockets[*i]->Send (reinterpret_cast<const uint8_t*>(filterInfo.GetString()), filterInfo.GetSize(), 0);
    m_peersSockets[*i]->Send(delimiter, 1, 0);
  }
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

  AdvertiseNewTransactionInv(InetSocketAddress::ConvertFrom(m_local).GetIpv4(), transactionId, 0);

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
            case FILTER:
            {
              uint32_t filter = d["filter"].GetInt();
              filters[from] = filter;
              break;
            }
            case MODE:
            {
              ModeType mode = ModeType(d["mode"].GetInt());
              peersMode[InetSocketAddress::ConvertFrom(from).GetIpv4()] = mode;
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
                peersKnowTx[parsedInv].push_back(from);
                if(std::find(knownTxHashes.begin(), knownTxHashes.end(), parsedInv) != knownTxHashes.end()) {
                  continue;
                } else {
                  SaveTxData(parsedInv);
                  if (m_mode == SPY) {
                    heardTotal++;
                    firstTimeHops[hopNumber]++;
                    std::cout << knownTxHashes.size() << ", Fraction heard from creator: " << (float)firstTimeHops[0] / heardTotal << std::endl;
                    m_nodeStats->firstSpySuccess = (float)firstTimeHops[0] / heardTotal;
                  } else {
                    AdvertiseNewTransactionInv(from, parsedInv, hopNumber + 1);
                  }
                }
                // requestTxs.push_back(parsedInv);
              }

              // Do not need to send getData
              // if (requestTxs.size() == 0)
              //   break;
              //
              // rapidjson::Value   value;
              // rapidjson::Value   array(rapidjson::kArrayType);
              // d.RemoveMember("inv");
              //
              // for (auto tx_it = requestTxs.begin(); tx_it < requestTxs.end(); tx_it++)
              // {
              //   value.SetString(tx_it->c_str(), tx_it->size(), d.GetAllocator());
              //   array.PushBack(value, d.GetAllocator());
              // }
              //
              // d.AddMember("transactions", array, d.GetAllocator());
              // // SendMessage(INV, GET_DATA, d, from);
              break;
            }
            // case GET_DATA:
            // {
            //   NS_LOG_INFO ("GET_DATA");
            //   m_nodeStats->getDataReceivedMessages += 1;
            //   SendMessage(GET_DATA, TX, d, from);
            //   break;
            // }
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
BitcoinNode::AdvertiseNewTransactionInv (Address from, const std::string transactionHash, int hopNumber)
{
  NS_LOG_FUNCTION (this);

  uint numberHash = std::hash<std::string>()(transactionHash);
  int count = 0;

  for (std::vector<Ipv4Address>::const_iterator i = m_peersAddresses.begin(); i != m_peersAddresses.end(); ++i)
  {

    // Does not matter for this simulation
    // if (m_protocol == FILTERS_ON_LINKS && filters[from] && (numberHash % 8) != filters[from]) {
    //   continue;
    // }
    // if (peersMode[*i] == BLOCKS_ONLY) {
    //   continue;
    // }

    if (*i != InetSocketAddress::ConvertFrom(from).GetIpv4())
    {
      auto delay = 0;
      if (peersMode[*i] == SPY) {
        delay = invIntervalSeconds * 100;
        for (int k = 0; k < m_netGroups; k++) {
          delay = std::min(delay, PoissonNextSendTo(invIntervalSeconds * m_r, k));
        }
        Simulator::Schedule (Seconds(delay), &BitcoinNode::SendInvToNode, this, *i, transactionHash, hopNumber);
      } else if (count < 8) {
        delay = PoissonNextSend(invIntervalSeconds);
        Simulator::Schedule (Seconds(delay), &BitcoinNode::SendInvToNode, this, *i, transactionHash, hopNumber);
        count++;
      } else {
        auto netGroup = MurmurHash3Mixer(i->Get()) % m_netGroups;
        delay = PoissonNextSendTo(invIntervalSeconds * m_r, netGroup);
        Simulator::Schedule (Seconds(delay), &BitcoinNode::SendInvToNode, this, *i, transactionHash, hopNumber);
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
