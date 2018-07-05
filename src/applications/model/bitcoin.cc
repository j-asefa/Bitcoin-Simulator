/**
 * This file contains the definitions of the functions declared in bitcoin.h
 */


#include "ns3/application.h"
#include "ns3/event-id.h"
#include "ns3/ptr.h"
#include "ns3/traced-callback.h"
#include "ns3/address.h"
#include "ns3/log.h"
#include "bitcoin.h"

namespace ns3 {


const char* getMessageName(enum Messages m)
{
  switch (m)
  {
    case INV: return "INV";
    case GET_DATA: return "GET_DATA";
    case TX: return "TX";
    case FILTER_REQUEST: return "FILTER_REQUEST";
    case MODE: return "MODE";
    case BLOCK: return "BLOCK";
    case UPDATE_FILTER_BEGIN: return "UPDATE_FILTER_BEGIN";
    case UPDATE_FILTER_END: return "UPDATE_FILTER_END";
  }
}


const char* getProtocolType(enum ProtocolType m)
{
  switch (m)
  {
    case STANDARD_PROTOCOL: return "STANDARD_PROTOCOL";
    case FILTERS_ON_INCOMING_LINKS: return "FILTERS_ON_INCOMING_LINKS";
    case PREFERRED_DESTINATIONS: return "PREFERRED_DESTINATIONS";
    case OUTGOING_FILTERS: return "OUTGOING_FILTERS";
    case DANDELION_LIKE: return "DANDELION_LIKE";
  }
}

}// Namespace ns3
