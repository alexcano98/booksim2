// $Id$

/*
Copyright (c) 2007-2011, Trustees of The Leland Stanford Junior University
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this 
list of conditions and the following disclaimer in the documentation and/or 
other materials provided with the distribution.
Neither the name of the Stanford University nor the names of its contributors 
may be used to endorse or promote products derived from this software without 
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _SYNTHETICTRAFFICMANAGER_HPP_
#define _SYNTHETICTRAFFICMANAGER_HPP_

#include <vector>

#include "trafficmanager.hpp"
#include "traffic.hpp"
#include "stats.hpp"

class SyntheticTrafficManager : public TrafficManager {

protected:

  vector<string> _traffic;
  vector<TrafficPattern *> _traffic_pattern;

  vector<int> _packet_size;

  vector<int> _reply_class;
  vector<int> _request_class;

  vector<vector<int> > _qtime;
  vector<vector<bool> > _qdrained;

  vector<Stats *> _tlat_stats;     
  vector<Stats *> _overall_min_tlat;  
  vector<Stats *> _overall_avg_tlat;  
  vector<Stats *> _overall_max_tlat;  

  vector<vector<Stats *> > _pair_tlat;

  virtual void _RetirePacket(Flit * head, Flit * tail, int dest);

  virtual void _Inject( );

  virtual bool _PacketsOutstanding( ) const;

  virtual void _ResetSim( );

  virtual void _ClearStats( );

  virtual void _UpdateOverallStats( );

  virtual string _OverallStatsCSV(int c) const;

  SyntheticTrafficManager( const Configuration &config, const vector<Network *> & net );

public:

  virtual ~SyntheticTrafficManager( );

  virtual void WriteClassStats(int c, ostream & os = cout) const;
  virtual void DisplayOverallStats(int c, ostream & os = cout) const;

};

#endif