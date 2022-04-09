// $Id$

/*
Copyright (c) 2007-2012, Trustees of The Leland Stanford Junior University
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

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

#ifndef _HYPERX_HPP_
#define _HYPERX_HPP_

#include "network.hpp"
#include "routefunc.hpp"


class Hyperx : public Network {

  bool _mesh;

  int _k;
  int _n;
  int _c;
  int _xr;

  void _ComputeSize( const Configuration &config );
  void _BuildNet( const Configuration &config );

  int _getChannel(int node, int dim, int offset);

public:
  Hyperx( const Configuration &config, const string & name, bool mesh );
  static void RegisterRoutingFunctions();

  int GetN( ) const;
  int GetK( ) const;


//  double Capacity( ) const;

  void InsertRandomFaults( const Configuration &config );
  int calculateRouter(int node);
  int calculateExitPort(int node);
  int calculateDOR_routers(int nodo_destino, int nodo_actual);
  int calculateDOR_ugal(int inyector_destino, int nodo_actual);
  int find_distance_hyperx (int src, int dest);
  int find_ran_intm_hyperx (int src, int dest);


};

void dor_hyperx( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject );
void adaptive_xyyx_hyperx( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject );
void adaptive_dor_exit_hyperx( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject );
void valiant_hyperx( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject );
void ugal_hyperx( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject );
void adaptive_escalera_hyperx( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject );
void omni_war_hyperx( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject );
void adaptive_escape_hyperx( const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject );


#endif
