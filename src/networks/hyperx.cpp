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

/*hyperx.cpp
*
*       Hyperx
*
*/

#include "booksim.hpp"
#include <vector>
#include <sstream>
#include "hyperx.hpp"

#include "random_utils.hpp"
#include "misc_utils.hpp"
#include "globals.hpp"


Hyperx::Hyperx( const Configuration &config, const string & name, bool mesh ) :Network( config, name )
{
  _mesh = mesh;

  _ComputeSize( config );
  _Alloc( );
  _BuildNet( config );
}


void Hyperx::_ComputeSize( const Configuration &config )
{
  _k = config.GetInt( "k" );
  _n = config.GetInt( "n" );

  gK = _k; gN = _n;
  _size     = powi( _k, _n );
  _channels = (_k)*_n*_size; //aqui puede se k-1, sin contar inyectores

  _nodes = _size;

  node_vectors = (int *) malloc(gN*_size*sizeof(int));
}

void Hyperx::RegisterRoutingFunctions() {


  gRoutingFunctionMap["dor_hyperx"] = &dor_hyperx;


}

void Hyperx::_BuildNet( const Configuration &config )
{

  int adj_nodes[_k-1];
  int input_node[_k-1];
  int output_node[_k-1];

  int latency = 1; //Esto igual se cambia en el futuro

  ostringstream router_name;

  //latency type, noc or conventional network
  bool use_noc_latency;
  use_noc_latency = (config.GetInt("use_noc_latency")==1);

  for ( int node = 0; node < _size; ++node ) {

    router_name << "router";

    //_routers[node] = Router::NewRouter( config, this, router_name.str( ),node, 2*_n + 1, 2*_n + 1 );

    _routers[node] = Router::NewRouter( config, this, router_name.str( ),node, (_k-1)*_n + 1, (_k-1)*_n + 1); //+1,pero en un futuro +c


    _timed_modules.push_back(_routers[node]);

    router_name.str("");

    for ( int dim = 0; dim < _n; ++dim ) {

      int salto  = powi( _k, dim );

      node_vectors[node * gN + dim] = ( node / salto ) % _k; // % _k//aqui va el powi

      int adj= node;
      int chan_input= dim * (gK-1) -1;
      int chan_output = dim * (gK-1) -1;

      for(int counter = 0; counter < gK-1; ++counter){ //bucle para todas las adyacencias en una dimension de un nodo

        if( ( (adj/ salto ) % _k ) == (_k-1) ){ //estamos en el borde

          //Aqui habria que darle la vuelta al adj
          adj = adj - (_k-1) * salto - salto; //el ultimo (- salto para dejarlo bien)
          //chan_input+= (-_k+1); //lo pongo apuntando al siguiente, puede ser negativo, va hacia atras
          //chan_output+= (-_k+1); //lo pongo apuntando al siguiente,  puede ser negativo, va hacia atras

        }
        adj = (adj + salto);
        chan_input+=1;
        chan_output+=1;


        adj_nodes[counter] = adj;
        input_node[counter] = chan_input;
        output_node[counter] = chan_output;

        int channel = _getChannel(node, dim, chan_input);

        //INPUT CHANNEL
        _routers[node]->AddInputChannel( _chan[channel], _chan_cred[channel] );

        if(use_noc_latency){
          _chan[channel]->SetLatency( latency );
          _chan_cred[channel]->SetLatency( latency );
        } else {
          _chan[channel]->SetLatency( 1 );
          _chan_cred[channel]->SetLatency( 1 );
        }

        //OUTPUT CHANNEL
        _routers[node]->AddOutputChannel( _chan[channel], _chan_cred[channel] );

        if(use_noc_latency){
          _chan[channel]->SetLatency( latency );
          _chan_cred[channel]->SetLatency( latency );
        } else {
          _chan[channel]->SetLatency( 1 );
          _chan_cred[channel]->SetLatency( 1 );
        }

      }

      /** HYPERX CONSTRUCCIÓN DE LA DIMENSION DE UN NODO HASTA AQUI. **/

      //debug

        printf("ID: %d, dimension: %d: ",_routers[node]->GetID(), dim);
        for(int i = 0; i< gK-1; i++){
          printf("(%d, %d, %d)  ",adj_nodes[i], output_node[i], _getChannel(node, dim, i));
        }
        printf("\n");
        fflush(stdout);


    }

    //injection and ejection channel, always 1 latency
    _routers[node]->AddInputChannel( _inject[node], _inject_cred[node] );
    _routers[node]->AddOutputChannel( _eject[node], _eject_cred[node] );
    _inject[node]->SetLatency( 1 );
    _eject[node]->SetLatency( 1 );

  }

  //Empezamos el debug:

  for ( int node = 0; node < _size; ++node ) {

    printf("ID: %d, inputs: %d, outputs: %d \n",_routers[node]->GetID(), _routers[node]->NumInputs(), _routers[node]->NumOutputs());

  }

  fflush(stdout);


}

/*
* node: escalar del nodo al que vamos a poner el canal de input
* dim: dimension en la que se encuentra el canal
* offset: el desplazamiento respecto a node del canal que le vamos a poner dentro de la dimension
*/
int Hyperx::_getChannel(int node, int dim, int offset)
{
  // The base channel -- nos colocamos en el nodo
  int base = (_k-1)*_n*node;

  // The offset -- nos colocamos en la dimension dentro del nodo + offset
  int off  = (_k-1)*dim + offset;

  return ( base + off );
}


int Hyperx::GetN() const
{
  return _n;
}

int Hyperx::GetK() const
{
  return _k;
}

/*legacy, not sure how this fits into the own scheme of things*/
void Hyperx::InsertRandomFaults(const Configuration &config)
{

}

double Hyperx::Capacity( ) const
{
  return (double)_k / ( _mesh ? 8.0 : 4.0 );
}


void dor_hyperx( const Router *r, const Flit *f, int in_channel,
  OutputSet *outputs, bool inject )
  {

    int vcBegin = 0, vcEnd = gNumVCs-1;

    int out_port;

    int nodo_actual = r->GetID();
    int nodo_destino = f->dest;
    if(inject) {

      out_port = -1;

    } else if(nodo_destino == nodo_actual){
      out_port = gN*(gK -1);

    }else{

      //printf("%d, end: %d \n", f->vc, vcEnd);
      int salida = 0;
      int i = 0;

      for(i = 0; i< gN && salida == 0; i++){

        salida = (node_vectors[nodo_destino * gN+i] - node_vectors[nodo_actual * gN+i] + gK) %gK;

      }
      i--; //cogemos el último

      //printf("salida: %d, i: %d \n", salida-1, i);
      assert(salida != 0); //esto se puede cambiar

      out_port = i *(gK-1) + salida - 1; //esto es lo normalizado.

    //  printf("actual %d, destino %d, outport %d \n", nodo_actual, nodo_destino, out_port);
      fflush(stdout);
    }


    outputs->Clear( );



    outputs->AddRange( out_port , vcBegin, vcEnd );
  }
