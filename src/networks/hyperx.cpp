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
  _c = config.GetInt( "c" );

  gC = _c;gK = _k; gN = _n;
  _size     = powi( _k, _n );
  _channels = (_k-1)*_n*_size; //sin contar inyectores

  _nodes = _size * _c;

  node_vectors = (int *) malloc(gN*_size*sizeof(int));
}

void Hyperx::RegisterRoutingFunctions() {

  gRoutingFunctionMap["dor_hyperx"] = &dor_hyperx;
  gRoutingFunctionMap["adaptive_xyyx_hyperx"] = &adaptive_xyyx_hyperx;
  gRoutingFunctionMap["adaptive_dor_escape_hyperx"] = &adaptive_dor_escape_hyperx;
  gRoutingFunctionMap["valiant_hyperx"] = &valiant_hyperx;

}

void Hyperx::_BuildNet( const Configuration &config )
{

  int latency = 1; //Esto igual se cambia en el futuro

  ostringstream router_name;

  //latency type, noc or conventional network
  bool use_noc_latency;
  use_noc_latency = (config.GetInt("use_noc_latency")==1);

  for ( int node = 0; node < _size; ++node ) {

    _routers[node] = Router::NewRouter( config, this, router_name.str( ),node, (_k-1)*_n + _c, (_k-1)*_n + _c); //+1,pero en un futuro +c
    _timed_modules.push_back(_routers[node]);
    router_name.str("");

    for ( int dim = 0; dim < _n; dim++ ) {

      int salto  = powi( _k, dim );

      node_vectors[node * gN + dim] = ( node / salto ) % _k; // % _k//aqui va el powi

      int adj_nodes[_k-1]; //debug
      int input_node_channel[_k-1]; //debug
      int output_node_channel[_k-1]; //debug

      int adj= node;
      int chan_input= dim * (_k-1) + (_k-2);
      int chan_output = dim * (_k-1);

      for(int counter = 0; counter < _k-1; ++counter){ //bucle para todas las adyacencias en una dimension de un nodo

        if( ( (adj/ salto ) % _k ) == (_k-1) ){ //estamos en el borde

          //Aqui habria que darle la vuelta al adj
          adj = adj - (_k-1) * salto - salto; //el ultimo (- salto para dejarlo bien)
          //chan_input+= (-_k+1); //lo pongo apuntando al siguiente, puede ser negativo, va hacia atras
          //chan_output+= (-_k+1); //lo pongo apuntando al siguiente,  puede ser negativo, va hacia atras
        }

        adj = (adj + salto);

        /*if(node == 3 && dim == 1){
        printf("%d, %d, %d \n",adj, _getChannel(adj, dim, chan_input), _getChannel(node, dim, chan_output) );
        fflush(stdout);
      }*/

      adj_nodes[counter] = adj;

      int channel_input = _getChannel(adj, dim, chan_input); //le sacamos el canal al nodo adyacente para llegar a nosotros
      input_node_channel[counter] = channel_input;

      //INPUT CHANNEL
      _routers[node]->AddInputChannel( _chan[channel_input], _chan_cred[channel_input] );

      if(use_noc_latency){
        _chan[channel_input]->SetLatency( latency );
        _chan_cred[channel_input]->SetLatency( latency );
      } else {
        _chan[channel_input]->SetLatency( 1 );
        _chan_cred[channel_input]->SetLatency( 1 );
      }


      int channel_output = _getChannel(node, dim, chan_output);
      output_node_channel[counter] = channel_output;

      //OUTPUT CHANNEL
      _routers[node]->AddOutputChannel( _chan[channel_output], _chan_cred[channel_output] );

      if(use_noc_latency){
        _chan[channel_output]->SetLatency( latency );
        _chan_cred[channel_output]->SetLatency( latency );
      } else {
        _chan[channel_output]->SetLatency( 1 );
        _chan_cred[channel_output]->SetLatency( 1 );
      }

      chan_input = chan_input - 1;
      chan_output = chan_output + 1;


    } //fin ady dim node

    /** HYPERX CONSTRUCCIÓN DE LA DIMENSION DE UN NODO HASTA AQUI. **/

    //debug

    printf("ID: %d, dimension: %d: ",_routers[node]->GetID(), dim);
    for(int i = 0; i < _k-1; i++){
      printf("(adj: %d, output: %d, input: %d)  ",adj_nodes[i], output_node_channel[i], input_node_channel[i]);
    }
    printf("\n");
    fflush(stdout);

  } //fin dim node

  for(int i = 0; i< _c; i++){
    int link = (node * _c) +i;
    //injection and ejection channel, always 1 latency
    _routers[node]->AddInputChannel( _inject[link], _inject_cred[link] );
    _routers[node]->AddOutputChannel( _eject[link], _eject_cred[link] );
    _inject[link]->SetLatency( 1 );
    _eject[link]->SetLatency( 1 );
  }

} //fin node

//Empezamos el debug:
printf("\n ================================================================== \n");
for ( int node = 0; node < _size; ++node ) {

  printf("ID: %d, inputs: %d, outputs: %d \n",_routers[node]->GetID(), _routers[node]->NumInputs(), _routers[node]->NumOutputs());

}
printf("\n");
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
  int base = (_k-1) * _n * node + offset;

  // The offset -- nos colocamos en la dimension dentro del nodo + offset
  //int off  = (_k-1) + offset;

  return base;
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

/*double Hyperx::Capacity( ) const
{
return (double)_k / ( _mesh ? 8.0 : 4.0 );
}*/


int calculateRouter(int node){
  return node/gC;
}

int calculateExitPort(int node){
  return node%gC;
}

int calculateDOR(int nodo_destino, int nodo_actual){

  int salida = 0;
  int i = 0;

  for(i = 0; i< gN ; i++){

    salida = (node_vectors[nodo_destino * gN+i] - node_vectors[nodo_actual * gN+i] + gK) %gK;
    if(salida != 0) break;

  }
  assert(salida != 0); //esto se puede cambiar, meter aqui la salida tmbn

  return i *(gK-1) + salida - 1; //esto es lo normalizado.

}


void dor_hyperx( const Router *r, const Flit *f, int in_channel,
  OutputSet *outputs, bool inject )
  {

    int vcBegin = 0, vcEnd = gNumVCs-1;

    int out_port;

    int nodo_actual = r->GetID();
    int nodo_destino = calculateRouter(f->dest);
    if(inject) {

      out_port = -1;

    } else if(nodo_destino == nodo_actual){

      /*printf("salgo\n" );
      fflush(stdout);*/
      out_port = gN * (gK -1) + calculateExitPort(f->dest);

    }else{

      //printf("%d, end: %d \n", f->vc, vcEnd);
      /*int salida = 0;
      int i = 0;

      for(i = 0; i< gN ; i++){

      salida = (node_vectors[nodo_destino * gN+i] - node_vectors[nodo_actual * gN+i] + gK) %gK;

      if(salida != 0) break;

    }


    //printf("salida: %d, i: %d \n", salida-1, i);
    assert(salida != 0); //esto se puede cambiar, meter aqui la salida tmbn

    out_port = i *(gK-1) + salida - 1; //esto es lo normalizado.*/

    out_port = calculateDOR(nodo_destino, nodo_actual);

    //  printf("actual %d, destino %d, outport %d \n", nodo_actual, nodo_destino, out_port);
    //fflush(stdout);
  }


  outputs->Clear( );



  outputs->AddRange( out_port , vcBegin, vcEnd );
}


void adaptive_xyyx_hyperx( const Router *r, const Flit *f, int in_channel,
  OutputSet *outputs, bool inject )
  {

    int vcBegin = 0, vcEnd = gNumVCs-1;

    int out_port = -1;

    int nodo_actual = r->GetID();
    int nodo_destino = calculateRouter(f->dest);

    if(inject) {

      out_port = -1;

    } else if(nodo_destino == nodo_actual){

      out_port = gN * (gK -1) + calculateExitPort(f->dest);

    }else{

      int dimension_salida = -1;
      int flits_disponibles = -1;
      vector<int> creditos = r->FreeCredits();

      for(int i = 0; i< gN ; i++){

        int salida = (node_vectors[nodo_destino * gN+i] - node_vectors[nodo_actual * gN+i] + gK) %gK;

        if(salida != 0){ //Si hay que recorrer esta salida...

          int sitios_libres = 0;
          for(int canal = 0; canal <gNumVCs; canal++){
            sitios_libres += creditos[i*gNumVCs + canal];
          }

          if(sitios_libres > flits_disponibles){
            flits_disponibles = sitios_libres;
            dimension_salida = i;
            out_port = i *(gK-1) + salida - 1; //esto es lo normalizado.
          }

        }

      }

      assert(out_port != -1); //no deberia ser...

      if(in_channel >= gN * (gK-1)){ //inyeccion
        vcEnd -= gNumVCs/2;
      }else{
        vcBegin += gNumVCs/2;
      }

    }


    outputs->Clear( );

    outputs->AddRange( out_port , vcBegin, vcEnd );
  }



  void adaptive_dor_escape_hyperx( const Router *r, const Flit *f, int in_channel,
    OutputSet *outputs, bool inject )
    {

      int vcBegin = 0, vcEnd = gNumVCs-2; //quitamos el último canal
      bool es_dor = f->vc == gNumVCs-1;

      int out_port = -1;

      int nodo_actual = r->GetID();
      int nodo_destino = calculateRouter(f->dest);

      if(inject) {

        out_port = -1;

      } else if(nodo_destino == nodo_actual){

        out_port = gN * (gK -1) + calculateExitPort(f->dest);

        if(es_dor){ //inyeccion
          vcEnd++;
          vcBegin = vcEnd;
        }

      }else{

        int dimension_salida = -1;
        int flits_disponibles = -1;
        vector<int> creditos = r->FreeCredits();

        for(int i = 0; i< gN ; i++){

          int salida = (node_vectors[nodo_destino * gN+i] - node_vectors[nodo_actual * gN+i] + gK) %gK;

          if(salida != 0){ //Si hay que recorrer esta salida...

            int sitios_libres = 0;
            for(int canal = 0; canal < (gNumVCs -1); canal++){
              sitios_libres += creditos[i*gNumVCs + canal];
            }

            if(sitios_libres > flits_disponibles){ //es DOR ante empates a 0 flits...
              flits_disponibles = sitios_libres;
              dimension_salida = i;
              out_port = i *(gK-1) + salida - 1; //esto es lo normalizado.

              if(es_dor) break; //Entramos a la primera solo
            }

          }

        }

        assert(out_port != -1); //no deberia ser...

        if(flits_disponibles == 0 || es_dor){ //inyeccion
          vcEnd++;
          vcBegin = vcEnd;
        }

      }


      outputs->Clear( );

      outputs->AddRange( out_port , vcBegin, vcEnd );
    }


    void valiant_hyperx( const Router *r, const Flit *f, int in_channel,
      OutputSet *outputs, bool inject )
      {
        // ( Traffic Class , Routing Order ) -> Virtual Channel Range
        int vcBegin = 0, vcEnd = gNumVCs-1;
        if ( f->type == Flit::READ_REQUEST ) {
          vcBegin = gReadReqBeginVC;
          vcEnd = gReadReqEndVC;
        } else if ( f->type == Flit::WRITE_REQUEST ) {
          vcBegin = gWriteReqBeginVC;
          vcEnd = gWriteReqEndVC;
        } else if ( f->type ==  Flit::READ_REPLY ) {
          vcBegin = gReadReplyBeginVC;
          vcEnd = gReadReplyEndVC;
        } else if ( f->type ==  Flit::WRITE_REPLY ) {
          vcBegin = gWriteReplyBeginVC;
          vcEnd = gWriteReplyEndVC;
        }
        assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

        int out_port;


        if(inject) {

          out_port = -1;

        } else {

          if ( in_channel >= (gK-1)*gN ){
            f->ph = 0;
            f->intm = RandomInt( powi( gK, gN )*gC-1);
          }

          int nodo_actual = r->GetID();
          int nodo_intermedio = calculateRouter(f->intm);
          int nodo_destino = calculateRouter(f->dest); //aqui sacamos todos los routers de la red.

          if(nodo_intermedio == nodo_actual|| nodo_destino== nodo_actual){
            f->ph = 1;
          }

          //each class must have at least 2 vcs assigned or else valiant valiant will deadlock
          int const available_vcs = (vcEnd - vcBegin + 1) / 2;
          assert(available_vcs > 0);

          if(nodo_destino != nodo_actual) {

            if(f->ph == 0) {
              out_port = calculateDOR(nodo_intermedio, nodo_actual);
              vcEnd -= available_vcs;
            } else {
              assert(f->ph == 1);
              out_port = calculateDOR(nodo_destino, nodo_actual);
              vcBegin += available_vcs;
            }


          }else{
            out_port = gN * (gK -1) + calculateExitPort(f->dest); //sacamos el destino exacto
            assert(f->ph == 1);
            vcBegin += available_vcs;
          }

        }

        outputs->Clear( );

        outputs->AddRange( out_port , vcBegin, vcEnd );
      }
