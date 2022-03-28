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
	gRoutingFunctionMap["adaptive_dor_exit_hyperx"] = &adaptive_dor_exit_hyperx;
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

	//printf("%d, %d\n", nodo_destino,nodo_actual);
	//fflush(stdout);
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

//FIXME: THIS IS A WRONG!
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
			int min_occupancy = INT_MAX;
			vector<int> creditos = r->FreeCredits();

			if(in_channel >= gN * (gK-1)){ //inyeccion

					for(int i = 0; i< gN ; i++){

					int salida = (node_vectors[nodo_destino * gN+i] - node_vectors[nodo_actual * gN+i] + gK) %gK;

					if(salida != 0){ //Si hay que recorrer esta salida...

						int occupancy = 0;
						int puerto = i *(gK-1) + salida - 1;

						occupancy = r->GetUsedCredit(puerto);

						if(occupancy < min_occupancy){
							min_occupancy = occupancy;
							dimension_salida = i;
							out_port = puerto; //esto es lo normalizado.
						}

					}

				}
			
			}


			assert(out_port != -1); //no deberia ser...


			if(in_channel >= gN * (gK-1)){ //inyeccion
				vcEnd -= gNumVCs/2;
			}else{
				vcBegin += gNumVCs/2;
			}

			//printf("entry: %d, start: %d, end: %d \n", f->vc, vcBegin, vcEnd);
		}


		outputs->Clear( );

		outputs->AddRange( out_port , vcBegin, vcEnd );
	}

void adaptive_escalera_hyperx( const Router *r, const Flit *f, int in_channel,
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
			int min_occupancy = INT_MAX;
			vector<int> creditos = r->FreeCredits();

			for(int i = 0; i< gN ; i++){

				int salida = (node_vectors[nodo_destino * gN+i] - node_vectors[nodo_actual * gN+i] + gK) %gK;

				if(salida != 0){ //Si hay que recorrer esta salida...

					int occupancy = 0;
					int puerto = i *(gK-1) + salida - 1;

					/*for(int canal = vcBegin; canal < vcEnd; canal++){
						sitios_libres += creditos[puerto*gNumVCs + canal];
					}*/

					occupancy = r->GetUsedCredit(puerto);

					if(occupancy < min_occupancy){
						min_occupancy = occupancy;
						dimension_salida = i;
						out_port = puerto; //esto es lo normalizado.
					}

				}

			}

			assert(out_port != -1); //no deberia ser...


			if(in_channel >= gN * (gK-1)){ //inyeccion
				vcEnd -= gNumVCs/2;
			}else{
				vcBegin += gNumVCs/2;
			}

			printf("entry: %d, start: %d, end: %d \n", f->vc, vcBegin, vcEnd);
		}


		outputs->Clear( );

		outputs->AddRange( out_port , vcBegin, vcEnd );
	}



	void adaptive_dor_exit_hyperx( const Router *r, const Flit *f, int in_channel,
		OutputSet *outputs, bool inject )
		{

			int vcBegin = 0, vcEnd = gNumVCs-1; //quitamos el último canal
			int out_port = -1;

			int nodo_actual = r->GetID();
			int nodo_destino = calculateRouter(f->dest);

			if(inject) {

				out_port = -1;

			} else if(nodo_destino == nodo_actual){

				out_port = gN * (gK -1) + calculateExitPort(f->dest);

			}else if(f->vc == vcEnd){

				out_port = calculateDOR(nodo_destino, nodo_actual);

				vcBegin = vcEnd;

			}else{

				vcEnd--; //quitamos el ultimo canal.
				int flits_disponibles_max = 0;

				for(int i = 0; i< gN ; i++){

					int salida = (node_vectors[nodo_destino * gN+i] - node_vectors[nodo_actual * gN+i] + gK) %gK;

					if(salida != 0){ //Si hay que recorrer esta salida...

						int sitios_libres = 0;
						int puerto = i *(gK-1) + salida - 1;
						vector<int> creditos = r->FreeCredits();

						for(int canal = 0; canal < (gNumVCs -1); canal++){
							sitios_libres += creditos[puerto * (gNumVCs) + canal];
						}

						if(sitios_libres > flits_disponibles_max){

							flits_disponibles_max = sitios_libres;
							out_port = puerto; //esto es lo normalizado.
							break; // no hace falta coger el maximo en principio....

						}

					}

				}

				if(flits_disponibles_max <= 1){ //si no hay flits...
					out_port = calculateDOR(nodo_destino, nodo_actual);
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
						f->intm = RandomInt( powi( gK, gN ) -1); //apuntamos a un router de la red... es mas comodo.
					}

					int nodo_actual = r->GetID();
					int nodo_intermedio = f->intm;
					int nodo_destino = calculateRouter(f->dest); //aqui sacamos todos los routers de la red.

					if(nodo_intermedio == nodo_actual || nodo_destino== nodo_actual){
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
						out_port = gN * (gK -1) + calculateExitPort(f->dest); //sacamos el outpor exacto
						assert(f->ph == 1);
						vcBegin += available_vcs;
					}

				}

				outputs->Clear( );

				outputs->AddRange( out_port , vcBegin, vcEnd );
			}


	/*    //same as ugal except uses xyyx routing
			void ugal_xyyx_hyperx( const Router *r, const Flit *f, int in_channel,
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

						int dest  = flatfly_transformation(f->dest);

						int rID =  r->GetID();
						int _concentration = gC;
						int found;
						int debug = 0;
						int tmp_out_port, _ran_intm;
						int _min_hop, _nonmin_hop, _min_queucnt, _nonmin_queucnt;
						int threshold = 2;


						if ( in_channel < gC ){ //ESTO ESTA MAL, ESTA VIEJO
							if(gTrace){
								cout<<"New Flit "<<f->src<<endl;
							}
							f->ph   = 0;
						}

						if(gTrace){
							int load = 0;
							cout<<"Router "<<rID<<endl;
							cout<<"Input Channel "<<in_channel<<endl;
							//need to modify router to report the buffere depth
							load +=r->GetBufferOccupancy(in_channel);
							cout<<"Rload "<<load<<endl;
						}

						if (debug){
							cout << " FLIT ID: " << f->id << " Router: " << rID << " routing from src : " << f->src <<  " to dest : " << dest << " f->ph: " <<f->ph << " intm: " << f->intm <<  endl;
						}
						// f->ph == 0  ==> make initial global adaptive decision
						// f->ph == 1  ==> route nonminimaly to random intermediate node
						// f->ph == 2  ==> route minimally to destination

						found = 0;

						if (f->ph == 1){
							dest = f->intm;
						}

						if (dest >= rID*_concentration && dest < (rID+1)*_concentration) { //LLEGA A UN NODO DESTINO, PERO LA CONDICION ES BASURA
							if (f->ph == 1) {
								f->ph = 2;
								dest = flatfly_transformation(f->dest);
								if (debug)   cout << "      done routing to intermediate ";
							}
							else  {
								found = 1;
								out_port = dest % gC;
								if (debug)   cout << "      final routing to destination ";
							}
						}

						if (!found) { //SI NO SE HA LLEGAO AL FINAL

							int const xy_available_vcs = (vcEnd - vcBegin + 1) / 2;
							assert(xy_available_vcs > 0);

							// randomly select dimension order at first hop
							bool x_then_y = ((in_channel < gC) ?
							(RandomInt(1) > 0) :
							(f->vc < (vcBegin + xy_available_vcs)));

							if (f->ph == 0) {
								//find the min port and min distance
								_min_hop = find_distance(flatfly_transformation(f->src),dest);
								if(x_then_y){
									tmp_out_port =  flatfly_outport(dest, rID);
								} else {
									tmp_out_port =  flatfly_outport_yx(dest, rID);
								}
								if (f->watch){
									cout << " MIN tmp_out_port: " << tmp_out_port;
								}
								//sum over all vcs of that port
								_min_queucnt =   r->GetUsedCredit(tmp_out_port);

								//find the nonmin router, nonmin port, nonmin count
								_ran_intm = find_ran_intm(flatfly_transformation(f->src), dest);
								_nonmin_hop = find_distance(flatfly_transformation(f->src),_ran_intm) +    find_distance(_ran_intm, dest);
								if(x_then_y){
									tmp_out_port =  flatfly_outport(_ran_intm, rID);
								} else {
									tmp_out_port =  flatfly_outport_yx(_ran_intm, rID);
								}

								if (f->watch){
									cout << " NONMIN tmp_out_port: " << tmp_out_port << endl;
								}
								if (_ran_intm >= rID*_concentration && _ran_intm < (rID+1)*_concentration) {
									_nonmin_queucnt = numeric_limits<int>::max();
								} else  {
									_nonmin_queucnt =   r->GetUsedCredit(tmp_out_port);
								}

								if (debug){
									cout << " _min_hop " << _min_hop << " _min_queucnt: " <<_min_queucnt << " _nonmin_hop: " << _nonmin_hop << " _nonmin_queucnt :" << _nonmin_queucnt <<  endl;
								}

								if (_min_hop * _min_queucnt   <= _nonmin_hop * _nonmin_queucnt +threshold) {

									if (debug) cout << " Route MINIMALLY " << endl;
									f->ph = 2;
								} else {
									// route non-minimally
									if (debug)  { cout << " Route NONMINIMALLY int node: " <<_ran_intm << endl; }
									f->ph = 1;
									f->intm = _ran_intm;
									dest = f->intm;
									if (dest >= rID*_concentration && dest < (rID+1)*_concentration) {
										f->ph = 2;
										dest = flatfly_transformation(f->dest);
									}
								}
							}

							//dest here should be == intm if ph==1, or dest == dest if ph == 2
							if(x_then_y){
								out_port =  flatfly_outport(dest, rID);
								if(out_port >= gC) {
									vcEnd -= xy_available_vcs;
								}
							} else {
								out_port =  flatfly_outport_yx(dest, rID);
								if(out_port >= gC) {
									vcBegin += xy_available_vcs;
								}
							}

							// if we haven't reached our destination, restrict VCs appropriately to avoid routing deadlock
							if(out_port >= gC) {

								int const ph_available_vcs = xy_available_vcs / 2;
								assert(ph_available_vcs > 0);

								if(f->ph == 1) {
									vcEnd -= ph_available_vcs;
								} else {
									assert(f->ph == 2);
									vcBegin += ph_available_vcs;
								}
							}

							found = 1;
						}

						if (!found) {
							cout << " ERROR: output not found in routing. " << endl;
							cout << *f; exit (-1);
						}

						if (out_port >= gN*(gK-1) + gC)  {
							cout << " ERROR: output port too big! " << endl;
							cout << " OUTPUT select: " << out_port << endl;
							cout << " router radix: " <<  gN*(gK-1) + gK << endl;
							exit (-1);
						}

						if (debug) cout << "        through output port : " << out_port << endl;
						if(gTrace){cout<<"Outport "<<out_port<<endl;cout<<"Stop Mark"<<endl;}

					}

					outputs->Clear( );

					outputs->AddRange( out_port , vcBegin, vcEnd );
				}*/
