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
#include <limits>
#include "hyperx.hpp"
#include <numeric>
#include <list>
#include <random>
#include <algorithm>

#include "random_utils.hpp"
#include "misc_utils.hpp"
#include "globals.hpp"
#include "outputset.hpp"

#define HBIAS 3
#define OBIAS 8
#define MAX_OUTPUTS 6
// cREO QUE Lo mejor esta entre 3 y 4 Output ports....

Hyperx::Hyperx(const Configuration &config, const string &name, bool mesh) : Network(config, name)
{
	_mesh = mesh;

	_ComputeSize(config);
	_Alloc();
	_BuildNet(config);
}

void Hyperx::_ComputeSize(const Configuration &config)
{
	_k = config.GetInt("k");
	_n = config.GetInt("n");
	_c = config.GetInt("c");
	_xr = config.GetInt("xr");
	buff_size = config.GetInt("vc_buf_size") * config.GetInt("num_vcs");

	assert(_xr == _c); // cosas para el trafico tornado entre otros....

	gC = _c;
	gK = _k;
	gN = _n;
	_size = powi(_k, _n);
	_channels = (_k - 1) * _n * _size; // sin contar inyectores

	_nodes = _size * _c;

	node_vectors = (int *)malloc(gN * _size * sizeof(int));
}

void Hyperx::RegisterRoutingFunctions()
{

	gRoutingFunctionMap["dor_hyperx"] = &dor_hyperx;
	gRoutingFunctionMap["O1_turn_hyperx"] = &O1_turn_hyperx;
	gRoutingFunctionMap["adaptive_dor_exit_hyperx"] = &adaptive_dor_exit_hyperx; // wrong manner to approach channel escape
	gRoutingFunctionMap["adaptive_escape_hyperx"] = &adaptive_escape_hyperx;
	gRoutingFunctionMap["valiant_hyperx"] = &valiant_hyperx;
	gRoutingFunctionMap["adaptive_escalera_hyperx"] = &adaptive_escalera_hyperx;
	gRoutingFunctionMap["ugal_hyperx"] = &ugal_hyperx;

	gRoutingFunctionMap["omni_war_hyperx"] = &omni_war_random_hyperx;			 //&omni_war_hyperx;  // &omni_war_priority_hyperx;
	gRoutingFunctionMap["omni_war_priority_hyperx"] = &omni_war_random_hyperx; //&omni_war_priority_hyperx; // &omni_war_priority_hyperx;
	gRoutingFunctionMap["omni_war_random_hyperx"] = &omni_war_random_hyperx;
	gRoutingFunctionMap["omni_dor_hyperx"] = &omni_dor_hyperx;

	gRoutingFunctionMap["dal_wormhole_hyperx"] = &dal_wormhole_hyperx;
	gRoutingFunctionMap["dal_vct_hyperx"] = &dal_vct_random_hyperx; // &dal_vct_hyperx;
	gRoutingFunctionMap["dal_vct_random_hyperx"] = &dal_vct_random_hyperx;
	gRoutingFunctionMap["dal_vct_random2_hyperx"] = &dal_vct_random2_hyperx;
	gRoutingFunctionMap["dal_vct_turned_hyperx"] = &dal_vct_turned_hyperx;

	gRoutingFunctionMap["minimal_turn_model_hyperx"] = &minimal_turn_model_hyperx;
	gRoutingFunctionMap["missrouting_turn_model_hyperx"] = &missrouting_turn_model_hyperx;
}

void Hyperx::_BuildNet(const Configuration &config)
{

	int latency = 1; // Esto igual se cambia en el futuro

	ostringstream router_name;

	// latency type, noc or conventional network
	bool use_noc_latency;
	use_noc_latency = (config.GetInt("use_noc_latency") == 1);

	for (int node = 0; node < _size; ++node)
	{

		_routers[node] = Router::NewRouter(config, this, router_name.str(), node, (_k - 1) * _n + _c, (_k - 1) * _n + _c); //+1,pero en un futuro +c
		_timed_modules.push_back(_routers[node]);
		router_name.str("");

		for (int dim = 0; dim < _n; dim++)
		{

			int salto = powi(_k, dim);

			node_vectors[node * gN + dim] = (node / salto) % _k; // % _k//aqui va el powi

			int adj_nodes[_k - 1];			 // debug
			int input_node_channel[_k - 1];	 // debug
			int output_node_channel[_k - 1]; // debug

			int adj = node;
			int chan_input = dim * (_k - 1) + (_k - 2);
			int chan_output = dim * (_k - 1);

			for (int counter = 0; counter < _k - 1; ++counter)
			{ // bucle para todas las adyacencias en una dimension de un nodo

				if (((adj / salto) % _k) == (_k - 1))
				{ // estamos en el borde

					// Aqui habria que darle la vuelta al adj
					adj = adj - (_k - 1) * salto - salto; // el ultimo (- salto para dejarlo bien)
														  // chan_input+= (-_k+1); //lo pongo apuntando al siguiente, puede ser negativo, va hacia atras
														  // chan_output+= (-_k+1); //lo pongo apuntando al siguiente,  puede ser negativo, va hacia atras
				}

				adj = (adj + salto);

				/*if(node == 3 && dim == 1){
				printf("%d, %d, %d \n",adj, _getChannel(adj, dim, chan_input), _getChannel(node, dim, chan_output) );
				fflush(stdout);
				}*/

				adj_nodes[counter] = adj;

				int channel_input = _getChannel(adj, dim, chan_input); // le sacamos el canal al nodo adyacente para llegar a nosotros
				input_node_channel[counter] = channel_input;

				// INPUT CHANNEL
				_routers[node]->AddInputChannel(_chan[channel_input], _chan_cred[channel_input]);

				if (use_noc_latency)
				{
					_chan[channel_input]->SetLatency(latency);
					_chan_cred[channel_input]->SetLatency(latency);
				}
				else
				{
					_chan[channel_input]->SetLatency(1);
					_chan_cred[channel_input]->SetLatency(1);
				}

				int channel_output = _getChannel(node, dim, chan_output);
				output_node_channel[counter] = channel_output;

				// OUTPUT CHANNEL
				_routers[node]->AddOutputChannel(_chan[channel_output], _chan_cred[channel_output]);

				if (use_noc_latency)
				{
					_chan[channel_output]->SetLatency(latency);
					_chan_cred[channel_output]->SetLatency(latency);
				}
				else
				{
					_chan[channel_output]->SetLatency(1);
					_chan_cred[channel_output]->SetLatency(1);
				}

				chan_input = chan_input - 1;
				chan_output = chan_output + 1;

			} // fin ady dim node

			/** HYPERX CONSTRUCCIÓN DE LA DIMENSION DE UN NODO HASTA AQUI. **/

			// debug

			/*printf("ID: %d, dimension: %d: ",_routers[node]->GetID(), dim);
			for(int i = 0; i < _k-1; i++){
				printf("(adj: %d, output: %d, input: %d)  ",adj_nodes[i], output_node_channel[i], input_node_channel[i]);
			}
			printf("\n");
			fflush(stdout);*/

		} // fin dim node

		for (int i = 0; i < _c; i++)
		{
			int link = (node * _c) + i;
			// injection and ejection channel, always 1 latency
			_routers[node]->AddInputChannel(_inject[link], _inject_cred[link]);
			_routers[node]->AddOutputChannel(_eject[link], _eject_cred[link]);
			_inject[link]->SetLatency(1);
			_eject[link]->SetLatency(1);
		}

	} // fin node

	// Empezamos el debug:
	/*printf("\n ================================================================== \n");
	for ( int node = 0; node < _size; ++node ) {

		printf("ID: %d, inputs: %d, outputs: %d \n",_routers[node]->GetID(), _routers[node]->NumInputs(), _routers[node]->NumOutputs());

	}
	printf("\n");
	fflush(stdout);*/

	// loop through _k

	/*int num_restrictions = _k / 2;

	std::vector<int> numeros(_k - 1);
	std::iota(numeros.begin(), numeros.end(), 0);
	extern map<int, vector<int>> node_restrictions_hyperx;
	// node_restrictions_hyperx
	for (size_t i = 0; i < _k; i++)
	{
		node_restrictions_hyperx[i].resize(num_restrictions);
		std::random_shuffle(numeros.begin(), numeros.end());

		for (size_t j = 0; j < num_restrictions; j++)
		{
			node_restrictions_hyperx[i][j] = numeros[j];
		}
	}*/


	std::vector<int> enlaces((_k-1) * _k );
	std::iota(enlaces.begin(), enlaces.end(), 0);
	std::shuffle(enlaces.begin(), enlaces.end(), std::mt19937(std::random_device()()));
	//print secuencia de enlaces
	/*printf("\n ================================================================== \n");
	for (int i = 0; i < _n *(_k-1) * _k; i++)
	{
		printf("%d, ",enlaces[i]);
	}
	printf("\n");*/

	//iter the _size
	for (int node = 0; node < _k; node++)
	{
		//iter the through k-1 links
		for (int adj = 0; adj < (_k - 1); adj++)
		{

			int link_actual = (node * (_k - 1)) + adj;
			int link_actual_order = enlaces[link_actual];
			int next_node = (node + adj +1) % _k;
			//print status
			/*printf("\n ================================================================== \n");
			printf("Link: %d, node: %d, next_node: %d, link_actual_order: %d \n", link_actual, node, next_node, link_actual_order);
			fflush(stdout);*/
			
			assert(next_node != node);
			
			for (int j = 0; j < (_k - 1); j++)
			{
				int link_next_node = (next_node * (_k - 1)) + j;
				int link_next_node_order = enlaces[link_next_node];
				int next_next_node = (next_node + j +1) % _k;

				if(link_next_node_order > link_actual_order && next_next_node != node){
					
					link_restrictions_hyperx[node][next_next_node].push_back(next_node);
				}	
			}
		}
	}
}

/*
 * node: escalar del nodo al que vamos a poner el canal de input
 * dim: dimension en la que se encuentra el canal
 * offset: el desplazamiento respecto a node del canal que le vamos a poner dentro de la dimension
 */
int Hyperx::_getChannel(int node, int dim, int offset)
{
	// The base channel -- nos colocamos en el nodo
	int base = (_k - 1) * _n * node + offset;

	// The offset -- nos colocamos en la dimension dentro del nodo + offset
	// int off  = (_k-1) + offset;

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

int calculateRouter(int node)
{
	return node / gC;
}

int calculateExitPort(int node)
{
	return node % gC;
}

int calculateDOR_routers(int nodo_destino, int nodo_actual)
{ // router y router

	int salida = 0;
	int i = 0;

	for (i = 0; i < gN; i++)
	{

		salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;
		if (salida != 0)
			break;
	}

	// printf("%d, %d\n", nodo_destino,nodo_actual);
	// fflush(stdout);
	assert(salida != 0); // esto se puede cambiar, meter aqui la salida tmbn

	return i * (gK - 1) + salida - 1; // esto es lo normalizado.
}

int calculateDORYX_routers(int nodo_destino, int nodo_actual)
{

	int salida = 0;
	int i = gN - 1;

	for (i = gN - 1; i >= 0; i--)
	{

		salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

		if (salida != 0)
		{
			break;
		}
	}

	// printf("%d, %d\n", nodo_destino,nodo_actual);
	// fflush(stdout);
	assert(i >= 0 && i < gN); // esto se puede cambiar, meter aqui la salida tmbn

	return i * (gK - 1) + salida - 1; // esto es lo normalizado.
}

int calculateDOR_ugal(int inyector_destino, int nodo_actual)
{

	int salida = 0;
	int i = 0;
	int nodo_destino = calculateRouter(inyector_destino);
	// int nodo_actual = calculateRouter(inyector_actual);

	for (i = 0; i < gN; i++)
	{

		salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;
		if (salida != 0)
			break;
	}

	if (salida == 0)
	{

		return gN * (gK - 1) + calculateExitPort(inyector_destino);
	}
	else
	{

		return i * (gK - 1) + salida - 1; // esto es lo normalizado.
	}
}

int find_distance_hyperx(int src, int dest)
{
	int dist = 0;
	int dim = gN;

	int src_tmp = (int)src / gC;
	int dest_tmp = (int)dest / gC;

	for (int d = 0; d < dim; d++)
	{

		int src_id = src_tmp % gK;
		int dest_id = dest_tmp % gK;
		if (src_id != dest_id)
			dist++;
		src_tmp = (int)(src_tmp / gK);
		dest_tmp = (int)(dest_tmp / gK);
	}

	return dist;
}

int getDataForVCRange(int vcBegin, int vcEnd, vector<int> &data)
{
	int suma = 0;
	for (int i = vcBegin; i <= vcEnd; i++)
	{
		suma += data[i];
	}
	return suma / (vcEnd - vcBegin + 1);
}

/*Para Ugal...*/
int find_ran_intm_hyperx(int src, int dest)
{
	int _dim = gN;
	int _dim_size;
	int _ran_dest = 0;
	int debug = 0;

	if (debug)
		cout << " INTM node for  src: " << src << " dest: " << dest << endl;

	src = (int)(src / gC);
	dest = (int)(dest / gC);

	_ran_dest = RandomInt(gC - 1);
	if (debug)
		cout << " ............ _ran_dest : " << _ran_dest << endl;

	for (int d = 0; d < _dim; d++)
	{
		_dim_size = powi(gK, d) * gC;
		if ((src % gK) == (dest % gK))
		{
			_ran_dest += (src % gK) * _dim_size;
			if (debug)
				cout << "    share same dimension : " << d << " int node : " << _ran_dest << " src ID : " << src % gK << endl;
		}
		else
		{
			// src and dest are in the same dimension "d" + 1
			// ==> thus generate a random destination within
			_ran_dest += RandomInt(gK - 1) * _dim_size;
			if (debug)
				cout << "    different  dimension : " << d << " int node : " << _ran_dest << " _dim_size: " << _dim_size << endl;
		}
		src = (int)(src / gK);
		dest = (int)(dest / gK);
	}

	if (debug)
		cout << " intermediate destination NODE: " << _ran_dest << endl;
	return _ran_dest;
}

int getPriority_dal(int buff_size, int occupancy, int distance_to_dest, int missrouting, int escape_channel)
{

	int occ_prio = (buff_size - occupancy) / gNumVCs + OBIAS;
	int hops_prio = (gN - distance_to_dest + HBIAS);

	// int prio_miss = (buff_size - r->GetUsedCredit(puerto_miss)) * (gN - distance_to_dest);

	if (escape_channel)
	{
		hops_prio = 1;

		if (missrouting)
		{
			hops_prio = 0.5;
		}
	}
	assert(occ_prio >= 0 && hops_prio >= 0);

	return occ_prio * hops_prio;
}

int getPriority_omni(int free_credits, int distance_to_dest, int missroute)
{ // int buff_size, int occupancy, int distance_to_dest

	int occ_bias = 8; // buffer size = 8...
	int hop_bias = 3; // hops = 3...
	return (free_credits + occ_bias) * (gN - distance_to_dest + hop_bias);

	// return occ_bias * 0.5;
}

void dor_hyperx(const Router *r, const Flit *f, int in_channel,
				OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	int out_port;

	int nodo_actual = r->GetID();
	int nodo_destino = calculateRouter(f->dest);
	if (inject)
	{

		out_port = -1;
	}
	else if (nodo_destino == nodo_actual)
	{

		/*printf("salgo\n" );
		fflush(stdout);*/
		out_port = gN * (gK - 1) + calculateExitPort(f->dest);
	}
	else
	{

		// printf("%d, end: %d \n", f->vc, vcEnd);
		/*int salida = 0;
		int i = 0;

		for(i = 0; i< gN ; i++){

		salida = (node_vectors[nodo_destino * gN+i] - node_vectors[nodo_actual * gN+i] + gK) %gK;

		if(salida != 0) break;

		}


		//printf("salida: %d, i: %d \n", salida-1, i);
		assert(salida != 0); //esto se puede cambiar, meter aqui la salida tmbn

		out_port = i *(gK-1) + salida - 1; //esto es lo normalizado.*/

		out_port = calculateDOR_routers(nodo_destino, nodo_actual);

		//  printf("actual %d, destino %d, outport %d \n", nodo_actual, nodo_destino, out_port);
		// fflush(stdout);
		// print used credits of the output port
		//printf("OCCUPANCY port %d: %d\n", out_port, r->GetUsedCredit(out_port));
	}

	outputs->Clear();

	outputs->AddRange(out_port, vcBegin, vcEnd);
}

// ZYX or XYZ routing...
void O1_turn_hyperx(const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	assert(gNumVCs >= 2);
	assert(gNumVCs >= gN);

	int out_port = -1;
	outputs->Clear();
	int nodo_actual = r->GetID();
	int nodo_destino = calculateRouter(f->dest);

	if (inject)
	{
		// out_port = -1;
		outputs->AddRange(-1, vcBegin, vcEnd);
	}
	else if (nodo_destino == nodo_actual)
	{

		out_port = gN * (gK - 1) + calculateExitPort(f->dest);
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{
		vector<int> free_credits = r->FreeCredits();
		int available_vcs = gNumVCs / 2;

		if (in_channel >= gN * (gK - 1)) // inyeccion
		{

			int out_xy = calculateDOR_routers(nodo_destino, nodo_actual);
			int out_yx = calculateDORYX_routers(nodo_destino, nodo_actual);

			int prio_xy = getDataForVCRange(out_xy * gNumVCs + vcBegin, out_xy * gNumVCs + available_vcs - 1, free_credits);
			int prio_yx = getDataForVCRange(out_yx * gNumVCs + available_vcs, out_yx * gNumVCs + vcEnd, free_credits);

			outputs->AddRange(out_xy, vcBegin, available_vcs - 1, prio_xy);
			outputs->AddRange(out_yx, available_vcs, vcEnd, prio_yx);
		}
		else
		{
			int dimension_salida = f->vc / available_vcs;
			if (dimension_salida == 1)
			{

				out_port = calculateDORYX_routers(nodo_destino, nodo_actual);
				vcBegin = available_vcs;
			}
			else
			{

				out_port = calculateDOR_routers(nodo_destino, nodo_actual);
				vcEnd = available_vcs - 1;
			}

			// printf("vcBegin: %d, vcEnd: %d, hops: %d", vcBegin, vcEnd, f->hops);
			// fflush(stdout);

			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
	}
}

void adaptive_escalera_hyperx(const Router *r, const Flit *f, int in_channel,
							  OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	int out_port = -1;

	outputs->Clear();

	int nodo_actual = r->GetID();
	int nodo_destino = calculateRouter(f->dest);

	if (inject)
	{

		// out_port = -1;
		outputs->AddRange(-1, vcBegin, vcEnd);
	}
	else if (nodo_destino == nodo_actual)
	{

		out_port = gN * (gK - 1) + calculateExitPort(f->dest);
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{

		int available_vcs = gNumVCs / gN;
		vector<int> free_credits = r->FreeCredits();

		if (in_channel >= gN * (gK - 1))
		{ // inyeccion
			vcBegin = 0;
			vcEnd = available_vcs - 1;
		}
		else
		{
			vcBegin = f->hops * available_vcs;
			vcEnd = vcBegin + available_vcs - 1;
		}

		// printf("vcBegin: %d, vcEnd: %d, hops: %d",vcBegin, vcEnd, f->hops);
		// fflush(stdout);

		for (int i = 0; i < gN; i++)
		{

			int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

			if (salida != 0)
			{ // Si hay que recorrer esta salida...

				int puerto = i * (gK - 1) + salida - 1;
				int prio_min = getDataForVCRange(puerto * gNumVCs + vcBegin, puerto * gNumVCs + vcEnd, free_credits);
				outputs->AddRange(puerto, vcBegin, vcEnd, prio_min);

				/*occupancy = r->GetUsedCredit(puerto);
				if (occupancy < min_occupancy)
				{
					min_occupancy = occupancy;
					dimension_salida = i;
					out_port = puerto;
				}*/
			}
		}
	}
}

int getCreditOutportVC_hyperx(int outport, int vc, vector<int> &creditos)
{
	return creditos[outport * gNumVCs + vc];
}

int maxCreditsVC_hyperx(int outport, int canales, vector<int> &creditos)
{

	int vc = -1;
	int max_creditos = -1;
	for (int i = 0; i < canales; i++)
	{
		int c = getCreditOutportVC_hyperx(outport, i, creditos);

		if (c > max_creditos)
		{
			max_creditos = c;
			vc = i;
		}
	}

	return vc;
}

void adaptive_dor_exit_hyperx(const Router *r, const Flit *f, int in_channel,
							  OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1; // quitamos el último canal

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

	int out_port = -1;

	int nodo_actual = r->GetID();
	int nodo_destino = calculateRouter(f->dest);

	if (inject)
	{

		out_port = -1;
	}
	else if (nodo_destino == nodo_actual)
	{

		out_port = gN * (gK - 1) + calculateExitPort(f->dest);
	}
	else if (f->vc == vcEnd)
	{

		out_port = calculateDOR_routers(nodo_destino, nodo_actual);

		vcBegin = vcEnd;
	}
	else
	{

		vcEnd--; // quitamos el ultimo canal.

		int vc_salida;
		int max_creditos_total = -1;

		for (int i = 0; i < gN; i++)
		{

			int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

			if (salida != 0)
			{ // Si hay que recorrer esta salida...

				int puerto = i * (gK - 1) + salida - 1;

				vector<int> creditos = r->FreeCredits();

				int canalvc = maxCreditsVC_hyperx(puerto, gNumVCs - 1, creditos);
				int max_cred = getCreditOutportVC_hyperx(puerto, canalvc, creditos);

				if (max_cred >= max_creditos_total)
				{
					max_creditos_total = max_cred;
					vc_salida = canalvc;
					out_port = puerto;
				}
			}
		}

		assert(max_creditos_total != -1);

		if (max_creditos_total <= 0)
		{ // si no hay flits...
			out_port = calculateDOR_routers(nodo_destino, nodo_actual);
			vcEnd++;
			vcBegin = vcEnd;
		}
		else
		{
			vcBegin = vc_salida;
			vcEnd = vc_salida;
		}
	}

	outputs->Clear();

	outputs->AddRange(out_port, vcBegin, vcEnd);
}

void valiant_hyperx(const Router *r, const Flit *f, int in_channel,
					OutputSet *outputs, bool inject)
{
	// ( Traffic Class , Routing Order ) -> Virtual Channel Range
	int vcBegin = 0, vcEnd = gNumVCs - 1;
	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

	int out_port;

	if (inject)
	{

		out_port = -1;
	}
	else
	{

		int nodo_actual = r->GetID();

		if (in_channel >= (gK - 1) * gN)
		{
			f->ph = 0;

			f->intm = RandomInt(powi(gK, gN) - 1); // apuntamos a un router de la red... es mas comodo.
			while(f->intm == nodo_actual)
			{
				f->intm = RandomInt(powi(gK, gN) - 1);
			}
		}
		
		int nodo_intermedio = f->intm;
		int nodo_destino = calculateRouter(f->dest); // aqui sacamos todos los routers de la red.

		if (nodo_intermedio == nodo_actual)
		{
			f->ph = 1;
		}

		// each class must have at least 1 vcs assigned or else valiant will deadlock
		int const available_vcs = (vcEnd - vcBegin + 1) / 2;
		assert(available_vcs > 0);

		if (nodo_destino == nodo_actual)
		{
			f->ph = 1;
			out_port = gN * (gK - 1) + calculateExitPort(f->dest); // sacamos el outpor exacto
		}
		else if (f->ph == 0)
		{
			out_port = calculateDOR_routers(nodo_intermedio, nodo_actual);
			vcBegin = available_vcs;
		}
		else
		{
			assert(f->ph == 1);
			out_port = calculateDOR_routers(nodo_destino, nodo_actual);
			vcEnd = available_vcs - 1;
		}
	}

	outputs->Clear();

	outputs->AddRange(out_port, vcBegin, vcEnd);
}

void ugal_hyperx(const Router *r, const Flit *f, int in_channel,
				 OutputSet *outputs, bool inject)
{
	// ( Traffic Class , Routing Order ) -> Virtual Channel Range
	int vcBegin = 0, vcEnd = gNumVCs - 1;
	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

	int out_port;

	if (inject)
	{

		out_port = -1;
	}
	else
	{

		int dest = f->dest; // flatfly_transformation(f->dest);

		int rID = r->GetID();
		// int _concentration = gC;
		int found;
		int debug = 0;
		int tmp_out_port, _ran_intm;
		int _min_hop, _nonmin_hop, _min_queucnt, _nonmin_queucnt;
		int threshold = 2;

		if (in_channel >= (gK - 1) * gN)
		{ // inyección
			if (gTrace)
			{
				cout << "New Flit " << f->src << endl;
			}
			f->ph = 0;
		}

		if (gTrace)
		{
			int load = 0;
			cout << "Router " << rID << endl;
			cout << "Input Channel " << in_channel << endl;
			// need to modify router to report the buffere depth
			load += r->GetBufferOccupancy(in_channel);
			cout << "Rload " << load << endl;
		}

		if (debug)
		{
			cout << " FLIT ID: " << f->id << " Router: " << rID << " routing from src : " << f->src << " to dest : " << dest << " f->ph: " << f->ph << " intm: " << f->intm << endl;
		}
		// f->ph == 0  ==> make initial global adaptive decision
		// f->ph == 1  ==> route nonminimaly to random intermediate node
		// f->ph == 2  ==> route minimally to destination

		found = 0;

		if (f->ph == 1)
		{
			dest = f->intm;
		}

		if (calculateRouter(dest) == rID)
		{ // if (dest >= rID*_concentration && dest < (rID+1)*_concentration) {

			if (f->ph == 1)
			{
				f->ph = 2;
				dest = f->dest; // flatfly_transformation(f->dest);
				if (debug)
					cout << "      done routing to intermediate ";
			}
			else
			{
				found = 1;
				out_port = gN * (gK - 1) + calculateExitPort(dest); // dest % gC;
				if (debug)
					cout << "      final routing to destination ";
			}
		}

		if (!found)
		{

			if (f->ph == 0)
			{
				_min_hop = find_distance_hyperx(f->src, dest);
				_ran_intm = find_ran_intm_hyperx(f->src, dest);
				tmp_out_port = calculateDOR_ugal(dest, rID); // flatfly_outport(dest, rID);
				if (f->watch)
				{
					*gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
							   << " MIN tmp_out_port: " << tmp_out_port;
				}

				_min_queucnt = r->GetUsedCredit(tmp_out_port);

				_nonmin_hop = find_distance_hyperx(f->src, _ran_intm) + find_distance_hyperx(_ran_intm, dest);
				tmp_out_port = calculateDOR_ugal(_ran_intm, rID); // flatfly_outport(_ran_intm, rID);

				if (f->watch)
				{
					*gWatchOut << GetSimTime() << " | " << r->FullName() << " | "
							   << " NONMIN tmp_out_port: " << tmp_out_port << endl;
				}

				if (calculateRouter(_ran_intm) == rID)
				{ // if (_ran_intm >= rID*_concentration && _ran_intm < (rID+1)*_concentration) { //si ya estamos en intermedio pues fijo vamos minimo
					_nonmin_queucnt = numeric_limits<int>::max();
				}
				else
				{
					_nonmin_queucnt = r->GetUsedCredit(tmp_out_port);
				}

				if (debug)
				{
					cout << " _min_hop " << _min_hop << " _min_queucnt: " << _min_queucnt << " _nonmin_hop: " << _nonmin_hop << " _nonmin_queucnt :" << _nonmin_queucnt << endl;
				}

				if (_min_hop * _min_queucnt <= _nonmin_hop * _nonmin_queucnt + threshold)
				{

					if (debug)
						cout << " Route MINIMALLY " << endl;
					f->ph = 2;
				}
				else
				{
					// route non-minimally
					if (debug)
					{
						cout << " Route NONMINIMALLY int node: " << _ran_intm << endl;
					}
					f->ph = 1;
					f->intm = _ran_intm;
					dest = f->intm;
					if (calculateRouter(dest) == rID)
					{ // if (dest >= rID*_concentration && dest < (rID+1)*_concentration) {
						f->ph = 2;
						dest = f->dest; // flatfly_transformation(f->dest);
					}
				}
			}

			// find minimal correct dimension to route through
			out_port = calculateDOR_ugal(dest, rID); // flatfly_outport(dest, rID);

			// if we haven't reached our destination, restrict VCs appropriately to avoid routing deadlock
			if (out_port < (gK - 1) * gN)
			{ //>= (gK-1)*gN

				int const available_vcs = (vcEnd - vcBegin + 1) / 2;
				assert(available_vcs > 0);
				if (f->ph == 1)
				{
					vcEnd -= available_vcs;
				}
				else
				{
					assert(f->ph == 2);
					vcBegin += available_vcs;
				}
			}

			found = 1;
		}

		if (!found)
		{
			cout << " ERROR: output not found in routing. " << endl;
			cout << *f;
			exit(-1);
		}

		if (out_port >= gN * (gK - 1) + gC)
		{
			cout << " ERROR: output port too big! " << endl;
			cout << " OUTPUT select: " << out_port << endl;
			cout << " router radix: " << gN * (gK - 1) + gK << endl;
			exit(-1);
		}

		if (debug)
			cout << "        through output port : " << out_port << endl;
		if (gTrace)
		{
			cout << "Outport " << out_port << endl;
			cout << "Stop Mark" << endl;
		}
	}

	outputs->Clear();

	outputs->AddRange(out_port, vcBegin, vcEnd);
}
/*
void omni_war_hyperx(const Router *r, const Flit *f, int in_channel,
					 OutputSet *outputs, bool inject)
{
	// ( Traffic Class , Routing Order ) -> Virtual Channel Range
	int vcBegin = 0, vcEnd = gNumVCs - 1;
	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

	int out_port = -1;

	int nodo_destino = calculateRouter(f->dest);
	int nodo_actual = r->GetID();

	if (inject)
	{

		out_port = -1;
	}
	else if (nodo_destino == nodo_actual)
	{

		out_port = gN * (gK - 1) + calculateExitPort(f->dest);
	}
	else
	{

		if (in_channel >= (gK - 1) * gN)
		{ // inyección
			vcBegin = 0;
			vcEnd = 0;
		}
		else
		{

			vcBegin = f->vc + 1;
			vcEnd = f->vc + 1;
		}
		// assert(false);

		int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC); // esto es un apaño... xD
		int available_vcs = gNumVCs - vcEnd;									//+1 -1

		int missroute = available_vcs - distance_to_dest;

		assert(missroute > -1);
		vector<int> creditos = r->FreeCredits();

		int vcCredits_max = -1;
		int dimension_salida = -1;

		for (int i = 0; i < gN; i++)
		{

			int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

			if (salida != 0)
			{ // Si hay que recorrer esta salida...

				int puerto = i * (gK - 1) + salida - 1;
				int vcfree = getCreditOutportVC_hyperx(puerto, vcEnd, creditos);

				if (vcfree >= vcCredits_max)
				{
					vcCredits_max = vcfree;
					dimension_salida = i;
					out_port = puerto;
				}
			}
		}

		assert(dimension_salida != -1);

		if (missroute)
		{

			int vcCredits_max_missroute = -1;
			int out_port_missroute = -1;

			int prohibido = -(in_channel % (gK - 1)) + gK - 2 + in_channel - in_channel % (gK - 1);
			// printf("in channel: %d prohibido port: %d\n",in_channel, prohibido);
			for (int dim = 0; dim < gN; dim++)
			{

				int salida = (node_vectors[nodo_destino * gN + dim] - node_vectors[nodo_actual * gN + dim] + gK) % gK;

				if (salida != 0)
				{ // Si hay que recorrer esta salida...

					for (int i = 0; i <(gK - 1); i++)
					{
						int puerto = dim * (gK - 1) + i;
						int vcfree = getCreditOutportVC_hyperx(puerto, vcEnd, creditos);

						if (puerto != prohibido && vcfree > vcCredits_max_missroute)
						{ // es 1 hop mas....
							vcCredits_max_missroute = vcfree;
							out_port_missroute = puerto;
						}
					}
				}
			}

			if (r->GetUsedCredit(out_port) * distance_to_dest > r->GetUsedCredit(out_port_missroute) * (distance_to_dest + 1))
			{ // si realmente renta.... puede ser que se coja
				out_port = out_port_missroute;
				// printf("out: %d %d \n", out_port, vcCredits_max_missroute);	//un puerto minimo, pero da igual
			}
		}
		assert(out_port != -1);
	}

	outputs->Clear();

	outputs->AddRange(out_port, vcBegin, vcEnd);
}

void omni_war_priority_hyperx_fake(const Router *r, const Flit *f, int in_channel,
							  OutputSet *outputs, bool inject)
{
	// ( Traffic Class , Routing Order ) -> Virtual Channel Range
	int vcBegin = 0, vcEnd = gNumVCs - 1;
	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

	outputs->Clear();
	int out_port = -1;

	int nodo_destino = calculateRouter(f->dest);
	int nodo_actual = r->GetID();

	if (inject)
	{

		out_port = -1;

		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else if (nodo_destino == nodo_actual)
	{

		out_port = gN * (gK - 1) + calculateExitPort(f->dest);
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{

		int prohibido = -1;

		if (in_channel >= (gK - 1) * gN)
		{ // inyección
			vcBegin = 0;
			vcEnd = 0;
		}
		else
		{
			vcBegin = f->vc + 1;
			vcEnd = f->vc + 1;
			prohibido = -(in_channel % (gK - 1)) + gK - 2 + in_channel - in_channel % (gK - 1);
		}
		// assert(false);

		int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC); // esto es un apaño... xD
		int available_vcs = gNumVCs - vcEnd;									//+1 -1

		int missroute = available_vcs - distance_to_dest;

		assert(missroute > -1);

		for (int i = 0; i < gN; i++)
		{

			int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

			if (salida != 0)
			{ // Si hay que recorrer esta salida...

				int puerto_min = i * (gK - 1) + salida - 1;
				int prio = (buff_size - r->GetUsedCredit(puerto_min)) * (gN - distance_to_dest + 1); //+1 es por ser minimo, un salto menos...
				outputs->AddRange(puerto_min, vcBegin, vcEnd, prio);

				if (missroute > 1)
				{
					for (int k_salida = 0; k_salida < (gK - 2); k_salida++)
					{
						if (k_salida != salida - 1)
						{
							int puerto_miss = i * (gK - 1) + k_salida; // salida - 1
							if (r->GetUsedCredit(puerto_min) * distance_to_dest > r->GetUsedCredit(puerto_miss) * (distance_to_dest + 1))
							{
								int prio_miss = (buff_size - r->GetUsedCredit(puerto_miss)) * (gN - distance_to_dest);
								outputs->AddRange(puerto_miss, vcBegin, vcEnd, prio_miss);
							}
						}
					}
				}
			}
		}
	}
}*/

void adaptive_escape_hyperx(const Router *r, const Flit *f, int in_channel,
							OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	assert(gNumVCs >= 2);

	outputs->Clear();

	int out_port;

	if (inject)
	{

		out_port = -1;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
		else
		{

			vector<int> free_credits = r->FreeCredits();
			out_port = calculateDOR_routers(targetr, r->GetID());
			outputs->AddRange(out_port, vcEnd, vcEnd, 0);

			int nodo_destino = targetr;
			int nodo_actual = r->GetID();

			for (int i = 0; i < gN; i++)
			{

				int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

				if (salida != 0)
				{
					int puerto = i * (gK - 1) + salida - 1;
					outputs->AddRange(puerto, vcBegin, vcEnd - 1, getDataForVCRange(puerto * gNumVCs + vcBegin, puerto * gNumVCs + vcEnd, free_credits));
				}
			}
		}
	}
}

/*
void dal_hyperx_fake(const Router *r, const Flit *f, int in_channel,
				OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	assert(gNumVCs >= 2);

	outputs->Clear();

	int out_port;

	if (inject)
	{

		out_port = -1;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
		else
		{
			int avaliable_vcs = (vcEnd - vcBegin + 1) / 2;
			int nodo_destino = targetr;
			int nodo_actual = r->GetID();

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);
			out_port = calculateDOR_routers(targetr, r->GetID());

			int prio_dor = buff_size - r->GetUsedCredit(out_port); // con menos prioridad para que se coja menos...
			outputs->AddRange(out_port, vcBegin, vcBegin + avaliable_vcs - 1, prio_dor);

			if (f->vc >= avaliable_vcs)
			{ // adaptativo

				for (int i = 0; i < gN; i++)
				{

					int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

					if (salida != 0)
					{ // Si hay que recorrer esta salida...

						int puerto_min = i * (gK - 1) + salida - 1;

						int prio = (buff_size - r->GetUsedCredit(puerto_min)) * (gN - distance_to_dest + 1); //+1 es por ser minimo, un salto menos...
						outputs->AddRange(puerto_min, avaliable_vcs, vcEnd, prio);												   // minimo

						for (int k_salida = 0; k_salida < (gK - 2); k_salida++)
						{
							if (k_salida != salida - 1)
							{

								int puerto_miss = i * (gK - 1) + k_salida; // salida - 1
								int nodo_fuente = f->src / gC;
								int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

								if (condicion == 0 && r->GetUsedCredit(puerto_min) * i > r->GetUsedCredit(puerto_miss) * (i + 1) )
								{
									int prio_miss = (buff_size - r->GetUsedCredit(puerto_miss)) * (gN - distance_to_dest);
									outputs->AddRange(puerto_miss, avaliable_vcs, vcEnd, prio_miss); // miss
								}																	 // no se ha desviado aun en esta dimension
							}
						}
					}
				}
			}
		}
	}
}
*/
int valid_output(int score, int vector[])
{

	int lowest = vector[0];
	int index = 0;
	for (int i = 0; i < MAX_OUTPUTS; i++)
	{
		if (vector[i] < lowest)
		{
			lowest = vector[i];
			index = i;
		}
	}

	if (score > lowest)
	{
		return index;
	}
	return -1;
}

void omni_war_priority_hyperx(const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject)
{
	// ( Traffic Class , Routing Order ) -> Virtual Channel Range
	int vcBegin = 0, vcEnd = gNumVCs - 1;
	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

	outputs->Clear();

	int out_port = -1;

	int nodo_destino = calculateRouter(f->dest);
	// int nodo_fuente = calculateRouter(f->src);
	int nodo_actual = r->GetID();

	if (inject)
	{

		out_port = -1; // outport means the port that will be used to inject the flit

		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else if (nodo_destino == nodo_actual)
	{

		out_port = gN * (gK - 1) + calculateExitPort(f->dest);
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{

		int canales_por_salto = gNumVCs / (2 * gN); // lo seteo asi para hacer un missrouting por dimension
		if (in_channel >= (gK - 1) * gN)
		{				 // inyección
			vcBegin = 0; // vcBegin is the first VC of the injection channel
			vcEnd = canales_por_salto - 1;
			assert(f->hops == 0);
		}
		else
		{
			vcBegin = (f->hops) * canales_por_salto;
			vcEnd = vcBegin + canales_por_salto - 1;
		}

		int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC); // esto es un apaño... xD
		int available_vcs = (gNumVCs - vcEnd - 1) / canales_por_salto + 1;

		vector<int> free_credits = r->FreeCredits();

		int missroute = available_vcs - distance_to_dest;

		assert(missroute > -1);
		// printf("missroute: %d\n", missroute);
		for (int i = 0; i < gN; i++)
		{

			int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

			if (salida != 0)
			{ // Si hay que recorrer esta salida...

				int puerto_min = i * (gK - 1) + salida - 1;
				int prio_min = getPriority_omni(getDataForVCRange(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, free_credits), distance_to_dest, 0);

				// if (in_channel != puerto_min){
				outputs->AddRange(puerto_min, vcBegin, vcEnd, prio_min);
				//}

				if (missroute > 0)
				{

					for (int k_salida = 0; k_salida < (gK - 1); k_salida++)
					{
						if (k_salida != salida - 1)
						{

							int puerto_miss = i * (gK - 1) + k_salida;
							int prio_miss = getPriority_omni(getDataForVCRange(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits), distance_to_dest + 1, 1);
							// printf("prio_miss: %d, prio_min:%d \n", prio_miss, prio_min);
							if (in_channel != puerto_miss && prio_miss > prio_min)
							{
								// assert(false);
								outputs->AddRange(puerto_miss, vcBegin, vcEnd, prio_miss);
							}
						}
					}
				}
			}
		}
	}
}

int getPriority2_omni(int free_credits, int distance_to_dest, int missroute)
{ // int buff_size, int occupancy, int distance_to_dest

	int occ_bias = 8;

	int hop_bias = 2 * gN;
	int cycles = 2 + RandomInt(1);

	return occ_bias * 0.5;
}

int getDataForVCRange2(int vcBegin, int vcEnd, vector<int> &data)
{
	int suma = 0;
	for (int i = vcBegin; i <= vcEnd; i++)
	{
		suma += data[i];
	}
	return suma;
}

int getMeanDataExcluding(int vcBegin, int vcEnd, vector<int> &data, int puerto_min, int dim)
{
	int suma = 0;
	for (int i = 0; i < gK - 1; i++)
	{
		if (puerto_min != (gK - 1) * dim + i)
			suma += getDataForVCRange2(((gK - 1) * dim + i) * gNumVCs + vcBegin, ((gK - 1) * dim + i) * gNumVCs + vcEnd, data);
	}
	return suma / (gK - 2);
}

void omni_war_random_hyperx(const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject)
{
	// ( Traffic Class , Routing Order ) -> Virtual Channel Range
	int vcBegin = 0, vcEnd = gNumVCs - 1;
	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

	outputs->Clear();

	int out_port = -1;

	int nodo_destino = calculateRouter(f->dest);
	int nodo_fuente = calculateRouter(f->src);
	int nodo_actual = r->GetID();

	if (inject)
	{

		out_port = -1; // outport means the port that will be used to inject the flit

		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else if (nodo_destino == nodo_actual)
	{

		out_port = gN * (gK - 1) + calculateExitPort(f->dest);
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{
		assert((2 * gN <= gNumVCs));
		int canales_por_salto = gNumVCs / (2 * gN); // lo seteo asi para hacer un missrouting por dimension
		if (in_channel >= (gK - 1) * gN)
		{				 // inyección
			vcBegin = 0; // vcBegin is the first VC of the injection channel
			vcEnd = canales_por_salto - 1;
			assert(f->hops == 0);
		}
		else
		{
			vcBegin = (f->hops) * canales_por_salto;
			vcEnd = vcBegin + canales_por_salto - 1;
		}

		int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC); // esto es un apaño... xD
		int available_vcs = (gNumVCs - vcEnd - 1) / canales_por_salto + 1;

		vector<int> free_credits = r->FreeCredits();
		vector<int> ocupancy = r->UsedCredits();

		int missroute = available_vcs - distance_to_dest;

		assert(missroute > -1);
		// printf("missroute: %d\n", missroute);
		for (int i = 0; i < gN; i++)
		{

			int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

			if (salida != 0)
			{ // Si hay que recorrer esta salida...

				int puerto_min = i * (gK - 1) + salida - 1;
				int free_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, free_credits);
				int free_mean = getMeanDataExcluding(vcBegin, vcEnd, free_credits, puerto_min, i);

				int occ_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, ocupancy);
				int occ_mean = getMeanDataExcluding(vcBegin, vcEnd, ocupancy, puerto_min, i);

				int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

				if (missroute == 0 || condicion != 0)
				{

					outputs->AddRange(puerto_min, vcBegin, vcEnd, free_min + free_mean);
				}
				else
				{

					int k_salida = RandomInt(gK - 2);
					int puerto_miss = i * (gK - 1) + k_salida;
					int occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
					int free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
					int k_repuesto = -1;

					for (int t = 0; t < gK - 1; t++)
					{
						if ((puerto_miss != puerto_min) && (puerto_miss != in_channel))
							k_repuesto = k_salida;

						if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss))
							break;

						k_salida = (k_salida + 1) % (gK - 1);
						puerto_miss = i * (gK - 1) + k_salida;
						occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
						free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						// print a summary of the iteration
					}
					assert(k_repuesto != -1);

					if (puerto_min == puerto_miss) //avoid not being caught
					{
						puerto_miss = i * (gK - 1) + k_repuesto;
						occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
						free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
					}


					if ((occ_min - occ_mean - 3 * canales_por_salto >= 0))
					{
						
						outputs->AddRange(puerto_miss, vcBegin, vcEnd, free_miss + free_mean);

					}else if(occ_mean + canales_por_salto >= occ_min){

						outputs->AddRange(puerto_min, vcBegin, vcEnd, free_min + free_mean + 3 * canales_por_salto);

					} else if (occ_mean >= 5 *canales_por_salto ) //(occ_min == 8 * canales_por_salto && occ_mean < 7 * canales_por_salto)
					{ 
						outputs->AddRange(puerto_min, vcBegin, vcEnd, free_min + free_mean + 3 * canales_por_salto);
						outputs->AddRange(puerto_miss, vcBegin, vcEnd, free_miss + free_mean);
					}
					else
					{
						outputs->AddRange(puerto_min, vcBegin, vcEnd, free_min + free_mean); //free_min DENERIAAA
						//assert(false);
					}
				}
			}
		}
	}
}

void omni_dor_hyperx(const Router *r, const Flit *f, int in_channel, OutputSet *outputs, bool inject)
{
	// ( Traffic Class , Routing Order ) -> Virtual Channel Range
	int vcBegin = 0, vcEnd = gNumVCs - 1;
	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

	outputs->Clear();

	int out_port = -1;

	int nodo_destino = calculateRouter(f->dest);
	int nodo_fuente = calculateRouter(f->src);
	int nodo_actual = r->GetID();

	int canales_por_salto = gNumVCs - 1;

	if (inject)
	{
		out_port = -1;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else if (nodo_destino == nodo_actual)
	{

		out_port = gN * (gK - 1) + calculateExitPort(f->dest);
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{
		int dim_entrada = in_channel / (gK - 1);
		int missroute = 0;
		int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);

		if (dim_entrada < gN && node_vectors[nodo_fuente * gN + dim_entrada] != node_vectors[nodo_actual * gN + dim_entrada] && node_vectors[nodo_destino * gN + dim_entrada] != node_vectors[nodo_actual * gN + dim_entrada])
		{
			// dateline
			vcBegin = canales_por_salto;
			vcEnd = gNumVCs - 1;
		}
		else
		{
			missroute = 1;
			vcBegin = 0;
			vcEnd = canales_por_salto - 1;
		}

		/*
		if (in_channel >= (gK - 1) * gN)
		{ // inyección
			vcBegin = 0;
			vcEnd = canales_por_salto - 1;
		}
		else
		{
			vcBegin = canales_por_salto;
			vcEnd = vcBegin + canales_por_salto - 1;
		}
		*/

		vector<int> free_credits = r->FreeCredits();
		assert(missroute > -1);

		for (int i = 0; i < gN; i++)
		{

			int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

			if (salida != 0)
			{ // Si hay que recorrer esta salida...

				int puerto_min = i * (gK - 1) + salida - 1;
				int prio_min = getPriority_omni(getDataForVCRange(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, free_credits), distance_to_dest, 0);

				outputs->AddRange(puerto_min, vcBegin, vcEnd, prio_min);

				if (missroute > 0)
				{
					// assert(false);
					for (int k_salida = 0; k_salida < (gK - 1); k_salida++)
					{
						if (k_salida != salida - 1)
						{

							int puerto_miss = i * (gK - 1) + k_salida;
							int prio_miss = getPriority_omni(getDataForVCRange(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits), distance_to_dest + 1, 1);

							if (in_channel != puerto_miss && prio_miss > prio_min)
							{
								outputs->AddRange(puerto_miss, vcBegin, vcEnd, prio_miss);
							}
						}
					}
				}
				break; // ES DOR, SOLO EL PRIMERO!!
			}
		}
	}
}

void dal_vct_hyperx(const Router *r, const Flit *f, int in_channel,
					OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	assert(gNumVCs >= 2);

	outputs->Clear();

	int out_port = -1;

	if (inject)
	{

		out_port = -1;
		f->ph = 0;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
		else
		{

			int escape_vcs = 1;
			int nodo_destino = targetr;
			int nodo_actual = r->GetID();

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);
			int out_port_dor = calculateDOR_routers(targetr, r->GetID());

			int prio_dor = getPriority_dal(buff_size, r->GetUsedCredit(out_port_dor), distance_to_dest, 0, 1);
			outputs->AddRange(out_port_dor, vcBegin, vcBegin + escape_vcs - 1, prio_dor);

			// adaptativo
			for (int i = 0; i < gN; i++)
			{
				int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

				if (salida != 0)
				{ // Si hay que recorrer esta salida...

					int puerto_min = i * (gK - 1) + salida - 1;

					int prio_min = getPriority_dal(buff_size, r->GetUsedCredit(puerto_min), distance_to_dest, 0, 0);

					outputs->AddRange(puerto_min, escape_vcs, vcEnd, prio_min); // minimo , prio_min

					for (int k_salida = 0; k_salida < (gK - 1); k_salida++)
					{
						if (k_salida != salida - 1)
						{
							int puerto_miss = i * (gK - 1) + k_salida; // salida - 1
							int nodo_fuente = f->src / gC;
							int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;
							int prio_miss = getPriority_dal(buff_size, r->GetUsedCredit(puerto_miss), distance_to_dest + 1, 0, 0);

							if (condicion == 0 && prio_miss > prio_min)
							{
								outputs->AddRange(puerto_miss, escape_vcs, vcEnd, prio_miss);
							}
						}
					}
				}
			}
		}
	}
}

void dal_vct_random_hyperx(const Router *r, const Flit *f, int in_channel,
						   OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;
	//int vcTrueBegin = vcBegin;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	assert(gNumVCs >= 2);

	outputs->Clear();

	int out_port = -1;

	if (inject)
	{

		out_port = -1;
		f->ph = 0;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
		else
		{

			int escape_vcs = 1;
			int num_canales = gNumVCs - escape_vcs;
			int nodo_destino = targetr;
			int nodo_actual = r->GetID();
			int nodo_fuente = calculateRouter(f->src);

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);

			vector<int> free_credits = r->FreeCredits();
			vector<int> ocupancy = r->UsedCredits();

			
			int dor = 0;

			// adaptativo
			vcBegin += escape_vcs;
			for (int i = 0; i < gN; i++)
			{
				int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

				if (salida != 0)
				{ // Si hay que recorrer esta salida...

					int puerto_min = i * (gK - 1) + salida - 1;

					
					int free_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, free_credits);
					int free_mean = getMeanDataExcluding(vcBegin, vcEnd, free_credits, puerto_min, i);

					int occ_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, ocupancy);
					int occ_mean = getMeanDataExcluding(vcBegin, vcEnd, ocupancy, puerto_min, i);


					int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

					//int escape_prio = 1; //getDataForVCRange2(puerto_min * gNumVCs , puerto_min * gNumVCs , free_credits); //primer canal solo
					int hops_bias = f->hops +1;
					
					if (condicion == 0)
					{
						int k_salida = RandomInt(gK - 2);
						int puerto_miss = i * (gK - 1) + k_salida;
						int occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
						int free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						//int k_repuesto = -1;
						
						std::vector<int> salidas(gK - 1);
						std::iota(salidas.begin(), salidas.end(), 0);
						//shuffle with a fixed seed
						std::mt19937 g(RandomInt(INT_MAX));
						std::shuffle(salidas.begin(), salidas.end(), std::mt19937(std::random_device()()));
						
						/*for (int t = 0; t < gK - 1; t++)
						{
							k_salida = salidas[t];
							puerto_miss = i * (gK - 1) + k_salida;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
							
							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel))
								break;
						}*/

						/*
						for (int t = 0; t < gK - 1; t++)
						{
							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel))
								k_repuesto = k_salida;

							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss))
								break;

							k_salida = (k_salida + 1) % (gK - 1);
							puerto_miss = i * (gK - 1) + k_salida;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						}
						assert(k_repuesto != -1);

						if (puerto_min == puerto_miss)
						{
							puerto_miss = i * (gK - 1) + k_repuesto;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						}*/

						//print a summary of the occupancies
						//printf("occ_min: %d, occ_mean: %d, occ_miss: %d, free_min: %d, free_mean: %d, free_miss: %d\n", occ_min, occ_mean, occ_miss, free_min, free_mean, free_miss);


						/*if (occ_min - occ_mean - 3 * num_canales >= 0)
						{
							for (int t = 0; t < gK - 1; t++)
							{
								k_salida = salidas[t];
								puerto_miss = i * (gK - 1) + k_salida;
								occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
								free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
								
								if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss)){
									
										outputs->AddRange(puerto_miss, vcBegin, vcEnd,  free_miss  +1);
										max_miss++;
										if(max_miss == 3) break;
								}
									
							}
							
						}*/
						int max_miss = gK;
						if ((occ_min - occ_mean - 3 * num_canales >= 0))
						{
							for (int t = 0; t < gK - 1; t++)
							{
								k_salida = salidas[t];
								puerto_miss = i * (gK - 1) + k_salida;
								occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
								free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
								
								if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss)){
										outputs->AddRange(puerto_miss, vcBegin, vcEnd,  2);
										max_miss--;
										if(max_miss == 0) break;
								}	
							}
							//outputs->AddRange(puerto_miss, vcBegin, vcEnd, free_miss + free_mean +1);

						} else if(free_min >= num_canales){ //hay hueco minimo
							
							outputs->AddRange(puerto_min, vcBegin, vcEnd, 5); //(free_min + 3 * num_canales) +1
							
							for (int t = 0; t < gK - 1; t++)
							{
								k_salida = salidas[t];
								puerto_miss = i * (gK - 1) + k_salida;
								occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
								free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
								
								if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss)){
										outputs->AddRange(puerto_miss, vcBegin, vcEnd,  1);
										max_miss--;
										if(max_miss == 0) break;
								}	
							}

						} else if ( free_min < num_canales && free_mean > num_canales) //no hay hueco minimo pero si que hay hueco en media
						{ 
							outputs->AddRange(puerto_min, vcBegin, vcEnd, 4);//free_min + num_canales +1
							//outputs->AddRange(puerto_miss, vcBegin, vcEnd, free_miss + free_mean +1);
							for (int t = 0; t < gK - 1; t++)
							{
								k_salida = salidas[t];
								puerto_miss = i * (gK - 1) + k_salida;
								occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
								free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
								
								if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss)){
										outputs->AddRange(puerto_miss, vcBegin, vcEnd,  3); //free_miss  +1
										max_miss--;
										if(max_miss == 0) break;
								}	
							}
						}
						else //if(free_mean != 0)
						{
							outputs->AddRange(puerto_min, vcBegin, vcEnd, 2); 
							for (int t = 0; t < gK - 1; t++)
							{
								k_salida = salidas[t];
								puerto_miss = i * (gK - 1) + k_salida;
								occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
								free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
								
								if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss)){
										outputs->AddRange(puerto_miss, vcBegin, vcEnd, 1);
										max_miss--;
										if(max_miss == (gK-1)) break;
								}	
							}
						}/*else{

							outputs->AddRange(puerto_min, vcBegin, vcEnd, 6);
						
						}*/
					}
					else
					{

						outputs->AddRange(puerto_min, vcBegin, vcEnd, 6);

					}

					if(dor == 0){
						dor = 1;
						outputs->AddRange(puerto_min, 0, 0, 0);
					}
				}
			}
			/*printf("\n ======================================================== \n");
			set<OutputSet::sSetElement>::const_iterator iter = outputs->GetSet().begin( );

			while(iter!= outputs->GetSet().end( )){
				//print the details of the output set
				printf("port: %d vc start: %d vc end: %d prio: %d\n", iter->output_port, iter->vc_start, iter->vc_end, iter->pri);
				iter++;
			}*/
		}
	}
}

/*
void dal_vct_random_hyperx(const Router *r, const Flit *f, int in_channel,
						   OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;
	//int vcTrueBegin = vcBegin;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	assert(gNumVCs >= 2);

	outputs->Clear();

	int out_port = -1;

	if (inject)
	{

		out_port = -1;
		f->ph = 0;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
		else
		{

			int escape_vcs = 1;
			int num_canales = gNumVCs - escape_vcs;
			int nodo_destino = targetr;
			int nodo_actual = r->GetID();
			int nodo_fuente = calculateRouter(f->src);

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);

			vector<int> free_credits = r->FreeCredits();
			vector<int> ocupancy = r->UsedCredits();

			
			int dor = 0;

			// adaptativo
			vcBegin += escape_vcs;
			for (int i = 0; i < gN; i++)
			{
				int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

				if (salida != 0)
				{ // Si hay que recorrer esta salida...

					int puerto_min = i * (gK - 1) + salida - 1;

					
					int free_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, free_credits);
					int free_mean = getMeanDataExcluding(vcBegin, vcEnd, free_credits, puerto_min, i);

					int occ_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, ocupancy);
					int occ_mean = getMeanDataExcluding(vcBegin, vcEnd, ocupancy, puerto_min, i);


					int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

					//int escape_prio = 1; //getDataForVCRange2(puerto_min * gNumVCs , puerto_min * gNumVCs , free_credits); //primer canal solo
					int hops_bias = f->hops +1;
					
					if (condicion == 0)
					{
						int k_salida = RandomInt(gK - 2);
						int puerto_miss = i * (gK - 1) + k_salida;
						int occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
						int free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						//int k_repuesto = -1;
						
						std::vector<int> salidas(gK - 1);
						std::iota(salidas.begin(), salidas.end(), 0);
						//shuffle with a fixed seed
						std::mt19937 g(RandomInt(INT_MAX));
						std::shuffle(salidas.begin(), salidas.end(), std::mt19937(std::random_device()()));
						
						/*for (int t = 0; t < gK - 1; t++)
						{
							k_salida = salidas[t];
							puerto_miss = i * (gK - 1) + k_salida;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
							
							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel))
								break;
						}*/

						/*
						for (int t = 0; t < gK - 1; t++)
						{
							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel))
								k_repuesto = k_salida;

							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss))
								break;

							k_salida = (k_salida + 1) % (gK - 1);
							puerto_miss = i * (gK - 1) + k_salida;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						}
						assert(k_repuesto != -1);

						if (puerto_min == puerto_miss)
						{
							puerto_miss = i * (gK - 1) + k_repuesto;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						} /*

						//print a summary of the occupancies
						//printf("occ_min: %d, occ_mean: %d, occ_miss: %d, free_min: %d, free_mean: %d, free_miss: %d\n", occ_min, occ_mean, occ_miss, free_min, free_mean, free_miss);
						
						int max_miss = 0;

						

						/*if (occ_min - occ_mean - 3 * num_canales >= 0)
						{
							for (int t = 0; t < gK - 1; t++)
							{
								k_salida = salidas[t];
								puerto_miss = i * (gK - 1) + k_salida;
								occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
								free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
								
								if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss)){
									
										outputs->AddRange(puerto_miss, vcBegin, vcEnd,  free_miss  +1);
										max_miss++;
										if(max_miss == 3) break;
								}
									
							}
							
						}*/ 
						/*if (free_mean >= num_canales)
						{ // 
						if(free_min > 0)
							outputs->AddRange(puerto_min, vcBegin, vcEnd, (free_min + 3 * num_canales) +1 );
						else
							outputs->AddRange(puerto_min, vcBegin, vcEnd, 1 );
						//outputs->AddRange(puerto_min, vcBegin, vcEnd, (free_min + num_canales) +1);
						//outputs->AddRange(puerto_miss, vcBegin, vcEnd,  (free_miss + free_mean)  +1);

						//if (free_mean > num_canales){
							for (int t = 0; t < gK - 1; t++)
							{
								k_salida = salidas[t];
								puerto_miss = i * (gK - 1) + k_salida;
								occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
								free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
								
								if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss)){
										outputs->AddRange(puerto_miss, vcBegin, vcEnd,  free_miss  +1);
										max_miss++;
										if(max_miss == gK) break;
								}	
							}
					//	}
						/*}else{
							outputs->AddRange(puerto_min, vcBegin, vcEnd, (free_min + 3 * num_canales) +1 );
						}/*
					}
					else
					{
						if(free_min > 0)
							outputs->AddRange(puerto_min, vcBegin, vcEnd, (free_min + 3 * num_canales) +1 );
						else
							outputs->AddRange(puerto_min, vcBegin, vcEnd, 2 );
					}

					if(dor == 0){
						dor = 1;
						outputs->AddRange(puerto_min, 0, 0, 0);
					}
				}
			}
			/*printf("\n ======================================================== \n");
			set<OutputSet::sSetElement>::const_iterator iter = outputs->GetSet().begin( );

			while(iter!= outputs->GetSet().end( )){
				//print the details of the output set
				printf("port: %d vc start: %d vc end: %d prio: %d\n", iter->output_port, iter->vc_start, iter->vc_end, iter->pri);
				iter++;
			}/*
		}
	}
}*/


void dal_vct_random2_hyperx(const Router *r, const Flit *f, int in_channel,
						   OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;
	//int vcTrueBegin = vcBegin;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	assert(gNumVCs >= 2);

	outputs->Clear();

	int out_port = -1;

	if (inject)
	{

		out_port = -1;
		f->ph = 0;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
		else
		{

			int escape_vcs = 1;
			int num_canales = gNumVCs - escape_vcs;
			int nodo_destino = targetr;
			int nodo_actual = r->GetID();
			int nodo_fuente = calculateRouter(f->src);

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);

			vector<int> free_credits = r->FreeCredits();
			vector<int> ocupancy = r->UsedCredits();

			
			int dor = 0;

			// adaptativo
			vcBegin += escape_vcs;
			int free_mean_total = 0;
			int free_min_total = 0;
			
			for (int i = 0; i < gN; i++)
			{
				int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;
				if (salida != 0)
				{
					int puerto_min = i * (gK - 1) + salida - 1;
					free_mean_total += getMeanDataExcluding(vcBegin, vcEnd, free_credits, -1, i);

					int free_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, free_credits);
					free_min_total += free_min;
				}
			}

			free_mean_total = free_mean_total / distance_to_dest;
			free_min_total = free_min_total / distance_to_dest;
			//add dor escpape route
			int out_port_dor = calculateDOR_routers(targetr, r->GetID());
			outputs->AddRange(out_port_dor, 0, 0, 0);

			
			for (int i = 0; i < gN; i++)
			{
				int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

				if (salida != 0)
				{ // Si hay que recorrer esta salida...

					int puerto_min = i * (gK - 1) + salida - 1;
					
					int free_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, free_credits);
					int free_mean = getMeanDataExcluding(vcBegin, vcEnd, free_credits, puerto_min, i);

					int occ_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, ocupancy);
					int occ_mean = getMeanDataExcluding(vcBegin, vcEnd, ocupancy, puerto_min, i);


					int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

					//int escape_prio = 1; //getDataForVCRange2(puerto_min * gNumVCs , puerto_min * gNumVCs , free_credits); //primer canal solo
					int hops_bias = f->hops +1;
					
					if (condicion == 0)
					{
						int k_salida = RandomInt(gK - 2);
						int puerto_miss = i * (gK - 1) + k_salida;
						int occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
						int free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						int k_repuesto = -1;

						for (int t = 0; t < gK - 1; t++)
						{
							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel))
								k_repuesto = k_salida;

							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss))
								break;

							k_salida = (k_salida + 1) % (gK - 1);
							puerto_miss = i * (gK - 1) + k_salida;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						}
						assert(k_repuesto != -1);

						if (puerto_min == puerto_miss)
						{
							puerto_miss = i * (gK - 1) + k_repuesto;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						}


						if ((occ_min - occ_mean - 3 * num_canales >= 0))
						{
							
							outputs->AddRange(puerto_miss, vcBegin, vcEnd, free_mean* hops_bias* (free_miss) +1);

						}else if(occ_mean + num_canales > occ_min)
						{

							outputs->AddRange(puerto_min, vcBegin, vcEnd, free_mean* hops_bias* (free_min + (3 * num_canales)) +1);
							outputs->AddRange(puerto_miss, vcBegin, vcEnd, free_mean* hops_bias* (free_miss) +1);
						} else if (occ_mean >= 5 *num_canales ) //(occ_min == 8 * canales_por_salto && occ_mean < 7 * canales_por_salto)
						{ 
							outputs->AddRange(puerto_min, vcBegin, vcEnd, free_mean* hops_bias* (free_min + ( 3 * num_canales)) +1);
							outputs->AddRange(puerto_miss, vcBegin, vcEnd, free_mean* hops_bias* (free_miss )  +1);
						}
						else
						{
							outputs->AddRange(puerto_min, vcBegin, vcEnd, free_mean* hops_bias* (free_min + (3 * num_canales)) +1); //free_min DENERIAAA
							//assert(false);
						}
							
					}
					else
					{
						outputs->AddRange(puerto_min, vcBegin, vcEnd, free_mean* hops_bias* (free_min + (3 * num_canales)) + 1);
					}

					if(dor == 0){
						dor = 1;
						outputs->AddRange(puerto_min, 0, 0, 0);
					}
				}
			}
		}
	}
}

void dal_wormhole_hyperx(const Router *r, const Flit *f, int in_channel,
						 OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	assert(gNumVCs >= 2);

	outputs->Clear();

	int out_port = -1;
	int escape_vcs = 1; // (vcEnd - vcBegin + 1) / 2; //

	if (inject)
	{

		out_port = -1;
		f->ph = 0;
		outputs->AddRange(out_port, vcBegin + escape_vcs, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, escape_vcs, vcEnd);
		}
		else
		{
			if (f->vc == vcBegin)
			{
				f->ph = 1; // Minimal to Destination
			}

			int nodo_destino = targetr;
			int nodo_actual = r->GetID();

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);
			int out_port_dor = calculateDOR_routers(targetr, r->GetID());

			int prio_dor = getPriority_dal(buff_size, r->GetUsedCredit(out_port_dor), distance_to_dest, 0, 1);
			outputs->AddRange(out_port_dor, vcBegin, vcBegin + escape_vcs - 1, prio_dor);

			for (int i = 0; i < gN; i++)
			{
				int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

				if (salida != 0)
				{ // Si hay que recorrer esta salida...

					int puerto_min = i * (gK - 1) + salida - 1;

					int prio_min = getPriority_dal(buff_size, r->GetUsedCredit(puerto_min), distance_to_dest, 0, 0);

					outputs->AddRange(puerto_min, escape_vcs, vcEnd, prio_min); // minimo , prio_min

					if (f->ph == 0)
					{
						for (int k_salida = 0; k_salida < (gK - 1); k_salida++)
						{
							if (k_salida != salida - 1)
							{
								int puerto_miss = i * (gK - 1) + k_salida; // salida - 1
								int nodo_fuente = f->src / gC;
								int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;
								int prio_miss = getPriority_dal(buff_size, r->GetUsedCredit(puerto_miss), distance_to_dest + 1, 0, 0);

								if (condicion == 0 && prio_miss > prio_min)
								{
									outputs->AddRange(puerto_miss, escape_vcs, vcEnd, prio_miss);
								}
							}
						}
					}
				}
			}
		}
	}
}

void minimal_turn_model_hyperx(const Router *r, const Flit *f, int in_channel,
							   OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));

	outputs->Clear();

	int out_port = -1;

	if (inject)
	{
		f->ph = 0;
		out_port = -1;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
		else
		{

			int nodo_destino = targetr;
			int nodo_actual = r->GetID();

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);
			int fase = f->ph;

			if (fase == 0)
			{
				fase = 1;
				for (int i = 0; i < gN; i++)
				{
					int salida = node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i];

					if (salida < 0) // Si hay que recorrer esta salida...
					{
						fase = 0;
						int puerto_min = i * (gK - 1) + (salida + gK) - 1;
						int prio_min = getPriority_dal(buff_size, r->GetUsedCredit(puerto_min), distance_to_dest, 0, 0);
						outputs->AddRange(puerto_min, vcBegin, vcEnd, prio_min);
					}
				}
			}

			f->ph = fase;

			if (fase == 1)
			{

				for (int i = 0; i < gN; i++)
				{
					int salida = node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i];

					if (salida > 0) // Si hay que recorrer esta salida...
					{
						fase = 0;
						int puerto_min = i * (gK - 1) + salida - 1;
						int prio_min = getPriority_dal(buff_size, r->GetUsedCredit(puerto_min), distance_to_dest, 0, 0);
						outputs->AddRange(puerto_min, vcBegin, vcEnd, prio_min);
					}
				}
			}
		}
	}
}

void missrouting_turn_model_hyperx(const Router *r, const Flit *f, int in_channel,
								   OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	// assert(gNumVCs >= 2);

	outputs->Clear();

	int out_port = -1;

	if (inject)
	{
		f->ph = 0;
		out_port = -1;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
		else
		{
			// assert(false);
			int nodo_destino = targetr;
			int nodo_actual = r->GetID();
			int nodo_fuente = calculateRouter(f->src);

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);
			int fase = f->ph;

			if (fase == 0)
			{
				fase = 1;
				for (int i = 0; i < gN; i++)
				{
					int salida = node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i];

					if (salida < 0) // Si actual es mayor que salida.....
					{
						fase = 0;
						int index = (salida + gK);
						int puerto_min = i * (gK - 1) + index - 1;
						int prio_min = getPriority_dal(buff_size, r->GetUsedCredit(puerto_min), distance_to_dest, 0, 0);
						outputs->AddRange(puerto_min, vcBegin, vcEnd, prio_min);

						int condicion = node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i];

						if (condicion == 0)
						{ // si ya se ha hecho missrouting....

							/*for(int j = index+1; j < (gK-1); j++){

								int puerto_miss = i * (gK - 1) + j;
								int prio_miss = getPriority_dal(buff_size, r->GetUsedCredit(puerto_miss), distance_to_dest+1, 0, 0);
								if(j != index- 1 && prio_miss > prio_min){
									outputs->AddRange(puerto_miss, vcBegin, vcEnd, prio_miss);
								}
							}

							for(int j = index-1; j >= (index + salida) && j >= 0; j--){ //salida is negative number

								int puerto_miss = i * (gK - 1) + j;
								int prio_miss = getPriority_dal(buff_size, r->GetUsedCredit(puerto_miss), distance_to_dest+1, 0, 0);
								if(j != index- 1 && prio_miss > prio_min){
									outputs->AddRange(puerto_miss, vcBegin, vcEnd, prio_miss);
								}
							}*/

							int start = gK - node_vectors[nodo_actual * gN + i] - 1; // el menos 1 es porque empezamos a contar en 0...
							// del enlace mínimo a gk-1 es el missrouting que te saca de la fase....
							for (int j = start; j < (gK - 1); j++)
							{
								int puerto_miss = i * (gK - 1) + j;
								int prio_miss = getPriority_dal(buff_size, r->GetUsedCredit(puerto_miss), distance_to_dest + 1, 0, 0);
								if (j != index - 1 && prio_miss > prio_min)
								{
									outputs->AddRange(puerto_miss, vcBegin, vcEnd, prio_miss);
								}
							}
						}
					}
				}
			}

			f->ph = fase;

			if (fase == 1)
			{

				int dim_entrada = in_channel / (gK - 1);
				int alineado_fuente = node_vectors[nodo_fuente * gN + dim_entrada] - node_vectors[nodo_actual * gN + dim_entrada];

				for (int i = 0; i < gN; i++)
				{
					int salida = node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i];

					if (salida > 0) // Si actual es menor que salida
					{
						fase = 0;
						int puerto_min = i * (gK - 1) + salida - 1;
						int prio_min = getPriority_dal(buff_size, r->GetUsedCredit(puerto_min), distance_to_dest, 0, 0);
						outputs->AddRange(puerto_min, vcBegin, vcEnd, prio_min);

						int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;
						int index = salida;

						if (condicion == 0)
						{

							for (int j = 0; j < index - 1; j++) // entre actual y salida
							{
								int puerto_miss = i * (gK - 1) + j;
								int prio_miss = getPriority_dal(buff_size, r->GetUsedCredit(puerto_miss), distance_to_dest + 1, 0, 0);
								if (prio_miss > prio_min)
								{
									outputs->AddRange(puerto_miss, vcBegin, vcEnd, prio_miss);
								}
							}
						}
					}
				}
			}
		}
	}
}

void add_escape_turn_model(OutputSet *outputs, const Router *r, const Flit *f, int nodo_destino,
						   int nodo_actual, int vcBegin, int vcEnd, int dim, int puerto_min, int prio_min, int missrouting)
{	

	vector<int> free_credits = r->FreeCredits();
	//int free_min = getDataForVCRange2(puerto_min + vcBegin, puerto_min + vcEnd,free_credits);
	

	if(missrouting && link_restrictions_hyperx[nodo_actual][nodo_destino].size() > 0){
		//get a random element of the list
		int random_element = RandomInt(link_restrictions_hyperx[nodo_actual][nodo_destino].size() -1);
		
		int nodo_intermedio = link_restrictions_hyperx[nodo_actual][nodo_destino][random_element];
		
		int salida_miss = (node_vectors[nodo_actual * gN + dim] - node_vectors[nodo_intermedio * gN + dim] + gK) % gK;
		assert(nodo_actual != nodo_intermedio);
		int puerto_miss = dim * (gK - 1) + salida_miss - 1;
		//int free_miss = getDataForVCRange2(puerto_miss + vcBegin, puerto_miss + vcEnd, free_credits);
		

		outputs->AddRange(puerto_miss, vcBegin, vcEnd, missrouting - 1); //max(free_miss - 3 * (vcEnd - vcBegin +1)
		assert(puerto_min != puerto_miss);
		if(prio_min != 0){
			outputs->AddRange(puerto_min, vcBegin, vcEnd, prio_min -1);
		}
		
		
	}else{
		outputs->AddRange(puerto_min, vcBegin, vcEnd, prio_min -1);
	}

}

void dal_vct_turned_hyperx(const Router *r, const Flit *f, int in_channel,
						   OutputSet *outputs, bool inject)
{

	
	int vcBegin = 0, vcEnd = gNumVCs - 1;
	int vcTrueBegin = vcBegin;

	if (f->type == Flit::READ_REQUEST)
	{
		vcBegin = gReadReqBeginVC;
		vcEnd = gReadReqEndVC;
	}
	else if (f->type == Flit::WRITE_REQUEST)
	{
		vcBegin = gWriteReqBeginVC;
		vcEnd = gWriteReqEndVC;
	}
	else if (f->type == Flit::READ_REPLY)
	{
		vcBegin = gReadReplyBeginVC;
		vcEnd = gReadReplyEndVC;
	}
	else if (f->type == Flit::WRITE_REPLY)
	{
		vcBegin = gWriteReplyBeginVC;
		vcEnd = gWriteReplyEndVC;
	}
	assert(((f->vc >= vcBegin) && (f->vc <= vcEnd)) || (inject && (f->vc < 0)));
	assert(gNumVCs >= 2);

	outputs->Clear();

	int out_port = -1;

	if (inject)
	{

		out_port = -1;
		f->ph = 0;
		outputs->AddRange(out_port, vcBegin, vcEnd);
	}
	else
	{ // si no se inyecta

		int dest = f->dest;
		int targetr = (int)(dest / gC);

		if (targetr == r->GetID())
		{ // if we are at the final router, yay, output to client
			out_port = gN * (gK - 1) + calculateExitPort(dest);
			outputs->AddRange(out_port, vcBegin, vcEnd);
		}
		else
		{

			int escape_vcs = 1;
			int num_canales = gNumVCs - escape_vcs;
			int nodo_destino = targetr;
			int nodo_actual = r->GetID();
			int nodo_fuente = calculateRouter(f->src);

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC);

			vector<int> free_credits = r->FreeCredits();
			vector<int> ocupancy = r->UsedCredits();

			
			int dor = 0;

			// adaptativo
			vcBegin += escape_vcs;
			for (int i = 0; i < gN; i++)
			{
				int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

				if (salida != 0)
				{ // Si hay que recorrer esta salida...

					int puerto_min = i * (gK - 1) + salida - 1;

					
					int free_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, free_credits);
					int free_mean = getMeanDataExcluding(vcBegin, vcEnd, free_credits, puerto_min, i);

					int occ_min = getDataForVCRange2(puerto_min * gNumVCs + vcBegin, puerto_min * gNumVCs + vcEnd, ocupancy);
					int occ_mean = getMeanDataExcluding(vcBegin, vcEnd, ocupancy, puerto_min, i);


					int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

					int escape_prio = getDataForVCRange2(puerto_min * gNumVCs , puerto_min * gNumVCs , free_credits); //primer canal solo
					int escape_miss = 0;

					if (condicion == 0)
					{
						int k_salida = RandomInt(gK - 2);
						int puerto_miss = i * (gK - 1) + k_salida;
						int occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
						int free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						int k_repuesto = -1;

						for (int t = 0; t < gK - 1; t++)
						{
							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel))
								k_repuesto = k_salida;

							if ((puerto_miss != puerto_min) && (puerto_miss != in_channel) && ((occ_mean) >= occ_miss))
								break;

							k_salida = (k_salida + 1) % (gK - 1);
							puerto_miss = i * (gK - 1) + k_salida;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						}
						assert(k_repuesto != -1);

						if (puerto_min == puerto_miss)
						{
							puerto_miss = i * (gK - 1) + k_repuesto;
							occ_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, ocupancy);
							free_miss = getDataForVCRange2(puerto_miss * gNumVCs + vcBegin, puerto_miss * gNumVCs + vcEnd, free_credits);
						}

						//=========== THE ROUTING COMPUTATION ===========

						if ((occ_min - occ_mean - 3 * num_canales >= 0))
						{
							outputs->AddRange(puerto_miss, vcBegin, vcEnd, free_miss + free_mean +1 ); //missrouting
							escape_prio= 0;
							escape_miss= free_miss + free_mean;
						}
						else if ((occ_min == 8 * num_canales && occ_mean < 7 * num_canales))
						{ 
							outputs->AddRange(puerto_min, vcBegin, vcEnd, free_mean + 3 * num_canales +1); //idk
							outputs->AddRange(puerto_miss, vcBegin, vcEnd, free_miss + free_mean +1);
							escape_prio= free_mean + 3 * num_canales;
							escape_miss =free_miss + free_mean;
						}
						else
						{
							outputs->AddRange(puerto_min, vcBegin, vcEnd, free_min + free_mean+ 1); //min
							escape_prio= free_min + free_mean;
						
						}
					}
					else
					{

						outputs->AddRange(puerto_min, escape_vcs, vcEnd, free_min + free_mean +1 ); //min
						escape_prio= free_min + free_mean;
					}

					if(dor == 0){
						dor = 1;
						//outputs->AddRange(puerto_min, vcBegin - escape_vcs, vcBegin - 1, escape_prio);
						add_escape_turn_model(outputs, r, f, nodo_destino, nodo_actual, 0, 0, i, puerto_min, escape_prio, escape_miss);
					}
				}

			}
		}
	}
}