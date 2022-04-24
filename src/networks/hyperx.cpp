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

#include "random_utils.hpp"
#include "misc_utils.hpp"
#include "globals.hpp"

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
	gRoutingFunctionMap["O1_turn_hyperx"] = &adaptive_xyyx_hyperx;
	gRoutingFunctionMap["adaptive_dor_exit_hyperx"] = &adaptive_dor_exit_hyperx; //wrong manner to approach channel escape
	gRoutingFunctionMap["adaptive_escape_hyperx"] = &adaptive_escape_hyperx;
	gRoutingFunctionMap["valiant_hyperx"] = &valiant_hyperx;
	gRoutingFunctionMap["adaptive_escalera_hyperx"] = &adaptive_escalera_hyperx;
	gRoutingFunctionMap["ugal_hyperx"] = &ugal_hyperx;

	gRoutingFunctionMap["omni_war_hyperx"] = &omni_war_hyperx;					 // &omni_war_priority_hyperx;
	gRoutingFunctionMap["omni_war_priority_hyperx"] = &omni_war_priority_hyperx; // &omni_war_priority_hyperx;
	gRoutingFunctionMap["dal_hyperx"] = &dal_hyperx;
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
	int i = 0;
	int dim_help = -1;

	for (i = 0; i < gN; i++)
	{

		int temp = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

		if (temp != 0)
		{
			salida = temp;
			dim_help = i;
		}
	}

	// printf("%d, %d\n", nodo_destino,nodo_actual);
	// fflush(stdout);
	assert(salida != 0 && dim_help != -1); // esto se puede cambiar, meter aqui la salida tmbn

	return dim_help * (gK - 1) + salida - 1; // esto es lo normalizado.
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
	int _dim = gN;

	int src_tmp = (int)src / gC;
	int dest_tmp = (int)dest / gC;

	for (int d = 0; d < _dim; d++)
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
	}

	outputs->Clear();

	outputs->AddRange(out_port, vcBegin, vcEnd);
}

// ZYX or XYZ routing...
void adaptive_xyyx_hyperx(const Router *r, const Flit *f, int in_channel,
						  OutputSet *outputs, bool inject)
{

	assert(gNumVCs>=gN);
	int vcBegin = 0, vcEnd = gNumVCs - 1;
	int available_vcs = gNumVCs / 2;

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
	else
	{

		int dimension_salida = -1;

		if (in_channel >= gN * (gK - 1)) // inyeccion
		{ 

			int out_xy = calculateDOR_routers(nodo_destino, nodo_actual);
			int out_yx = calculateDORYX_routers(nodo_destino, nodo_actual);

			if(r->GetUsedCredit(out_xy) <= r->GetUsedCredit(out_yx) ){

				dimension_salida = 0;
			
			}else{
				
				dimension_salida = 1;
			
			}
			
		}
		else
		{ 
			dimension_salida = f->vc / available_vcs;
			//printf("dimension_salida: %d \n", dimension_salida);
		}

		if (dimension_salida == 1)
		{
			
			out_port = calculateDORYX_routers(nodo_destino, nodo_actual);
			vcBegin = available_vcs;

		}
		else
		{
			
			out_port = calculateDOR_routers(nodo_destino, nodo_actual);
			vcEnd = available_vcs -1;


		}

		
		assert(out_port != -1); // no deberia ser...

	}

	outputs->Clear();

	outputs->AddRange(out_port, vcBegin, vcEnd);
}

void adaptive_escalera_hyperx(const Router *r, const Flit *f, int in_channel,
							  OutputSet *outputs, bool inject)
{

	int vcBegin = 0, vcEnd = gNumVCs - 1;

	int out_port = -1;
	
	assert(gNumVCs>=gN);

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
	else
	{

		int dimension_salida = -1;
		int min_occupancy = INT_MAX;
		int available_vcs = gNumVCs / gN;

		for (int i = 0; i < gN; i++)
		{

			int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

			if (salida != 0)
			{ // Si hay que recorrer esta salida...

				int occupancy = 0;
				int puerto = i * (gK - 1) + salida - 1;

				occupancy = r->GetUsedCredit(puerto);

				if (occupancy < min_occupancy)
				{
					min_occupancy = occupancy;
					dimension_salida = i;
					out_port = puerto; // esto es lo normalizado.
				}
			}
		}

		assert(out_port != -1); // no deberia ser...

		if (in_channel >= gN * (gK - 1))
		{ // inyeccion
			vcBegin = 0;
			vcEnd = available_vcs-1;
		}
		else
		{
			int d = f->vc / available_vcs;
			vcBegin = d*available_vcs;
			vcEnd = vcBegin + available_vcs -1;
		}
	}

	outputs->Clear();

	outputs->AddRange(out_port, vcBegin, vcEnd);
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

		if (in_channel >= (gK - 1) * gN)
		{
			f->ph = 0;
			f->intm = RandomInt(powi(gK, gN) - 1); // apuntamos a un router de la red... es mas comodo.
		}

		int nodo_actual = r->GetID();
		int nodo_intermedio = f->intm;
		int nodo_destino = calculateRouter(f->dest); // aqui sacamos todos los routers de la red.

		if (nodo_intermedio == nodo_actual || nodo_destino == nodo_actual)
		{
			f->ph = 1;
		}

		// each class must have at least 2 vcs assigned or else valiant valiant will deadlock
		int const available_vcs = (vcEnd - vcBegin + 1) / 2;
		assert(available_vcs > 0);

		if (nodo_destino != nodo_actual)
		{

			if (f->ph == 0)
			{
				out_port = calculateDOR_routers(nodo_intermedio, nodo_actual);
				vcEnd -= available_vcs;
			}
			else
			{
				assert(f->ph == 1);
				out_port = calculateDOR_routers(nodo_destino, nodo_actual);
				vcBegin += available_vcs;
			}
		}
		else
		{
			out_port = gN * (gK - 1) + calculateExitPort(f->dest); // sacamos el outpor exacto
			assert(f->ph == 1);
			vcBegin += available_vcs;
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
		int _concentration = gC;
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
		// printf("gNumVCs: %d\n", gNumVCs);
		// printf("vEnd: %d\n", vcEnd);
		// printf("avaliable_vcs: %d\n", available_vcs);
		// printf("distance_to_dest: %d\n", distance_to_dest);
		int missroute = available_vcs - distance_to_dest;
		// printf("missroute: %d\n", missroute);
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

					for (int i = 0; i < (gK - 2); i++)
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
			// vcCredits_max *(1.5) < vcCredits_max_missroute
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

void omni_war_priority_hyperx(const Router *r, const Flit *f, int in_channel,
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

		int prohibido = -(in_channel % (gK - 1)) + gK - 2 + in_channel - in_channel % (gK - 1);

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
}

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
			out_port = calculateDOR_routers(targetr, r->GetID());
			outputs->AddRange(out_port, vcEnd, vcEnd, 0);

			// if (f->vc != vcEnd)
			//dsf
			//{ // adaptativo

			int nodo_destino = targetr;
			int nodo_actual = r->GetID();

			for (int i = 0; i < gN; i++)
			{

				int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

				if (salida != 0)
				{ // Si hay que recorrer esta salida...

					int puerto = i * (gK - 1) + salida - 1;
					outputs->AddRange(puerto, vcBegin, vcEnd - 1, 1);
				}
			}
			//}
		}
	}
}

void dal_hyperx(const Router *r, const Flit *f, int in_channel,
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
		outputs->AddRange(out_port, vcBegin + 1, vcEnd);
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
			int avaliable_vcs = (vcEnd - vcBegin + 1)/2;
			int nodo_destino = targetr;
			int nodo_actual = r->GetID();

			int distance_to_dest = find_distance_hyperx(f->dest, nodo_actual * gC); 
			out_port = calculateDOR_routers(targetr, r->GetID());

			int prio_dor = buff_size - r->GetUsedCredit(out_port); // con menos prioridad para que se coja menos...
			outputs->AddRange(out_port, vcBegin, vcBegin+ avaliable_vcs -1, prio_dor);

			if (f->vc >= avaliable_vcs)
			{ // adaptativo

			
				for (int i = 0; i < gN; i++)
				{

					int salida = (node_vectors[nodo_destino * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

					if (salida != 0)
					{ // Si hay que recorrer esta salida...

						int puerto_min = i * (gK - 1) + salida - 1;
						
						int prio = (buff_size - r->GetUsedCredit(puerto_min)) * (gN - distance_to_dest + 1); //+1 es por ser minimo, un salto menos...
						outputs->AddRange(puerto_min, avaliable_vcs, vcEnd, prio); // minimo

						for (int k_salida = 0; k_salida < (gK - 2); k_salida++)
						{
							if (k_salida != salida - 1)
							{

								int puerto_miss = i * (gK - 1) + k_salida; // salida - 1
								int nodo_fuente = f->src / gC;
								int condicion = (node_vectors[nodo_fuente * gN + i] - node_vectors[nodo_actual * gN + i] + gK) % gK;

								if (condicion == 0 && r->GetUsedCredit(puerto_min) * i > r->GetUsedCredit(puerto_miss) * (i + 1))
								{
									int prio_miss = (buff_size - r->GetUsedCredit(puerto_miss)) * (gN - distance_to_dest);
									outputs->AddRange(puerto_miss, avaliable_vcs, vcEnd, prio_miss); // miss
								}													  // no se ha desviado aun en esta dimension
							}
						}
					}
				}
			}
		}
	}
}