#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include "erros_code.hpp"

using namespace std;

class Graph_Cycle
{
		int n, m;
		vector<vector<int> > cost, flow, cap;
		vector<int> b;

		const int INF;

		bool has_negative_edge(){
				return false;
		}

		bool has_parallel_edge(){
				return false;
		}

		bool is_valid_inbalance(){
				int sum = 0;
				for (int i = 0; i < n; ++i)
						sum += b[i];
				return sum == 0;		
		}

		vector<int> find_cycle(int start, const vector<int>& prev){
				vector<int> cycle;
				vector<bool> mark(n+1, false);

				int u = prev[start];
				while (u != start){
						mark[u] = true;
						assert(prev[u] != -1);
						if (mark[prev[u]]) break;
						u = prev[u];
				}

				start = u;
				u = prev[u];
				while (u != start){
						cycle.push_back(u);
						u = prev[u];	
				}
				cycle.push_back(start);
				reverse(cycle.begin(), cycle.end());
				return cycle;
		}

		int get_edge_cap(int u, int v)
		{
			return cap[u][v] - flow[u][v];
		}

		int find_bottleneck_cycle(const vector<int>& cycle)
		{
				int bottleneck = INF;
				int ncycle = cycle.size();
				for (int i = 0; i < ncycle; ++i)
						bottleneck = min(bottleneck, get_edge_cap(cycle[i], cycle[(i+1)%ncycle]));
				return bottleneck;
		}


		// trys to improve flow by finding a negative 
		// cost cycle and passing flow through it.
		bool improve_flow(){
				// creates a source vertice with a
				// zero cost and zero capacity
				// edge for every other vertex.
				// This way we can run bellman just one time
				int source = n;
				for (int i = 0; i < n; ++i)
						set_edge(source, i, 0, 0);

				vector<int> dist(n + 1, INF);
				vector<int> prev(n + 1, -1);

				dist[source] = 0;

				// bellman-ford
				for (int step = 0; step < n + 1 - 1;  ++step)
						for (int i = 0; i < n + 1; ++i)
								for (int j = 0; j < n + 1; ++j) if (i != j && get_edge_cap(i, j) > 0)
										if (dist[j] > dist[i] + cost[i][j]){
												dist[j] = dist[i] + cost[i][j];
												prev[j] = i;
										}
				// detect negative cycle						
				for (int i = 0; i < n + 1; ++i)
						for (int j = 0; j < n + 1; ++j) if (i != j && get_edge_cap(i,j) > 0)
								if (dist[j] > dist[i] + cost[i][j]){ // negative cycle found
										vector<int> cycle = find_cycle(j, prev);
										int bottleneck = find_bottleneck_cycle(cycle);
										// update flow
										int ncycle = cycle.size();
										for (int i = 0; i < ncycle; ++i){
												int u = cycle[i], v = cycle[(i+1)%ncycle];
												flow[u][v] += bottleneck;
												flow[v][u] -= bottleneck;
										}
										return true;		
								}

				return false;
		}

		/*
		   A DFS to try to find a capacited path from where to sink.
		   If there exist a path like that, send flow through it.
		 */
		int send_flow(int where, int sink, vector<bool>& mark, int f)
		{
				if (where == sink) 
						return f;

				mark[where] = true;
				for (int next = 0; next < n; ++next) if (!mark[next]){
						int cap_edge = cap[where][next] - flow[where][next];
						if (cap_edge <= 0) continue;
						int sent_flow = send_flow(next, sink, mark, min(f, cap_edge));
						if (sent_flow > 0) {
								flow[where][next] += sent_flow;
								flow[next][where] -= sent_flow;
								return sent_flow;
						}
				}
				return 0;	
		}



		public:
		/* 
		   Ford-Fulkerson max-flow algorithm 
		 */
		int max_flow(int source, int sink)
		{
				int f = 0;
				vector<bool> mark(n, false);

				while (1){
						int sent = send_flow(source, sink, mark, INF);
						if (!sent) break;	
						f += sent;
						fill(mark.begin(), mark.end(), false);
					
				}	
				return f;
		}

		bool feaseable_solution(){
				int source = n;
				int sink = n + 1;
				for (int i = 0; i < n; ++i)
						if (b[i] > 0) set_edge(source, i, 0, b[i]);
						else if (b[i] < 0) set_edge(i, sink, 0, -b[i]);
				n += 2;
				int f = max_flow(source, sink);
				// removing source and sink edges
				n -= 2;
				for (int i = 0; i < n; ++i)
						if (b[i] > 0) del_edge(source, i);
						else if (b[i] < 0) del_edge(i, sink);
				return f != 0;		
		}


		Graph_Cycle(int n): INF(0x3f3f3f3f){
				this->n = n;
				// create two extra spots for sink and source
				cost.resize(n + 2);
				flow.resize(n + 2);
				cap.resize(n + 2);
				b.resize(n + 2);
				for (int i = 0; i < n + 2; ++i){
						flow[i].resize(n + 2);
						cost[i].resize(n + 2);
						cap[i].resize(n + 2);
				}
		}

		// sets a edge v-->u with capacity e_cap and cost e_cost
		void set_edge(int v, int u, int e_cost, int e_cap){
				cost[v][u] = e_cost;
				cost[u][v] = -e_cost;
				flow[v][u] = flow[u][v] = 0;
				cap[v][u] = e_cap;
				cap[u][v] = 0;
		}

		bool is_edge_set(int u, int v){
				assert (0 <= u && u < n);
				assert (0 <= v && v < n);
				return cap[u][v] != 0;
		}
		
		// sets a edge v-->u with capacity e_cap and cost e_cost
		void del_edge(int v, int u){
				cost[v][u] = 0;
				cost[u][v] = 0;
				flow[v][u] = flow[u][v] = 0;
				cap[v][u] = 0;
				cap[u][v] = 0;
		}

		
		// sets inbalance of node v to a
		void set_inbalance(int v, int a){
				assert(0 <= v && v < n);
				b[v] = a;
		}

		ERROR_CODE min_cost_flow(){
				if (has_negative_edge())
						return MCF_ERROR_NEGATIVE_EDGE_COST;
				if (has_parallel_edge())
						return MCF_ERROR_PARALLEL_EDGE;
				if (!is_valid_inbalance())
						return MCF_ERROR_BAD_INBALANCE;
				if (!feaseable_solution())
						return MCF_ERROR_NO_FEASEABLE_SOLUTION;

				while (true)
					if (!improve_flow()) break;
				return MCF_SUCCESS;			
		}

		// basically a debug helper function
		bool is_valid_flow(){
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j) if (flow[i][j] > 0)
							if (flow[i][j] > cap[i][j]) return false;
				if (improve_flow()) return false;
				return true;
		}
	
		int calc_flow_cost(){
				int f = 0;
				for (int i = 0; i < n; ++i)
						for (int j = 0; j < n; ++j) if (flow[i][j] > 0)
								f += flow[i][j] * cost[i][j];
				return f;
		}
};

