#include <ctime>
#include "../cycle_canceling_mean.cpp"
#include "../cycle_canceling.cpp"

class UnitTests
{
public:
	double total_time_mean;
	double total_time_cycle;

	UnitTests(){
		total_time_mean = 0;
		total_time_cycle = 0;
	}
/*
	void test_maxflow()
	{
		Graph grafo(3);
		grafo.set_edge(0, 1, 2, 10);
		grafo.set_edge(1, 2, 2, 10);
		assert(grafo.max_flow(0, 2) == 10);
		cerr<< "test_maxflow: SUCCESS!"<<'\n';
	}

	void test_maxflow2(){
		Graph grafo(4);
		grafo.set_edge(0, 1, 2, 10);
		grafo.set_edge(0, 2, 2, 10);
		grafo.set_edge(1, 3, 2, 10);
		grafo.set_edge(2, 3, 2, 10);
		grafo.set_edge(1, 2, 2, 1);
		assert(grafo.max_flow(0, 3) == 20);
		cerr<< "test_maxflow2: SUCCESS!"<<'\n';
	}

	void test_mincostflow()
	{
		Graph grafo(4);
		grafo.set_edge(0, 1, 2, 10);
		grafo.set_edge(0, 2, 2, 10);
		grafo.set_edge(1, 3, 2, 10);
		grafo.set_edge(2, 3, 2, 10);
		grafo.set_edge(1, 2, 2, 1);
		grafo.set_inbalance(0, 20);
		grafo.set_inbalance(3, -20);
		assert(grafo.min_cost_flow() == MCF_SUCCESS);
		assert(grafo.is_valid_flow());
		assert(grafo.calc_flow_cost() == 80);
		cerr<< "test_mincostflow: SUCCESS!"<<'\n';
	}

	void test_mincostflow2()
	{
		Graph grafo(4);

		//se_edgee(v, u, cost, cap)
		grafo.set_edge(0, 1, 10, 2);
		grafo.set_edge(0, 2, 10, 1);
		grafo.set_edge(1, 3, 100, 3);
		grafo.set_edge(2, 3, 10, 2);
		grafo.set_edge(1, 2, 10, 1);
		grafo.set_inbalance(0, 2);
		grafo.set_inbalance(3, -2);
		assert(grafo.min_cost_flow() == MCF_SUCCESS);
		assert(grafo.is_valid_flow());
		assert(grafo.calc_flow_cost() == 50);
		cerr<< "test_mincostflow2: SUCCESS!"<<'\n';
	}

	void test_mincostflow3()
	{
		Graph grafo(4);

		//se_edgee(v, u, cost, cap)
		grafo.set_edge(0, 1, 10, 2);
		grafo.set_edge(0, 2, 10, 1);
		grafo.set_edge(1, 3, 100, 3);
		grafo.set_edge(2, 3, 10, 2);
		grafo.set_edge(1, 2, 10, 1);
		grafo.set_inbalance(0, 2);
		grafo.set_inbalance(1, -1);
		grafo.set_inbalance(3, -1);
		assert(grafo.min_cost_flow() == MCF_SUCCESS);
		assert(grafo.is_valid_flow());
		assert(grafo.calc_flow_cost() == 30);
		cerr<< "test_mincostflow3: SUCCESS!"<<'\n';
	}

	void test_mincostflow4()
	{
		Graph grafo(4);

		//se_edgee(v, u, cost, cap)
		grafo.set_edge(0, 1, 10, 2);
		grafo.set_edge(0, 2, 10, 1);
		grafo.set_edge(1, 3, 100, 3);
		grafo.set_edge(2, 3, 10, 2);
		grafo.set_edge(1, 2, 10, 1);
		grafo.set_inbalance(0, 3);
		grafo.set_inbalance(1, -1);
		grafo.set_inbalance(2, -1);
		grafo.set_inbalance(3, -1);
		assert(grafo.min_cost_flow() == MCF_SUCCESS);
		assert(grafo.is_valid_flow());
		assert(grafo.calc_flow_cost() == 50);
		cerr<< "test_mincostflow4: SUCCESS!"<<'\n';
	}

	void test_mincostflow5()
	{
		Graph grafo(3);

		//se_edgee(v, u, cost, cap)
		grafo.set_edge(0, 1, 4, 4);
		grafo.set_edge(0, 2, 4, 4);
		grafo.set_edge(1, 2, 4, 4);
		grafo.set_inbalance(0, 6);
		grafo.set_inbalance(1, -1);
		grafo.set_inbalance(2, -5);
		assert(grafo.min_cost_flow() == MCF_SUCCESS);
		assert(grafo.is_valid_flow());
		assert(grafo.calc_flow_cost() == 28);
		cerr<< "test_mincostflow5: SUCCESS!"<<'\n';
	}
*/
	int floyd_warshall(int src, int sink, vector<vector<int> >& graph, vector<vector<int> >& cost)
	{
		#define INF 0x3f3f3f3f
		int n = graph.size();
		vector<vector<int> > dist(n);
		for (int i = 0; i < n; ++i) dist[i].resize(n);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j) 
				if (graph[i][j]) 
					dist[i][j] = cost[i][j];
				else
					dist[i][j] = INF;
					
		for (int k = 0; k < n; ++k)
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j) 
					dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
		return dist[src][sink];
	}

		

	/*
		Formulates the minimum path problem
		as a minimum cost flow. Validate the
		found solution comparing it to Floyd-Warshall.
		
		The instances are random graphs with at most
		50 vertex and at most 25 * 50 non-parallel edges 	
	*/
	void test_minimum_path(int nv, int ne)
	{
		srand(time(NULL));

		// random number of vertexs
		//int n = rand()%100 + 2;
		int n = nv;
		Graph_Mean grafo_mean(n);
		Graph_Cycle grafo_cycle(n);

		// random number of edges
		//int m = rand()%(n/2 * n);
		int m = ne;

		// force the graph to have a path
		// between 0 and n - 1
		//grafo_mean.set_edge(0, n - 1, 1000 * 1000, 1);
		//grafo_cycle.set_edge(0, n - 1, 1000 * 1000, 1);


		// generate the remaining edges randomly
		for (int i = 0; i < m - 1; ++i)
		{
			int u, v, c;
			// generate a new edge, check if it or its inverse were
			// already generated. Parallel edges are not allowed.
			do
			{
				u = rand()%n;
				v = rand()%n;
				ca = (rand() + 1) %100000; // no zero capacity
				co = (rand() + 1) %100000; // no zero cost
			}while (grafo_mean.is_edge_set(u, v) || grafo_mean.is_edge_set(v, u));

			grafo_mean.set_edge(u, v, co, ca);
			grafo_cycle.set_edge(u, v, co, ca);

			mp_graph[u][v] = 1;
			mp_cost[u][v] = c;
		}
	
		grafo_mean.set_inbalance(0, 1);
		grafo_cycle.set_inbalance(0, 1);

		grafo_mean.set_inbalance(n-1, -1);
		grafo_cycle.set_inbalance(n-1, -1);

		int expected_value = floyd_warshall(0, n - 1, mp_graph, mp_cost);

		// Mean
		time_t tstart, tend;
		time(&tstart);
		assert(grafo_mean.min_cost_flow() == MCF_SUCCESS);
		
		//assert(grafo.is_valid_flow());
		int value = grafo_mean.calc_flow_cost();
		time(&tend);
		total_time_mean += difftime(tend, tstart);
		assert(value == expected_value);

		// Cycle
		time(&tstart);
		assert(grafo_cycle.min_cost_flow() == MCF_SUCCESS);
		
		//assert(grafo.is_valid_flow());
		 value = grafo_cycle.calc_flow_cost();
		time(&tend);
		total_time_cycle += difftime(tend, tstart);
		assert(value == expected_value);

		cerr<< "test_minimum_path: SUCCESS!"<<'\n';
	}

};

int main(int argc, char* argv[])
{
	UnitTests tests;
/*	tests.test_maxflow();
	tests.test_maxflow2();
	tests.test_mincostflow();
	tests.test_mincostflow2();
	tests.test_mincostflow3();
	tests.test_mincostflow4();
	tests.test_mincostflow5();*/
	int nv, ne;
	assert(argc == 3);
	nv = atoi(argv[1]);
	ne = atoi(argv[2]);

	int ninterations = 100;
	for (int i = 0; i < ninterations; ++i){
		try	
		{
			tests.test_minimum_path(nv, ne);
		}
		catch (exception& e)
		{
			cerr<<"Exception raised!!! "<<e.what()<<endl;
		}
	 }
	 cout<<"Average time mean = "<<tests.total_time_mean/ninterations<<"s"<<endl;
	 cout<<"Average time cycle = "<<tests.total_time_cycle/ninterations<<"s"<<endl;
	return 0;
}
