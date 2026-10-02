#include <iostream>
#include <global/global.h>
#include <global/SCAN.h>
vector<int> BOTBIN::rnk = vector<int>();
vector<int> SCAN_Forest_Index::rnk=vector<int>();
vector<int> SCAN_PPT_Index::rnk=vector<int>();
using namespace std;
int main(int argc, char * argv[])
{
	string dataset = argv[1];
	string method = argv[2];
	string pattern = argv[3];
	vector<string> parameters;
	if(pattern=="query")
	{
		float epsilon = atof(argv[4]);
		int mu = atoi(argv[5]);
		for(int i=6;i<argc;i++)
		{
			string parameter=argv[i];
			parameters.push_back(parameter);
		}
		SCAN worker(dataset,method,parameters);
		worker.single_cluster(epsilon,mu);
	}
	else if(pattern=="construct")
	{
		for(int i=4;i<argc;i++)
		{
			string parameter=argv[i];
			parameters.push_back(parameter);
		}
		SCAN worker(dataset,method,parameters);
		worker.index_construct();
	}
	else if(pattern=="exp-query")
	{
		// vector<float> epsilons={0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1};
		// vector<float> epsilons={
		// 	// 0.90,
		// 	// 0.89,0.88,0.87,0.86,
		// 	// 0.85,
		// 	// 0.84,0.83,0.82,0.81,
		// 	0.80,
		// 	0.79,0.78,0.77,0.76,
		// 	0.75,
		// 	0.74,0.73,0.72,0.71,
		// 	0.70,
		// 	0.69,0.68,0.67,0.66,
		// 	0.65,
		// 	0.64,0.63,0.62,0.61,
		// 	0.60,
		// 	0.59,0.58,0.57,0.56,
		// 	0.55,
		// 	0.54,0.53,0.52,0.51,
		// 	0.50,
		// 	0.49,0.48,0.47,0.46,
		// 	0.45,
		// 	0.44,0.43,0.42,0.41,
		// 	0.40,
		// 	0.39,0.38,0.37,0.36,
		// 	0.35,
		// 	0.34,0.33,0.32,0.31,
		// 	0.30,
		// 	0.29,0.28,0.27,0.26,
		// 	0.25,
		// 	0.24,0.23,0.22,0.21,
		// 	0.20,
		// 	// 0.19,0.18,0.17,0.16,
		// 	// 0.15,
		// 	// 0.14,0.13,0.12,0.11,
		// 	// 0.10
		// };
		vector<float> epsilons = {
			0.8,
			0.7,
			0.6,
			0.5,
			0.4,
			0.3,
			0.2
		};
		vector<int> mus={2,3,4,5,6,7,8,9,10,11,12,13,14,15};
		for(int i=4;i<argc;i++)
		{
			string parameter=argv[i];
			parameters.push_back(parameter);
		}
		SCAN worker(dataset,method,parameters);
		bool output_labels=false;
		// if(method=="GS-Index")
		// {
		// 	output_labels=true;
		// }
		worker.batch_cluster(epsilons,mus,false,output_labels);
	}
	else if(pattern=="test")
	{
		for(int i=4;i<argc;i++)
		{
			string parameter=argv[i];
			parameters.push_back(parameter);
		}
		SCAN worker(dataset,method,parameters);
		worker.index_test();
	}
	else if(pattern=="update")
	{
		string update_type=argv[4];
		int u=atoi(argv[5]);
		int v=atoi(argv[6]);
		int update_way=atoi(argv[7]);
		for(int i=8;i<argc;i++)
		{
			string parameter=argv[i];
			parameters.push_back(parameter);
		}
		SCAN worker(dataset,method,parameters);
		worker.update(update_type,u,v,update_way);
	}
	else if(pattern=="exp-update")
	{
		string update_type=argv[4];
		int update_way=atoi(argv[5]);
		for(int i=6;i<argc;i++)
		{
			string parameter=argv[i];
			parameters.push_back(parameter);
		}
		SCAN worker(dataset,method,parameters);
		worker.batch_update(update_type, update_way);
	}
	else
	{
		cout<<"error!!! no such pattern of parameters"<<endl;
	}
	return 0;
}