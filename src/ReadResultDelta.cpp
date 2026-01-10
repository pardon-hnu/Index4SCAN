#include <bits/stdc++.h>
using namespace std;

struct Result {
    double ari  = -1;
    double time = -1;
};
vector<string> datasets = {
	//real graphs
	// "amazon",
    // "youtube",
    // "skitter",
    // "pokec",
    // "topcats",
    // "livejournal1",
    "orkut2",
    // "indochina",
    "uk2002",
    // "uk2005",
    // "it2004",
    // "friendster"
    // "sk2005"
};
void read_approximate_result(string file_path,map<string, Result> & results)
{
    ifstream fin(file_path);
    string line;
    string current_param = "";
    double current_time = -1;
    double current_ari  = -1;
    while (getline(fin, line)) {
        /* ---------- 1. epsilon & mu ---------- */
        if (line.find("epsilon:") != string::npos &&
            line.find("mu:") != string::npos) {

            double eps;
            int mu;
            sscanf(line.c_str(), "epsilon: %lf mu: %d", &eps, &mu);

            ostringstream oss;
            oss << "epsilon=" << eps << ",mu=" << mu;
            current_param = oss.str();

            current_time = -1;
            current_ari  = -1;
        }
        /* ---------- 2. TIME COST ---------- */
        else if (line.find("[TIME COST]") != string::npos) {
            sscanf(line.c_str(), "[TIME COST] %lf s.", &current_time);
        }
        /* ---------- 3. ARI ---------- */
        else if (line.find("ari index score:") != string::npos) {
            sscanf(line.c_str(),
                   "[EVALUATION] ari index score: %lf",
                   &current_ari);
        }
        /* ---------- 4. commit result ---------- */
        if (!current_param.empty() &&
            current_time >= 0 &&
            current_ari  >= 0) {
            results[current_param] = {current_ari, current_time};
            // reset to avoid duplicate insert
            current_param.clear();
            current_time = -1;
            current_ari  = -1;
        }
    }
    fin.close();
}
int main() {
    // vector<string> epsilons = {
    //     "0.8",
    //     "0.75",
    //     "0.7",
    //     "0.65",
    //     "0.6",
    //     "0.55",
    //     "0.5",
    //     "0.45",
    //     "0.4",
    //     "0.35",
    //     "0.3",
    //     "0.25",
    //     "0.2",
    // };
    vector<string> epsilons = {
        "0.8",
        "0.79","0.78","0.77","0.76",
        "0.75",
        "0.74","0.73","0.72","0.71",
        "0.7",
        "0.69","0.68","0.67","0.66",
        "0.65",
        "0.64","0.63","0.62","0.61",
        "0.6",
        "0.59","0.58","0.57","0.56",
        "0.55",
        "0.54","0.53","0.52","0.51",
        "0.5",
        "0.49","0.48","0.47","0.46",
        "0.45",
        "0.44","0.43","0.42","0.41",
        "0.4",
        "0.39","0.38","0.37","0.36",
        "0.35",
        "0.34","0.33","0.32","0.31",
        "0.3",
        "0.29","0.28","0.27","0.26",
        "0.25",
        "0.24","0.23","0.22","0.21",
        "0.2",
    };
    vector<string> mus={
        "2","3","4","5","6","7","8","9","10","11","12","13","14","15"
    };
    unordered_set<string> need_keys;
    for(auto eps:epsilons)
    {
        for(auto mu:mus)
        {
            string key="epsilon="+eps+",mu="+mu;
            need_keys.insert(key);
        }
    }
    vector<string> deltas={
        "5","10","20","50","100"
    };
    for(auto dataset:datasets)
    {
        /* ---------- 输出结果 ---------- */
        cout << "==== Parsed Results of "<<dataset<<" ====\n";
        for(auto delta:deltas)
        {
            cout<<"-------------------"<<"delta: "<<delta<<"-------------------\n";
            string forest_result_file="res/query/"+dataset+"/forest_"+delta+".txt";
            string ppt_result_file="res/query/"+dataset+"/ppt_"+delta+".txt";
            map<string, Result> ppt_results,forest_results;
            read_approximate_result(forest_result_file,forest_results);
            read_approximate_result(ppt_result_file,ppt_results);

            double total_forest_ari=0.0;
            double total_forest_time=0.0;
            double count=0.0;
            for (auto &kv : forest_results) {
                // cout << kv.first
                //      << "  ARI=" << kv.second.ari
                //      << "  Time=" << kv.second.time << " s\n";
                if(need_keys.find(kv.first)!=need_keys.end())
                {
                    total_forest_ari+=kv.second.ari;
                    total_forest_time+=kv.second.time;
                    count++;
                }
            }

            cout<<" average forest ari: "<<total_forest_ari/count<<endl;
            cout<<" average forest time: "<<total_forest_time/count<<endl;

            double total_ppt_ari=0.0;
            double total_ppt_time=0.0;
            count=0.0;
            for (auto &kv : ppt_results) {
                // cout << kv.first
                //      << "  ARI=" << kv.second.ari
                //      << "  Time=" << kv.second.time << " s\n";
                if(need_keys.find(kv.first)!=need_keys.end())
                {
                    total_ppt_ari+=kv.second.ari;
                    total_ppt_time+=kv.second.time;
                    count++;
                }
            }
            cout<<" average ppt ari: "<<total_ppt_ari/count<<endl;
            cout<<" average ppt time: "<<total_ppt_time/count<<endl;
        }
    }
    return 0;
}