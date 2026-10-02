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
    // "friendster",
    // "sk2005"
};
void read_exact_result(string file_path,map<string, Result> & results)
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
        /* ---------- 4. commit result ---------- */
        if (!current_param.empty() &&
            current_time >= 0) {
            results[current_param] = {current_ari, current_time};
            // reset to avoid duplicate insert
            current_param.clear();
            current_time = -1;
        }
    }
    fin.close();
}
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
    //     vector<string> epsilons = {
    //     "0.8","0.79","0.78","0.77","0.76","0.75","0.74","0.73","0.72","0.71",
    //     "0.7","0.69","0.68","0.67","0.66","0.65","0.64","0.63","0.62","0.61",
    //     "0.6","0.59","0.58","0.57","0.56","0.55","0.54","0.53","0.52","0.51",
    //     "0.5","0.49","0.48","0.47","0.46","0.45","0.44","0.43","0.42","0.41",
    //     "0.4","0.39","0.38","0.37","0.36","0.35","0.34","0.33","0.32","0.31",
    //     "0.3","0.29","0.28","0.27","0.26","0.25","0.24","0.23","0.22","0.21",
    //     "0.2",
    // };
    // vector<string> epsilons = {
    //     "0.8",
    //     "0.7",
    //     "0.6",
    //     "0.5",
    //     "0.4",
    //     "0.3",
    //     "0.2"
    // };
    vector<string> epsilons = {
        "0.6"
    };
    // vector<string> mus={
    //     "5"
    // };
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
    for(auto dataset:datasets)
    {
        string gs_result_file="res/query/"+dataset+"/gs-index.txt";
        string botbin_result_file="res/query/"+dataset+"/botbin.txt";
        string forest_result_file="res/query/"+dataset+"/forest.txt";
        string ppt_result_file="res/query/"+dataset+"/ppt.txt";
        string vdstar_result_file="res/query/"+dataset+"/vd.txt";
        string vdstar_not_result_file="res/query/"+dataset+"/vdnot.txt";
        // if(dataset=="uk2002")
        // {
        //     gs_result_file="res/varying_parameter_uk2002/debug_uk2002_gs.res";
        //     botbin_result_file="res/varying_parameter_uk2002/debug_uk2002_botbin.res";
        //     forest_result_file="res/varying_parameter_uk2002/uk2002_forest_mu.res";
        //     ppt_result_file="res/varying_parameter_uk2002/debug_uk2002_ppt.res";
        // }
        map<string, Result> gs_results,botbin_results,ppt_results,forest_results,vdstar_results,vdstar_not_results;
        read_exact_result(gs_result_file,gs_results);
        read_approximate_result(botbin_result_file,botbin_results);
        read_approximate_result(forest_result_file,forest_results);
        read_approximate_result(ppt_result_file,ppt_results);
        read_exact_result(vdstar_result_file,vdstar_results);
        read_exact_result(vdstar_not_result_file,vdstar_not_results);
        /* ---------- 输出结果 ---------- */
        cout << "==== Parsed Results of "<<dataset<<" ====\n";
        double total_vdstar_time=0.0;
        double count=0.0;
        for(auto &kv:vdstar_results)
        {
            if(need_keys.find(kv.first)!=need_keys.end())
            {
                // cout << kv.first<< "  Time=" << kv.second.time << " s\n";
                total_vdstar_time+=kv.second.time;
                count++;
            }
        }
        cout<<" average vdstar time: "<<total_vdstar_time/count<<endl;

        double total_vdstar_not_time=0.0;
        count=0.0;
        for(auto &kv:vdstar_not_results)
        {
            if(need_keys.find(kv.first)!=need_keys.end())
            {
                cout << kv.first<< "  Time=" << kv.second.time << " s\n";
                total_vdstar_not_time+=kv.second.time;
                count++;
            }
        }
        cout<<" average vdstar not time: "<<total_vdstar_not_time/count<<endl;

        double total_gs_time=0.0;
        count=0.0;
        for (auto &kv : gs_results) {
            if(need_keys.find(kv.first)!=need_keys.end())
            {
                // cout << kv.first<< "  Time=" << kv.second.time << " s\n";
                total_gs_time+=kv.second.time;
                count++;
            }
        }
        cout<<" average gs time: "<<total_gs_time/count<<endl;

        double total_botbin_ari=0.0;
        double total_botbin_time=0.0;
        count=0.0;
        for (auto &kv : botbin_results) {
            if(need_keys.find(kv.first)!=need_keys.end())
            {
                // cout << kv.first
                //  << "  ARI=" << kv.second.ari
                //  << "  Time=" << kv.second.time << " s\n";
                total_botbin_ari+=kv.second.ari;
                total_botbin_time+=kv.second.time;
                count++;
            }
        }
        cout<<" average botbin ari: "<<total_botbin_ari/count<<endl;
        cout<<" average botbin time: "<<total_botbin_time/count<<endl;

        double total_forest_ari=0.0;
        double total_forest_time=0.0;
        double forest_max_speed_up_over_gs=0.0;
        double forest_max_speed_up_over_botbin=0.0;
        double forest_max_speed_up_over_vdstar=0.0;
        count=0.0;
        for (auto &kv : forest_results) {
            if(need_keys.find(kv.first)!=need_keys.end())
            {
                // cout << kv.first
                //  << "  ARI=" << kv.second.ari
                //  << "  Time=" << kv.second.time << " s\n";
                double temp_forest_time=kv.second.time;
                double temp_gs_time=gs_results[kv.first].time;
                double temp_botbin_time=botbin_results[kv.first].time;
                double temp_vdstar_time=vdstar_results[kv.first].time;
                double temp_forest_speed_up_over_gs=temp_gs_time/temp_forest_time;
                double temp_forest_speed_up_over_botbin=temp_botbin_time/temp_forest_time;
                double temp_forest_time_over_vdstar=temp_vdstar_time/temp_forest_time;
                if(forest_max_speed_up_over_gs<temp_forest_speed_up_over_gs)
                {
                    forest_max_speed_up_over_gs=temp_forest_speed_up_over_gs;
                }
                if(forest_max_speed_up_over_botbin<temp_forest_speed_up_over_botbin)
                {
                    forest_max_speed_up_over_botbin=temp_forest_speed_up_over_botbin;
                }
                if(forest_max_speed_up_over_vdstar<temp_forest_time_over_vdstar)
                {
                    forest_max_speed_up_over_vdstar=temp_forest_time_over_vdstar;
                }
                total_forest_ari+=kv.second.ari;
                total_forest_time+=kv.second.time;
                count++;
            }
        }

        cout<<" average forest ari: "<<total_forest_ari/count<<endl;
        cout<<" average forest time: "<<total_forest_time/count<<endl;
        cout<<" average speedup over gs: "<<total_gs_time/total_forest_time<<endl;
        cout<<" max speedup over gs: "<<forest_max_speed_up_over_gs<<endl;
        cout<<" average speedup over botbin: "<<total_botbin_time/total_forest_time<<endl;
        cout<<" max speedup over botbin: "<<forest_max_speed_up_over_botbin<<endl;
        cout<<" average speedup over vdstar: "<<total_vdstar_time/total_forest_time<<endl;
        cout<<" max speedup over vdstar: "<<forest_max_speed_up_over_vdstar<<endl;
        double total_ppt_ari=0.0;
        double total_ppt_time=0.0;
        double ppt_max_speed_up_over_gs=0.0;
        double ppt_max_speed_up_over_botbin=0.0;
        double ppt_max_speed_up_over_vdstar=0.0;
        count=0.0;
        for (auto &kv : ppt_results) {
            if(need_keys.find(kv.first)!=need_keys.end())
            {
                // cout << kv.first
                //  << "  ARI=" << kv.second.ari
                //  << "  Time=" << kv.second.time << " s\n";
                double temp_ppt_time=kv.second.time;
                double temp_gs_time=gs_results[kv.first].time;
                double temp_botbin_time=botbin_results[kv.first].time;
                double temp_vdstar_time=vdstar_results[kv.first].time;
                double temp_ppt_speed_up_over_gs=temp_gs_time/temp_ppt_time;
                double temp_ppt_speed_up_over_botbin=temp_botbin_time/temp_ppt_time;
                double temp_ppt_speed_up_over_vdstar=temp_vdstar_time/temp_ppt_time;
                if(ppt_max_speed_up_over_gs<temp_ppt_speed_up_over_gs)
                {
                    ppt_max_speed_up_over_gs=temp_ppt_speed_up_over_gs;
                }
                if(ppt_max_speed_up_over_botbin<temp_ppt_speed_up_over_botbin)
                {
                    ppt_max_speed_up_over_botbin=temp_ppt_speed_up_over_botbin;
                }
                if(ppt_max_speed_up_over_vdstar<temp_ppt_speed_up_over_vdstar)
                {
                    ppt_max_speed_up_over_vdstar=temp_ppt_speed_up_over_vdstar;
                }
                total_ppt_ari+=kv.second.ari;
                total_ppt_time+=kv.second.time;
                count++;
            }
        }
        cout<<" average ppt ari: "<<total_ppt_ari/count<<endl;
        cout<<" average ppt time: "<<total_ppt_time/count<<endl;
        cout<<" average speedup over gs: "<<total_gs_time/total_ppt_time<<endl;  
        cout<<" max speedup over gs: "<<ppt_max_speed_up_over_gs<<endl;
        cout<<" average speedup over botbin: "<<total_botbin_time/total_ppt_time<<endl;
        cout<<" max speedup over botbin: "<<ppt_max_speed_up_over_botbin<<endl;
        cout<<" average speedup over vdstar: "<<total_vdstar_time/total_ppt_time<<endl;
        cout<<" max speedup over vdstar: "<<ppt_max_speed_up_over_vdstar<<endl;
        ofstream fout;   
        string output_file="./res/time/"+dataset+".txt";
        fout.open(output_file);
        for (auto &kv : botbin_results) {
            if(need_keys.find(kv.first)!=need_keys.end())
            {
                auto key=kv.first;
                fout<<ppt_results[key].time<<" "<<forest_results[key].time<<" "<<gs_results[key].time<<" "<<botbin_results[key].time<<" "<<vdstar_results[key].time<<'\n';
        
            }
        }
        fout.close();

    }
    return 0;
}