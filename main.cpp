#include <bits/stdc++.h>
//#include "MT.h"

using namespace std;
float inf = numeric_limits<float>::infinity();
using ll = long long;
const double radius = 100;
const double GAMMA = pow(10, 6.0 / 10.0);
const double end_time = 10000;
constexpr double PI = 3.14159265358979323846264338;
const double theta = pow(10, 0.4);
//const double delta = 2.0 / pass_loss_exponent;
const double lambda_BS = 1.;
//const double POWER = 1.0;

//出力ファイルを開く
ofstream outputfile1;
ofstream outputfile2;
ofstream outputfile;

struct pos {
    double x;
    double y;
};

struct terminal {
    int id;
    pair<double, double> pos;
    //double to_BS[num_BS];
    double dst_to_BS;
    double dst_to_FBS;
    int nearest_BS_id;
    double level;
    double power;
    double IF_to_FBS;
    bool state; //True: in , False: out
    bool operator<( const terminal& right ) const {
            return power == right.power ? true : power < right.power;
    }
};


double urand(){ //0-1
    double m, a;
    m = RAND_MAX + 1.0;
    a = (rand() + 0.5)/m;
    a = (rand() + a)/m;
    return (rand() + a)/m;
}

//乱数発生
double my_rand(double Min, double Max) {
    mt19937 mt{ random_device{}() };
    //random_device rd;
    //default_random_engine eng(rd());
    uniform_real_distribution<double> distr(Min, Max);
    return distr(mt);
}

//円形領域内の座標をランダムに取得
pair<double, double> coordinate() {
//    double x, y;
//    while (true) {
//        x = radius * urand();
//        y = radius * urand();
//        double r_sq = x * x + y * y;
//        if (r_sq < radius * radius) break;
//    }
    
    
    double r = radius * sqrt(urand()); //[0, 1)
    double theta = PI * (2. * urand() - 1.); //(-1, 1)
//    while (theta < -PI || PI < theta) {
//        theta = PI * (2. * urand() - 1.);
//    }
    double x = r * cos(theta);
    double y = r * sin(theta);
    pair<double, double> pos = make_pair(x, y);
    return pos;
}

//2点間の距離を計算
double cal_dst(pair<double, double> pos1, pair<double, double> pos2) {
    return sqrt((pos1.first - pos2.first)*(pos1.first - pos2.first) + (pos1.second - pos2.second)*(pos1.second - pos2.second));
}

double distance(pos pos1, pos pos2) {
    return sqrt((pos1.x - pos2.x) * (pos1.x - pos2.x) + (pos1.y - pos2.y) * (pos1.y - pos2.y));
}

//正規分布乱数発生
double normrand(){
    return sqrt(-2*log(urand()))*cos(2*M_PI*urand());
}


double gauss_rand(double mu, double sig) {
    //mt19937 mt{ std::random_device{}() };
    random_device rd;
    default_random_engine eng(rd());
    normal_distribution<> dist(mu, sig);
    double g = dist(eng);
    return g;
}


double exp_dist(double lambda) {
    //double g = genrand_real3();
    double g = urand();
    double tau = - log(1 - g) / lambda;
    return tau;
}

double poisson_dist(double lambda) {
    double sum = 0, k;
    for (k = 0; sum < 1; k++) {
        sum += exp_dist(lambda);
    }
    return  k;
}

//ポアソン乱数
double poisson(double lambda) {
    random_device rd;
    default_random_engine eng(rd());
    poisson_distribution<> dist(lambda);
    return dist(eng);
}



//送信成功確率 分布作成
void success_dist(int lambda, int L, int C, string s) {
    //データ読み込み
    string filename0 = to_string(lambda) + "_NOMA_" + to_string(L) + "_" + to_string(C) + s + ".txt";
    //string filename0 = "5_ALOHA_dst.txt";
    vector<pair<double, double>> suc_dist;
    ifstream readingfile;
    readingfile.open(filename0);
    ifstream ifs(filename0);
        if (ifs.fail()) {
           cerr << "Cannot open file\n";
           exit(0);
        }
    double x, y;
    while (ifs >> x >> y) {
        suc_dist.push_back(make_pair(y, x));
        //cout << x << " " << y << endl;
    }
    ifs.close();
    
    //first 距離 second 成功失敗
    unsigned long num_dev = suc_dist.size();
    sort(suc_dist.begin(), suc_dist.end());
    vector<double> total(50, 0);
    vector<double> suc_sum(50, 0);
    double j = 1.;
    for (int i = 0; i < num_dev; i++) {
        if (suc_dist.at(i).first < 0.05 * j) {
            total.at(j-1)++;
            if (suc_dist.at(i).second == 1) suc_sum.at(j-1)++;
        } else {
            j++;
            total.at(j-1)++;
            if (suc_dist.at(i).second == 1) suc_sum.at(j-1)++;
        }
    }
    
    for (int i = 1; i <= 20; i++) {
        cout << 0.05 * i << " " << suc_sum.at(i-1) / total.at(i-1) << " " << suc_sum.at(i-1) << " " << total.at(i-1) << endl;
    }
    
}




void ALOHA_pr(double lambda_IoT, double alpha, double noise) {
    //double success = 0;
    cout << "path loss exp : " << alpha << endl;
    //outputfile << "pass loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int i = 0; i < end_time; i++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
        //BSの座標を設定
        pair<double, double> origin = make_pair(0, 0);
        pair<double, double> BS_pos;
        pair<double, double> nearest_pos;
        double nd = inf;
        for (int i = 0; i < num_BS; i++) {
            BS_pos = coordinate();
            double dst = cal_dst(origin, BS_pos);
            if (dst < nd) {
                nd = dst;
                nearest_pos = BS_pos;
            }
        }
        
        double SI = 0, coef = 0.0;
        double H = exp_dist(1.0);
        coef = H / pow(nd, alpha);
        //端末情報を初期化
        for (int i = 1; i < num_IoT; i++) {
            pair<double, double> IoT_pos;
            IoT_pos = coordinate();
            //端末-基地局間のフェージング係数を設定
            double H = exp_dist(1.0);//gauss_rand(0, 1);
            double dst2 = cal_dst(IoT_pos, nearest_pos);
            SI += H / pow(dst2, alpha);
        }
        
        double SINR;
        int flag = 0;
        SINR = coef / (SI + noise);
        if (SINR > theta) flag = 1;
        
        outputfile << flag << " " << nd << endl;
        
        
        //進捗状況を表示
        double progress = i / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
}



void SIC_pr(double lambda_IoT, double alpha, double noise) {
    //double success = 0;
    cout << "path loss exp : " << alpha << endl;
    //outputfile << "pass loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int i = 0; i < end_time; i++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
        //BSの座標を設定
        pair<double, double> origin = make_pair(0, 0);
        pair<double, double> BS_pos;
        pair<double, double> nearest_pos;
        double nd = inf;
        for (int i = 0; i < num_BS; i++) {
            BS_pos = coordinate();
            double dst = cal_dst(origin, BS_pos);
            if (dst < nd) {
                nd = dst;
                nearest_pos = BS_pos;
            }
        }
        
        double SI = 0, coef = 0.0, greatest = 0.0;
        double H = exp_dist(1.0);
        coef = H / pow(nd, alpha);
        greatest = coef;
        //端末情報を初期化
        for (int i = 1; i < num_IoT; i++) {
            pair<double, double> IoT_pos;
            IoT_pos = coordinate();
            
            //端末-基地局間のフェージング係数を設定
            double H = exp_dist(1.0);//gauss_rand(0, 1);
            double dst2 = cal_dst(IoT_pos, nearest_pos);
            double tmp = H / pow(dst2, alpha);
            if (tmp > greatest) greatest = tmp;
            SI += tmp;
        }
        
        double SINR = coef / (SI + noise);
        int flag = 0;
        if (SINR > theta) {
            flag = 1;
        } else {
            SINR = greatest / (SI - greatest + coef + noise);
            if (SINR > theta) {
                SINR = coef / (SI - greatest + noise);
                if (SINR > theta) flag = 1;
            }
        }
        
        outputfile << flag << " " << nd << endl;
        
        
        //進捗状況を表示
        double progress = i / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
}


//Power allocation Probability by distance
//原点にデバイスを一台置いて，そのデバイスが送信に成功した時の位置を記録
//あるスロットの全デバイスを調査した方がいいかも
void PA_NOMA_pr_dst(double lambda_IoT, double alpha, double noise, double L, int C, string s, double D) {
    double success = 0;
    cout << "path loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    vector<double> coef(L);
    if (D == 0) { //電力配分最適化
        if (L == 2) coef = {0, 0.47164};
        else if (L == 3) coef = {0, 0.36497, 0.595226};
        else if (L == 4) coef = {0, 0.175, 0.38608, 0.67117};
        else if (L == 5) coef = {0, 0.267434, 0.4051755, 0.545601, 0.720014};
        else if (L == 6) coef = {0, 0.1073693, 0.2292, 0.364867, 0.532356, 0.76407};
        else if (L == 7) coef = {0, 0.090711, 0.189755, 0.299831, 0.428612, 0.581466, 0.798057};
        else if (L == 8) coef = {0, 0.07778, 0.1625, 0.2547469, 0.3570889, 0.475923, 0.61921, 0.828777};
        else if (L == 9) coef = {0, 0.068333, 0.1428356, 0.22055, 0.306667, 0.736885, 0.516842, 0.65526, 0.8447569};
        else if (L == 10) coef = {0, 0.06166, 0.125667, 0.195537, 0.27069, 0.35214, 0.443844, 0.551359, 0.68066, 0.8680};
    } else {
        for (int i = 0; i < L; i++) {
            coef.at(i) = D * sqrt(i / L);
        }
    }

    cout << "Progress is 0%";
    for (int t = 0; t < end_time; t++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
        //BSの座標を設定
        pair<double, double> origin = make_pair(0, 0);
        vector<pair<double, double>> BS_pos(num_BS);
        vector<pair<double, double>> near_BS_pos; //原点に近いセルだけ電力制御することにする．近い座標だけ保存
        pair<double, double> FBS_pos; //Focused Base station（原点のデバイスから最も近い基地局）
        double nearest = inf;
        for (int i = 0; i < num_BS; i++) {
            pair<double, double> bspos = coordinate();
            BS_pos.at(i) = bspos;
            double dst = cal_dst(origin, bspos);
            if (dst < 20.0) {
                near_BS_pos.push_back(bspos);
                if (dst < nearest) {
                    nearest = dst;
                    FBS_pos = bspos;
                }
            }
        }
        int num_near_BS = (int)near_BS_pos.size();  //参照避けによる高速化
        
        //ノイズを入れるならノイズは送信電力の比にする
        double SI = 0, power = 0.0;
        vector<double> ACL;
        
        //原点デバイス
        double dst_0 = cal_dst(origin, FBS_pos), level_d0;
        for (int l = L - 1; ;l--) { //電力レベル決定
            if (coef.at(l) < dst_0) { //Radius.at(l) < dst_to_NBS
                level_d0 = l;
                break;
            }
        }
        //if (L==1) level_d0 = 0;
        double P = theta * pow(theta + 1, L - level_d0 - 1) * exp_dist(1.0);
        ACL.emplace_back(P); SI += P;
        int cha_d0 = rand() % C;
        
        //端末情報を初期化
        /* デバイスがすること（原点セルの座標だけ可視化する場合）
         1.原点セルの基地局との距離を計算し，離れてたらただ干渉信号として足す
         〜近い場合
         1.（比較的近い基地局リストの）全基地局との距離を計算して，一番近い基地局を探す．
         2.基地局までの距離から，割り当てられる送信電力を計算する
         3.
         〜原点セルにデバイスがある場合　割り当てられた電力で信号を足す（フェージングは推定できないとする）
         〜原点セル以外の場合　原点セルの基地局への干渉信号として電力を足す
         */
        for (int i = 1; i < num_IoT; i++) {
            pair<double, double> device_pos = coordinate();
            double dst_to_FBS = cal_dst(device_pos, FBS_pos);
            int channel = rand() % C;
            if (channel != cha_d0) continue; //別のチャネルは無視して良い
            
            double H = exp_dist(1.0);
            if (dst_to_FBS > 20.0) {
                SI += theta * pow(theta + 1, L - rand() % (int)L - 1) * H * pow(urand() / dst_to_FBS, alpha);
            } else {
                double dst_to_NBS = dst_to_FBS;
                int cell_num = 0, level = 0;
                for (int k = 0; k < num_near_BS; k++) {
                    double tmp_dst = cal_dst(device_pos, near_BS_pos.at(k));
                    if (tmp_dst < dst_to_NBS) {
                        dst_to_NBS = tmp_dst;
                        cell_num = -1;
                    }
                }
                for (int l = L - 1; ;l--) { //電力レベル決定
                    if (sqrt((double)l / L) < dst_to_NBS) {
                        level = l;
                        break;
                    }
                }
                
                power = theta * pow(theta + 1, L - level - 1) * H * pow(dst_to_NBS / dst_to_FBS, alpha);
                SI += power;
                if (cell_num == 0) { //原点セルにデバイスがある場合
                    ACL.emplace_back(power);
                }
            }
            
        }
        
        int flag = 0;
        sort(ACL.rbegin(), ACL.rend());
        for (int i = 0; i < ACL.size(); i++) {
            SI -= ACL.at(i);
            double SINR = P / (SI + noise);
            //cout << P << " " << ACL.at(i) << " " << SINR << " " << dst_0 << " " << i << endl;
            if (SINR > theta) {
                if (ACL.at(i) == P) {
                    flag = 1; success++;
                    break;
                }
            } else break;
        }
        
        
        //cout << s << endl;

        outputfile << flag << " " << dst_0 << endl;
        
        //進捗状況を表示
        double progress = t / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    cout << success / end_time << endl;
    
    cout << lambda_IoT << " " << success / end_time << endl;
    cout << endl;
    cout << "Distribution" << endl;
    success_dist(lambda_IoT, L, C, s);
    
}


//Throughput of Power allocation
void PA_NOMA_thp(double lambda_IoT, double alpha, double noise, double L, int C, double D) {
    double success = 0;
    cout << "path loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    vector<double> coef(L);
    if (D == 0) {
        if (L == 2) coef = {0, 0.47164};
        else if (L == 3) coef = {0, 0.36497, 0.595226};
        else if (L == 4) coef = {0, 0.175, 0.38608, 0.67117};
        else if (L == 5) coef = {0, 0.267434, 0.4051755, 0.545601, 0.720014};
        else if (L == 6) coef = {0, 0.1073693, 0.2292, 0.364867, 0.532356, 0.76407};
        else if (L == 7) coef = {0, 0.090711, 0.189755, 0.299831, 0.428612, 0.581466, 0.798057};
        else if (L == 8) coef = {0, 0.07778, 0.1625, 0.2547469, 0.3570889, 0.475923, 0.61921, 0.828777};
        else if (L == 9) coef = {0, 0.068333, 0.1428356, 0.22055, 0.306667, 0.736885, 0.516842, 0.65526, 0.8447569};
        else if (L == 10) coef = {0, 0.06166, 0.125667, 0.195537, 0.27069, 0.35214, 0.443844, 0.551359, 0.68066, 0.8680};
    } else {
        for (int i = 0; i < L; i++) {
            coef.at(i) = D * sqrt(i / L / PI);
        }
    }
    
    cout << "Progress is 0%";
    for (int t = 0; t < end_time; t++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
        //BSの座標を設定
        pair<double, double> origin = make_pair(0, 0);
        vector<pair<double, double>> BS_pos(num_BS);
        vector<pair<double, double>> near_BS_pos; //原点に近いセルだけ電力制御することにする．近い座標だけ保存
        //Focused Base station（原点のデバイスから最も近い基地局）
        near_BS_pos.push_back(origin);
        for (int i = 0; i < num_BS; i++) {
            BS_pos.at(i) = coordinate();
            double dst = cal_dst(origin, BS_pos.at(i));
            if (dst < 20.0) {
                near_BS_pos.push_back(BS_pos.at(i));
            }
        }
        int num_near_BS = (int)near_BS_pos.size();  //参照避け
        
        //ノイズを入れるならノイズは送信電力の比にする
        double power = 0.0;
        vector<double> SI(C, 0);
        vector<vector<double>> ACL(C, vector<double>());
        
        //端末情報を初期化
        /* デバイスがすること（原点セルの座標だけ可視化する場合）
         1.原点セルの基地局との距離を計算し，離れてたらただ干渉信号として足す
         〜近い場合
         1.（比較的近い基地局リストの）全基地局との距離を計算して，一番近い基地局を探す．
         2.基地局までの距離から，割り当てられる送信電力を計算する
         3.
         〜原点セルにデバイスがある場合　割り当てられた電力で信号を送信（フェージングは推定できないとする）
         〜原点セル以外の場合　原点セルの基地局への干渉信号として電力を足す
         */
        for (int i = 1; i < num_IoT; i++) {
            pair<double, double> device_pos = coordinate();
            double dst_to_FBS = cal_dst(device_pos, origin);
            int channel = rand() % C;
            
            double H = exp_dist(1.0);
            if (dst_to_FBS > 5.0) {
                SI.at(channel) += theta * pow(theta + 1, L - rand() % (int)L - 1) * H * pow(urand() / dst_to_FBS, alpha);
            } else {
                double dst_to_NBS = dst_to_FBS;
                int cell_num = 0, level = 0; //デバイスから一番近い基地局を探す．リストを抜けたら原点が最近接
                for (int k = 1; k < num_near_BS; k++) {
                    double tmp_dst = cal_dst(device_pos, near_BS_pos.at(k));
                    if (tmp_dst < dst_to_NBS) {
                        dst_to_NBS = tmp_dst;
                        cell_num = k;
                    }
                }
                for (int l = L - 1; ;l--) { //電力レベル決定
                    if (coef.at(l) < dst_to_NBS) {
                        level = l;
                        break;
                    }
                }
                
                power = theta * pow(theta + 1, L - level - 1) * H * pow(dst_to_NBS / dst_to_FBS, alpha);
                //干渉信号和．自分の信号含む
                SI.at(channel) += power;
                if (cell_num == 0) {
                    ACL.at(channel).emplace_back(power);
                }
            }
        }
        
        
        for (int i = 0; i < C; i++) {
            sort(ACL.at(i).rbegin(), ACL.at(i).rend());
            for (int j = 0; j < ACL.at(i).size(); j++) {
                double P = ACL.at(i).at(j);
                SI.at(i) -= P;
                double SINR = P / (SI.at(i) + noise);
                if (SINR > theta) success++;
                else break;
            }
        }
        
        
        //outputfile << success / end_time << endl;
        
        
        //進捗状況を表示
        double progress = t / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
    cout << lambda_IoT << " " << success / end_time << endl;
}




//Point visualization Power allocation
void PA_NOMA_point_visualize(double lambda_IoT, double alpha, double noise, double L, string s, double D) {
    //データ読み込み
    string filename0 = "BS_pos2.txt";
    //vector<pair<double, double>> bspos;
    ifstream readingfile;
    readingfile.open(filename0);
    ifstream ifs(filename0);
        if (ifs.fail()) {
           cerr << "Cannot open file\n";
           exit(0);
        }
    
    //int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
    //BSの座標を設定
    pair<double, double> origin = make_pair(0, 0);
    
    vector<pair<double, double>> near_BS_pos; //原点に近いセルだけ電力制御することにする．近い座標だけ保存
    near_BS_pos.push_back(origin);
    
    double x, y;
    while (ifs >> x >> y) {
        pair<double, double> BS_pos = make_pair(x, y);
        double dst = cal_dst(origin, BS_pos);
        if (dst < 5.0) {
            near_BS_pos.push_back(BS_pos);
        }
    }
    ifs.close();
    int num_near_BS = (int)near_BS_pos.size();  //参照避け
    
    vector<double> coef(L);
    if (D == 0) {
        if (L == 2) coef = {0, 0.47164};
        else if (L == 3) coef = {0, 0.36497, 0.595226};
        else if (L == 4) coef = {0, 0.175, 0.38608, 0.67117};
        else if (L == 5) coef = {0, 0.267434, 0.4051755, 0.545601, 0.720014};
        else if (L == 6) coef = {0, 0.1073693, 0.2292, 0.364867, 0.532356, 0.76407};
        else if (L == 7) coef = {0, 0.090711, 0.189755, 0.299831, 0.428612, 0.581466, 0.798057};
        else if (L == 8) coef = {0, 0.07778, 0.1625, 0.2547469, 0.3570889, 0.475923, 0.61921, 0.828777};
        else if (L == 9) coef = {0, 0.068333, 0.1428356, 0.22055, 0.306667, 0.736885, 0.516842, 0.65526, 0.8447569};
        else if (L == 10) coef = {0, 0.06166, 0.125667, 0.195537, 0.27069, 0.35214, 0.443844, 0.551359, 0.68066, 0.8680};
    } else {
        for (int i = 0; i < L; i++) {
            coef.at(i) = D * sqrt(i / L);
        }
    }
    
    
    double success = 0;
    cout << "path loss exp : " << alpha << endl;
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int t = 0; t < end_time; t++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        vector<terminal> device(num_IoT);
        
        //ノイズを入れるならノイズは送信電力の比にする
        double power = 0.0;
        vector<double> SI(num_near_BS, 0);
        //vector<vector<double>> ACL(num_near_BS, vector<double>());
        vector<vector<pair<double, double>>> ACL(num_near_BS, vector<pair<double, double>>());
        
        //端末情報を初期化
        for (int i = 0; i < num_IoT; i++) {
            pair<double, double> device_pos = coordinate();
            double dst_to_O = cal_dst(origin, device_pos);
            
            double H = exp_dist(1.0);
            if (dst_to_O > 5.0) {
                for (int j = 0; j < num_near_BS; j++) {
                    SI.at(j) += theta * pow(theta + 1, L - rand() % (int)L - 1) * H * pow(urand() / dst_to_O, alpha);
                }
            } else {
                double dst_to_NBS = dst_to_O;
                int cell_num = 0, level = 0; //デバイスから一番近い基地局を探す．リストを抜けたら原点が最近接
                for (int k = 1; k < num_near_BS; k++) {
                    double tmp_dst = cal_dst(device_pos, near_BS_pos.at(k));
                    if (tmp_dst < dst_to_NBS) {
                        dst_to_NBS = tmp_dst;
                        cell_num = k;
                    }
                }
                device.at(i).nearest_BS_id = cell_num;
                for (int l = L - 1; ;l--) { //電力レベル決定
                    if (coef.at(l) < dst_to_NBS) {
                        level = l;
                        break;
                    }
                }
                
                for (int j = 0; j < num_near_BS; j++) {
                    double dst_to_FBS = cal_dst(device_pos, near_BS_pos.at(j));
                    power = theta * pow(theta + 1, L - level - 1) * H * pow(dst_to_NBS / dst_to_FBS, alpha);
                    //干渉信号和．自分の信号含む
                    SI.at(j) += power;
                    ACL.at(j).emplace_back(make_pair(power, i));
                    device.at(i).pos = device_pos;
                }
            }
        }
        
        
        for (int i = 0; i < num_near_BS; i++) {
            sort(ACL.at(i).rbegin(), ACL.at(i).rend());
            for (int j = 0; j < ACL.at(i).size(); j++) {
                int id = ACL.at(i).at(j).second;
                if (device.at(id).nearest_BS_id != i) continue;
                double P = ACL.at(i).at(j).first;
                SI.at(i) -= P;
                double SINR = P / (SI.at(i) + noise);
                if (SINR > theta) {
                    success++;
                    outputfile << device.at(id).pos.first << " " << device.at(id).pos.second << endl;
                } else break;
            }
        }
        
        
        //outputfile << success / end_time << endl;
        
        
        //進捗状況を表示
        double progress = t / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
    cout << lambda_IoT << " " << success / end_time << endl;
}



void SIC_pos(double lambda_IoT, double alpha, double noise) {
    //データ読み込み
    string filename0 = "BS_pos2.txt";
    //vector<pair<double, double>> bspos;
    ifstream readingfile;
    readingfile.open(filename0);
    ifstream ifs(filename0);
        if (ifs.fail()) {
           cerr << "Cannot open file\n";
           exit(0);
        }
    
    //BSの座標を設定
    pair<double, double> origin = make_pair(0, 0);
    
    vector<pair<double, double>> near_BS_pos; //原点に近いセルだけ電力制御することにする．近い座標だけ保存
    near_BS_pos.push_back(origin);
    
    double x, y;
    while (ifs >> x >> y) {
        pair<double, double> BS_pos = make_pair(x, y);
        double dst = cal_dst(origin, BS_pos);
        if (dst < 5.0) {
            near_BS_pos.push_back(BS_pos);
        }
    }
    ifs.close();
    int num_BS = (int)near_BS_pos.size();
    
    double success = 0;
    cout << "path loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int i = 0; i < end_time; i++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        vector<terminal> device(num_IoT);

        vector<double> SI(num_BS, 0);
        vector<vector<pair<double, double>>> ACL(num_BS, vector<pair<double, double>>());
        
        //端末情報を初期化
        for (int i = 1; i < num_IoT; i++) {
            pair<double, double> IoT_pos;
            IoT_pos = coordinate();
            double dst = cal_dst(origin, IoT_pos);
            if (dst > 5.0) continue;
            
            double nd = dst;
            int nd_id = 0;
            for (int j = 0; j < num_BS; j++) {
                //端末-基地局間のフェージング係数を設定
                double H = exp_dist(1.0);//gauss_rand(0, 1);
                double dst2 = cal_dst(IoT_pos, near_BS_pos.at(j));
                if (dst2 < nd) {
                    nd = dst2;
                    nd_id = j;
                }
                double power = H / pow(dst2, alpha);
                ACL.at(j).emplace_back(make_pair(power, i));
                SI.at(j) += power;
                device.at(i).pos = IoT_pos;
            }
            device.at(i).nearest_BS_id = nd_id;
        }
        
        for (int i = 0; i < num_BS; i++) {
            sort(ACL.at(i).rbegin(), ACL.at(i).rend());
            for (int j = 0; j < 2; j++) {
                int id = ACL.at(i).at(j).second;
                if (device.at(id).nearest_BS_id != i) continue; //一番近い基地局へのアクセスだけ見る(重複カウント回避)
                double P = ACL.at(i).at(j).first;
                SI.at(i) -= P;
                double SINR = P / (SI.at(i) + noise + 0.001);
                if (SINR > theta) {
                    success++;
                    if (j == 0) outputfile << device.at(id).pos.first << " " << device.at(id).pos.second << endl;
                    else outputfile1 << device.at(id).pos.first << " " << device.at(id).pos.second << endl;
                } else break;
            }
        }
        
        
        //進捗状況を表示
        double progress = i / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    cout << lambda_IoT << " " << success / end_time << endl;
    
}



//Power allocation , number of device of each level
void av_num_level(double lambda_IoT, double L, double D) {
    vector<double> av_level(L, 0);
    
    bool p_flag[4] = {true, true, true, true};
    vector<double> coef(L);
    if (D == 0) {
        if (L == 2) coef = {0, 0.47164};
        else if (L == 3) coef = {0, 0.36497, 0.595226};
        else if (L == 4) coef = {0, 0.175, 0.38608, 0.67117};
        else if (L == 5) coef = {0, 0.267434, 0.4051755, 0.545601, 0.720014};
        else if (L == 6) coef = {0, 0.1073693, 0.2292, 0.364867, 0.532356, 0.76407};
        else if (L == 7) coef = {0, 0.090711, 0.189755, 0.299831, 0.428612, 0.581466, 0.798057};
        else if (L == 8) coef = {0, 0.07778, 0.1625, 0.2547469, 0.3570889, 0.475923, 0.61921, 0.828777};
        else if (L == 9) coef = {0, 0.068333, 0.1428356, 0.22055, 0.306667, 0.736885, 0.516842, 0.65526, 0.8447569};
        else if (L == 10) coef = {0, 0.06166, 0.125667, 0.195537, 0.27069, 0.35214, 0.443844, 0.551359, 0.68066, 0.8680};
    } else {
        for (int i = 0; i < L; i++) {
            coef.at(i) = D * sqrt(i / L);
        }
    }
    
    cout << "Progress is 0%";
    for (int t = 0; t < end_time; t++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
        //BSの座標を設定
        pair<double, double> origin = make_pair(0, 0);
        vector<pair<double, double>> BS_pos(num_BS);
        vector<pair<double, double>> near_BS_pos; //原点に近いセルだけ電力制御することにする．近い座標だけ保存
        //Focused Base station（原点のデバイスから最も近い基地局）
        near_BS_pos.push_back(origin);
        for (int i = 0; i < num_BS; i++) {
            BS_pos.at(i) = coordinate();
            double dst = cal_dst(origin, BS_pos.at(i));
            if (dst < 20.0) {
                near_BS_pos.push_back(BS_pos.at(i));
            }
        }
        int num_near_BS = (int)near_BS_pos.size();  //参照避け
        vector<int> sum_level(L, 0);
        
        for (int i = 1; i < num_IoT; i++) {
            pair<double, double> device_pos = coordinate();
            double dst_to_FBS = cal_dst(device_pos, origin);
            
            
            if (dst_to_FBS > 20.0) {
                continue;
            } else {
                double dst_to_NBS = dst_to_FBS;
                int level = 0; //デバイスから一番近い基地局を探す．リストを抜けたら原点が最近接
                for (int k = 1; k < num_near_BS; k++) {
                    double tmp_dst = cal_dst(device_pos, near_BS_pos.at(k));
                    if (tmp_dst < dst_to_NBS) {
                        dst_to_NBS = tmp_dst;
                    }
                }
                for (int l = L - 1; ;l--) { //電力レベル決定
                    if (coef.at(l) < dst_to_NBS) {
                        level = l;
                        sum_level.at(l)++;
                        break;
                    }
                }
            }
        }
        
        
        for (int i = 0; i < L; i++) {
            av_level.at(i) += sum_level.at(i) / (double)num_near_BS;
        }
        //outputfile << success / end_time << endl;
        
        
        //進捗状況を表示
        double progress = t / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
    for (int i = 0; i < L; i++) {
        cout << i+1 << " " << av_level.at(i) / end_time << endl;
    }
}



void gen_BS_pos() {
    int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
    string filename = "BS_pos3.txt";
    outputfile.open(filename);
    for (int i = 0; i < num_BS; i++) {
        pair<double, double> pos = coordinate();
        outputfile << pos.first << " " << pos.second << endl;
    }
}




int main() {
    srand((unsigned)time(NULL));
    //success_dist();
    //PA_NOMA_thp(1, 4.5, 0, 5);
    
    cout << "0:ALOHA, 1:SIC, 2:Power control, 3:Throughput (PA), 4:Distribution, " << endl;
    cout << "5:PA pos, 6:SIC pos, 7:AV of level, 8:Gen BS pos : ";
    double key; cin >> key; cout << endl;
    
    cout << "Put start lambda IoT : ";
    double k; cin >> k;

    string filename;

    double alpha = 4.5;
    if (key == 0) {
        filename = to_string((int)k) + "_ALOHA_dst.txt";
        outputfile.open(filename);
        ALOHA_pr(k, alpha, 0);
    } else if (key == 1) {
        filename = to_string((int)k) + "_SIC_dst.txt";
        outputfile.open(filename);
        SIC_pr(k, alpha, 0);
    } else if (key == 2) {
        double L; cout << "Power level : "; cin >> L;
        int C; cout << "Channel num : "; cin >> C;
        string s; cout << "Name : "; cin >> s;
        double D; cout << "Voronoi Radius :"; cin >> D; cout << endl;
        filename = to_string((int)k) + "_NOMA_" + to_string((int)L) + "_" + to_string((int)C) + s + ".txt";
        outputfile.open(filename);
        PA_NOMA_pr_dst(k, alpha, 0, L, C, s, D);
    } else if (key == 3) { //Throughput
        double L; cout << "Power level : "; cin >> L;
        int C; cout << "Channel num : "; cin >> C;
        double D; cout << "Voronoi Radius :"; cin >> D; cout << endl;
        filename = to_string((int)k) + "_NOMA_" + to_string((int)L) + "_" + to_string(D) + "_thp.txt";
        outputfile.open(filename);
        PA_NOMA_thp(k, alpha, 0, L, C, D);
    } else if (key == 4) {
        double L; cout << "Power level : "; cin >> L;
        int C; cout << "Channel num : "; cin >> C;
        string s; cout << "string : "; cin >> s; cout << endl;
        success_dist(k, L, C, s);
    } else if (key == 5) {
        double L; cout << "Power level : "; cin >> L;
        double D; cout << "Voronoi Radius :"; cin >> D;
        string s; cout << "Name : "; cin >> s; cout << endl;
        filename = to_string((int)k) + "_PAPoint_" + to_string((int)L) + "_" + s + ".txt";
        outputfile.open(filename);
        PA_NOMA_point_visualize(k, alpha, 0, L, s, D);
    } else if (key == 6) {
        string s; cout << "Name : "; cin >> s; cout << endl;
        filename = to_string((int)k) + "_SICPoint_" + s;
        outputfile.open(filename+".txt");
        outputfile1.open(filename+"second.txt");
        SIC_pos(k, alpha, 0);
        outputfile1.close();
    } else if (key == 7) {
        double L; cout << "Power level : "; cin >> L;
        double D; cout << "Voronoi Radius :"; cin >> D; cout << endl;
        av_num_level(k, L, D);
    } else if (key == 8) {
        gen_BS_pos();
    }


    outputfile.close();
    
}

