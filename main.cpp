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
    pair<double, double> pos;
    //double to_BS[num_BS];
    double dst_to_BS;
    double dst_to_FBS;
    int nearest_BS;
    double level;
    double coef;
    double IF_to_FBS;
    bool state; //True: in , False: out
    bool operator<( const terminal& right ) const {
            return coef == right.coef ? true : coef < right.coef;
    }
};


double urand(){
    double m, a;
    m = RAND_MAX + 1.0;
    a = (rand() + 0.5)/m;
    a = (rand() + a)/m;
    return (rand() + a)/m;
}

//乱数発生
//double my_rand(double Min, double Max) {
//    mt19937 mt{ random_device{}() };
//    //random_device rd;
//    //default_random_engine eng(rd());
//    uniform_real_distribution<double> distr(Min, Max);
//    return distr(mt);
//}

//円形領域内の座標をランダムに取得
pair<double, double> coordinate() {
    //init_genrand(i);
    double r = radius * sqrt(urand()); //[0, 1)
    double theta = PI * (2. * urand() - 1.); //(-1, 1)
    while (theta < -PI || PI < theta) {
        theta = PI * (2. * urand() - 1.);
    }
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


void ALOHA_pr(double lambda_IoT, double alpha, double noise) {
    //double success = 0;
    cout << "pass loss exp : " << alpha << endl;
    outputfile << "pass loss exp : " << alpha << endl;
    
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
        //端末情報を初期化
        for (int i = 0; i < num_IoT; i++) {
            pair<double, double> IoT_pos;
            if (i > 0) {
                IoT_pos = coordinate();
            }
            //端末-基地局間のフェージング係数を設定
            double H = exp_dist(1.0);//gauss_rand(0, 1);
            if (i == 0) {
                coef = H / pow(nd, alpha);
            } else {
                double dst2 = cal_dst(IoT_pos, nearest_pos);
                SI += H / pow(dst2, alpha);
            }
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
    cout << "pass loss exp : " << alpha << endl;
    outputfile << "pass loss exp : " << alpha << endl;
    
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
        //端末情報を初期化
        for (int i = 0; i < num_IoT; i++) {
            pair<double, double> IoT_pos;
            if (i > 0) {
                IoT_pos = coordinate();
            }
            //double dst = cal_dst(origin, IoT_pos);
            //if (radius < dst) continue;
            //端末-基地局間のフェージング係数を設定
            double H = exp_dist(1.0);//gauss_rand(0, 1);
            if (i == 0) {
                coef = H / pow(nd, alpha);
                greatest = coef;
            } else {
                double dst2 = cal_dst(IoT_pos, nearest_pos);
                double tmp = H / pow(dst2, alpha);
                if (tmp > greatest) greatest = tmp;
                SI += tmp;
            }
        }
        
        double SINR;
        int flag = 0;
        if (coef >= greatest) {
            SINR = coef / (SI + noise);
            if (SINR > theta) flag = 1;
        } else {
            SINR = greatest / (SI - greatest + coef + noise);
            if (SINR > theta) {
                SINR = coef / (SI + noise);
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


//Power allocation
void PA_NOMA_pr(double lambda_IoT, double alpha, double noise, double L) {
    double success = 0;
    cout << "pass loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    double R_L3[] = {0, 0.517848, 0.780638};
    double R_L4[] = {0, 0.448469, 0.659178, 0.839081};
    double R_L5[] = {0, 0.401123, 0.581634, 0.730387, 0.864729};
    vector<double> Radius;
    for (int i = 0; i <= L; i++) {
        Radius.push_back(sqrt((double)i / L));
//        if (L == 3) {
//            Radius.push_back(R_L3[i]);
//        } else if (L == 4) {
//            Radius.push_back(R_L4[i]);
//        } else if (L == 5) {
//            Radius.push_back(R_L5[i]);
//        } else {
//            cout << "Set power level 3 ~ 5" << endl;
//            return;
//        }
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
        double nd = inf;
        int index = 0;
        for (int i = 0; i < num_BS; i++) {
            BS_pos.at(i) = coordinate();
            double dst = cal_dst(origin, BS_pos.at(i));
            if (dst < 3.0) {
                near_BS_pos.push_back(BS_pos.at(i));
                if (dst < nd) {
                    nd = dst;
                    FBS_pos = BS_pos.at(i);
                    index = i;
                }
            }
        }
        BS_pos.at(index) = BS_pos.at(0);
        BS_pos.at(0) = FBS_pos;
        int num_near_BS = (int)near_BS_pos.size();  //参照避けによる高速化
        
        //ノイズを入れるならノイズは送信電力の比にする
        double SI = 0, coef = 0.0;
        vector<double> LV_SI(L, 0);  //レベル毎に干渉信号を足して，SIC適用後の計算を可能にする
        vector<int> acl; //原点セル基地局へのアクセスリスト
        vector<terminal> device(num_IoT);
        device.at(0).pos = origin;
        //端末情報を初期化
        for (int i = 1; i < num_IoT; i++) {
            pair<double, double> device_pos = coordinate();
            device.at(i).pos = device_pos;
            double dst_to_FBS = cal_dst(device.at(i).pos, FBS_pos);
            device.at(i).dst_to_FBS = dst_to_FBS;
            
            /* デバイスがすること（原点セルの座標だけ可視化する場合）
             1.原点セルの基地局との距離を計算し，離れてたらただ干渉信号として足す
             〜近い場合
             1.（比較的近い基地局リストの）全基地局との距離を計算して，一番近い基地局を探す．
             2.基地局までの距離から，割り当てられる送信電力を計算する
             3.
             〜原点セルにデバイスがある場合　割り当てられた電力で信号を足す（フェージングは推定できないとする）
             〜原点セル以外の場合　原点セルの基地局への干渉信号として電力を足す
             */
            
            double H = exp_dist(1.0);
            if (dst_to_FBS > 5.0) {
                SI += H / pow(dst_to_FBS, alpha);
            } else {
                double dst_to_NBS = inf;
                int ind = 0, level = 0;
                for (int k = 0; k < num_near_BS; k++) {
                    double tmp_dst = cal_dst(device_pos, near_BS_pos.at(k));
                    if (tmp_dst < dst_to_NBS) {
                        dst_to_NBS = tmp_dst;
                        ind = k;
                    }
                }
                for (int l = L - 1; ;l--) { //電力レベル決定
                    if (Radius.at(l) < dst_to_NBS) {
                        level = l;
                        break;
                    }
                }
                if (ind == 0) { //原点セルにデバイスがある場合
                    acl.push_back(i);
                    coef = theta * pow(theta + 1, L - level - 1) * H;
                    LV_SI.at(level) += coef;
                    device.at(i).coef = coef;
                    device.at(i).level = level;
                } else {
                    SI += theta * pow(theta + 1, L - level - 1) * H / pow(dst_to_FBS, alpha);
                }
            }
            
        }
        
        
        int s = (int)acl.size();
        for (int i = 0; i < s; i++) {
            int f = acl.at(i);
            double P = device.at(f).coef;
            int l = device.at(f).level;
            double si = 0;
            for (int j=l; j < L; j++) {   //j=l; j < L; j++ じゃない？
                si += LV_SI.at(j);
            }
            double SINR = P / (SI + si - P + noise);
            if (SINR > theta) success++;
            else break;
        }
        //cout << s << endl;
        
        //outputfile << nd << " " << flag << endl;
        
        
        //進捗状況を表示
        double progress = t / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    cout << success / end_time << endl;
    
}


int main() {
    cout << "ALOHA 0, SIC 1, Power control 2 :";
    double key; cin >> key; cout << endl;
//    cout << "Power allocate 0, Random any :";
//    double op; cin >> op; cout << endl;
    cout << "Put start lambda IoT : ";
    double k; cin >> k; cout << endl;

    string filename;
    if (key == 0) filename = to_string(k) + "ALOHA_dst.txt";
    else if (key == 1) filename = to_string(k) + "SIC_dst.txt";
    else {
        string filename2 = "Second_decode_pos.txt";
        string filename1 = "First_decode_pos.txt";
        outputfile2.open(filename2);
        outputfile1.open(filename1);
    }

    outputfile.open(filename);
    
    for (double L = 3.0; L <= 5; L++) {
        double alpha = 4.5;
        if (key == 0) ALOHA_pr(k, alpha, 0);
        else if (key == 1) SIC_pr(k, alpha, 0);
        else if (key == 2) PA_NOMA_pr(k, alpha, 0, L);
    }


    outputfile1.close();
    outputfile2.close();
    outputfile.close();
    
}

