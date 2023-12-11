#include <bits/stdc++.h>
//#include "MT.h"

using namespace std;
float inf = numeric_limits<float>::infinity();
using ll = long long;
const double radius = 100;
const double GAMMA = pow(10, 6.0 / 10.0);
const double end_time = 1000;
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
double my_rand(double Min, double Max) {
    mt19937 mt{ random_device{}() };
    //random_device rd;
    //default_random_engine eng(rd());
    uniform_real_distribution<double> distr(Min, Max);
    return distr(mt);
}

//円形領域内の座標をランダムに取得
pair<double, double> coordinate() {
    double x, y;
    while (true) {
        x = urand();
        y = urand();
        double r_sq = x * x + y * y;
        if (r_sq < 1.) break;
    }
    x *= radius;
    y *= radius;
    
//    double r = radius * sqrt(urand()); //[0, 1)
//    double theta = PI * (2. * urand() - 1.); //(-1, 1)
//    while (theta < -PI || PI < theta) {
//        theta = PI * (2. * urand() - 1.);
//    }
//    double x = r * cos(theta);
//    double y = r * sin(theta);
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
    cout << "pass loss exp : " << alpha << endl;
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


//Power allocation
void PA_NOMA_pr(double lambda_IoT, double alpha, double noise, double L, int C) {
    double success = 0;
    cout << "pass loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
//    double R_L3[] = {0, 0.517848, 0.780638};
//    double R_L4[] = {0, 0.448469, 0.659178, 0.839081};
//    double R_L5[] = {0, 0.401123, 0.581634, 0.730387, 0.864729};
//    vector<double> Radius;
//    for (int i = 0; i <= L; i++) {
//        Radius.push_back(sqrt((double)i / L));
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
    //}
    
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
        double SI = 0, si = 0, coef = 0.0;
        vector<double> acl(L, 0); //原点セル基地局へのアクセスリスト．0以外の値が既に格納されていれば衝突
        vector<vector<double>> ACL(C, vector<double>(L, 0));
        bool fail_flag = false;
        
        vector<double> IF_in_cell(C, 0);
        vector<int> FAIL_FLAG(C, L);
        
        //原点デバイス
        double dst_0 = cal_dst(origin, FBS_pos), level_d0;
        for (int l = L - 1; ;l--) { //電力レベル決定
            if (sqrt((double)l / L) < dst_0) { //Radius.at(l) < dst_to_NBS
                level_d0 = l;
                break;
            }
        }
        double P = theta * pow(theta + 1, L - level_d0 - 1) * exp_dist(1.0);
        int cha_d0 = rand() % C;
        ACL.at(cha_d0).at(level_d0) = P;
        
        //端末情報を初期化
        for (int i = 1; i < num_IoT; i++) {
            pair<double, double> device_pos = coordinate();
            double dst_to_FBS = cal_dst(device_pos, FBS_pos);
            int channel = rand() % C;
            
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
                int cell_num = 0, level = 0;
                for (int k = 0; k < num_near_BS; k++) {
                    double tmp_dst = cal_dst(device_pos, near_BS_pos.at(k));
                    if (tmp_dst < dst_to_NBS) {
                        dst_to_NBS = tmp_dst;
                        cell_num = k;
                    }
                }
                for (int l = L - 1; ;l--) { //電力レベル決定
                    if (sqrt((double)l / L) < dst_to_NBS) { //Radius.at(l) < dst_to_NBS
                        level = l;
                        break;
                    }
                }
                if (cell_num == 0 && channel == cha_d0) { //原点セルにデバイスがある，かつチャネルが同じ場合
                    coef = theta * pow(theta + 1, L - level - 1) * H;
                    if (level < level_d0) { //原点デバイスより高いレベルの信号がある場合
                        if (ACL.at(channel).at(level) == 0) ACL.at(channel).at(level) = coef; //(acl.at(level) == 0) acl.at(level) = coef;
                        else {
                            fail_flag = true; //高レベルで衝突，送信失敗
                            break;
                        }
                    } else si += coef;  //原点セル，下位電力レベルの干渉信号和
                } else {
                    SI += theta * pow(theta + 1, L - level - 1) * H * pow(dst_to_NBS, alpha) / pow(dst_to_FBS, alpha);
                }
            }
            
        }
        
        int flag = 0;
        if (!fail_flag) {
            double SINR = P / (SI / (double)C + si + noise); //-Pしてた...
            if (SINR > theta) {
                flag = 1;
                success++;
            }
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
    
}




//Throughput of Power allocation
void PA_NOMA_thp(double lambda_IoT, double alpha, double noise, double L, int C) {
    double success = 0;
    cout << "pass loss exp : " << alpha << endl;
    
    bool p_flag[4] = {true, true, true, true};
    
    cout << "Progress is 0%";
    for (int t = 0; t < end_time; t++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
        //BSの座標を設定
        pair<double, double> origin = make_pair(0, 0);
        vector<pair<double, double>> BS_pos(num_BS);
        vector<pair<double, double>> near_BS_pos; //原点に近いセルだけ電力制御することにする．近い座標だけ保存
        //pair<double, double> FBS_pos; //Focused Base station（原点のデバイスから最も近い基地局）
        near_BS_pos.push_back(origin);
        for (int i = 0; i < num_BS; i++) {
            BS_pos.at(i) = coordinate();
            double dst = cal_dst(origin, BS_pos.at(i));
            if (dst < 3.0) {
                near_BS_pos.push_back(BS_pos.at(i));
            }
        }
        int num_near_BS = (int)near_BS_pos.size();  //参照避けによる高速化
        
        //ノイズを入れるならノイズは送信電力の比にする
        double SI = 0, coef = 0.0;
        vector<double> IF_in_cell(C, 0);
        vector<vector<vector<double>>> ACL(C, vector<vector<double>>(L, vector<double>(1000, 0)));
        vector<int> FAIL_FLAG(C, L);
        //vector<double> acl(L, 0); //原点セル基地局へのアクセスリスト．0以外の値が既に格納されていれば衝突
        //int fail_flag = L;
        
        //vector<double> LV_SI(L, 0);  //レベル毎に干渉信号を足して，SIC適用後の計算を可能にする
        
        //端末情報を初期化
        for (int i = 1; i < num_IoT; i++) {
            pair<double, double> device_pos = coordinate();
            double dst_to_FBS = cal_dst(device_pos, origin);
            int channel = rand() % C;
            
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
                SI += H / pow(0.5 / dst_to_FBS, alpha);
            } else {
                double dst_to_NBS = dst_to_FBS;
                int cell_num = 0, level = 0;
                for (int k = 1; k < num_near_BS; k++) {
                    double tmp_dst = cal_dst(device_pos, near_BS_pos.at(k));
                    if (tmp_dst < dst_to_NBS) {
                        dst_to_NBS = tmp_dst;
                        cell_num = k;
                        break;
                    }
                }
                for (int l = L - 1; ;l--) { //電力レベル決定
                    if (sqrt((double)l / L) < dst_to_NBS) { //Radius.at(l) < dst_to_NBS
                        level = l;
                        break;
                    }
                }
                
                if (cell_num == 0) { //原点セルにデバイスがある場合
                    coef = theta * pow(theta + 1, L - level - 1) * H;
                    IF_in_cell.at(channel) += coef;  //原点セルの干渉信号和．自分の信号含む
                    int m = 0;
                    while (ACL.at(channel).at(level).at(m) != 0) m++;
                    ACL.at(channel).at(level).at(m) = coef;
                    
//                    if (ACL.at(channel).at(level).at(0) == 0) ACL.at(channel).at(level).at(0) = coef; //(acl.at(level) == 0) acl.at(level) = coef;
//                    else {
//                        FAIL_FLAG.at(channel) = level;
//                        //fail_flag = level; //このレベルの衝突フラグを立てる
//                    }
                } else {
                    SI += theta * pow(theta + 1, L - level - 1) * H * pow(dst_to_NBS, alpha) / pow(dst_to_FBS, alpha);
                }
            }
            
        }
        
        
        for (int i = 0; i < C; i++) {
            for (int j = 0; j < L; j++) {
                double sub = 0, tmp = 0;
                for (int k = 0; ; k++) {
                    if (ACL.at(i).at(j).at(k) == 0) {
                        IF_in_cell.at(i) -= sub;
                        break;
                    }
                    double P = ACL.at(i).at(j).at(k);
                    double SINR = P / (SI / (double)C + IF_in_cell.at(i) - P + noise);
                    if (SINR > theta) {
                        success++;
                        if (tmp < P) tmp = P;
                        sub = tmp;
                    }
                }
            }
        }
        
        
//        for (int j = 0; j < C; j++) {
//            int FF = FAIL_FLAG.at(j);
//            for (int i = 0; i < FF; i++) {
//                double P  = ACL.at(j).at(i).at(0); //チャネルj，レベルi の電力 //acl.at(i);
//                IF_in_cell.at(j) -= P;//si -= P;
//                double SINR = P / (SI / (double)C + IF_in_cell.at(j) + noise);
//                //cout << SINR << endl;
//                if (SINR > theta) success++;
//                else break;
//            }
//        }
        
        
        //cout << s << endl;
        
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


//送信成功確率 分布作成
void success_dist() {
    //データ読み込み
    string filename0 = "3_NOMA_52.txt";
    vector<pair<double, double>> suc_dist;
    ifstream readingfile;
    readingfile.open(filename0);
    ifstream ifs("3_NOMA_52.txt");
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
        cout << 0.05 * i << " " << suc_sum.at(i-1) / total.at(i-1) << endl;
    }
    
}



int main() {
    srand((unsigned)time(NULL));
    //success_dist();
    //PA_NOMA_thp(1, 4.5, 0, 5);
    
//    for (int i = 1; i <= 5; i+=2) {
//        string filename;
//        filename = to_string(i) + "_NOMA_L5_dst.txt";
//        outputfile.open(filename);
//        PA_NOMA_pr((double)i, 4.5, 0, 5);
//        outputfile.close();
//    }
//    for (int i = 1; i <= 5; i+=2) {
//        string filename;
//        filename = to_string(i) + "_NOMA_L1_dst.txt";
//        outputfile.open(filename);
//        PA_NOMA_pr((double)i, 4.5, 0, 1);
//        outputfile.close();
//    }
//    for (int i = 1; i <= 5; i+=2) {
//        string filename;
//        filename = to_string(i) + "_NOMA_L3_dst.txt";
//        outputfile.open(filename);
//        PA_NOMA_pr((double)i, 4.5, 0, 3);
//        outputfile.close();
//    }

    

    
    cout << "0:ALOHA, 1:SIC, 2:Power control, 3:Throughput (PA), 4:Distribution :";
    double key; cin >> key; cout << endl;
//    cout << "Power allocate 0, Random any :";
//    double op; cin >> op; cout << endl;
    if (key == 4) {
        success_dist();
        return 0;
    }
    cout << "Put start lambda IoT : ";
    double k; cin >> k; cout << endl;

    string filename;

    outputfile.open(filename);
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
        double L; cout << "Power level : "; cin >> L; cout << endl;
        int C; cout << "Channel num : "; cin >> C; cout << endl;
        filename = to_string((int)k) + "_NOMA_" + to_string((int)L) + ".txt";
        outputfile.open(filename);
        PA_NOMA_pr(k, alpha, 0, L, C);
    } else if (key == 3) { //Throughput
        double L; cout << "Power level : "; cin >> L; cout << endl;
        int C; cout << "Channel num : "; cin >> C; cout << endl;
        //filename = to_string((int)k) + "_NOMA_" + to_string((int)L) + "thp.txt";
        //outputfile.open(filename);
        PA_NOMA_thp(k, alpha, 0, L, C);
    }


    outputfile.close();
    
}

