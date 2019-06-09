#include <iostream>
#include <sstream>
#include <string>
#include "quadprog++.hh"
#include <fstream>
#include <cmath>

int main (int argc, char *const argv[]) {
	double G[MATRIX_DIM][MATRIX_DIM], g0[MATRIX_DIM], 
		CE[MATRIX_DIM][MATRIX_DIM], ce0[MATRIX_DIM], 
		CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM], 
		alpha[MATRIX_DIM];
	int n, m, p;
    int N;//データの次元数
    int data_n;//データ分割数
    int total_data_size = 0;//総データ数。学習用データと評価用データに別れる。
    int estimate_data_size = 0;//評価用データ数
    int training_data_size = 0;//学習用データ数
    int kernel = 0;//kernelの選択。0:内積 1:多項式カーネル 2:ガウスカーネル
    std::string file_name;//データファイルの名前
    double kernel_result(double* x, double* y, int kernel, int degree, double sigma, int d);//カーネルの計算結果を出力する関数

    double sigma = 10;//ガウスカーネルで用いるシグマ定数
    int d = 2;//多項式カーネルで用いるd定数
    
    double data[MATRIX_DIM][MATRIX_DIM] = {};    //データ格納先
    double label[MATRIX_DIM] = {};               //ラベル格納先
    double theta = 0;//閾値

    //ファイルを指定する
    std::cout << "データファイルを指定してください : ";
    std::cin >> file_name;
    //次元数を指定する
    std::cout << "データの次元数を指定してください : ";
    std::cin >> N;
    //カーネルを指定する
    std::cout << "使用するカーネルを指定してください(0:カーネルトリックなし　1:多項式カーネル　2:ガウスカーネル) : ";
    std::cin >> kernel;
    if(!(kernel == 0 || kernel == 1 || kernel == 2)){
        std::cout << "正しくカーネルを指定してください(0:カーネルトリックなし　1:多項式カーネル　2:ガウスカーネル)" << std::endl;
        exit(0);
    }
    if(kernel == 2){
        std::cout << "ガウスカーネルで使用するシグマ定数を指定してください : ";
        std::cin >> sigma;
    }
    //交差検定を行う際に何分割するか決定する
    std::cout << "交差検定を行う際のデータ分割数を指定してください : ";
    std::cin >> data_n;
    if(data_n < 0){
        std::cout << "データ分割数は正整数で指定してください" << std::endl;
        exit(0);
    }

    //ファイルを読み取る。
    std::ifstream ifs(file_name);
    if(ifs.fail()){
        std::cout << "ファイルが見つかりません : " << file_name << std::endl;
        exit(0);
    }

    //データを読み取り、dataとlabelに入れていく。
    {
    int index = 0;//今扱っているdataの番号
    std::string str;
    double num;
    char c;
    while (getline(ifs, str)) {
        std::istringstream iss(str);
        for(int i = 0; i < N; i++){
            iss >> num >> c;
            data[index][i] = num;
        }
        iss >> label[index];
        index++;
        //データの重複を排除する
        for(int i = 0; i < index; i++){
            if(data[index] == data[i]){
                index--;
                break;
            }
        }
        total_data_size++;
    }
    }
    
    estimate_data_size = total_data_size / data_n;//評価データの要素数
    training_data_size = total_data_size - estimate_data_size;//学習データの要素数

    n = training_data_size; //データの数
    m = n;
    p = 1;

    //g0の設定(全部-1)
    {
    for(int i = 0; i < n; i++){
        g0[i] = -1;
    }
    }

    //CEの設定（(label[0] label[1]...)）
    {
    for(int i = 0; i < n; i++){
        for(int j = 0; j < p; j++){
	        CE[i][j] = label[i];
        }
    }
    }

    //ce0の設定（0）
    {
    for (int i = 0; i < p; i++){
      ce0[i] = 0;
    }
    }

    //CIの設定（単位行列）
    {
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
	        if(i == j){
	            CI[i][j] = 1;
	        }else{
	            CI[i][j] = 0;
            }
        }
	}
    }

    //ci0の設定（全部0)
    {
    for (int i = 0; i < m; i++){
      ci0[i] = 0;
    }
    }

    //分割数に応じてdataをブロックごとにわけて、どのブロックを学習、どのブロックを評価という風にわけていく
    {
    int block_size = total_data_size / data_n;//ブロックに要素がいくつ入っているか

    //data配列のラベルの次の要素にブロック番号を表す配列を追加する。
    //日付順に並んでいるデータなど、順番に意味がある場合にはランダムにわりつけるべきだが、今回はそうでないため前から順に割り付けている。
    for(int i = 0; i < total_data_size; i++){
        int block_num = i / block_size;//ブロック番号
        data[i][N+1] = block_num;
    }
    }

    //交差検証開始
    {
    double accuracy[data_n];    //正解率のこと。予測結果全体と、答えがどれぐらい一致しているかを判断する指標。
    double precision[data_n];   //適合率のこと。予測を正と判断した中で、答えも正のもの。
    double recall[data_n];      //再現率のこと。答えが正の中で、予測が正とされたもの。
    double F_measure[data_n];   //F値のこと。予測精度の評価指標。PresicionとRecallの調和平均。
    for(int estimate_data_num = 0; estimate_data_num < data_n; estimate_data_num++){
        //estimate_data_num : 評価に使うブロック番号

        double training_data[MATRIX_DIM][MATRIX_DIM];   //このサイクルで使う学習データ
        double training_data_label[MATRIX_DIM];         //このサイクルで使う学習データのラベル
        double estimate_data[MATRIX_DIM][MATRIX_DIM];   //このサイクルで使う評価データ
        double estimate_data_label[MATRIX_DIM];         //このサイクルで使う評価データのラベル

        //データを実際に使用するtraining_dataに移す
        int training_now = 0;   //今移しているtraining_dataの番号
        int estimate_now = 0;   //今写しているestimate_dataの番号
        for(int i = 0; i < total_data_size; i++){
            if(data[i][N+1] == estimate_data_num){//データのブロック割付が評価ブロックならば
                for(int j = 0; j < N; j++){
                    estimate_data[estimate_now][j] = data[i][j];
                    estimate_data_label[estimate_now] = label[i];
                }
                estimate_now++;
            }else{//データのブロック割付が学習ブロックならば
                for(int j = 0; j < N; j++){
                    training_data[training_now][j] = data[i][j];
                    training_data_label[training_now] = label[i];
                }
                training_now++;
            }
        }

        //行列Gを計算する。その際行列を半正定値行列にするために対角線上成分に微小値を加える。
        //カーネルによって微小値を変化させている。
        {
        double delta = 0;//微小値
        if(kernel == 1){    //ガウスカーネル
            delta = 1.0e-7;
        }else{              //それ以外
            delta = 1.0e-9;
        }

        //training_dataを用いてGを計算
        for(int i = 0; i < training_data_size; i++){
            for(int j = 0; j < training_data_size; j++){
                G[i][j] = training_data_label[i] * training_data_label[j] * kernel_result(training_data[i], training_data[j], kernel, N, sigma, d);
                if(i == j){
                    G[i][j] += delta;
                }
                
            }
        }
        }

        //solve_quadprogの使用
        solve_quadprog(G, g0, n, CE, ce0, p, CI, ci0, m, alpha);

        double maxA = alpha[0];//alphaの最大値
        int maxAindex = 0;//alphaの最大値のインデックス
        for(int i = 0; i < training_data_size; i++){
            //std::cout << "alpha[" << i << "] : " << alpha[i] << std::endl;
            if(alpha[i] > maxA){
                maxA = alpha[i];
                maxAindex = i;
            }
        }

        
        // 重みを計算する
        double weight[MATRIX_DIM] = {};
        for(int i = 0; i < N; i++){
            for(int j = 0; j < training_data_size; j++){
                weight[i] += alpha[j] * training_data_label[j] * training_data[j][i];
            }
        }

        //閾値を計算する
        for(int i = 0; i < training_data_size; i++){
            theta += alpha[i] * training_data_label[i] * kernel_result(training_data[i], training_data[maxAindex], kernel, N, sigma, d);
        }
        theta -= training_data_label[maxAindex];

        //重みと閾値を用いて評価データのデータを評価する。
        double true_positive = 0;   //真陽性
        double true_negative = 0;   //真陰性
        double false_positive = 0;  //偽陽性
        double false_negative = 0;  //偽陰性
        //評価データをひとつずつ評価していく。
        for(int i = 0; i < estimate_data_size; i++){
            double sum = 0;
            bool result;
            for(int k = 0; k < training_data_size; k++){
                sum += alpha[k] * training_data_label[k] * kernel_result(training_data[k], estimate_data[i], kernel, N, sigma, d);
            }
            sum -= theta;
            if(sum > 0){
                result = true;
            }else{
                result = false;
            }

            if(result){
                if(estimate_data_label[i] == 1.0){
                    true_positive++;
                }else{
                    false_positive++;
                }
            }else{
                if(estimate_data_label[i] == 1.0){
                    false_negative++;
                }else{
                    true_negative++;
                }
            }
        }

        //結果を格納
        if((true_positive + false_positive) != 0){
            accuracy[estimate_data_num] = (true_negative + true_positive) / (true_negative + true_positive + false_negative + false_positive);
            precision[estimate_data_num] = true_positive / (true_positive + false_positive);
            recall[estimate_data_num] = true_positive / (true_positive + false_negative);
            F_measure[estimate_data_num] = 2 * recall[estimate_data_num] * precision[estimate_data_num] / (recall[estimate_data_num] + precision[estimate_data_num]);
        }else{
            accuracy[estimate_data_num] = -1;
            precision[estimate_data_num] = -1;
            recall[estimate_data_num] = -1;
            F_measure[estimate_data_num] = -1;
        }
        
        
    }

    //結果の平均を計算
    double accuracy_ave = 0;
    double precision_ave = 0;
    double recall_ave = 0;
    double F_measure_ave = 0;
    double true_data_num = data_n;//正しく評価値が計算できたブロックの数
    for(int i = 0; i < data_n; i++){
        if(accuracy[i] == -1){
            true_data_num--;
        }else{
            accuracy_ave += accuracy[i];
            precision_ave += precision[i];
            recall_ave += recall[i];
            F_measure_ave += F_measure[i];
        }
    }
    accuracy_ave /= true_data_num;
    precision_ave /= true_data_num;
    recall_ave /= true_data_num;
    F_measure_ave /= true_data_num;

    //結果を表示
    for(int i = 0; i < data_n; i++){
        std::cout << "estimate_block : " << i << std::endl;
        std::cout << "accuracy : " << accuracy[i] << std::endl;
        std::cout << "precision : " << precision[i] << std::endl;
        std::cout << "recall : " << recall[i] << std::endl;
        std::cout << "F_measure : " << F_measure[i] << std::endl;
        std::cout << "......................................" << std::endl;
    }
    std::cout << "Average" << std::endl;
    std::cout << "accuracy : " << accuracy_ave << std::endl;
    std::cout << "precision : " << precision_ave << std::endl;
    std::cout << "recall : " << recall_ave << std::endl;
    std::cout << "F_measure : " << F_measure_ave << std::endl;
    }

    return 0;
}

//カーネルを計算する関数
/*
double* x,y :計算する二つの配列
int kernel  :使用するカーネルを表す変数0:内積 1:多項式カーネル 2:ガウスカーネル)
int degree  :データの次元数。大域変数ではNにあたる
double sigma:ガウスカーネル計算時に使用するシグマ定数
*/
double kernel_result(double* x, double* y, int kernel, int degree, double sigma, int d){
    double get_norm(double* x, double* y, int degree);

    double result = 0;//結果を格納し、この変数をリターンする
    if(kernel == 0){//カーネルトリックなし（内積）
        for(int i = 0; i < degree; i++){
            result += x[i] * y[i];
        }
    }else if(kernel == 1){//多項式カーネル
        double inner_product = 0;
        for(int i = 0; i < degree; i++){
            inner_product += x[i] * y[i];
        }
        result = pow(1 + inner_product, d);
    }else if(kernel == 2){//ガウスカーネル
        double norm = get_norm(x, y, degree);
        result = exp(0 - norm / (2 * pow(sigma, 2)));
    }else{
        exit(1);
    }
    return result;
}

//二乗ノルムを得るための関数
/*
double* x, y    :計算する二つの配列
int degree      :データの次元数
*/
double get_norm(double* x, double* y, int degree){
    double result = 0;
    for(int i = 0; i < degree; i++){
        result += pow(x[i]-y[i], 2);
    }
    return result;
}