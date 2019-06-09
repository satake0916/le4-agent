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
    int kernel = 0;//kernelの選択。0:内積 1:多項式カーネル 2:ガウスカーネル
    std::string file_name;//データファイルの名前
    double kernel_result(double* x, double* y, int kernel, int degree, double sigma);//カーネルの計算結果を出力する関数

    double sigma = 10;//ガウスカーネルで用いるシグマ定数
    double maxX = -100000000;//データのx変数の最低値。データが二次元のときのみ可視化に用いる。
    double maxY = -100000000;//データのy変数の最低値。データが二次元のときのみ可視化に用いる。
    double minX = 1000000000;//データのx変数の最高値。データが二次元のときのみ可視化に用いる。
    double minY = 1000000000;//データのy変数の最高値。データが二次元のときのみ可視化に用いる。
    
    double data[MATRIX_DIM][MATRIX_DIM] = {};    //データ格納先
    double label[MATRIX_DIM] = {};               //ラベル格納先
    int index = 0;//データ数
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
    if(kernel < 0 || kernel > 2){
        std::cout << "正しくカーネルを指定してください(0:カーネルトリックなし　1:多項式カーネル　2:ガウスカーネル)" << std::endl;
        exit(0);
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

    //データを読み取り、dataとlabelに入れていく。maxX,maxY,minX,minYの設定もおこなっている。
    std::string str;
    double num;
    char c;
    while (getline(ifs, str)) {
        std::istringstream iss(str);
        for(int i = 0; i < N; i++){
            iss >> num >> c;
            data[index][i] = num;
            if(i == 0){
                if(maxX < num){
                    maxX = num;
                }else if(minX > num){
                    minX = num;
                }
            }else if(i == 1){
                if(maxY < num){
                    maxY = num;
                }else if(minY > num){
                    minY = num;
                }
            }
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
    }
    

    int estimate_data_size = n / data_n;//評価データの要素数
    int training_data_size = n - estimate_data_size;//学習データの要素数

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
    int block_size = index / data_n;//ブロックに要素がいくつ入っているか
    //double data_block[MATRIX_DIM][MATRIX_DIM][MATRIX_DIM];//data_block[いくつめのブロックか][上からなんこめのデータか][データの中の何次元目の情報か]
    //double label_block[MATRIX_DIM][MATRIX_DIM];//label_block[いくつめのブロックか][上からなんこめのデータのラベルか]
    //segmentation error
  
    //ブロックにデータを格納する
    for(int i = 0; i< index; n++){
        int block_num = i / block_size;//ブロック割つけ
        int number = i % block_size;//ブロック内での番号
        for(int j = 0; j < N; j++){
            data_block[block_num][number][j] = data[i][j];
        }

        //ラベルブロックにわりつけ
        label_block[block_num][number] = label[i];
    }

    for(int estimate_data_num = 0; estimate_data_num < data_n; estimate_data_num++){//交差検証開始
        //estimate_data_num : 評価に使うブロック番号

        double training_data[MATRIX_DIM][MATRIX_DIM];//このサイクルで使う学習データ
        double training_data_label[MATRIX_DIM];//このサイクルで使う学習データのラベル

        //ブロックにあるデータをtraining_dataに移す
        
        int now = 0;//今移しているtraining_dataの番号
        for(int now_block_num = 0; now_block_num < data_n; now_block_num++){
            if(estimate_data_num != now_block_num){//ブロックごとに学習用か評価用か判定し、学習用であればtraining_dataに詰めていく
                for(int i = 0; i < block_size; i++){//ブロック内での番号
                    for(int j = 0; j < N; j++){//データを詰める
                        training_data[now][j] = data_block[now_block_num][i][j];
                    }
                    training_data_label[now] = label_block[now_block_num][i];//ラベルを詰める
                    now++;
                }
            }
        }
        
    

        //行列Gを計算する。その際行列を半正定値行列にするために対角線上成分に微小値を加える。
        //カーネルによって微小値を変化させている。
        double delta = 0;//微小値
        if(kernel == 1){    //ガウスカーネル
            delta = 1.0e-7;
        }else{              //それ以外
            delta = 1.0e-9;
        }

        //training_dataを用いてGを計算
        for(int i = 0; i < training_data_size; i++){
            for(int j = 0; j < training_data_size; j++){
                G[i][j] = training_data_label[i] * training_data_label[j] * kernel_result(training_data[i], training_data[j], kernel, N, sigma);
                if(i == j){
                    G[i][j] += delta;
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
        double weight[N];
        for(int i = 0; i < N; i++){
            for(int j = 0; j < training_data_size; j++){
                weight[i] += alpha[j] * training_data_label[j] * training_data[j][i];
            }
        }
        std::cout << "w: " << weight[0];
        for(int i = 1; i < N; i++){
            std::cout << " " << weight[i];
        }
        std::cout << std::endl;

        //閾値を計算する
        for(int i = 0; i < training_data_size; i++){
            theta += alpha[i] * training_data_label[i] * kernel_result(training_data[i], training_data[maxAindex], kernel, N, sigma);
        }
        theta -= training_data_label[maxAindex];
        std::cout << "θ: " << theta << std::endl;

        //重みと閾値を用いて評価データのデータを評価する。
        int true_positive = 0;
        int true_negative = 0;
        int false_positive = 0;
        int false_negative = 0;
        for(int i = 0; i < estimate_data_size; i++){
            double sum = 0;
            int result;
            for(int k = 0; k < training_data_size; k++){
                sum += alpha[k] * training_data_label[k] * kernel_result(training_data[k], data_block[estimate_data_num][i], kernel, N, sigma);
            }
            sum -= theta;
            if(sum > 0){
                result = 1;
            }else{
                result = -1;
            }

            if(result == 1 && label_block[estimate_data_num][i] == 1){
                true_positive++;
            }else if(result == 1 && label_block[estimate_data_num][i] == -1){
                false_positive++;
            }else if(result == -1 && label_block[estimate_data_num][i] == 1){
                false_negative++;
            }else{
                true_negative++;
            }
        }

        //結果を表示
        std::cout << "estimate_block : " << estimate_data_num << std::endl;
        std::cout << "precision : " << true_positive / (true_positive + false_positive) << std::endl;
        std::cout << "recall : " << true_positive / (true_positive + false_negative) << std::endl;
    }
*/
    return 0;
}

//カーネルを計算する関数
/*
double* x,y :計算する二つの配列
int kernel  :使用するカーネルを表す変数0:内積 1:多項式カーネル 2:ガウスカーネル)
int degree  :データの次元数。大域変数ではNにあたる
double sigma:ガウスカーネル計算時に使用するシグマ定数
*/
double kernel_result(double* x, double* y, int kernel, int degree, double sigma){
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
        result = pow(1 + inner_product, 2);
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