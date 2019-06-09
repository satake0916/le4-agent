#include <iostream>
#include <sstream>
#include <string>
#include "QuadProg++.hh"
#include <fstream>
#include <cmath>
#include <chrono>

int main (int argc, char *const argv[]) {
	double G[MATRIX_DIM][MATRIX_DIM], g0[MATRIX_DIM], 
		CE[MATRIX_DIM][MATRIX_DIM], ce0[MATRIX_DIM], 
		CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM], 
		alpha[MATRIX_DIM];
	int n, m, p;
    int N;//データの次元数
    int kernel = 0;//kernelの選択。0:内積 1:多項式カーネル 2:ガウスカーネル
    std::string file_name;//データファイルの名前
    double kernel_result(double* x, double* y, int kernel, int degree, double sigma);//カーネルの計算結果を出力する関数

    double sigma =  10;//ガウスカーネルで用いるシグマ定数
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


    //ファイルを読み取る。
    {
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
    }

    n = index; //データの数
    m = n;
    p = 1;

    //行列Gを計算する。その際行列を半正定値行列にするために対角線上成分に微小値を加える。
    {   //カーネルによって微小値を変化させている。
    double delta = 0;//微小値
    if(kernel == 1){    //ガウスカーネル
        delta = 1.0e-7;
    }else{              //それ以外
        delta = 1.0e-9;
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            G[i][j] = label[i] * label[j] * kernel_result(data[i], data[j], kernel, N, sigma);
            if(i == j){
                G[i][j] += delta;
            }
        }
    }
    }

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

    /*時間の計測
    auto start = std::chrono::system_clock::now();
    */

    //solve_quadprogの使用
    solve_quadprog(G, g0, n, CE, ce0, p, CI, ci0, m, alpha);

    /*時間の計測
    auto end = std::chrono::system_clock::now();       
    auto dur = end - start;        
    auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << msec << " ms" << std::endl;
    */

    // solve_quadprogの結果出て来たalphaを出力し、その最大値をmaxAに、インデックスをmaxAindexに格納する。
    {
    double maxA = alpha[0];//alphaの最大値
    int maxAindex = 0;//alphaの最大値のインデックス
    for(int i = 0; i < n; i++){
        std::cout << "alpha[" << i << "] : " << alpha[i] << std::endl;
        if(alpha[i] > maxA){
            maxA = alpha[i];
            maxAindex = i;
        }
    }

    // 重みを計算する
    double weight[N];
    for(int i = 0; i < N; i++){
        for(int j = 0; j < index; j++){
            weight[i] += alpha[j] * label[j] * data[j][i];
        }
    }
    std::cout << "w: " << weight[0];
    for(int i = 1; i < N; i++){
        std::cout << " " << weight[i];
    }
    std::cout << std::endl;

    //閾値を計算する
    
    for(int i = 0; i < n; i++){
        theta += alpha[i] * label[i] * kernel_result(data[i], data[maxAindex], kernel, N, sigma);
    }
    theta -= label[maxAindex];
    std::cout << "θ: " << theta << std::endl;
    }

    //データの次元数が2のときのみgnuplotを用いて可視化する
    if(N == 2){
        FILE *gp;
        gp = popen("gnuplot -persist","w");
        fprintf(gp, "set multiplot\n");
        fprintf(gp, "unset key\n");
        fprintf(gp, "set xrange [%f:%f]\n", minX, maxX);
        fprintf(gp, "set yrange [%f:%f]\n", minY, maxY);
        fprintf(gp, "set lmargin at screen 0.1\n");
        fprintf(gp, "set rmargin at screen 0.9\n");
        fprintf(gp, "set bmargin at screen 0.1\n");
        fprintf(gp, "set tmargin at screen 0.9\n");
        fprintf(gp, "set isosample 200, 200\n");
        fprintf(gp, "set contour base\n");
        fprintf(gp, "set cntrparam levels discrete 0.0\n");
        fprintf(gp, "set nosurface\n");
        fprintf(gp, "set zeroaxis ls -1\n");
        fprintf(gp, "set view 0, 0\n");

        //カーネルごとに線かマージンをプロットする
        if(kernel == 0){
            //カーネルトリックなし
            fprintf(gp, "f(x, y) = 0");
            for(int i = 0; i < n; i++){
                fprintf(gp, " + %f * %f * (%f * x + %f * y)", alpha[i], label[i], data[i][0], data[i][1]);
            }
            fprintf(gp, " - %f\n", theta);
            fprintf(gp, "splot f(x, y) \n");
        }else if(kernel == 1){
            //多項式カーネル
            fprintf(gp, "f(x, y) = 0");
            for(int i = 0; i < n; i++){
                fprintf(gp, " + %f * %f * (1 + %f * x + %f * y) ** 2", alpha[i], label[i], data[i][0], data[i][1]);
            }
            fprintf(gp, " - %f\n", theta);
            fprintf(gp, "splot f(x, y)\n");
        }else{
            //ガウスカーネル
            //マージンをプロットする。データの範囲を100*100に区分けし、そのうちのマージンの部分を黒く表示する
            double dx = (maxX - minX) / 100;
            double dy = (maxY - minY) / 100;
            double epsilon = 0.4;//マージンとする範囲
            double nowX, nowY;
            for(int i = 0; i < 100; i++){
                for(int j = 0; j < 100; j++){
                    nowX = minX + dx * i;
                    nowY = minY + dy * j;
                    double nows[] = {nowX, nowY};
                    double sign = 0;
                    for(int k = 0; k < n; k++){
                        sign += alpha[k] * label[k] * kernel_result(data[k], nows, kernel, N, sigma);
                    }
                    sign -= theta;
                    if(sign > - epsilon && sign < epsilon){
                        fprintf(gp, "set object %d rect from %f, %f to %f, %f back linewidth 0 fillcolor rgb 'black'\n", i * 100 + j, nowX, nowY, nowX + dx, nowY + dy);
                    }
                }
            }
        }

    //データのプロット
    fprintf(gp, "plot '%s' using 1:($3==1?$2:1/0) w p pt 7 lc rgb 'blue' ps 2 title 'true';",file_name.c_str());
    fprintf(gp, "replot '%s' using 1:($3==-1?$2:1/0) w p pt 2 lc rgb 'red' ps 2 title 'false';",file_name.c_str());

    fprintf(gp, "unset multiplot\n");
    pclose(gp);
    }
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