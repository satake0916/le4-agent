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
		CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM];
	int n, m, p;
    int N;//データの次元数
    int data_n = 5;//データ分割数
    int total_data_size = 100;//総データ数。学習用データと評価用データに別れる。
    int estimate_data_size = 0;//評価用データ数
    int training_data_size = 0;//学習用データ数
    int kernel = 0;//kernelの選択。0:内積 1:多項式カーネル 2:ガウスカーネル
    std::string file_name;//データファイルの名前
    double kernel_result(double* x, double* y, int kernel, int degree, double sigma);//カーネルの計算結果を出力する関数

    double sigma = pow(5.0, 0.5);//ガウスカーネルで用いるシグマ定数
    double data[total_data_size][total_data_size];    //データ格納先
    double label[total_data_size];               //ラベル格納先
    int index = 0;//データ数
    double theta = 0;//閾値
    double epsilon = 0.1;
    double C = 1000;
    bool estimate_bool;
    std::string using_labels[MATRIX_DIM];//学習に使用する属性



    {
    //ファイルを指定する
    std::cout << "データファイルを指定してください : ";
    std::cin >> file_name;
    //次元数を指定する
    std::cout << "データの次元数を指定してください : ";
    std::cin >> N;
    //使う属性を選択する
    std::cout << "使用する属性を指定してください : ";
    for(int i = 0; i < N ; i++){
        std::string use_name;//属性名
        std::cin >> use_name;
        using_labels[i] = use_name;
    }
    //Cを指定する
    std::cout << "Cを指定してください : ";
    std::cin >> C;
    //カーネルを指定する
    std::cout << "使用するカーネルを指定してください(0:カーネルトリックなし　1:多項式カーネル　2:ガウスカーネル) : ";
    std::cin >> kernel;
    if(kernel < 0 || kernel > 2){
        std::cout << "正しくカーネルを指定してください(0:カーネルトリックなし　1:多項式カーネル　2:ガウスカーネル)" << std::endl;
        exit(0);
    }
    char estimate_char;
    std::cout << "交差検証を行いますか？(y/n)：";
    std::cin >> estimate_char;
    if(estimate_char == 'y'){
        estimate_bool = true;
    }else if(estimate_char == 'n'){
        estimate_bool = false;
    }else{
        std::cout << "yまたはnを入力してください";
        exit(0);
    }
    }  

    file_name = "SFlisting.csv";
    N = 4;
    using_labels[0] = "accommodates";
    using_labels[1] = "bedrooms";
    using_labels[2] = "beds";
    using_labels[3] = "minimum_nights";

    //csvファイルを読み取る。csv_dataにstring型で格納する。
    {
    std::ifstream ifs(file_name);
    if(ifs.fail()){
        std::cout << "ファイルが見つかりません : " << file_name << std::endl;
        exit(0);
    }

    std::string line;
    std::string csv_data[200][100];
    const std::string delim = ",";
    int row = 0;
    int now = 0;
    int label_num = 0;
    int col;
    while(now <= total_data_size){
        if(getline(ifs, line)){
            col = 0;
            //delimを区切り文字として切り分け、intに変換してdata[][]に格納する
            for(std::string::size_type spos, epos = 0;(spos = line.find_first_of(delim, epos)) != std::string::npos;){
                spos += 1;
                std::string token = line.substr(spos,(epos = line.find_first_of(delim, spos))-spos);
                if(token == ""){
                    token = "0";
                }
            csv_data[row][col++] = token;
            }
        row++;
        }
    now++;
    }

    //csv_dataからdata[][]に移す
    int using_labels_numbers[N];
    for(int i = 0; i < N; i++){
        using_labels_numbers[i] = -1;
    }
    int price_number;
    int s = 0;
    for(int i = 0; i < col; i++){
        for(std::string using_label: using_labels){
            if(using_label == csv_data[0][i]){
                using_labels_numbers[s] = i;
                s++;
                break;
            }
        }
        if("price" == csv_data[0][i]){
            price_number = i;
        }
    }

    //入力で指定された引数が存在するかどうか
    for(int i = 0; i < N; i++){
        if(using_labels_numbers[i] == -1){
            std::cout << "指定された属性が存在しません" << std::endl;
            exit(0);
        }
    }


    double max_numbers[N];
    double min_numbers[N];
    for(int i = 0; i < N; i++){
        max_numbers[i] = 0;
        min_numbers[i] = 10000;
    }
    double max_price = 0;
    double min_price = 10000;
    for(int i = 0; i < total_data_size; i++){
        for(int j = 0; j < N ; j++){
            data[i][j] = std::stod(csv_data[i+1][using_labels_numbers[j]]);
            //std::cout << csv_data[i+1][using_labels_numbers[j]] << " ";
            if(data[i][j] > max_numbers[j]){
                max_numbers[j] = data[i][j];
            }else if(data[i][j] < min_numbers[j]){
                min_numbers[j] = data[i][j];
            }
        }
        label[i] = std::stod(csv_data[i+1][price_number].erase(0,1));
        if(label[i] > max_price){
            max_price = label[i];
        }else if(label[i] < min_price){
            min_price = label[i];
        }
        //std::cout << label[i] << std::endl;
    }

    //正規化
    for(int i = 0; i < total_data_size; i++){
        for(int j = 0; j < N; j++){
            data[i][j] = (data[i][j] - min_numbers[j]) / (max_numbers[j] - min_numbers[j]);
            //std::cout << data[i][j] << " ";
        }
        //label[i] = (label[i] - min_price) / (max_price - min_price);
        //std::cout << label[i] << std::endl;
    }
    }
    
    if(estimate_bool){
        estimate_data_size = total_data_size / data_n;//評価データの要素数
        training_data_size = total_data_size - estimate_data_size;//学習データの要素数

        n = 2 * training_data_size;
        index = training_data_size;
    }else{
        n = 2 * total_data_size;
        index = total_data_size;
    }
    m = 2 * n;
    p = 1;

    //g0の設定
    {
    for(int i = 0; i < index; i++){
        g0[i] = -1 * label[i];
    }
    for(int i = index; i < n; i++){
        g0[i] = epsilon;
    }
    }

    //CEの設定（全部1）
    
    {
    for(int i = 0; i < index; i++){
        for(int j = 0; j < p; j++){
	        CE[i][j] = 1;
        }
        for(int i = index; i < n; i++){
            for(int j = 0; j < p; j++){
                CE[i][j] = 0;
            }
        }
    }
    }

    //ce0の設定（0）
    {
    for (int i = 0; i < p; i++){
      ce0[i] = 0;
    }
    }


    // 0 <= alpha, alpha* <= C
    {
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            CI[i][j] = 0;
        }
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < index; j++){
            if(i % index == j){
                CI[i][j] = 1;
            }
        }
        for(int j = index; j < n; j++){
            if(i % index == (j - index)){
                CI[i][j] = -1;
            }
        }
        for(int j = n; j < n + index; j++){
            if(i == (j - n)){
                CI[i][j] = -1;
            }else if((i - index) == (j - n)){
                CI[i][j] = 1;
            }
        }
        for(int j = (n + index); j < m; j++){
            if(i == (j - n - index)){
                CI[i][j] = 1;
            }else if((i - index) == (j - n - index)){
                CI[i][j] = -1;
            }
        }
	}
    }

    //ci0の設定（前半で以上、後半で以下を指定する)
    {
    for(int i = 0; i < index; i++){
        ci0[i] = 0;
    }
    for(int i = index; i < n; i++){
        ci0[i] = 2 * C;
    }
    for(int i = n; i < (n + index); i++){
        ci0[i] = 0;
    }
    for(int i = (n + index); i < m; i++){
        ci0[i] = 2 * C;
    }
    }
    
    if(estimate_bool){//交差検証をするとき
        //分割数に応じてdataをブロックごとにわけて、どのブロックを学習、どのブロックを評価という風にわけていく
        double result[data_n];//交差検証の評価値
        {
        int block_size = total_data_size / data_n;//ブロックに要素がいくつ入っているか

        //data配列のラベルの次の要素にブロック番号を表す配列を追加する。
        for(int i = 0; i < total_data_size; i++){
            int block_num = i / block_size;//ブロック番号
            data[i][N+1] = block_num;
        }
        }
        

        //交差検証開始
        {
        for(int estimate_data_num = 0; estimate_data_num < data_n; estimate_data_num++){
            std::cout << estimate_data_num << std::endl;
            //estimate_data_num : 評価に使うブロック番号
            double pre_alpha[total_data_size * 2];//[alpha1 - alpha*1, alpha2 - alpha*2, ..., alphan - alpha*n, alpha1 + alpha1, ..., alphan + alpha*n]
            double alpha[total_data_size];//[alpha1, alpha2, ..., alphan]
            double alpha_star[total_data_size];//[alpha*1, alpha*2, ..., alpha*n]

            double training_data[total_data_size][total_data_size];   //このサイクルで使う学習データ
            double training_data_label[total_data_size];         //このサイクルで使う学習データのラベル
            double estimate_data[total_data_size][total_data_size];   //このサイクルで使う評価データ
            double estimate_data_label[total_data_size];         //このサイクルで使う評価データのラベル

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
            {   //カーネルによって微小値を変化させている。
            double delta = 0;//微小値
            if(kernel == 2){    //ガウスカーネル
                delta = 1.0e-7;
            }else{              //それ以外
                delta = 1.0e-9;
            }

            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    G[i][j] = 0;
                    if(i == j){
                        G[i][j] += delta;
                    }
                }
            }

            for(int i = 0; i < index; i++){
                for(int j = 0; j < index; j++){
                    G[i][j] += kernel_result(training_data[i], training_data[j], kernel, N, sigma);
                }
            }
            }


            //solve_quadprogの使用
            solve_quadprog(G, g0, n, CE, ce0, p, CI, ci0, m, pre_alpha);

            //[alpha1 - alpha*1, alpha2 - alpha*2, ..., alphan - alpha*n, alpha1 + alpha1, ..., alphan + alpha*n]を出す
            // solve_quadprogの結果出て来たalphaを出力し、その最大値をmaxAに、インデックスをmaxAindexに格納する。

            for(int i = 0; i < index; i++){
                alpha[i] = (pre_alpha[i] + pre_alpha[index + i]) / 2;
                alpha_star[i] = (-1 * pre_alpha[i] + pre_alpha[index + i]) / 2;
            }

            double maxA = alpha[0];//alphaの最大値
            int maxAindex = 0;//alphaの最大値のインデックス
            for(int i = 0; i < training_data_size; i++){
                if(alpha[i] > maxA){
                    maxA = alpha[i];
                    maxAindex = i;
                }
            }

            
            // 重みを計算する
            double weight[100] = {};
            for(int i = 0; i < N; i++){
                for(int j = 0; j < training_data_size; j++){
                    weight[i] += (alpha[j] - alpha_star[j]) * training_data[j][i];
                }
                std::cout << "weight[" << i << "] : " << weight[i] << std::endl;
            }

            //閾値を計算する
            theta = 0;
            for(int i = 0; i < training_data_size; i++){
                theta += (alpha[i] - alpha_star[i]) * kernel_result(training_data[i], training_data[maxAindex], kernel, N, sigma);
            }
            theta += (epsilon - training_data_label[maxAindex]);
            std::cout << "theta : " << theta << std::endl;

            //重みと閾値を用いて評価データのデータを評価する。
            //評価データをひとつずつ評価していく。
            double f[estimate_data_num];//回帰式による評価値
            for(int i = 0; i < estimate_data_size; i++){
                f[i] = kernel_result(weight, estimate_data[i], kernel, N, sigma) - theta;
            }

            //結果を平均二乗誤差で評価
            double result_sum = 0;
            for(int i = 0; i < estimate_data_size; i++){
                result_sum += pow((estimate_data_label[i] - f[i]), 2);
            }
            result[estimate_data_num] = result_sum / estimate_data_size;
            std::cout << "平均二乗誤差 : " << result[estimate_data_num] << std::endl;
            std::cout << "............................" << std::endl;
            
            }
        }//交差検証おわり

        double sum = 0;
        for(int i = 0; i < data_n; i++){
            sum += result[i];
        }
        std::cout << "平均 : " << (sum / data_n) << std::endl;


    }else{
        double pre_alpha[total_data_size * 2];
        double alpha[total_data_size];
        double alpha_star[total_data_size];
        //行列Gを計算する。その際行列を半正定値行列にするために対角線上成分に微小値を加える。
        {   //カーネルによって微小値を変化させている。
        double delta = 0;//微小値
        if(kernel == 2){    //ガウスカーネル
            delta = 1.0e-7;
        }else{              //それ以外
            delta = 1.0e-9;
        }

        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                G[i][j] = 0;
                if(i == j){
                    G[i][j] += delta;
                }
            }
        }

        for(int i = 0; i < index; i++){
            for(int j = 0; j < index; j++){
                G[i][j] += kernel_result(data[i], data[j], kernel, N, sigma);
            }
        }
        }

        //solve_quadprogの使用
        solve_quadprog(G, g0, n, CE, ce0, p, CI, ci0, m, pre_alpha);

        //[alpha1 - alpha*1, alpha2 - alpha*2, ..., alphan - alpha*n, alpha1 + alpha1, ..., alphan + alpha*n]を出す
        // solve_quadprogの結果出て来たalphaを出力し、その最大値をmaxAに、インデックスをmaxAindexに格納する。

        for(int i = 0; i < index; i++){
            alpha[i] = (pre_alpha[i] + pre_alpha[index + i]) / 2;
            alpha_star[i] = (-1 * pre_alpha[i] + pre_alpha[index + i]) / 2;
        }

        double maxA = alpha[0];//alphaの最大値
        int maxAindex = 0;//alphaの最大値のインデックス
        for(int i = 0; i < total_data_size; i++){
            if(alpha[i] > maxA){
                maxA = alpha[i];
                maxAindex = i;
            }
        }

            
        // 重みを計算する
        double weight[100] = {}; 
        for(int i = 0; i < N; i++){
            for(int j = 0; j < total_data_size; j++){
                weight[i] += (alpha[j] - alpha_star[j]) * data[j][i];
            }
            std::cout << "weight[" << i << "] : " << weight[i] << std::endl;
        }

        //閾値を計算する
        for(int i = 0; i < total_data_size; i++){
            theta += (alpha[i] - alpha_star[i]) * kernel_result(data[i], data[maxAindex], kernel, N, sigma);
        }
        theta += (epsilon - label[maxAindex]);
        std::cout << "theta : " << theta << std::endl;

        // gnuplotを使ってプロットする
        if(N == 2){
            FILE *gp;
            gp = popen("gnuplot -persist", "w");
            fprintf(gp, "set multiplot\n");
            fprintf(gp, "set xrange [0:1]\n");
            fprintf(gp, "set yrange [0:1]\n");
            fprintf(gp, "set zrange [0:2]\n");
            fprintf(gp, "set ticslevel 0\n");

            //関数をプロット
            if(kernel == 0){
                //カーネルトリックなし
                fprintf(gp, "splot (%f * x + %f * y) - %f\n", weight[0], weight[1], theta);
            }else if(kernel == 1){
                //多項式カーネル
                fprintf(gp, "splot (1 + %f * x + %f * y) ** 2 - %f\n", weight[0], weight[1], theta);
            }else{
                //ガウスカーネル
                fprintf(gp, "splot exp(0 - ((x - %f) ** 2 + (y - %f) ** 2) / (2 * (%f ** 2))) - %f\n", weight[0], weight[1], sigma, theta);
            }

            //データ点をプロット
            fprintf(gp, "splot '-' \n");
            for (int i = 0; i < total_data_size; i++) {
                fprintf(gp, "%f\t%f\t%f\n", data[i][0], data[i][1], label[i]);
            }
       
            fprintf(gp, "unset multiplot\n");
            pclose(gp);      
        }
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
        result = pow((1 + inner_product), 2);
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