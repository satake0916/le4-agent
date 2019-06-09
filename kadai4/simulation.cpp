#include <iostream>
#include <sstream>
#include <string>
#include "QuadProg++.hh"
#include <fstream>
#include <cmath>
#include <chrono>
#include <time.h>
#include <float.h>

/*
価格推薦戦略プログラム

課題 4-1
SVR での予測が 100% 正しいとすると予測よりわずかに小さい価格を設定すれば収入を最大化することが できます。
しかし実際に SVR で 100% 正しい予測をすることはできません。
そのため収入を最大化しようと するあまり予測との差を小さくしすぎると SVR の誤差で競合施設の価格を大きく予測してしまい、
結果自己 の施設の価格の方が高くなってしまう可能性があります。
したがってまずは作成した SVR の誤差を確認し、 予測よりどれほど価格を設定すれば安全であるかを考えます。
相関係数を調べた結果から学習に使用する属 性を accomodates, bedrooms, beds, minimum nights の 4 属性で予測を行なったところ、
平均二乗誤差は 18486.5 となりました (カーネルトリックなし、C=1000)。すなわち実際の誤差は約 135 と考えることができ ます。
以上より SVR の誤差を帳消しするために予測値から 135 を引いた後、それより少し小さな値を取るた めに 0.9 倍した価格に設定する戦略を取ります。

課題 4-2
他の価格推薦エージェントも同様に予測値の少し下という価格推薦を行うと、価格競争が発生して双方の利 潤が失われてしまいます。
今回は利潤の合計を大きくすることを目標とし、価格推薦エージェントが自己の他 にもう 1 つ存在すると考えるとそれぞれの利潤がちょうど半分になることが理想です。
そこで毎回市場に参加 するのではなく、確率的に市場に参加することを考えます。こうすることで費用を肥大化させることなく利 潤合計を維持することが可能になります。
条件より 1/2 の確率で市場に参加するエージェントを作成します。 またとにかく入札することを防ぐために市場に参加するには市場価格の半分の費用がかかると仮定し、
それぞ れのエージェントにかかる費用も計算します。
*/

int main (int argc, char *const argv[]) {
	double G[MATRIX_DIM][MATRIX_DIM], g0[MATRIX_DIM],
		CE[MATRIX_DIM][MATRIX_DIM], ce0[MATRIX_DIM], 
		CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM];
	int n, m, p;
    int N;                      //データの次元数
    int data_n = 5;             //データ分割数。前半を訓練、後半を評価に用いる。
    int total_data_size = 100;  //総データ数。学習用データと評価用データに別れる。
    int estimate_data_size = 0; //評価用データ数
    int training_data_size = 0; //学習用データ数
    int kernel = 0;             //kernelの選択。0:内積 1:多項式カーネル 2:ガウスカーネル
    std::string file_name;      //データファイルの名前
    double kernel_result(double* x, double* y, int kernel, int degree, double sigma);//カーネルの計算結果を出力する関数

    double sigma = pow(5.0, 0.5);                       //ガウスカーネルで用いるシグマ定数
    double data[total_data_size][total_data_size];      //データ格納先
    double label[total_data_size];                      //ラベル格納先
    int index = 0;              //データ数
    double theta = 0;           //閾値
    double epsilon = 0.1;       //誤差許容範囲。本プログラムでは0.1に固定
    double C = 1000;            //判定ミスコスト。本プログラムでは1000に固定
    std::string using_labels[MATRIX_DIM];               //学習に使用する属性
    int simulation = 0;         //指定されたシミュレーション番号。1なら課題4-1、課題4-1用のシミュレートを行う。

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
        std::string use_name;   //属性名
        std::cin >> use_name;
        using_labels[i] = use_name;
    }
    //カーネルを指定する
    std::cout << "使用するカーネルを指定してください(0:カーネルトリックなし　1:多項式カーネル　2:ガウスカーネル) : ";
    std::cin >> kernel;
    if(kernel < 0 || kernel > 2){
        std::cout << "正しくカーネルを指定してください(0:カーネルトリックなし　1:多項式カーネル　2:ガウスカーネル)" << std::endl;
        exit(0);
    }
    //シミュレーション方法を指定する（1:課題4-1; 2:課題4-2）
    std::cout << "シミュレーション方法を指定してください（1:課題4-1; 2:課題4-2）：";
    std::cin >> simulation;
    } 
    if(simulation < 1 || simulation > 2){
        std::cout << "正しくシミュレーション方法を指定してください（1:課題4-1; 2:課題4-2）" << std::endl;
        exit(0);
    }

    
    //csvファイルを読み取る。csv_dataにstring型で格納する。
    {
    //ファイルをifstreamを用いて読み取る
    std::ifstream ifs(file_name);
    if(ifs.fail()){
        std::cout << "ファイルが見つかりません : " << file_name << std::endl;
        exit(0);
    }

    std::string line;               //1行を保存するローカル変数
    std::string csv_data[total_data_size * 2][total_data_size]; //csvデータをそのまま保存するローカル配列
    const std::string delim = ","; 
    int row = 0;        //列数を保持する
    int now = 0;        //データ数がtotal_data_sizeを超えないようにする
    int col;            //行数を保持する
    while(now <= total_data_size){
        if(getline(ifs, line)){
            col = 0;    //何行目を扱っているかを管理する
            //delimを区切り文字として切り分け、intに変換してdata[][]に格納する
            for(std::string::size_type spos, epos = 0;(spos = line.find_first_of(delim, epos)) != std::string::npos;){
                spos += 1;
                std::string token = line.substr(spos,(epos = line.find_first_of(delim, spos))-spos);
                if(token == ""){     //空白セルは0として扱う
                    token = "0";
                }
            csv_data[row][col++] = token;   //csv_dataに格納する
            }
        row++;
        }
    now++;
    }

    //属性の列番号を調べる。using_labels_numbersとprice_numbersは属性の列番号を表す。
    int using_labels_numbers[N];
    for(int i = 0; i < N; i++){
        using_labels_numbers[i] = -1;
    }
    int price_number = -1;
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
        if(price_number == -1){
            std::cout << "価格が存在しません" << std::endl;
            exit(0);
        }
    }


    double max_numbers[N];  //各属性の最大値を保持するローカル配列。
    double min_numbers[N];  //各属性の最小値を保持するローカル配列。
    for(int i = 0; i < N; i++){
        max_numbers[i] = 0;
        min_numbers[i] = DBL_MAX;
    }
    //csv_dataからdata[][]とlabel[]に移す。
    for(int i = 0; i < total_data_size; i++){
        for(int j = 0; j < N ; j++){
            data[i][j] = std::stod(csv_data[i+1][using_labels_numbers[j]]);
            if(data[i][j] > max_numbers[j]){
                max_numbers[j] = data[i][j];
            }else if(data[i][j] < min_numbers[j]){
                min_numbers[j] = data[i][j];
            }
        }
        label[i] = std::stod(csv_data[i+1][price_number].erase(0,1));   //文字列の先頭に$マークが付いているので外す。
    }

    //各属性の値を正規化する。
    for(int i = 0; i < total_data_size; i++){
        for(int j = 0; j < N; j++){
            data[i][j] = (data[i][j] - min_numbers[j]) / (max_numbers[j] - min_numbers[j]);
        }
    }

    //相関係数を計算してみる
    /*
    double price_average = 0;
    for(int j = 0; j < total_data_size; j++){
        price_average += label[j];
    }
    price_average /= total_data_size;
    for(int i = 0; i < N; i++){
        double sxy = 0;
        double sx = 0;
        double sy = 0;
        double average = 0;
        for(int j = 0; j < total_data_size; j++){
            average += data[j][i];
        }
        average /= total_data_size;
        for(int j = 0; j < total_data_size; j++){
            sxy += ( data[j][i] - average ) * ( label[j] - price_average );
            sx += ( data[j][i] - average ) * ( data[j][i] - average );
            sy += ( label[j] - price_average ) * ( label[j] - price_average );
        }
        sxy /= total_data_size;
        sx = pow((sx / total_data_size), 0.5);
        sy = pow((sy / total_data_size), 0.5);
        double result = sxy / ( sx * sy );
        std::cout << result << std::endl;
    }
    */
    }
    
    estimate_data_size = total_data_size / data_n;//評価データの要素数
    training_data_size = total_data_size - estimate_data_size;//学習データの要素数
    n = 2 * training_data_size;
    index = training_data_size;
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
        

    //データ前半8割を訓練用としてSVRで価格予測器を作成し、シミュレートする
        {
            //estimate_data_num : 評価に使うブロック番号
            double pre_alpha[total_data_size * 2];//[alpha1 - alpha*1, alpha2 - alpha*2, ..., alphan - alpha*n, alpha1 + alpha1, ..., alphan + alpha*n]
            double alpha[total_data_size];//[alpha1, alpha2, ..., alphan]
            double alpha_star[total_data_size];//[alpha*1, alpha*2, ..., alpha*n]

            double training_data[total_data_size][total_data_size];     //学習データ
            double training_data_label[total_data_size];                //学習データのラベル
            double estimate_data[total_data_size][total_data_size];     //評価データ
            double estimate_data_label[total_data_size];                //評価データのラベル

            //データを実際に使用するtraining_dataに移す
            int training_now = 0;   //今移しているtraining_dataの番号
            int estimate_now = 0;   //今写しているestimate_dataの番号
            for(int i = 0; i < total_data_size; i++){
                if(i < training_data_size){//データ前半をtraining_dataに格納する
                    for(int j = 0; j < N; j++){
                        training_data[training_now][j] = data[i][j];
                    }
                    training_data_label[training_now] = label[i];
                    training_now++;
                }else{//データ後半をestimate_dataに格納する
                    for(int j = 0; j < N; j++){
                        estimate_data[estimate_now][j] = data[i][j];
                    }
                    estimate_data_label[estimate_now] = label[i];
                    estimate_now++;
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
            double weight[total_data_size];
            for(int i = 0; i < total_data_size; i++){
                weight[i] = 0;
            }
            for(int i = 0; i < N; i++){
                for(int j = 0; j < training_data_size; j++){
                    weight[i] += (alpha[j] - alpha_star[j]) * training_data[j][i];
                }
            }

            //閾値を計算する
            theta = 0;
            for(int i = 0; i < training_data_size; i++){
                theta += (alpha[i] - alpha_star[i]) * kernel_result(training_data[i], training_data[maxAindex], kernel, N, sigma);
            }
            theta += (epsilon - training_data_label[maxAindex]);

            double f[estimate_data_size];//回帰式による評価値
            //SVRを用いた価格設定を冒頭の戦略通りに行う。
            for(int i = 0; i < estimate_data_size; i++){
                    f[i] = kernel_result(weight, estimate_data[i], kernel, N, sigma) - theta;
                    if(f[i] > 135){
                        f[i] -= 135;
                    }
                    f[i] *= 0.9;
            }

            if(simulation == 1){    //課題4-1
                //重みと閾値を用いて評価データの価格を予測し、その0.9倍の価格を市場に出した時の純利益を計算する
                //同時にすべての物件の平均価格で市場に出す単純な価格設定戦略での純利益も計算し、SVRを用いた価格推薦戦略の性能を評価する。
                double simple_f = 0;        //単純
                double svr_revenue = 0;     //SVRを用いた価格推薦戦略の純利益
                double simple_revenue = 0;  //単純な価格設定戦略の純利益
                double ideal_revenue = 0;   //評価データの実際の市場価格の合計
                for(int i = 0; i < training_data_size; i++){
                    simple_f += training_data_label[i];
                }
                simple_f /= training_data_size;
                //実際の市場の価格と比べてみる
                for(int i = 0; i < estimate_data_size; i++){
                    ideal_revenue += estimate_data_label[i];
                    std::cout << "f[" << i << "]:" << f[i] << " label[" << i << "]:" << estimate_data_label[i] << std::endl;
                    if(f[i] < estimate_data_label[i]){
                        svr_revenue += f[i];
                    }
                    if(simple_f < estimate_data_label[i]){
                        simple_revenue += simple_f;
                    }
                }
                //結果を表示
                std::cout << "理論上の最大収益(実際の市場価格の和) : ¥" << ideal_revenue << std::endl;
                std::cout << "SVRを用いて価格設定を行なったときの収益 : ¥" << svr_revenue << std::endl;
                std::cout << "単純な価格設定を行なった時の収益 : ¥" << simple_revenue << std::endl;
            }else if(simulation == 2){//課題4-2
                //民泊市場にはリスティングデータに従う貸し手と，価格推薦エージェント２体の計３つが存在するという条件のもとで，考案した入札戦略の性能を評価する．
                //エージェントAとBがそれぞれ1/2の確率で市場に参加するとする。
                //同時に同じ価格で出店した場合はランダムに選ばれるとし、出店する際には実際の市場の価格の0.5倍のコストがかかるとする。
                double agentA_revenue = 0;  //エージェントAの収入
                double agentB_revenue = 0;  //エージェントBの収入
                double agentA_cost = 0;     //エージェントAの費用
                double agentB_cost = 0;     //エージェントBの費用
                double ideal_revenue = 0;   //評価データの実際の市場価格の合計
                srand((unsigned int)time(NULL));    //乱数の初期値を現在時刻で設定する
                for(int i = 0; i < estimate_data_size; i++){
                    //エージェントA,Bが出店するかどうかをランダムで決定する
                    bool agentA_auction = true;
                    bool agentB_auction = true;
                    if(rand() % 2 == 0){    // 1/2の確率で出店する
                        agentA_auction = false;
                    }
                    if(rand() % 2 == 0){    // 1/2の確率で出店する
                        agentB_auction = false;
                    }
                    //費用の計算。市場価格の半額が必要費用と仮定する。
                    if(agentA_auction){
                        agentA_cost += (estimate_data_label[i] / 2);
                    }
                    if(agentB_auction){
                        agentB_cost += (estimate_data_label[i] / 2);
                    }
                    //収入を計算
                    if(f[i] < estimate_data_label[i]){
                        ideal_revenue += estimate_data_label[i];
                        if(agentA_auction && agentB_auction){
                            if(rand() % 2 == 0){
                                agentA_revenue += f[i];
                            }else{
                                agentB_revenue += f[i];
                            }
                        }else if(agentA_auction == true && agentB_auction == false){
                            agentA_revenue += f[i];
                        }else if(agentA_auction == false && agentB_auction == true){
                            agentB_revenue += f[i];
                        }else{
                            //何もなし
                        }
                    }
                }
                //結果を表示
                std::cout << "理論上の最大収益(実際の市場価格の和) : ¥" << ideal_revenue << std::endl;
                std::cout << "エージェントAの収益 : ¥" << agentA_revenue << std::endl;
                std::cout << "エージェントBの収益 : ¥" << agentB_revenue << std::endl;
                std::cout << "エージェントAの費用 : ¥" << agentA_cost << std::endl;
                std::cout << "エージェントBの費用 : ¥" << agentB_cost << std::endl;
            }else{
                exit(1);
            }
        }//シミュレーションおわり
    
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