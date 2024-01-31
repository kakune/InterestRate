# InterestRate
## 概要
何も考えずInterestRateというリポジトリ名にしたが、金利に限らず金融工学に関する様々な数値計算ライブラリの実装を目指す。


## 実行環境
- Intel C++ Compiler
    - 他のCompilerでも良いが、cpp_calculator/CMakeLists.txtを書き換える必要がある。
- Python 3.11.5  
    - `python_wrapper/.venv`にvenv環境を作成することを推奨。
    - 必要なライブラリは`python_wrapper/requirements.txt`。
- Visual Studio Code
    - このリポジトリをcloneして、ワークスペースに読み込めば、Linux環境であればそのまま動作すると思われる。


## 主要な構成
InterestRate  
├── cpp_calculator `C++による計算エンジン`  
├── python_wrapper `pythonによるC++実行およびグラフ描画`  
├── output `出力ファイル`  
├── parameters `パラメータiniファイル`  
└── scripts `デバッグのためのスクリプトなど`  


## 実行用コードの説明
実行用コードは、cppファイルと、それに対応する同名のpython wrapperファイルによって成る。

### SABR_implied_volatility

$$
\mathrm{d} F_{t}(T) = \sigma_{t}(F_{t}(T))^{\beta} \mathrm{d} W_{t}^{1} \\
\mathrm{d} \sigma_{t} = \nu \sigma_{t} \mathrm{d} W_{t}^{2} \\
\mathrm{d} W_{t}^{1}\mathrm{d} W_{t}^{2} = \rho \mathrm{d}t
$$
に従うフォワード価格$F_{t}(T)$のimplied volatilityを計算する。
ただし、$\log F_{t}(T)$および$\log \sigma_{t}$についてのSDEにより実装されている。  
パラメータは以下の通り。

|  パラメータ名  |  意味  |
| ---- | ---- |
|  NPath  |  サンプルパスの数  |
|  NTerms  |  時間の分割数  |
|  TimeMaturity  |  満期までの時間  |
|  MinStrike  |  ストライクの最小値  |
|  MaxStrike  |  ストライクの最大値  |
|  DStrike  |  ストライクの刻み幅  |
|  InitPrice  |  最初のフォワード価格 $f$  |
|  InitVol  |  最初のVolatility $\alpha$  |
|  Volvol  |  $\nu$  |
|  Exponent  |  $\beta$  |
|  Corr  |  $\rho$  |

- C++ファイル
    - 入力
        - parameterのiniファイルの絶対パス
        - parameterのsection名
        - outputのcsvの絶対パス
    - 出力
        - 入力で与えられたパスにcsvファイルを生成する。
            - 一列目: strike
            - 二列目: implied volatility

- Pythonファイル
    - C++をリリースモードでbuildする。
    - gPathParam, gNameSection, gPathOutput を入力に渡して実行する。
    - 出力されたcsvファイルのimplied volatilityに以下の近似解を重ねて描画し、gPathGraphに出力する。

$$
    \sigma_{t}^{\mathrm{IB}} \simeq \frac{\alpha}{(fK)^{(1-\beta)/2}} \frac{1 + \left( \frac{(1-\beta)^{2}\alpha^{2}}{24(1-\beta)^{2}} + \frac{\rho\beta\nu\alpha}{4(fK)^{(1-\beta)/2}} + \frac{(2-3\rho^{2})\nu^{2}}{24}\right)}{1 + \frac{(1-\beta)^{2}}{24}\left(\log \frac{f}{K}\right)^{2} + \frac{(1-\beta)^{4}}{1920}\left(\log \frac{f}{K}\right)^{4}} \frac{\zeta}{\chi(\zeta)} \\
    \zeta = \frac{\nu}{\alpha} (fK)^{(1-\beta)/2} \log \frac{f}{K} \\
    \chi(\zeta) = \log \left( \frac{\sqrt{1-2\rho \zeta + \zeta^{2}} + \zeta - \rho}{1 - \rho}\right)
$$



## 主要なライブラリの説明  

大体は説明をヘッダファイルに書いているので、それを読めば良い。
ヘッダファイルのコメントはDoxygenの文法に従う。

### C++エンジン  

#### 全般的なこと  
- ファイル名はsnake_case。
- `lib/foo`というディレクトリは、`LibFoo`という名前でCMakeのLibraryになっている。
- いずれかのCMakeのライブラリをリンクすれば、コンパイル時に`lib`フォルダがインクルードパスに入る。
- namespace名は、ディレクトリ名およびヘッダファイル名をUpperCamelCaseにして定義される。
- class名はUpperCamelCase。
- 変数名はUpperCamelCaseの語頭に次のような接頭辞を上から順につけて定義される。
    - ローカル変数の場合`l`、グローバル変数の場合`g`、引数の場合`in`、メンバ変数の場合`m`、ループにおける一時変数の場合`i`
    - 生ポインタの場合`p`、shared_ptrの場合`s`、unique_ptrの場合`u`
- template変数名はUpperCamelCalseの接尾辞に`_`をつけて定義される。


#### analytical  
- Black_Scholes
    - BS modelでのヨーロピアンオプション価格。

#### math  
- special_functions
    - 様々な特殊関数。
- findroot_1d
    - 1変数関数の根を求める。
- interpolate_1d
    - 独立変数と従属変数の組を与えると、その内挿関数を構成する。
    - NewtonSpline
        - Newton補間により、任意の次数の多項式により内挿関数を作る。

#### process  
- random
    - ブラウン運動に従うpathを生成する。
- asset
    - SDEに従う原資産モデル。
- short_rate
    - 金利のshort rateモデル。

#### utils  
- parameters
    - パラメータのiniファイルを読み込み、mapで管理する。


### Python wrapper

#### 全般的なこと
- 基本的な命名規則はエンジンと同じ。
- 簡単なものしか作成しないので、コメントをつけていない。

#### cpp
- cmake
    - cmakeによるbuildとそれによって生成されたファイルを実行する関数。

#### finance
- SABR
    - SABRの近似解。

#### plot
- graph
    - グラフ描画。

#### utils
- parameters
    - iniファイルをconfigparserにより読み込む。



