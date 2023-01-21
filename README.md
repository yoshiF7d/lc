# lc.py

https://user-images.githubusercontent.com/31358974/213721453-248d3e7c-7fbf-4894-94cf-74ec46ad2ee7.mp4

## 概要
HPLCの解析プログラムです。このプログラムは、ピークの高さや面積ではなく、形から濃度を推定します。このため、濃度によるピーク位置や形状の変化に対して耐性があります。

## 入力ファイル
ShimazuのLabSolutionsのASCII出力ファイルを扱います。違うデータ形式を扱う場合は、LCData.pyを書き換える必要があります。

## インストール
### pythonのインストール
powershellを立ち上げ、pythonと打つ。

### python module のインストール
powershellで以下のコマンドを入力する。
v
### ソースファイルの配置
ソースファイルを powershellを始めに立ち上げて、pwdと入力したときに表示されるフォルダに置きます。PATHの通し方などわかる人はその必要はありません。

## プログラムの使用
解析の流れは以下のようになります。
まず標準試料のデータ`STD.csv`から、ピーク形状と濃度の関係モデルを作成し、パラメータファイル`Params.csv`を出力します。そのあと、このモデルを用いて未知試料のデータに対してフィッティングをかけ、化合物の濃度をもとめます。

ファイル名の入力はドラッグアンドドロップを使ってもできます。
<details><summary>ドラッグアンドドロップによるファイル名の入力方法</summary>

https://user-images.githubusercontent.com/31358974/213726092-c3425de8-421a-4ee7-993b-e3e9744586ea.mp4
  
</details>


### 標準サンプルの情報が入ったcsvファイルを作る
以下の例のようなCSVファイルを作ります。グラフの左から順番にあらわれるピークに対応する化合物の情報を表の上から順番に書きます。`STD1.txt` などの標準試料のデータは `STD.csv` と同じフォルダに入れていおいてください。

STD.csv:
|   | STD1.txt | STD2.txt | STD3.txt | STD4.txt | STD5.txt |
|:-:|:-:|:-:|:-:|:-:|:-:|
|Maleic Acid|50|25|12.5|6.25|3.125|
|Glucose|100|50|25|12.5|6.25|
|Xylose|50|25|12.5|6.25|3.125|
|Glycerol|50|25|12.5|6.25|3.125|
|Acetic Acid|50|25|12.5|6.25|3.125|
|Ethanol|30|15|7.5|3.75|1.875|
|n-Buthanol|30|15|7.5|3.75|1.875|

またこのファイルにピーク位置の情報を含めることができます。この場合、入力されたピーク位置に近いピークを検出します。

STD.csv (ピーク位置の情報を含める場合):
|   | STD1.txt | STD2.txt | STD3.txt | STD4.txt | STD5.txt |Peak position (.txt で終わらない文字列。空欄でもよい)| 
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|Maleic Acid|50|25|12.5|6.25|3.125|9.3|
|Glucose|100|50|25|12.5|6.25|10|
|Xylose|50|25|12.5|6.25|3.125|10.6|
|Glycerol|50|25|12.5|6.25|3.125|14|
|Acetic Acid|50|25|12.5|6.25|3.125|15.05|
|Ethanol|30|15|7.5|3.75|1.875|21.63|
|n-Buthanol|30|15|7.5|3.75|1.875|33.8|

### パラメータファイルの作成
パラメータファイルの名前を`Params.csv`とします。以下のコマンドを実行します。`--checkParams` オプションをつけると、フィッティングの結果が表示されます。

```
python lc.py -s 'STD.csv' -p 'Params.csv' --checkParams
```

### シングルデータの分析
```
python lc.py -p 'Params.csv' 'data.txt'
```
### 複数データの分析
`data` が複数のデータファイルが入っているフォルダとします。
```
python lc.py -p 'Params.csv' 'data'
```
#### フィッティングの結果がみたい場合
`--plotDir` オプションに出力フォルダ名を指定してください。
```
python lc.py -p 'Params.csv' 'data' --plotDir 'plot'
```
#### 結果をファイルに出力したい場合
```
python lc.py -p 'Params.csv' 'data' -o 'out.csv'
```

## その他のオプション

### 使用するディテクターを指定する
デフォルトでは`LC Chromatogram(Detector B-Ch1)`のデータを使います。
ディテクターA1やA2を使いたい場合、以下のようにしてください。
```
python lc.py -p 'Params.csv' 'data' --header 'A1'
```
### ベースライン推定
デフォルトではベースライン推定がオンになっており、下図のようにベースラインの変動が除去されてからデータ処理されます。

<img src="https://user-images.githubusercontent.com/31358974/213746956-54c0d5dd-c634-4714-a49c-2e275324562d.png" alt="baseline" width="600"/>

これをオフにする場合`--noBaselineCorrection`オプションをつけてください。

### ピーク検出のパラメータ
* `--peakProminence` は検出するピークの高さの閾値を設定します。デフォルトは 1.0 です。 
* `--peakWidth` は検出するピークの幅の閾値を設定します。デフォルトは 15 です。
* `--peakInclude` は検出するピーク位置を含む範囲を設定します。デフォルトはすべての範囲を検出します。
範囲 [tmin1,tmax1],[tmin2,tmax2],...を含めるとき `--peakInclude '[[tmin1,tmax1],[tmin2,tmax2],...]'` のように入力します。
* `--peakExclude` は検出するピーク位置から除外する範囲を設定します。

### グラフ描画のオプション
* `--xlim` は描画するグラフのx範囲を設定します。`--xlim '[xmin,xmax]'` のように入力します。
* `--ylim` は描画するグラフのy範囲を設定します。

### モデル関数とパラメータファイルについて
このプログラムは歪ガウシアン関数を用いてピークをモデル化します。歪ガウシアン関数は以下の式で表せます。
```math
f(x) = \textrm{ph}\cdot\left(1+\textrm{erf}\left(\textrm{sk}\cdot\frac{x-x_0}{w}\right)\right)\cdot e^{-(x-x_0)^2/w^2}
```
ここで各変数の意味は以下のようになります。

|||
|---|---|
|x0|ピーク位置|
|ph|ピーク高さ|
|w|ピーク幅|
|sk|ピーク歪み|

パラメータファイルにはこれらの情報が保存されます。

## Overview
This is a program for analyzing HPLC. This program estimates the concentration not by peak height or area, but by shape. Therefore, it is resistant to changes in peak position and shape due to concentration.

## File format
This program handles ASCII output files from Shimazu's LabSolutions. If you want to handle a different data format, you will need to modify LCData.py.

## Installation
### Install Python
Open powershell and type 'python'
### Install Python modules
In powershell, input the following commands:
```
pip install numpy
pip install matplotlib
pip install scipy
pip install tabulate
```
### Place source files
Place the source files in the folder diplayed when you enter 'pwd' in powershell at the start. If you know how to set the PATH, this is not necessary.

## Using the Program
The analysis process is as follows:
First, create a peak shape and concentration relationship model from the standard sample data `STD.csv`, and output a parameter file `Params.csv`. Then, use this model to fit the unknown sample data and determine the concentration of the compounds.

### Create a CSV file with information about the standard sample
Create a CSV file like the following example. Write the information about the compounds corresponding to the peaks that appear in order from left to right on the graph, in order from top to bottom of the table. Please put the standard sample data such as `STD1.txt` in the same folder as `STD.csv`

STD.csv:
|   | STD1.txt | STD2.txt | STD3.txt | STD4.txt | STD5.txt |
|:-:|:-:|:-:|:-:|:-:|:-:|
|Maleic Acid|50|25|12.5|6.25|3.125|
|Glucose|100|50|25|12.5|6.25|
|Xylose|50|25|12.5|6.25|3.125|
|Glycerol|50|25|12.5|6.25|3.125|
|Acetic Acid|50|25|12.5|6.25|3.125|
|Ethanol|30|15|7.5|3.75|1.875|
|n-Buthanol|30|15|7.5|3.75|1.875|

You can also include peak position information in this file. In this case, it detects peaks closest to the input peak position.

STD.csv (When including peak position information):

|   | STD1.txt | STD2.txt | STD3.txt | STD4.txt | STD5.txt |Peak position (strings that do not end with .txt, can be blank)| 
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|Maleic Acid|50|25|12.5|6.25|3.125|9.3|
|Glucose|100|50|25|12.5|6.25|10|
|Xylose|50|25|12.5|6.25|3.125|10.6|
|Glycerol|50|25|12.5|6.25|3.125|14|
|Acetic Acid|50|25|12.5|6.25|3.125|15.05|
|Ethanol|30|15|7.5|3.75|1.875|21.63|
|n-Buthanol|30|15|7.5|3.75|1.875|33.8|

### Create a parameter file
The parameter file name is set to `Params.csv` in this example. Run the following command. If you add the option `--checkParams`, the fitting results will be displayed.
```
python lc.py -s 'STD.csv' -p 'Params.csv' --checkParams
```

### Analyze a single data
'data' is a folder with multiple data files.
```
python lc.py -p 'Params.csv' 'data'
```
If you want to see the fitting results. Specify the output folder name with the '--plotDir' option
```
python lc.py -p 'Params.csv' 'data' --plotDIr 'plot'
```
If you want to output the results to a file
```
python lc.py -p 'Params.csv' 'data' -o 'out.csv'
```

## Other options

### Specifying the detector to use
By default, data from the `LC Chromatogram (Detector B-Ch1)` is used. If you would like to use Detector A1 or A2, please do so as follows.
```
python lc.py -p 'Params.csv' 'data' --header 'A1'
```
### Baseline estimation
By default, baseline estimation is turned on, and the data is processed after the variation of the baseline is removed. To turn this off, please use the `--noBaselineCorrection` option.

### Peak detection parameters
* `--peakProminence` sets the threshold for the height of the peak to be detecte. The default is 1.0.
* `--peakWidth` sets the threshold for the width of the peak to be detected. The default is 15.
*  `--peakInclude` sets the range of peak position to be included. The default is to detect all ranges.
To include ranges [tmin1,tmax1],[tmin2,tmax2],..., input `--peakInclude '[[tmin1,tmax1],[tmin2,tmax2],...]'`.
*  `--peakExclude` sets the range of peak position to be excluded.

### Graph drawing options
* `--xlim` sets the x range of the graph to be drawn. Input as  `--xlim '[xmin,xmax]'`.
* `--ylim` sets the y range of the graph to be drawn.
* 

### About model function and parameter file
This program models peaks using a skewed Gaussian function. The skewed Gaussian function can be represented by the following quation.
```math
f(x) = \textrm{ph}\cdot\left(1+\textrm{erf}\left(\textrm{sk}\cdot\frac{x-x_0}{w}\right)\right)\cdot e^{-(x-x_0)^2/w^2}
```
Here, the meaning of each variable is as follows:

|||
|---|---|
|x0|Peak position|
|ph|Peak height|
|w|Peak width|
|sk|Peak skew|

This information is saved in the parameter file.
