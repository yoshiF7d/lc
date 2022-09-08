1:python をインストールする
powershellを立ち上げ、pythonと打つ

2:moduleをインストールする
pip install numpy
pip install matplotlib
pip install scipy
pip install tabulate

3:lc.pyを powershellを始めに立ち上げて、pwdと入力したときに表示されるフォルダに置く

4：標準サンプルの情報が入ったcsvファイルを作る
以下のように作る。STD1.txtなどは、STD.csvと同じフォルダにおく。左から順番にあらわれるピークをに対応する化合物の情報を上から順番に書く。
STD.csv:
,STD1.txt,STD2.txt,STD3.txt,STD4.txt,STD5.txt
Maleic Acid,10,5,2.5,1.25,0.625
Glucose,100,50,25,12.5,6.25
Xylose,50,25,12.5,6.25,3.125

5-7 これからファイル名をpowershellに入力するときは、エクスプローラーでファイルをクリックして、powershell上にドラッグアンドドロップするとファイル名のフルネームが入力されるので便利（このとき""は必要)

5：パラメータファイルを作成する。パラメータファイルを　Params.csv　とする。
python lc.py -s "STD.csv" -p "Params.csv" --checkParams

6-1：1データを分析する。data.txt　はLabSolutionsからASCII出力として出力したファイル
python lc.py -p "Params.csv" "data.txt"

6-2：複数データを分析する。data はデータが入ってるディレクトリ
python lc.py -p "Params.csv" "data"

6-3：複数データを分析する。フォルダ　plot　にフィッティングの結果をplotする
python lc.py -p "Params.csv" "data" --plotDir "plot"

6-3：複数データを分析する。ファイル　out.csv　に結果を出力する
python lc.py -p "Params.csv" "data" "out.csv"

デフォルトではディテクターBのデータを取得する。Aを使うときはlc.pyの
HEADER = "LC Chromatogram(Detector B-Ch1)"
を
HEADER = "LC Chromatogram(Detector A-Ch2)"
などにする。

#ベースライン推定がおかしい場合 
--noBLE オプションをつける。

#ピークの位置、面積、高さから検量線を作成する場合
以下のようなパラメータでParam.csvファイルを編集する。
x0 : ピーク位置
ph : ピーク高さ
w : 面積 / (ph * √π)
sk : 0 


######################################################################################################################################################

1:install python
launch powershell type 'python'

2:insall modules
pip install numpy
pip install matplotlib
pip install scipy
pip install tabulate

3:place lc.py on the home directory which is a folder shown by typing 'pwd' when you initially launch powershell

4:make a csv file that specifies information of standard samples.
place STD1.txt STD2.txt etc on the same folder as where STD.csv is.
The peaks of chemical species appear from left to right on the plot. That order correspond to the order in the STD.csv from top to the
bottom.

example:

STD.csv:
,STD1.txt,STD2.txt,STD3.txt,STD4.txt,STD5.txt
Maleic Acid,10,5,2.5,1.25,0.625
Glucose,100,50,25,12.5,6.25
Xylose,50,25,12.5,6.25,3.125

5 when you type filename on powershell, you can instead find that file in explorer, then drag and drop that file on the powershell.
(you still need to type "")

5：Make parameter file.
python lc.py -s "STD.csv" -p "Params.csv" --checkParams

6-A：analyse single data. "data.txt" is a ASCII file exported from LabSolutions 
python lc.py -p "Params.csv" "data.txt"

6-B：analyse multiple data. "data" is a folder that files are in
python lc.py -p "Params.csv" "data"

6-C：analyse multile data. Plot the result of fitting on folder "plot"
python lc.py -p "Params.csv" "data" --plotDir "plot"

6-D：analyse multiple data. Output result on csv file
python lc.py -p "Params.csv" "data" "out.csv"

By default this program use detector B signal. when you use A
change the line on lc.py
HEADER = "LC Chromatogram(Detector B-Ch1)"
to
HEADER = "LC Chromatogram(Detector A-Ch2)"

#when baseline estimation is erroneous.
use --noBLE option

#when you want to add a chemical specie whose peak position,peak height,peak area are known.
Add the chemical specie to the Param.csv file as below
x0 : peak position
ph : peak height
w : (peak area) / (ph * sqrt(pi))
sk : 0 

