# 浸透・地表流解析に関するコード

Green-Amptモデルによる浸透解析と，Diffusion Waveモデルによる地表流解析を実行するfortranコード

具体的な式については，[Engineering Geologyの論文](https://doi.org/10.1016/j.enggeo.2022.106808)に記載

### (※)注意
　論文中の式(3), (5), (10), (11) に誤植あり．以下が正しい式．
  ([Dolojan et al. (2021)](https://doi.org/10.1007/S10346-020-01617-X)や
  [Sayama et al. (2012)](https://doi.org/10.1080/02626667.2011.644245)も参照)   

  <img src="./img/EG_eq(3).png" width="70%">

  <img src="./img/EG_eq(5).png" width="70%">
  
  <img src="./img/EG_eq(10).png" width="60%">
  
  <img src="./img/EG_eq(11).png" width="60%">


## 入力ファイル
 
### input/: 解析領域の座標データと解析パラメータ

 - coordinate.txt       
  座標データ: 1列目にx座標，2列目にy座標，3列目に標高値   
  座標はグリッドを想定
  以下の図に示す順番でデータを格納    

  <img src="./img/coordinate.png" width="40%">

 - num_node.txt  
  x,y方向の節点数

 - inp_par.txt   
  解析パラメータのデータを格納（空間的に一様なものとして仮定）   
  1行目：飽和透水係数(m/s)，Suction head at Wetting front(m)   
  2行目：初期体積含水率，飽和体積含水率     
  3行目：マニングの粗度係数(s/m^(1/3))，地表流解析の実施有無(1:地表流計算なし, 2:あり)     



### raindata/: 降雨データ

 - rain.txt  
  降水量の時空間分布データ(行方向：時間，列方向：空間，単位：m/s)       
  空間的に一様な場合は1列のデータを作成    

 - time_condition.txt      
  降雨の時間ステップ数(=rain.txtの行数)と1時間ステップの時間間隔(s)     
  設定した時間間隔で結果も出力

 - coordinate_rain.txt   
  降水量の座標データ


## コード内の主な変数

 - 座標・時間変数     

 | 変数名   | 次元    | 説明                                      |   
 | :------: | :----:  | :-------------                            |   
 | nx, ny   | -       | x,y方向の節点数，途中から要素数を表現     |   
 | xx (yy)  | nx (ny) | x,y方向の座標(m)                          |   
 | gl       | nx×ny  | 標高値の行列(m)                           |   
 | dx, dy   | -       | x,y方向の解像度(m)                        |   
 | nti      | -       | 計算の総時間ステップ数                    |   
 | tstep    | -       | 時間ステップ数                            |    
 | dti      | -       | 降雨データの時間間隔(出力値の時間間隔) (s)|    
 | dt       | -       | 解析の時間間隔 (s)                        |    
 | slp      | nx×ny  | 斜面勾配                                  |    


 - 降水量データ      

 | 変数名   | 次元    | 説明                                 |    
 | :------: | :-----: | :-------------                       |   
 | nr       | -       | 降水量データの空間次元数             |   
 | rain     | nti×nr | 降水量の時空間分布(m/s)              |   
 | raint    | nr      | ある時間での降水量(m/s)              |    
 | rid      | nx×ny  | 各計算点が参照する降水量データの番号 |   


 - 浸透解析     

 | 変数名   | 次元    | 説明                                 |    
 | :------: | :-----: | :-------------                       |   
 | k0       | -       | 飽和透水係数(m/s)                    |    
 | thi, ths | -       | 初期(飽和)体積含水率(-)              |    
 | psi      | -       | Suction Head at Wetting front(m)     |    
 | fn       | -       | ある時間での累積浸透量(m)            |    
 | zzn      | -       | 鉛直方向の浸透深さ(m)                |    


 - 地表流解析     

 | 変数名   | 次元    | 説明                                 |    
 | :------: | :-----: | :-------------                       |   
 | qsx,qsy  | 2       | x(y)方向の単位幅流量(m^2/s)          |    
 | nn       | -       | マニングの粗度係数(s/m^(1/3))        |    
 | hsn      | nx×ny  | 地表水位 (m)                         |    
 | hmx      | nx×ny  | 地表水位の最大値(m)                  |    



## 計算方法

 予めoutputフォルダを作成しておく必要あり
 
 - コンパイル
  
 シングル計算    

 > gfortran infil_sflow.f90 -o run    


 並列化(OpenMP)    

 > gfortran -fopenmp infil_sflow.f90 -o run 

 runの部分は別の名前でもOK     



## 出力データ
 
 - water-depth-level_****.txt    
  出力結果：各要素の値として出力（=num_nodeで設定した数より1ずつ小さい空間次元）   
  ****は時間ステップ番号(0001,0002,...)    
  1列目：x座標   
  2列目：y座標   
  3列目：Wetting front depth (m)   
  4列目：地表水位 (m)    

 - max_s_water_level.txt   
   地表水位の最大値のデータ    
   地表流計算を実施した時のみ出力



## 参考文献
 - [K.Tozato et al.(2022)](https://doi.org/10.1016/j.enggeo.2022.106808), Limit equilibrium method-based 3D slope stability analysis for wide area considering influence of rainfall, Engineering Geology, Vol. 308, pp. 106808.      
 - [N. L. J. Dolojan et al. (2021)](https://doi.org/10.1007/S10346-020-01617-X), Mapping method of rainfall-induced landslide hazards by infiltration and slope stability analysis, Landslides, Vol. 6, pp. 2039–2057. 
 - [Sayama et al. (2012)](https://doi.org/10.1080/02626667.2011.644245), Rainfall-runoffinundation analysis of the 2010 Pakistan flood in the Kabul river basin. Hydrol. Sci. J., Vol. 57, pp. 298–312. 




