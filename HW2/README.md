# HW1
## Introduction
資料夾中總共有1份程式碼，arithmetic.py這個檔案主要可執行256-symbol/binary Arithmetic Coding，也可以選擇使用fixed probability model和PPM。 

## How to setup and execute
### <font color="#ff0">注意</font>
<font color="#ff0">若想驗證下列指令，強烈建議同時多開terminal執行，應該總共有8種演算法要驗證，每筆花費20分鐘至5小時不等，強烈建議開8個terminal，6小時內就完事了。</font>

### arithmetic.py
環境請使用python3，python可能會跑不動。  
因為library會使用到
```

import six
import argparse 

```
通常python都會內建argparse，six我不確定是不是都內建好，沒內建好的話請使用
``` pip install six```

接下來講解如何使用arithmetic.py。  

Command:  
```python3 arithmetic.py infilename outfilename [--algorithm ALGORITHM] [--order ORDER]```   

其中infilename必填，infilename就是填要壓縮前的檔案名稱，包含副檔名也要填。  
outfilename也是選填，outfilename就是填要壓縮後的檔案名稱，包含副檔名也要填。  
而--algorithm是選填，主要分成下列四種模式:
1. 256_AC 
2. binary_AC
3. 256_PPM
4. binary_PPM

256_AC代表著256-symbol AC with fixed probability model，同時也是default algorithm。  
binary_AC代表著binary AC with fixed probability model。  
256_PPM代表著256-symbol AC with PPM。
binary_PPM代表著binary AC with PPM。

最後--order是在使用--algorithm=256_PPM和--algorithm=binary_PPM才能使用，可以填入任何非負整數，default integer=2。

<font color="#ff0">懶人包:  </font>
1. 想對alexnet.pth進行256-symbol AC with fixed probability model，並輸出256AC.cps，請打:   
```python3 arithmetic.py alexnet.pth 256AC.cps --algorithm=256_AC```
2. 想對alexnet.pth進行binary AC with fixed probability model，並輸出binaryAC.cps，請打:   
```python3 arithmetic.py alexnet.pth binaryAC.cps --algorithm=binary_AC```
3. 想對alexnet.pth進行256-symbol AC with PPM(order=2)，並輸出256PPMod2.cps，請打:   
```python3 arithmetic.py alexnet.pth 256PPMod2.cps --algorithm=256_PPM --order=2```
4. 想對alexnet.pth進行256-symbol AC with PPM(order=1)，並輸出256PPMod1.cps，請打:   
```python3 arithmetic.py alexnet.pth 256PPMod1.cps --algorithm=256_PPM --order=1```
5. 想對alexnet.pth進行256-symbol AC with PPM(order=0)，並輸出256PPMod0.cps，請打:   
```python3 arithmetic.py alexnet.pth 256PPMod0.cps --algorithm=256_PPM --order=0```
6. 想對alexnet.pth進行binary AC with PPM(order=2)，並輸出binaryPPMod2.cps，請打:   
```python3 arithmetic.py alexnet.pth binaryPPMod2.cps --algorithm=binary_PPM --order=2```
7. 想對alexnet.pth進行binary AC with PPM(order=1)，並輸出binaryPPMod1.cps，請打:   
```python3 arithmetic.py alexnet.pth binaryPPMod1.cps --algorithm=binary_PPM --order=1```
8. 想對alexnet.pth進行binary AC with PPM(order=0)，並輸出binaryPPMod0.cps，請打:   
```python3 arithmetic.py alexnet.pth binaryPPMod0.cps --algorithm=binary_PPM --order=0```
