# HW1
## Introduction
資料夾中總共有2份程式碼，huffman.py這個檔案主要可執行basic和extended huffman coding(只要symbol的bits數為8的倍數)，而adaptive_huffman.cpp這個檔案則是可執行8bits adaptive huffman coding。  
一開始我也想把這兩個演算法都用python實現，可惜python的效能太糟了，所以就變成了兩份程式碼。
## How to setup and execute
### huffman.py
環境請使用python3，python可能會跑不動。  
因為library會使用到
```

import six
import argparse 
import os
```
通常python都會內建os和argparse，six我不確定是不是都內建好，沒內建好的話請使用
``` pip install six```

接下來講解如何使用huffman.py。  

Command:  
```python3 huffman.py mode infilename [--outfilename OUTFILENAME] [--symbol_length SYMBOL_LENGTH]```   
其中mode為必填，mode包含了compress，decompress和check_pmf。   
infilename也是必填，infilename就是填要壓縮/解壓縮的檔案名稱，包含副檔名也要填。  
而--outfilename是選填，--outfilename就是填壓縮/解壓縮的輸出檔案名稱，包含副檔名也要填。  
而--symbol_length也是選填，代表著每個symbol的byte長度，default為1byte(8bits)，而這個flag只支援compress和decompress mode，check_pmf我寫死了8bits不支援。

懶人包:  
想對alexnet.pth執行8 bits Huffman coding，並得到輸出檔名為alexnet.pth.hcps1:   
```python3 huffman.py compress alexnet.pth --symbol_length=1```   
想對alexnet.pth.hcps1執行8 bits Huffman decoding，並將輸出檔案取名為alexnet_copy.pth:   
```python3 huffman.py decompress alexnet.pth.hcps1 --outfilename=alexnet_copy.pth --symbol_length=1``` 

想對alexnet.pth執行32 bits Huffman coding，並得到輸出檔名為alexnet.pth.hcps4:   
```python3 huffman.py compress alexnet.pth --symbol_length=4```   
想對alexnet.pth.hcps4執行32 bits Huffman decoding，並將輸出檔案取名為alexnet_copy2.pth:   
```python3 huffman.py decompress alexnet.pth.hcps4 --outfilename=alexnet_copy2.pth --symbol_length=4``` 

想對alexnet.pth執行16 bits Huffman coding，並得到輸出檔名為alexnet.pth.hcps2:   
```python3 huffman.py compress alexnet.pth --symbol_length=2```   
想對alexnet.pth.hcps2執行16 bits Huffman decoding，並將輸出檔案取名為alexnet_copy3.pth:   
```python3 huffman.py decompress alexnet.pth.hcps2 --outfilename=alexnet_copy3.pth --symbol_length=2``` 

想看alexnet.pth在不同區間下，機率分布前10名高的symbol:    
```python3 huffman.py check_pmf alexnet.pth``` 

### adaptive_huffman.cpp
環境:要有c++11   
編譯:
``` g++ -std=c++11 adaptive_huffman.cpp -o adaptive_huffman```

Command:   
```./adaptive_huffman [mode] [inputfilename] [outputfilename]```

其中mode為必填，必須填compress或decompress。  
其中inputfilename也是必填，inputfilename就是填要壓縮/解壓縮的檔案名稱，包含副檔名也要填。
而outputfilename是選填，outputfilename就是填壓縮/解壓縮的輸出檔案名稱，包含副檔名也要填。  

懶人包:   
想對alexnet.pth執行8 bits Adaptive Huffman coding，並得到輸出檔名為alexnet.pth.ahcps:    
```./adaptive_huffman compress alexnet.pth```    
想對alexnet.pth.ahcps執行8 bits Huffman decoding，並將輸出檔案取名為alexnet_copy2.pth:   
```./adaptive_huffman decompress alexnet.pth.ahcps alexnet_copy2.pth```    