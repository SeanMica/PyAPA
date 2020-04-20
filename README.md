# PyAPA
Python based Auto-Picking and Associator

## Reference 

Chang, Y.-H., S.-H. Hung, and Y.-L. Chen, 2019: A fast algorithm for automatic 
phase picker and event location: Application to the 2018 Hualien earthquake sequences. 
Terr. Atmos. Ocean. Sci., 30, 435-448, doi: 10.3319/TAO.2018.12.23.01

## About

The PyAPA is modified based on PhasePApy[https://github.com/austinholland/PhasePApy].
I modified some of algorithms of PhasePApy and improved the efficiency so that it can
use on larger seiemometer arrays.

## Compile

To compile aicd.f in phasepicker/
```bash
cd phasepicker
f2py -m aicd -c aicd.f
```

## Notification

Because the current earthquake location program is too old, it is impossible to find the source 
code owner, so it cannot be included in this source code. 
If you want to use code, please consider to find the available earthquake location programs and modify the input and output parts in __associator/search_hypo.py__.