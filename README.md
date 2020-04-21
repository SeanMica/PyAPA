# PyAPA
Python-based algorithm for phase auto-picking and associator

## Citation

Chang, Y.-H., S.-H. Hung, and Y.-L. Chen, 2019: A fast algorithm for automatic 
phase picker and event location: Application to the 2018 Hualien earthquake sequences. 
Terr. Atmos. Ocean. Sci., 30, 435-448, doi: 10.3319/TAO.2018.12.23.01

## About

The PyAPA is modified from [PhasePApy](https://github.com/austinholland/PhasePApy) so that 
it can be applied to automatically picking P and S phases recorded by large-scale seismic
networks or large numbers of array stations with much improved efficiency.

## Compile

To compile __aicd.f__ in __phasepicker/__
```bash
cd phasepicker
f2py -m aicd -c aicd.f
```

## Notification

The earthquake location program presented in our paper (see Citation) is owned by the
[Central Weather Bureau](https://www.cwb.gov.tw/) (CWB) in Taiwan. The CWB has used the
program (perhaps coded and modified by different people) for decades to routinely report
the origin times and hypocenters of local earthquakes in the vicinity of Taiwan. As the 
exact source information of the program is unknown, it is not included here. If you
consider to use our auto-picking and associator algorithm for your own earthquake location
problems, you can replace the location program with any available open-source programs by
simply modifying the input and output parts in __associator/search_hypo.py__.
