function [LCp,UCp,Cpxfoil,Cfxfoil]=naca0012_exp_data

% Lower surface
% x/c    y/c      C_p (M=0.3,Re=1.85E6,alpha=0)
LCp=[    
0.0000   0.0       0.9980  
0.0197  -0.02334  -0.1115 
0.0296  -0.02823  -0.2491
0.0590  -0.03813  -0.3847
0.0697  -0.04084  -0.4068
0.0989  -0.04670  -0.4336
0.1195  -0.04986  -0.4188
0.1801  -0.05613  -0.4162
0.1997  -0.05742  -0.4087
0.2200  -0.05841  -0.3987
0.2400  -0.05910  -0.3818
0.2902  -0.05996  -0.3630
0.3200  -0.05996  -0.3495
0.3495  -0.05953  -0.3392
0.4103  -0.05750  -0.3106
0.4396  -0.05622  -0.2261
0.4698  -0.05470  -0.2331
0.5001  -0.05292  -0.2116
0.5302  -0.05092  -0.1948
0.5899  -0.04647  -0.1776
0.6201  -0.04398  -0.1514
0.6496  -0.04139  -0.1455
0.7099  -0.03564  -0.1027
0.7393  -0.03263  -0.0797
0.7699  -0.02941  -0.0746
0.7998  -0.02609  -0.0394
0.8298  -0.02276  -0.0247
0.8598  -0.01924   0.0064
0.8895  -0.01569   0.0240
0.9195  -0.01193   0.0686
0.9497  -0.00804   0.0996];

% Upper surface
% x/c    y/c      C_p (M=0.3,Re=1.85E6,alpha=0)
UCp=[    
0.9701  0.00547   0.1172
0.9102  0.01323   0.0513
0.8797  0.01702   0.0150
0.8503  0.02056  -0.0089
0.8198  0.02408  -0.0345
0.7900  0.02738  -0.0772
0.7304  0.03367  -0.0965
0.6994  0.03679  -0.1157
0.6696  0.03961  -0.1413
0.6396  0.04225  -0.1515
0.6095  0.04486  -0.1664
0.5792  0.04727  -0.1709
0.5494  0.04953  -0.2036
0.5192  0.05171  -0.2039
0.4893  0.05364  -0.2408
0.4591  0.05538  -0.2368
0.3994  0.05804  -0.2935
0.3691  0.05900  -0.3266
0.3393  0.05979  -0.3367
0.3102  0.06005  -0.3575
0.2793  0.05996  -0.3651
0.2492  0.05944  -0.3939
0.2293  0.05879  -0.4029
0.2096  0.05791  -0.4108
0.1895  0.05680  -0.4189
0.1697  0.05536  -0.4326
0.1496  0.05352  -0.4405
0.1293  0.05121  -0.4263
0.1098  0.04848  -0.4361
0.0995  0.04681  -0.4284
0.0893  0.04497  -0.4101
0.0803  0.04319  -0.4129
0.0505  0.03577  -0.3431
0.0298  0.02841  -0.2848
0.0191  0.02314  -0.1336];

%(M=0.3,Re=1.85E6,alpha=0)
Cpxfoil=[1.00000  0.22268
   0.99276  0.20162
   0.97987  0.16184
   0.96465  0.12677
   0.94764  0.09541
   0.92954  0.06845
   0.91089  0.04493
   0.89197  0.02448
   0.87294  0.00614
   0.85385 -0.01042
   0.83472 -0.02574
   0.81558 -0.03989
   0.79643 -0.05317
   0.77727 -0.06578
   0.75811 -0.07772
   0.73895 -0.08928
   0.71978 -0.10003
   0.70061 -0.11050
   0.68145 -0.12024
   0.66229 -0.12939
   0.64314 -0.13753
   0.62400 -0.14518
   0.60487 -0.15433
   0.58575 -0.17531
   0.56665 -0.21640
   0.54757 -0.22028
   0.52850 -0.22750
   0.50946 -0.23636
   0.49045 -0.24605
   0.47146 -0.25632
   0.45250 -0.26695
   0.43358 -0.27788
   0.41470 -0.28897
   0.39587 -0.30022
   0.37708 -0.31154
   0.35834 -0.32298
   0.33966 -0.33429
   0.32106 -0.34559
   0.30252 -0.35686
   0.28407 -0.36767
   0.26572 -0.37851
   0.24748 -0.38885
   0.22937 -0.39867
   0.21141 -0.40800
   0.19363 -0.41645
   0.17607 -0.42405
   0.15879 -0.43038
   0.14187 -0.43509
   0.12542 -0.43798
   0.10962 -0.43808
   0.09469 -0.43515
   0.08088 -0.42833
   0.06844 -0.41662
   0.05753 -0.39964
   0.04813 -0.37807
   0.04014 -0.34887
   0.03336 -0.31389
   0.02761 -0.27106
   0.02270 -0.22001
   0.01849 -0.15729
   0.01486 -0.08459
   0.01172  0.00460
   0.00902  0.11011
   0.00670  0.23346
   0.00475  0.37422
   0.00314  0.53322
   0.00186  0.69486
   0.00093  0.84224
   0.00032  0.95698
   0.00003  1.01631];

%(M=0.3,Re=1.85E6,alpha=0)
% s        x        y     Ue/Vinf    Dstar     Theta      Cf
 Cfxfoil=[0.00000  1.00000  0.00126  0.88229  0.003272  0.001998  0.001664
   0.00731  0.99276  0.00227  0.89403  0.003059  0.001898  0.001806
   0.02032  0.97987  0.00405  0.91583  0.002711  0.001729  0.002085
   0.03568  0.96465  0.00611  0.93466  0.002451  0.001592  0.002322
   0.05284  0.94764  0.00838  0.95121  0.002242  0.001476  0.002532
   0.07110  0.92954  0.01074  0.96522  0.002076  0.001379  0.002709
   0.08990  0.91089  0.01311  0.97730  0.001936  0.001294  0.002865
   0.10896  0.89197  0.01547  0.98769  0.001817  0.001220  0.003002
   0.12814  0.87294  0.01779  0.99693  0.001711  0.001153  0.003128
   0.14736  0.85385  0.02007  1.00520  0.001615  0.001091  0.003247
   0.16662  0.83472  0.02230  1.01280  0.001525  0.001032  0.003361
   0.18588  0.81558  0.02449  1.01977  0.001441  0.000977  0.003473
   0.20515  0.79643  0.02663  1.02627  0.001360  0.000924  0.003586
   0.22442  0.77727  0.02871  1.03241  0.001282  0.000873  0.003701
   0.24369  0.75811  0.03075  1.03820  0.001207  0.000823  0.003818
   0.26296  0.73895  0.03275  1.04377  0.001134  0.000774  0.003940
   0.28223  0.71978  0.03469  1.04893  0.001063  0.000726  0.004060
   0.30148  0.70061  0.03658  1.05393  0.000994  0.000679  0.004176
   0.32074  0.68145  0.03842  1.05857  0.000929  0.000633  0.004271
   0.33998  0.66229  0.04021  1.06291  0.000869  0.000587  0.004310
   0.35920  0.64314  0.04194  1.06675  0.000823  0.000544  0.004212
   0.37842  0.62400  0.04361  1.07035  0.000802  0.000504  0.003817
   0.39762  0.60487  0.04523  1.07465  0.000835  0.000469  0.002878
   0.41680  0.58575  0.04679  1.08444  0.000970  0.000436  0.001431
   0.43596  0.56665  0.04828  1.10338  0.001188  0.000397  0.000379
   0.45510  0.54757  0.04970  1.10516  0.001214  0.000391  0.000304
   0.47421  0.52850  0.05105  1.10845  0.001200  0.000383  0.000294
   0.49329  0.50946  0.05233  1.11248  0.001169  0.000374  0.000306
   0.51235  0.49045  0.05353  1.11687  0.001132  0.000364  0.000327
   0.53136  0.47146  0.05465  1.12152  0.001091  0.000354  0.000357
   0.55035  0.45250  0.05568  1.12630  0.001049  0.000344  0.000391
   0.56929  0.43358  0.05661  1.13120  0.001007  0.000333  0.000430
   0.58819  0.41470  0.05745  1.13616  0.000966  0.000323  0.000471
   0.60704  0.39587  0.05818  1.14116  0.000925  0.000312  0.000516
   0.62584  0.37708  0.05880  1.14618  0.000886  0.000302  0.000564
   0.64458  0.35834  0.05930  1.15123  0.000847  0.000291  0.000617
   0.66326  0.33966  0.05968  1.15620  0.000809  0.000281  0.000672
   0.68187  0.32106  0.05992  1.16115  0.000772  0.000270  0.000731
   0.70041  0.30252  0.06002  1.16607  0.000735  0.000260  0.000797
   0.71886  0.28407  0.05996  1.17078  0.000701  0.000249  0.000863
   0.73721  0.26572  0.05974  1.17547  0.000665  0.000239  0.000939
   0.75546  0.24748  0.05935  1.17994  0.000631  0.000228  0.001021
   0.77358  0.22937  0.05877  1.18416  0.000597  0.000217  0.001108
   0.79156  0.21141  0.05799  1.18817  0.000563  0.000207  0.001206
   0.80936  0.19363  0.05699  1.19178  0.000530  0.000196  0.001313
   0.82696  0.17607  0.05577  1.19503  0.000497  0.000185  0.001434
   0.84430  0.15879  0.05431  1.19773  0.000464  0.000174  0.001569
   0.86131  0.14187  0.05259  1.19973  0.000432  0.000162  0.001720
   0.87788  0.12542  0.05061  1.20096  0.000399  0.000151  0.001900
   0.89383  0.10962  0.04837  1.20100  0.000367  0.000140  0.002100
   0.90897  0.09469  0.04591  1.19975  0.000336  0.000128  0.002335
   0.92303  0.08088  0.04326  1.19685  0.000305  0.000117  0.002612
   0.93577  0.06844  0.04050  1.19186  0.000276  0.000107  0.002927
   0.94704  0.05753  0.03771  1.18458  0.000249  0.000097  0.003258
   0.95683  0.04813  0.03497  1.17528  0.000224  0.000088  0.003688
   0.96525  0.04014  0.03233  1.16259  0.000202  0.000080  0.004102
   0.97248  0.03336  0.02979  1.14721  0.000181  0.000072  0.004590
   0.97873  0.02761  0.02736  1.12814  0.000163  0.000066  0.005085
   0.98416  0.02270  0.02502  1.10503  0.000146  0.000059  0.005701
   0.98895  0.01849  0.02275  1.07603  0.000132  0.000054  0.006216
   0.99320  0.01486  0.02054  1.04152  0.000117  0.000049  0.006913
   0.99701  0.01172  0.01837  0.99770  0.000105  0.000044  0.007521
   1.00046  0.00902  0.01622  0.94349  0.000094  0.000040  0.008136
   1.00362  0.00670  0.01407  0.87622  0.000084  0.000036  0.008545
   1.00654  0.00475  0.01191  0.79305  0.000074  0.000032  0.008898
   1.00924  0.00314  0.00974  0.68788  0.000067  0.000029  0.008764
   1.01177  0.00186  0.00755  0.56214  0.000061  0.000027  0.007860
   1.01416  0.00093  0.00536  0.41680  0.000056  0.000025  0.006453
   1.01643  0.00032  0.00317  0.25224  0.000053  0.000024  0.004138
   1.01858  0.00003  0.00104  0.08322  0.000053  0.000024  0.001386];

   