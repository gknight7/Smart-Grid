%% Comparison linearized load flow with conventional load flow
% for balanced operation
clc
clear all
%        Node1 - Node2 - rkm -   xkm -    bkm/2 - P - Q - alpha    
Feeder = [  1    2     0.0025    0.0026    0.0001    0.0000    0.0000    0
            2    3     0.0034    0.0036    0.0001    0.0000    0.0000    0
            3    4     0.0040    0.0023    0.0000    0.0850    0.0400    2
            4    5     0.0013    0.0008    0.0000    0.0000    0.0000    0
            5    6     0.0040    0.0023    0.0000    0.0850    0.0400    2
            5    7     0.0022    0.0012    0.0000    0.0000    0.0000    0
            7    8     0.0042    0.0013    0.0000    0.0420    0.0210    0
            7    9     0.0022    0.0012    0.0000    0.0850    0.0400    1
            9    10    0.0038    0.0021    0.0000    0.0420    0.0210    0
            10   11    0.0043    0.0025    0.0000    0.1400    0.0700    1
            11   12    0.0027    0.0015    0.0000    0.1260    0.0620    0
            12   13    0.0027    0.0015    0.0000    0.0000    0.0000    0
            13   14    0.0026    0.0008    0.0000    0.0850    0.0400    0
            13   15    0.0027    0.0015    0.0000    0.0420    0.0210    1
            10   16    0.0068    0.0022    0.0000    0.0000    0.0000    0
            16   17    0.0167    0.0054    0.0001    0.0420    0.0210    2
            16   18    0.0026    0.0008    0.0000    0.0850    0.0400    0
             3   19    0.0031    0.0010    0.0000    0.0420    0.0210    0
            19   20    0.0019    0.0011    0.0000    0.0420    0.0210    0
            20   21    0.0026    0.0008    0.0000    0.1260    0.0630    0
            20   22    0.0037    0.0012    0.0000    0.0420    0.0210    1
             2   23    0.0024    0.0014    0.0000    0.0850    0.0400    0
            23   24    0.0035    0.0020    0.0000    0.0000    0.0000    0
            24   25    0.0010    0.0003    0.0000    0.0380    0.0180    1
            25   26    0.0068    0.0022    0.0000    0.0850    0.0400    2
            24   27    0.0054    0.0031    0.0000    0.0850    0.0400    0
            27   28    0.0040    0.0023    0.0000    0.0000    0.0000    0
            28   29    0.0037    0.0012    0.0000    0.0420    0.0210    0
            27   30    0.0120    0.0039    0.0000    0.0000    0.0000    0
            30   31    0.0016    0.0005    0.0000    0.1610    0.0800    1
            30   32    0.0099    0.0032    0.0000    0.0420    0.0210    2
             2   33    0.0052    0.0017    0.0000    0.0000    0.0000    0
            33   34    0.0031    0.0010    0.0000    0.0850    0.0400    0
            33   35    0.0042    0.0013    0.0000    0.0930    0.0440    2];
% This function calculates the load flow using the back-forward sweep algorithm        
res = Radial_Load_Flow(Feeder);
% This function calculates the load flow using the linear method
resap = Linear_Load_Flow(Feeder);
% display a comparison
figure(1)
bar(abs(res.Vs-resap.Vs));
grid on
xlabel('NODE');
ylabel('\epsilon_{k}');
title('Single phase case')
format short
disp('    V(pu)    Ang(deg)   V(aprox) Ang(aprox)')
disp([abs(res.Vs),angle(res.Vs*180/pi),abs(resap.Vs),angle(resap.Vs*180/pi)])        