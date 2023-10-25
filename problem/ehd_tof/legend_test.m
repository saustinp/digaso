h(1) = plot(1:10,sin(1:10),'r');
hold on
h(2) = plot(1:10,cos(1:10),'r');
h(3) = plot(1:10,cos(1:10) + sin(1:10),'b');
legend(h([1,3]),{'data1','data2'})