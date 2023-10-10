
for j=2:622
    fn1 = ['mv' num2str(j) '.png'];
    fn2 = sprintf('%s%05d.png','mv',j-1);    
    movefile(fn1,fn2);
end

