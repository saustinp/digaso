
etaw = linspace(0,2*pi,50); 
x0 = cos(etaw)*rw+xw;
y0 = sin(etaw)*rw+yw;

nu = 1e-3;
ds = loginc(linspace(0.002,0.1,500),2);
Ua(:,1,:) = UDG(:,1,:);
Ua(:,2,:) = UDG(:,3,:);
Ua(:,3,:) = UDG(:,5,:);    
[Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);        
figure(2);clf;
plot(etaw,Peek1,etaw,R*ones(size(etaw)));
pause(1);

if max(Peek1)>R
    disp('max(Peek1)>R');
    kappa1 = kappa;
    kappa2 = kappa;
    while (max(Peek1)>R)
        kappa2 = kappa2*1.2;
        app.arg = {1,nu,eps0,kappa2,xw,yw,t1,t2,Wind,tau};     
        [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
        Ua(:,1,:) = UDG(:,1,:);
        Ua(:,2,:) = UDG(:,3,:);
        Ua(:,3,:) = UDG(:,5,:);    
        [Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);      
        figure(2);clf; plot(etaw,Peek1,etaw,R*ones(size(etaw))); pause(1);
        [kappa1 kappa2]
        if max(Peek1)>R
            kappa1 = kappa2;
        end
    end
else
    disp('max(Peek1)<R');
    kappa1 = kappa;
    kappa2 = kappa;
    while (max(Peek1)<R)
        kappa1 = kappa1/1.2;
        app.arg = {1,nu,eps0,kappa1,xw,yw,t1,t2,Wind,tau};     
        [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
        Ua(:,1,:) = UDG(:,1,:);
        Ua(:,2,:) = UDG(:,3,:);
        Ua(:,3,:) = UDG(:,5,:);    
        [Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);     
        figure(2);clf; plot(etaw,Peek1,etaw,R*ones(size(etaw))); pause(1);
        [kappa1 kappa2]
        if max(Peek1)<R
            kappa2 = kappa1;
        end
    end    
end

while abs(max(Peek1)-R)>0.001    
    kappa = (kappa1+kappa2)/2;    
    app.arg = {1,nu,eps0,kappa,xw,yw,t1,t2,Wind,tau};     
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

    Ua(:,1,:) = UDG(:,1,:);
    Ua(:,2,:) = UDG(:,3,:);
    Ua(:,3,:) = UDG(:,5,:);    
    [Peek1,x1,y1,s1,Ex1,Ey1] = PeekIntegral(mesh, 0, 1, 0*Ua, Ua, x0, y0, ds, 2);
    figure(2);clf; plot(etaw,Peek1,etaw,R*ones(size(etaw))); pause(1);
    if max(Peek1)>R
        kappa1 = kappa;
    else
        kappa2 = kappa;
    end
    [kappa1 kappa2 abs(max(Peek1)-R)]
end

