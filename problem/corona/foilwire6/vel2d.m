setapplicationpath('FM/poi');

porder = 2;
elemtype = 0;
nodetype = 1;
hybrid = 'hdg';

kappa = 1;
c = [0,0]; 
tau = 1;

app1.uqpk=0;
app1.source = 'source';
app1.flux = 'flux';
app1.fbou = 'fbou';
app1.fhat = 'fhat';
app1.localsolve=1;
app1.arg = {kappa,tau};
app1.bcm = [3;4;3];
app1.bcs = [0;0;0];
app1.bcd = [];
app1.bcv = [];

app1.hybrid = hybrid;
app1.tdep = false;
app1.wave = false;
app1.alag = false;
app1.flg_q = 1;
app1.flg_p = 0;
app1.flg_g = 0;

app1.fc_q = 1;
app1.fc_u = 0;
app1.fc_p = 0;

app1.nd = 2;
app1.nch  = 1;                       % Number of componets of UH
app1.nc   = app1.nch*(app1.nd+1);    % Number of componeents of UDG
app1.ncu = 1;

app1.time = [];
app1.dtfc = [];
app1.alpha = [];

% [mesh1,pva,pvw] = mkmesh_foilwire2(porder,3,xp,the,R);
% master = mkmaster(mesh1,2*porder);
% [master,mesh1] = preprocess(master,mesh1,hybrid);

UDGp = initu(mesh1,{0;0;0});
UHp=inituhat(master,mesh1.elcon,UDGp,1);

% HDG solver
[UDGp,UHp] = hdg_solve(master,mesh1,app1,UDGp,UHp,0*UDGp);

figure(1); clf; scaplot(mesh1,UDGp(:,1,:),[],1,0); 
axis equal; axis tight; colormap jet;
hold on; 
plot(pva(:,1),pva(:,2),'-k','LineWidth',1.5);
plot(pvw(:,1),pvw(:,2),'-k','LineWidth',1.5);

figure(2); clf; scaplot(mesh1,-UDGp(:,2,:),[0 1.5],1);
axis equal; axis tight; colormap jet;
hold on; 
plot(pva(:,1),pva(:,2),'-k','LineWidth',1.5);
plot(pvw(:,1),pvw(:,2),'-k','LineWidth',1.5);

figure(3); clf; scaplot(mesh1,-UDGp(:,3,:),[-.25 .25],1);
axis equal; axis tight; colormap jet;
hold on; 
plot(pva(:,1),pva(:,2),'-k','LineWidth',1.5);
plot(pvw(:,1),pvw(:,2),'-k','LineWidth',1.5);

mesh1.dgnodes(:,3,:) = -UDGp(:,2,:);
mesh1.dgnodes(:,4,:) = -UDGp(:,3,:);

return;


