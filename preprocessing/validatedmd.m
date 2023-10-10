hybrid='hdg';
ngrid=9;
porder = 4;
morder = 4;
nproc = 4;
elemtype=0;
nodetype=1;
preconditioner = 0;
overlappinglevel = 1;
preconditionerlevel=1;

kappa = 1; 
tau = 1;
app.nfile=2;
app.appname = 'poisson';
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.localsolve=1;
app.arg = {kappa,tau};
app.bcm = [5];
app.bcs = [0];
app.bcd = [];
app.bcv = [];

app.overlappinglevel=overlappinglevel;
app.preconditionerlevel=preconditionerlevel;
app.preconditioner=preconditioner;
app.hybrid = hybrid;
app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 0;

app.morder = [morder morder];
app.porder = [porder porder];
app.nodetype = nodetype;
app.pgauss = 2*[porder porder];
app.pgaussR = 2*[porder porder];
app.quadtype = [0 0 0];
app.dt = 0;

app.nd = 2;
app.ncd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 1;
app.ncq = app.nd;
app.ncp = 0;
app.nco = 0;

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh = mkmesh_square(ngrid,ngrid,morder,1,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*(porder+1));
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{1;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

check = 1;
elementtype = elemtype*ones(mesh.ne,1);
appmpi=digasopre(app,'poi',mesh.p,mesh.t'-1,mesh.dgnodes,UDG,UH,[],elementtype,{'true'},nproc,0,check);

