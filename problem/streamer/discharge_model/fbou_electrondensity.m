function [fh,fh_udg,fh_uh] = fbou_electrondensity(ib,ui,nl,p,udg,uh,param,time)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = nc/3;     % 3 components to UDG: (U,QX,QY) -> assuming 2D. For example, in 2D if U has 4 equations, nc will be 12

%%%%%%%%%       Copied from fluxgen_electrondensity.m
% Read in values from the p vector -> this normally only contains the (r,z) coordinates but we can also pass in additional fields as appended elements
% r = p(:,1);
r = 1;
tau = param{end};

switch ib
    case 1 % Electrode: species set to homogeneous nuemann and potential set to dirichlet
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);
        
        % ORIGINAL BC: 0 NEUMANN FOR SPECIES
        % [fh_tmp,fh_udg_tmp,fh_uh_tmp] = fhat_axisymmetric(nl,p,udg,uh,param,time);

        % fh(:,[1,2]) = fh_tmp(:,[1,2]);
        % fh_udg(:,[1,2],:) = fh_udg_tmp(:,[1,2],:);
        % fh_uh(:,[1,2],:) = fh_uh_tmp(:,[1,2],:);


        % NEW BC: INFLOW/OUTFLOW FOR SPECIES
        Er = udg(:,6);
        Ez = udg(:,9);
        E = [Er Ez];
        Edotn = dot(E,nl, 2);

        % For electrons, (-mue*E)dot(n)<0, or (dphi/dz dot n)<0 corresponds to inflow because the electrons flow in the opposite direction of the potential field.
        %
        % NOTE: the inflow/outflow convention doesn't make sense for the positive ions, which are stationary in this model. However, this will be important for future models that include ion drift.
        % The inflow/outflow will be of the opposite sign for positive ions

        % Determine inflow/outflow nodes. This is general and accounts for nodes on the same boundary being of different inflow/outflow types.
        nodes_e_in=find(-Edotn < 0);     % Electron inflow
        nodes_e_out = [find(-Edotn > 0); find(-Edotn == 0)];     % Electron outflow, captures edge case where Edotn == 0 (unlikely)

        % If E pointing out of the domain: electrons flow in -> inflow
        % Dirichlet condition, ne=ni= 1e13 (nondimensional=>10), same as background density
        fh(nodes_e_in,1) = r.*tau .*(ui(nodes_e_in,1)-uh(nodes_e_in,1));
        fh_uh(nodes_e_in,1,1) = -r.*tau;

        fh(nodes_e_out,2) = r.*tau .*(ui(nodes_e_out,2)-uh(nodes_e_out,2));  % The electron outflow nodes are inflow for positive ions. Again, this doesn't really apply for the streamer case, but will in future models.
        fh_uh(nodes_e_out,2,2) = -r.*tau;

        % If E pointing into the domain: electrons flow out -> outflow
        % Extrapolate, u=uhat
        fh(nodes_e_out,1) = r.*tau.*(udg(nodes_e_out,1)-uh(nodes_e_out,1));
        fh_udg(nodes_e_out,1,1) = r.*tau;
        fh_uh(nodes_e_out,1,1) = -r.*tau;

        % Opposite for positive ions
        fh(nodes_e_in,2) = r.*tau.*(udg(nodes_e_in,2)-uh(nodes_e_in,2));
        fh_udg(nodes_e_in,2,2) = r.*tau;
        fh_uh(nodes_e_in,2,2) = -r.*tau;

        % Potential: dirichlet (same as original)
        fh(:,3) = r.*tau .*(ui(:,3)-uh(:,3));
        fh_uh(:,3,3) = -r.*tau;

        % error;
        
    case 2  % Right "farfield"
        % Species + potential all have homogeneous neumann
        [fh,fh_udg,fh_uh] = fhat_axisymmetric(nl,p,udg,uh,param,time);

    case 3  % Symmetry boundary: extrapolate m=u or u=uhat
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);

        fh(:,1) = r.*tau.*(udg(:,1)-uh(:,1));
        fh_udg(:,1,1) = r.*tau;
        fh_uh(:,1,1) = -r.*tau;

        fh(:,2) = r.*tau.*(udg(:,2)-uh(:,2));
        fh_udg(:,2,2) = r.*tau;
        fh_uh(:,2,2) = -r.*tau;

        fh(:,3) = r.*tau.*(udg(:,3)-uh(:,3));
        fh_udg(:,3,3) = r.*tau;
        fh_uh(:,3,3) = -r.*tau;

    otherwise
        error('unknown boundary type');
end
