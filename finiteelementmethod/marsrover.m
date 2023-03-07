%% Pre-processor

%%%%%%% constants %%%%%%%%%
thickness = 0.01;
kTi = 17;
kGlass = 0.8;
rhoTi = 4620;
rhoGlass = 3860;
cTi = 523;
cGlass = 670;
EGlass = 67e9;
ETi = 110e9;
nyGlass = 0.2;
nyTi = 0.34;
alphaGlass = 7e-6;
alphaTi = 9.4e-6;
TinfInitial = 40;       % initial outer conditions for transient problem
Tinf = 40;              % outer temperature during transient problem
Tc = 20;                % temperature at the right boundary
Tstressfree = 20;       % temperature when structure stress free
alphaC = 100;

DTi = kTi * eye(2);     % temperature const. matrix titanium
DGlass = kGlass * eye(2); % temperature const. matrix glass

%mech const. matrix glass
DmGlass = EGlass/((1+nyGlass)*(1-2*nyGlass)) * [1-nyGlass nyGlass 0;
                                                nyGlass 1-nyGlass 0;
                                                0 0 1/2*(1-2*nyGlass)];
                                            
%mech const. matrix titanium

DmTi = ETi/((1+nyTi)*(1-2*nyTi)) * [1-nyTi nyTi 0;
                                    nyTi 1-nyTi 0;
                                    0 0 1/2*(1-2*nyTi)];
eq = 0; %load

ptype = 2; %plain strain


%%%%%%%%%%%% geometry %%%%%%%%%

coord = p' / 100; % divide by 100 -> unit [m]

enod=t(1:3,:)';                     % nodes of elements
nelm=size(enod,1);                  % number of elements
nnod=size(coord,1);                 % number of nodes
dof=(1:nnod)';                      % dof number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number

edof = zeros(nelm, 4); %edof matrix for temperature
edof_S = zeros(nelm, 7); %edof matrix for mechanical

for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:)=[ie,enod(ie,:)];
end

% Check which segments that should have convections and displacement boundary conditions


er = e([1 2 5],:); % Reduced e
conv_segments = [1 11 22 31 8 9 10]; % Choosen boundary segments for convection
static_segments = [20 21 8 9 10 17 30 18 19]; % Boundary segments for displacement conditions
edges_stat = [];
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];
    end
    if ismember(er(3,i), static_segments)
        edges_stat = [edges_stat er(1:2,i)];
    end
end

[ex,ey]=coordxtr(edof,coord,dof,3);


%%%%%%%%Element calculations, assembly and solutions to lin. systems and stresses

Kt = zeros(nnod);       %Temperature K matrix
C = zeros(nnod);        %Transient C matrix
Km = zeros(nnod * 2);   %Mech K matrix

% loop through all elements and calculate element matrices -> assemble
                                 
for elnr = 1:nelm
    if t(4, elnr) == 1 %Check if in titanium
        Kte = flw2te(ex(elnr, :), ey(elnr,:), thickness, DTi, eq);
        Ce = plantml(ex(elnr, :), ey(elnr,:), rhoTi * cTi * thickness);
        Kme = plante(ex(elnr, :), ey(elnr,:), [ptype thickness], DmTi);
    else                %in glass
        Kte = flw2te(ex(elnr, :), ey(elnr,:), thickness, DGlass, eq);
        Ce = plantml(ex(elnr, :), ey(elnr,:), rhoGlass * cGlass * thickness);
        Kme = plante(ex(elnr, :), ey(elnr,:), [ptype thickness], DmGlass);
    end
    
   indx = edof(elnr,2:end);
   indx_S = edof_S(elnr,2:end);
   Kt(indx,indx) = Kt(indx,indx)+Kte;
   C(indx, indx) = C(indx, indx) + Ce;
   Km(indx_S, indx_S) = Km(indx_S,indx_S) + Kme;
   
end

Kc = zeros(nnod);           % K matrix from convection
f = zeros(nnod, 1);         % f vector from convection with Tinf
fInitial = zeros(nnod, 1);  % f vector from convection with Tinf_initial

%go through all node pairs along convection border to calc Kc and f_conv

for pairnr = 1:length(edges_conv)
    nodes = edges_conv(:, pairnr);
    coord1 = coord(nodes(1), :);
    coord2 = coord(nodes(2), :);
    L = norm(coord1 - coord2); %calc length between nodes
    
    Kce = L/6 * [2 1; 1 2]; % element stiffness matrix
    fe = L/2 * [1; 1];      % element f vector
    
    if coord1(1) > max(coord(:,1)) / 2 %if at right side we use Tc
        feInitial = fe * Tc;
        fe = fe * Tc;
       
    else                            % else we use Tinf/Tinfinitial
        feInitial = fe * TinfInitial;
        fe = fe * Tinf;
    end
    
    Kc(nodes,nodes) = Kc(nodes,nodes) + Kce;
    f(nodes) = f(nodes) + fe;
    fInitial(nodes) = fInitial(nodes) + feInitial;
    
end

% multiply by thickness and alphaC
fInitial = alphaC * thickness * fInitial;
f = alphaC * thickness * f;
Kc = alphaC * thickness * Kc;


[Tstart, Qv] = solveq(Kt + Kc, fInitial); %solve system
eT = extract(edof, Tstart); %element temps

f0 = zeros(nnod*2,1); %f vector for mech. problem from temperature distribution

%loop through all elements, calc all element vectors f0e and assemble

for elnr = 1:nelm
    deltaT = mean(eT(elnr,:)) - Tstressfree; %for deltaT, use mean of the temperatures in the nodes at that element
    if t(4, elnr) == 1              %if in titanium
        epsilon0 = (1+nyTi) * alphaTi * deltaT * [1 1 0]';
        f0e = plantf(ex(elnr,:), ey(elnr, :), [ptype thickness], (DmTi * epsilon0)');
    else                            %in glass
        epsilon0 = (1+nyGlass) * alphaGlass * deltaT * [1 1 0]';
        f0e = plantf(ex(elnr,:), ey(elnr, :), [ptype thickness], (DmGlass * epsilon0)');
    end 
    
    index_S = edof_S(elnr, 2:end);
    f0(index_S) = f0(index_S) + f0e;
end


%find out what node numbers have displacement dirichlet conditions

bc = [];


for edge = 1:length(edges_stat(1,:))
    nodes = edges_stat(:, edge);
    if abs(coord(nodes(1),1) - coord(nodes(2),1)) < 1e-9 %if at a vertical border
        bc = [bc; nodes(1) 0; nodes(2) 0];               %ux = 0
    else                                                %else at a horizontal border
        bc = [bc; dof_S(nodes(1),2) 0; dof_S(nodes(2),2) 0]; %uy = 0
   
    end
end

bc = unique(bc,'rows'); %get rid of duplicates

[u, tb] = solveq(Km, f0, bc); %solve system

eu = extract(edof_S, u); %element displacements



stress = zeros(nelm, 4); % [sigx, sigy, sigz, tauxy] for all elements


%loop through all elements to calc element stresses
for elnr = 1:nelm
    deltaT = mean(eT(elnr,:)) - Tstressfree;
    if t(4, elnr) == 1 %in titanium
        %es is the element stress incl. temperature
        [es,~] = plants(ex(elnr,:), ey(elnr,:), [ptype thickness], DmTi, eu(elnr,:));
        tempstress = DmTi * (1+nyTi) * alphaTi * deltaT * [1 1 0]'; %stress from temp (D*eps0)
        sigma = es - tempstress'; %subtract temp. stress
        sigmazz = nyTi * (sigma(1) + sigma(2)) - alphaTi * ETi * deltaT;
    else                %in glass
        [es,~] = plants(ex(elnr,:), ey(elnr,:), [ptype thickness], DmGlass, eu(elnr,:));
        tempstress =  DmGlass * (1+nyGlass) * alphaGlass * deltaT * [1 1 0]';
        sigma = es - tempstress';
        sigmazz = nyGlass * (sigma(1) + sigma(2)) - alphaGlass * EGlass * deltaT;
    end
    
    
    stress(elnr, [1 2 4]) = sigma; 
    stress(elnr, 3) = sigmazz;
end

%calculate the von Mises stress at all elements
vonMisesEl = sqrt(stress(:,1).^2 + stress(:,2).^2 + stress(:,3).^2 - ...
            stress(:,1) .* stress(:,2) - stress(:,1) .* stress(:,3) - ...
            stress(:,2) .* stress(:,3) + 3*stress(:,4).^2);

vonMisesNode = zeros(nnod,1);

%Calculate von Mises at each node (average of connected elements)
for i=1:size(coord,1)
    [c0,c1]=find(edof(:,2:4)==i);
    vonMisesNode(i,1)=sum(vonMisesEl(c0))/size(c0,1);
end

%% calc dispmag

T = zeros(nnod * 2); % allocate T matrix

for elnr = 1:nelm
    if t(4, elnr) == 3                                      %check if in lens
        Te = plantmlmod(ex(elnr,:), ey(elnr,:), thickness); %calculate int(N^T*N) with modified plantml
        index = edof_S(elnr, 2:end);
        T(index,index) = T(index,index) + Te;               %insert
    end 
end

dispMag = u' * T * u;

%% plot disp

% Calculate displaced coordinates
mag = 100; % Magnification (due to small deformations)
exd = ex + mag*eu(:,1:2:end);
eyd = ey + mag*eu(:,2:2:end);
figure()
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(ex',-ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3) % plot the bottom half
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3)
patch(exd',-eyd',[0 0 0],'FaceAlpha',0.3)
axis equal
title('Displacement field [Magnitude enhancement 100]')

%% plot stress
figure()
eStress = extract(edof,vonMisesNode);
patch(ex',ey',eStress','EdgeColor','none')
hold on
patch(ex',-ey',eStress', 'EdgeColor','none')
colormap(jet);
caxis([0.6 14]*1e7)
c = colorbar;
c.Label.String = '[N/m^2]';
title('Effective von Mises stressfield [N/m^2]')
xlabel('x-position [m]')
ylabel('y-position [m]')

%% time
totaltime = 500;
dt = 5;
nofsteps = totaltime/dt;

Tresult = zeros(nnod, nofsteps + 1);
Tresult(:, 1) = Tstart;

for step = 1:nofsteps
    Tnext = (C + dt * (Kt + Kc))\(C * Tresult(:, step) + dt * f); %solve lin. eq. system for each time step
    Tresult(:, step+1) = Tnext; %save result in matrix
end

figure()

%plot snapshots
for k = 4:-1:0   
    eT = extract(edof, Tresult(:,1 + k*6));
    subplot(2,3,k+1)
    patch(ex',ey',eT','EdgeColor','none')
    hold on
    patch(ex',-ey',eT','EdgeColor','none')
    title(['At ' num2str(dt*(k*6)) ' s'])
    colormap(hot);
    c = colorbar;
    c.Label.String = '[C]';
    xlabel('x-position [m]')
    ylabel('y-position [m]')

end

%plot at stationary
subplot(2,3,6)
eT = extract(edof,Tresult(:,end));
patch(ex',ey',eT','EdgeColor','none')
hold on
patch(ex',-ey',eT','EdgeColor','none')
title(['At ' num2str(totaltime) ' s'])
colormap(hot);
c = colorbar;
c.Label.String = '[C]';
xlabel('x-position [m]')
ylabel('y-position [m]')

suptitle('Temperature distribution evolution [C]')

    
    