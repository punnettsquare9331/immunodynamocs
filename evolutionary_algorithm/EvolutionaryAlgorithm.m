function [sol]=EvolutionaryAlgorithm

% The copyright belongs to Kang-Ling Liao, Assistant Professor, 
% Departments of Mathematics and Biological Sciences, University of Manitoba.

clear all
close all
warning('off', 'all')


%1. Load the data 
% =================================================
load('noDox.mat');
load('Dox.mat');
tt=13; % tt is the total duration of Data 
tspan=2:0.25:tt;  % unit is min
 
%2. Set Parameters setting
% =================================================
%UCP2 Model 2
NOP=0;  % number of old para
NNP=27;  % number of new para
TPara=NOP+NNP;
VN=5; % variable number 
AllNum=[NOP,NNP,TPara,VN];

% %UCP2 Model 3
% NOP=0;  % number of old para
% NNP=18;  % number of new para
% TPara=NOP+NNP;
% VN=3; % variable number 
% AllNum=[NOP,NNP,TPara,VN];

M=10000; % Generation times
PN=100; % Parent number
CN=50; % Children number
Parents=zeros(PN,TPara);
Children=zeros(PN*CN,TPara);

mu=0.4;     % percentage of Parent to have macromutation
lambda=0.5; % max possible percent change of a parameter
beta=2;     % effective energy

AllParents=zeros(M,PN,TPara);

%3. Define the first generation
% =================================================
%UCP2 Model 2
Para=zeros(1,27); % parameter set 
Para(1)= (0.1+9.9*rand)*0.5436; %lambda_C
Para(2)= (0.1+9.9*rand)*0.9; %tildeC
Para(3)= (0.1+9.9*rand)*24; %eta_T
Para(4)= (0.1+9.9*rand)*0.17; %d_C
Para(5)= (0.1+9.9*rand)*8*10^(-6); %tildelambda_DC
Para(6)= (0.1+9.9*rand)*0.4; %K_C
Para(7)= (0.1+9.9*rand)*2.4328*10^(-5); %lambda_DU
Para(8)= (0.1+9.9*rand)*4.425*10^(-4); %K_U
Para(9)= (0.1+9.9*rand)*0.1; %d_D
Para(10)= (0.1+9.9*rand)*1.0016*10^(-6); %tildelambda_M_2C
Para(11)= (0.1+9.9*rand)*0.0069; %lambda_M_2Q
Para(12)= (0.1+9.9*rand)*5.8494*10^(-20); %K_Q
Para(13)= (0.1+9.9*rand)*0.015; %d_M_2
Para(14)= (0.1+9.9*rand)*48.6*10^(-4); %tildelambda_T
Para(15)= (0.1+9.9*rand)*4*10^-5; %K_D
Para(16)= (0.1+9.9*rand)*2.87*10^(-5); %barK_M_2
Para(17)= (0.1+9.9*rand)*2.9247*10^(-20); %barK_Q
Para(18)= (0.1+9.9*rand)*97.2; %lambda_TU
Para(19)= (0.1+9.9*rand)*0.18; %d_T
Para(20)= (0.1+9.9*rand)*0.075; %lambda_UC
Para(21)= (0.1+9.9*rand)*0.0075; %lambda_UM_2
Para(22)= (0.1+9.9*rand)*5.73*10^(-5); %K_M_2
Para(23)= (0.1+9.9*rand)*33.271; %d_U
Para(24)= (0.1+9.9*rand)*2.49*10^(-7); %rho_P
Para(25)= (0.1+9.9*rand)*5.22*10^(-7); %rho_L
Para(26)= (0.1+9.9*rand)*10^(-3); %epsilon
Para(27)= (0.1+9.9*rand)*2030.385; %lambda_LU

% %UCP2 Model 3
% Para=zeros(1,18); % parameter set 
% Para(1)=; %lambda_C
% Para(2)=; %tildeC
% Para(3)=; %eta_T
% Para(4)=; %d_C
% Para(5)=; %tildelambda_T
% Para(6)=; %K_C
% Para(7)=; %barK_Q
% Para(8)=; %lambda_TU
% Para(9)=; %K_U
% Para(10)=; %d_T
% Para(11)=; %lambda_UC
% Para(12)=; %lambda_UQ
% Para(13)=; %K_Q
% Para(14)=; %d_U
% Para(15)=; %rho_P
% Para(16)=; %rho_L
% Para(17)=; %epsilon
% Para(18)=; %lambda_LU

for pn=1:PN  % Give parameters
    %Parents(pn,1:NOP)=Para(1:NOP); %if old parameters 
    for Set=NOP+1:NOP+NNP          
        Parents(pn,Set)=1000*rand(1); % Initial guesses
    end    
end


% Give IC for (C, D, M2, T, U)

init_IC=[0.1, 1.6164*10^(-5), 5.73*10^(-6), 1.929*10^-4, 4.425*10^(-5)];
%4. Start Generation
% =================================================
for Gen=1:M
fprintf('Gen=%d \n', Gen)

  % Macromutation for the 1st Parents
    PopuIndex=[1:PN];
    SamInd=randsample(PopuIndex, mu*PN);
    for k=1:mu*PN
        MutUpRate=1.5;
        MutLowRate=0.5; 
        for Set=NOP+1:NOP+NNP 
            Parents(SamInd(k),Set)=rand(1)*(Parents(SamInd(k),Set)*MutUpRate-Parents(SamInd(k),Set)*MutLowRate)+Parents(SamInd(k),Set)*MutLowRate;
        end 
    end

    ChildrenIndex = zeros(1, PN*CN);
    ChilSco       = zeros(1, PN*CN);
    theta         = zeros(1, PN*CN);

    parfor idx = 1 : PN*CN          % *** parfor variable IS the row index ***
    	pn = ceil(idx / CN);        % parent number (1…PN)

    	% ---- build one child ------------------------------------------------
    	child = Parents(pn,:);      % start from its parent
    	for Set = NOP+1 : NOP+NNP
        	delta       = rand;
        	child(Set)  = child(Set) * (1 + lambda*(delta - 0.5));
    	end

    	% ---- evaluate & store ----------------------------------------------
    	C_sol = UCP2Model2Odes_min(AllNum, child, init_IC, tspan);  % model 2
    	% C_sol = UCP2Model3Odes_min(AllNum, child, init_IC, tspan);  % model 3

    	Children(idx,:)    = child;             % unique row = idx
    	ChildrenIndex(idx) = idx;
    	ChilSco(idx)       = C_sol;
    	theta(idx)         = exp(-beta * C_sol);
    end
    
    PareSco = zeros(1, PN);    
    % Score & reorder Parents 
    parfor pn=1:PN 
        P_sol = UCP2Model2Odes_min(AllNum,Parents(pn,:),init_IC,tspan); %UCP2 model 2
        %P_sol = UCP2Model3Odes_min(AllNum,Parents(pn,:),init_IC,tspan); %UCP2 model 3
        PareSco(pn)=P_sol;
    end
    [Order_PareSco,P_INDEX]=sort(PareSco);

    for pn=1:PN
        ReOderParents(pn,:)=Parents(P_INDEX(pn),:);     
    end

    % Select the Parents for next generation
     % 10% are chosen from the old Parents with top score
     for pn=1:PN*0.1 
        Parents(pn,:)=ReOderParents(pn,:);
        ParentsScore(pn)=Order_PareSco(pn);
     end
     
     % 90% are chosen from children by theta
     New_theta=theta;
     for pn=1:PN*0.9
         if any(New_theta~=zeros(size(New_theta))) %LINE CHANGED
            FindIndex(pn)=randsample(ChildrenIndex, 1, true, New_theta);
         else
            FindIndex(pn)=randsample(ChildrenIndex, 1, true);
         end
         New_theta(FindIndex(pn))=0; % remove from the current index
     end
      
     for pn=1:PN*0.9 
        Parents(pn+PN*0.1,:)=Children(FindIndex(pn),:);
        ParentsScore(pn+PN*0.1)=ChilSco(FindIndex(pn));
     end
    
     % Rescore the min score at each generation
     % ParentsScore
     MinScore(Gen)=min(ParentsScore)
     %pause
    
     save('EvoAlgDataUCP2Model2.mat','Gen','Parents','ParentsScore','MinScore') %UCP2 model 2
     %save('EvoAlgDataUCP2Model3.mat','Gen','Parents','ParentsScore','MinScore') % UCP2 model 3

end
save('EvoAlgDataUCP2Model2.mat','Gen','Parents','ParentsScore','MinScore') %UCP2 model 2
%save('EvoAlgDataUCP2Model3.mat','Gen','Parents','ParentsScore','MinScore') % UCP2 model 3












