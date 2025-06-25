function AllError = UCP2Model3Odes_min(AllNum,ParaStart,IC,tspan)

% The copyright belongs to Kang-Ling Liao, Assistant Professor, 
% Departments of Mathematics and Biological Sciences, University of Manitoba.

NOP  = AllNum(1); % number of old para
NNP  = AllNum(2);  % number of new para
TPara= AllNum(3);
VN   = AllNum(4);

params = ParaStart;
load('noDox.mat');
load('Dox.mat');

% ======  Part I fitting =============
options = odeset('RelTol',1e-6);% SET THE TOLERENCE OF THE SOLVER
duration=tspan(end);               
sol = ode15s(UCP2Model3Odes_fit(params),tspan,IC,options);% SOLVE THE ODE

Error=zeros(1,VN);
AllError=0;  
New=zeros(VN,length(tspan));
for j=1:VN
    New(j,1:length(tspan))=spline(sol.x(:), sol.y(j,:),tspan);   
    Error(j)=0;
    for q=1:length(noDox(1,:))
        Error(j)=Error(j)+(New(j,noDox(1,q))-noDox(j+1,q))^2;
    end
    % for q=1:length(Dox(1,:))
    %     Error(j)=Error(j)+(New(j,Dox(1,q))-Dox(j+1,q))^2;
    % end  
    
    AllError=AllError+Error(j);
end

        
            
    


