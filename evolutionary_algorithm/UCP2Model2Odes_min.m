function AllError = UCP2Model2Odes_min(AllNum,ParaStart,IC,tspan)

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
options = odeset('RelTol',1e-8);% SET THE TOLERENCE OF THE SOLVER
duration=tspan(end);               
solDox = ode15s(UCP2Model2Odes_fit(params,1),tspan,IC,options); % SOLVE THE ODE with dox
IC(5) = 0;
solNoDox = ode15s(UCP2Model2Odes_fit(params,0),tspan,IC,options); % SOLVE THE ODE without dox

Error=zeros(1,VN);
AllError=0;  
NewDox=zeros(VN,length(tspan));
NewNoDox=zeros(VN,length(tspan));
for j=1:VN
    NewDox(j,1:length(tspan))=spline(solDox.x(:), solDox.y(j,:),tspan);   
    NewNoDox(j,1:length(tspan))=spline(solNoDox.x(:), solNoDox.y(j,:),tspan);  
    Error(j)=0;
    for q=1:length(noDox(1,:))
        Error(j)=Error(j)+(NewNoDox(j,noDox(1,q))-noDox(j+1,q))^2;
    end
    for q=1:length(Dox(1,:))
        Error(j)=Error(j)+(NewDox(j,Dox(1,q))-Dox(j+1,q))^2;
    end    
    
    AllError=AllError+Error(j);
end

        
            
    


