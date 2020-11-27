%Author(s) : Hrithik Agarwal, Hrishikesh Kotwal
%This code solve Q1 of Group Assignment 1 of ME608 course

m=2; % It is in form of CmHn, where C is carbon and H is hydrogen
n=6;
Af=4.338*(10^8); % constant from mass flow rate of fuel
Et=15098; % constant from mass flow rate of fuel (activation energy/universal gas constant)
x=0.1; % constant from mass flow rate of fuel
y=1.65; % constant from mass flow rate of fuel
del_H=4*(10^7); % Heat of combustion
cp=0.122*(10^4); % specific heat at constant pressure
phi=0.6; % equivalence ratio
dif=0.8; % percentage of reaction is complete
m_r=0.2; % mass flow rate
dia=0.03; % diameter
area=pi*dia*dia/4; % area
astio=(m+(n/4)); % stoichiometric A/F
mf=(12*m)+n; % mass of fuel
mo=32; % mass of o2
mpr= (m*44 + n*18/2)/(m+n/2); % mass of product i.e mCO2 + n/2 H2O
mn=28; % mass of N2
mh2o=18; % mass of H2O
ru =8315; % universal gas constant
P=101325 ; % Pressure in Pascal
m_u=mo*astio/mf; %


% initial conditions or at stating of flow x=0

mmixi= mf + (astio/phi)*(mo+3.76*mn); % mass of mixture
yfi=(mf)/mmixi; % mass fraction of fuel
yoi=(astio*mo/phi)/mmixi; % mass fraction of oxygen
ypi=0; % mass fraction of product
yni=1-yfi-yoi-ypi; % mass fraction of nitrogen

% conditon at total completion final conditions

mmixl=dif*((m+n/2)*mpr + (astio/phi-(m+n/4))*mo + (3.76*astio*mn/phi)) + (1-dif)*(mf + astio* (mo+3.76*mn)/phi ); % mass of mixture

yfl=(1-dif)*(mf)/mmixl;                               % mass fraction of fuel
yol=(astio/phi-dif*(m+n/2))*mo/mmixl;                  % mass fraction of oxygen
ypl=dif*(m+n/2)*mpr/mmixl;                              % mass fraction of product
ynl=1-yfl-yol-ypl;                                      % mass fraction of nitrogen

% initial guess / place takers

yf=yfi;           
yo=yoi;
yp=ypi;
yn=yni;


X=0;                   % length
X_prev=0;              % for previous iteration
del_x=10^-5;           % del X
Ti=1000;                % inital Temprature
T=1000;                  % initial guess/ place taker
i=0;                     
while i==0
  mmix=1/(yf/mf + yo/mo + yp/mpr +yn/mn);       % mass of mixture
  rho=P*mmix/(ru*T);                              % density 
  mftd=Af*(rho^(x+y))*exp(-Et/T)*(yf^x)*(yo^y);        % euation for mass flow rate of fuel 

  Tn=T+del_x*(area/m_r)*(del_H/cp)*mftd;               % Temprature equation 
  yf=cp*(Ti-Tn)/del_H + yfi;                           % mass fraction of fuel        
  yo=-(yfi-yf)*m_u + yoi;                               % mass fraction of oxygen  
  yp=(yfi-yf)*(m_u+1) + ypi;                           % mass fraction of product
  yn=1-yf-yo-yp;                                       % mass fraction of nitrogen
  X=X+del_x;                                           % Distance update 
  T=Tn;                                                 % Temprature update
  if yp > ypl
      L = X_prev+((ypl-yp_prev)*(X-X_prev)/(yp-yp_prev));       % to get anser using linear interpolation
    break;
  end
  if yp < 0
    print(-1);
    break;
  end
  yp_prev=(yfi-yf)*(m_u+1) + ypi;                          % updating product from previous iteration
  X_prev=X_prev+del_x;                                      % updating X for previous iteration
end

%disp(X)
disp(L)                                                       %display length
%disp(T)
%disp (mmixi)
%disp (mmixl)