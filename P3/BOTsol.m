function BOTsol;


n=50;

% Read data sensor 1
fileID=fopen('sensorA.txt','r');
yB=fscanf(fileID,'%f');
fclose(fileID);
% Read data sensor 2
fileID=fopen('sensorB.txt','r');
yA=fscanf(fileID,'%f');
fclose(fileID);

y=[yB' ; yA'];
figure(2); 
clf; subplot(1,2,1), plot(yB,'ok'); title('Sensor 1 angle'); xlabel('Time');
subplot(1,2,2), plot(yA,'ok'); title('Sensor 2 angle'); xlabel('Time');

A=eye(4);
delD=1/60;
A(1,3)=delD;
A(2,4)=delD;
% Initial pdfs
mprior=[10;30;15;-15];
sigprior=diag([5^2;5^2;3^2;3^2]);
SigPP=diag([0.1^2;0.1^2;0.7^2;0.7^2]);
Sigeps=diag([0.087^2;0.087^2]);
mprop=zeros(4,n);
mprop(1:4,1)=mprior;
mupd=zeros(4,n);
sigpred=zeros(4,4,n);
sigpred(1:4,1:4,1)=sigprior;
sigupd=zeros(4,4,n);
sserr=zeros(2,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EKF PART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recursive estimation
for i=1:n,
  
  
  % Linearize
  gm=zeros(2,4);
  gm(1,1:4)=(1/(1+(mprop(1,i)/mprop(2,i))^2))*[1/mprop(2,i) -mprop(1,i)/mprop(2,i)^2  0 0]; 
  gm(2,1:4)=(1/(1+((40-mprop(2,i))/(40-mprop(1,i)))^2))*[(40-mprop(2,i))/(40-mprop(1,i))^2 -1/(40-mprop(1,i)) 0 0]; 
  % Update
  Smat=gm*sigpred(1:4,1:4,i)*gm'+Sigeps;
  Kmat=sigpred(1:4,1:4,i)*gm'*inv(Smat);
  mupd(1:4,i)=mprop(1:4,i)+Kmat*(y(1:2,i)-[atan(mprop(1,i)/mprop(2,i));atan((40-mprop(2,i))/(40-mprop(1,i)))]);
  sigupd(1:4,1:4,i)=(eye(4)-Kmat*gm)*sigpred(1:4,1:4,i);
  % Sort some variables for plotting
  sserr(1,i)=sqrt(sigupd(1,1,i));
  sserr(2,i)=sqrt(sigupd(2,2,i));
  mekf(1:4,i)=mupd(1:4,i);
  % Predict (according to the reference system)
  if (i<n)
    mprop(1:4,i+1)=A*mupd(1:4,i);
    sigpred(1:4,1:4,i+1)=A*sigupd(1:4,1:4,i)*A'+SigPP;
  end;
  
  % Plotting
  figure(3)
  clf;
  sp=plot(mupd(1,i),mupd(2,i),'ko');
  hold;
  title('Extended Kalman filter');
  set(sp,'LineWidth',2);
  st=text(2500,2300,'ship');
  set(st,'FontSize',14);
  plot(mupd(1,1:i),mupd(2,1:i),'r');
  plot(mupd(1,1:i)+2*sserr(1,1:i),mupd(2,1:i)+2*sserr(2,1:i),'b-.');
  plot(mupd(1,1:i)-2*sserr(1,1:i),mupd(2,1:i)-2*sserr(2,1:i),'b-.');
  
  if (i>1)
      axis([-5 45 -5 45])
    pause(0.2);
    use=0;
  else
      axis([-5 45 -5 45])
    pause(1);
    use=0;
  end;  
end;

%% Particle Filter

% Sample initial ensembles
B=1000;
L=chol(sigprior)';
xB=zeros(4,B);
xB=mprior*ones(1,B)+L*randn(4,B);
xBupd=zeros(4,B);
mm=zeros(4,n);
mlow=zeros(4,n);
mhigh=zeros(4,n);

% Recursive estimation
for i=1:n,

  yPB=[atan(xB(1,:)./xB(2,:));atan((40-xB(2,:))./(40-xB(1,:)))];
  l=-0.5*inv(Sigeps(1,1))*(y(1:2,i)*ones(1,B)-yPB).*(y(1:2,i)*ones(1,B)-yPB); 
  w=sum(l,1);
  ew=exp(w);
  nw=ew/sum(ew);
  Fw=cumsum(nw);
  xBupd=zeros(4,B);
  for bb=1:B,
      Ur=rand;
      [iiS,indS]=find(Fw>Ur);
      ind=min(indS);
      xBupd(1:4,bb)=xB(1:4,ind);
  end;
  % Sort some variables for plotting
  mm(1:4,i)=mean(xBupd')';
  for ik=1:4,
      ms=sort(xBupd(ik,:));
      mlow(ik,i)=ms(1,ceil(0.05*B));
      mhigh(ik,i)=ms(1,floor(0.95*B));
  end
  % Predict (according to the reference system)
  if (i<n)
    xB=A*xBupd+chol(SigPP)'*randn(4,B);
  end;
  
  % Plotting
  figure(4)
  clf;
  sp=plot(mm(1,1),mm(2,1),'ko');
  hold;
  title('Particle filter');
  set(sp,'LineWidth',2);
  st=text(2500,2300,'ship');
  set(st,'FontSize',14);
  plot(mm(1,1:i),mm(2,1:i),'r');
  plot(mlow(1,1:i),mlow(2,1:i),'b-.');
  plot(mhigh(1,1:i),mhigh(2,1:i),'b-.');
  
  if (i>1)
      axis([-5 45 -5 45])
    pause(0.2);
    use=0;
  else
      axis([-5 45 -5 45])
    pause(1);
    use=0;
  end;  
end;


%% EnKF

% Sample initial ensembles
B=1000;
L=chol(sigprior)';
xB=zeros(4,B);
xB=mprior*ones(1,B)+L*randn(4,B);
xBupd=zeros(4,B);
mm=zeros(4,n);
mlow=zeros(4,n);
mhigh=zeros(4,n);

% Recursive estimation
for i=1:n,

  yPB=[atan(xB(1,:)./xB(2,:));atan((40-xB(2,:))./(40-xB(1,:)))]+chol(Sigeps)'*randn(2,B);
  Syhat=inv(B)*(yPB-mean(yPB')'*ones(1,B))*(yPB-mean(yPB')'*ones(1,B))';
  Sxyhat=inv(B)*(xB-mean(xB')'*ones(1,B))*(yPB-mean(yPB')'*ones(1,B))';
  Kmat=Sxyhat/Syhat;
  xBupd(1:4,:)=xB(1:4,:)+Kmat*(y(1:2,i)*ones(1,B)-yPB);
  % Sort some variables for plotting
  mm(1:4,i)=mean(xBupd')';
  for ik=1:4,
      ms=sort(xBupd(ik,:));
      mlow(ik,i)=ms(1,ceil(0.05*B));
      mhigh(ik,i)=ms(1,floor(0.95*B));
  end
  % Predict (according to the reference system)
  if (i<n)
    xB=A*xBupd+chol(SigPP)'*randn(4,B);
  end;
  
  % Plotting
  figure(5)
  clf;
  sp=plot(mm(1,1),mm(2,1),'ko');
  hold;
  title('Ensemble Kalman filter');
  set(sp,'LineWidth',2);
  st=text(2500,2300,'ship');
  set(st,'FontSize',14);
  plot(mm(1,1:i),mm(2,1:i),'r');
  plot(mlow(1,1:i),mlow(2,1:i),'b-.');
  plot(mhigh(1,1:i),mhigh(2,1:i),'b-.');
  
  if (i>1)
      axis([-5 45 -5 45])
    pause(0.2);
    use=0;
  else
      axis([-5 45 -5 45])
    pause(1);
    use=0;
  end;  
end;

