function SS = SS_over_rh(x,chain,data);

chain1=chain(1000:50:end,1); %thin the chains a bit 
chain2=chain(1000:50:end,2);
chain3=chain(1000:50:end,3);

[A B C] = Body_Water_Model_rh(x,median(chain1(:,1)), median(chain2(:,1)), median(chain3(:,1)));
ydata = data.ydata; 
yerr   = data.yerr;
sigma2 = yerr.^2;

ind = knnsearch(B,data.xdata(:,2));

E = B(ind);
D = A(ind);
SS = sum(((D - ydata).^2)./sigma2);
clear A B C chain1 chain2 chain3

