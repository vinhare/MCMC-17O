function ss = TripleOxgenss(theta,data)

sample = data.xdata(:,1)';
ydata  = data.ydata;
xdata  = data.xdata;
yerr   = data.yerr;

ymodel = BWM_Emu_st(sample,theta,xdata);
sigma2 = yerr.^2;
ss = sum(((ymodel - ydata).^2)./sigma2);
