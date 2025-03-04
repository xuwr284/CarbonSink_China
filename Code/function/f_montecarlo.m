function [fitAGB,q25,q75] = f_montecarlo(p,se,xall)
fitAGB=[];
frequency=10000;%iter; % 100000 
for ai=1:length(xall)
    agei=xall(ai);
    for i=1:frequency
        a=normrnd(p(1),se(1));
        b=normrnd(p(2),se(2));
        c=normrnd(p(3),se(3));
        tp=a*(1-exp(-1/c*agei)).^b;
        fitAGB(i,ai)=real(tp);
    end
end
fitAGB(fitAGB==Inf|fitAGB==-Inf)=nan;
q25=quantile(fitAGB,0.25);
q75=quantile(fitAGB,0.75);
end
