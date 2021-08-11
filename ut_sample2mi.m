function mi=ut_sample2mi(data)
% computes mutual information mi between a binary random variable
% with symmetric PDF f(x) given by the samples <data>
% (histogram with bin centers <c> and <h> samples per bin)
% mi=integral f(x) log_2 (2*f(x)/(f(x)+f(-x))) dx
c=linspace(-20,+20,301); h=hist(data,c);
mi=trapz(c,h.*log(2*h./(h+fliplr(h)+1e-30)+1e-30))/trapz(c,h)/log(2);