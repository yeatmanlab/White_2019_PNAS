%% function dprime = AgToDprime(Ag) 
% Converts between area under the ROC curve (Ag) to dprime. 
% Treats Ag like unbiased p(correct) in a 2AFC task. 
% by Alex L. White at the University of Washington 
function dprime = AgToDprime(Ag) 

dprime = sqrt(2)*norminv(Ag);


%% What's the relationship betwen dprime and area under the ROC curve? 
% %First let's test out computing area under the curve, Ag, for a range of d'
% %values. 
% % 
% %classic SDT model, assuming a value of d', the separation of two distributions with unit variance 
% %values of dprime: 
% ds =  0:0.1:5; 
% 
% %then we vary criterion and compute hit and false alarm rates along the way
% cs = 10:(-0.01):-10; 
% 
% %false alarm rates (don't depend on d')
% frs = 1-normcdf(cs,0,1); 
% 
% %Compute Ag corresponding to each value of d'
% Ags = NaN(1,length(ds));
% for di = 1:length(ds)
%     %hit rates given this value of dprime   
%     hrs = 1-normcdf(cs,ds(di),1);
%     %the compute area under the curve by summing up area of the trapezoids:
%     Ags(di) = computeAROC(hrs,frs);
% end
% 
% %Now, from basic SDT math, we have a way to convert from p(corr) to
% %d-prime, assuming neutral criterion. 
% %The formula is:
% %   d' = sqrt(2)*z(pc)
% % and  applies for 2-interval forced-choice experiments. The sqrt(2) factor basically corrects for the fact that 
% % the observer has 2x as much information. z(x) is the same as matlab's
% % norminv(x), inverse of normCDF. 
% % 
% %Can we treat Ag (area under curve) just like 2AFC p(corr), to translate from 
% % from Ag to dprime? 
% predDs = sqrt(2)*norminv(Ags); %convert to d' (assuming Ag is like p(corr))
% predAgs = normcdf(ds/sqrt(2)); %convert dprime to Ag also
% 
% figure;  hold on;
% plot(ds,Ags,'b-','LineWidth',4);
% plot(predDs,predAgs,'g-','LineWidth',2);
% 
% xlabel('d-prime'); 
% ylabel('Ag');

%So, yes, we can summarize the relationship between Ag (area under ROC
%curve) and dprime as follows: 
%% dprime = sqrt(2)*norminv(Ag)
%% Ag = normcdf(dprime/sqrt(2))