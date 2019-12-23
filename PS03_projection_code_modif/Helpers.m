classdef Helpers


%--------------------------------------------------------------------------
% Static methods (miscellaenous functions)
%--------------------------------------------------------------------------
   methods (Static)
       
       function f=logopen(logfname,perm)
          if nargin<2
              perm='a';
          end
          if isempty(logfname)
                f=1;
            else
                f=fopen(logfname,perm);
          end 
       end
       
       function logclose(f)
          if f>1
              fclose(f);
          end
       end
       
       function st=vecToStruct(vec,names)
           st=num2cell(vec);
           st=cell2struct(st,names,1);
       end
       
       function vec=structToVec(st)
            vec=struct2cell(st);
            vec=cell2mat(vec);
       end
       
       function [P_Rouw, z_Rouw] = rouwen(rho_Rouw, mu_uncond, sig_uncond, S_Rouw, n_R)
           %ROUWEN   Rouwenhorst's method (1995) to approximate an AR(1) process using
           %   a  finite state Markov process.
           %
           %   For details, see Rouwenhorst, G., 1995: Asset pricing  implications of
           %   equilibrium business cycle models, in Thomas Cooley (ed.), Frontiers of
           %   Business Cycle Research, Princeton University Press, Princeton, NJ.
           %
           %   Suppose we need to approximate the following AR(1) process:
           %
           %                   y'=rho_Rouw*y+e
           %
           %   where abs(rho_Rouw)<1, sig_uncond=std(e)/sqrt(1-rho_Rouw^2) and
           %   mu_uncond denotes E(y), the unconditional mean of y. Let n_R be the
           %   number of grid points. n_R must be a positive integer greater than one.
           %
           %   [P_Rouw, z_Rouw] = rouwen(rho_Rouw, mu_uncond, sig_uncond, n_R) returns
           %   the discrete state space of n_R grid points for y, z_Rouw, and
           %   the centrosymmetric transition matrix P_Rouw. Note that
           %
           %       1. z_Rouw is a column vector of n_R real numbers.
           %       2. The (i,j)-th element of P_Rouw is the conditional probability
           %          Prob(y'=z_Rouw(i)|y=z_Rouw(j)), i.e.
           %
           %                 P_Rouw(i,j)=Prob(y'=z_Rouw(i)|y=z_Rouw(j))
           %
           %           where z_i is the i-th element of vector z_Rouw. Therefore
           %
           %           P_Rouw(1,j)+P_Rouw(2,j)+ ... +P_Rouw(n,j)=1 for all j.
           %
           %   See also HITM_Z and HITM_S on how to simulate a Markov processes using
           %   a transition matrix and the grids.
           %
           %   Damba Lkhagvasuren, June 2005
           
           % CHECK IF abs(rho)<=1
           if abs(rho_Rouw)>1
               error('The persistence parameter, rho, must be less than one in absolute value.');
           end
           
           % CHECK IF n_R IS AN INTEGER GREATER THAN ONE.
           if n_R <1.50001 %| mod(n_R,1)~=0
               error('For the method to work, the number of grid points (n_R) must be an integer greater than one.');
           end
           
           % CHECK IF n_R IS AN INTEGER.
           if mod(n_R,1)~=0
               warning('the number of the grid points passed to ROUWEN is not an integer. The method rounded n_R to its nearest integer.')
               n_R=round(n_R);
               disp('n_R=');
               disp(n_R);
           end
           
           skewfun=@(p)(S_Rouw - (2*p-(1+rho_Rouw))/sqrt((n_R-1)*(1-p)*(p-rho_Rouw)));
           options=optimset('Display','none');
           [p,~,exit]=fsolve(skewfun,(rho_Rouw+1)/2,options);
           if exit<=0
              error('Could not solve for p given skewness'); 
           else
              q=rho_Rouw + 1 - p;
              if q<0 || q>1
                  error('Level of skewness not attainable given number of nodes'); 
              end
           end

           % GRIDS
           step_R = sig_uncond*sqrt((n_R-1)*(2-p-q)^2/(4*(1-p)*(1-q)));
           M=mu_uncond-step_R*(q-p)/(2-p-q);
           z_Rouw=[-1:2/(n_R-1):1]';
           z_Rouw=M+step_R*z_Rouw;
           
           % CONSTRUCTION OF THE TRANSITION PROBABILITY MATRIX           
           P_Rouw=[ p  (1-p);
               (1-q) q];
           
           for i_R=2:n_R-1
               a1R=[P_Rouw zeros(i_R, 1); zeros(1, i_R+1)];
               a2R=[zeros(i_R, 1) P_Rouw; zeros(1, i_R+1)];
               a3R=[zeros(1,i_R+1); P_Rouw zeros(i_R,1)];
               a4R=[zeros(1,i_R+1); zeros(i_R,1) P_Rouw];
               P_Rouw=p*a1R+(1-p)*a2R+(1-q)*a3R+q*a4R;
               P_Rouw(2:i_R, :) = P_Rouw(2:i_R, :)/2;
           end
           
           P_Rouw=P_Rouw';
           
           for i_R = 1:n_R
               P_Rouw(:,i_R) = P_Rouw(:,i_R)/sum(P_Rouw(:,i_R));
           end
       end
           
       
       function des = hitm_s(PPP,a1)
           %HITM_S simulates a time series for states (not actual values) of a finite
           %   state Markov chain using its transition matrix, PPP, and a column
           %   vector of N random numbers drawn from the standard uniform distribution
           %   on the interval(0,1), a1.
           %
           %   The length of the simulated time series is given by the size of a1,
           %   i.e. if a1 is an Nx1 column vector of the random numbers, the size of
           %   des will be Nx1.
           %
           %                   des = hitm_s(PPP,a1)
           %
           %   where
           %
           %       PPP - the transition matrix,
           %       a1 - a column vector of random numbers drawn from the standard
           %            uniform distribution on the interval(0,1), and
           %       des - the simulated time series for the states. If PPP is an MxM
           %             matrix, each period des will take one of the following M
           %             integers, {1,2,...,M}.
           %
           %   Note that the method assumes that (i,j)-th element of PPP is the
           %   conditional probability Prob(state(t+1)=i|state(t)=j), i.e.
           %
           %                 PPP(i,j)=Prob(state(t+1)=i|state(t)=j)
           %
           %   where i and j denote the numbers of the states. Therefore
           %
           %           PPP(1,j)+PPP(2,j)+ ... +PPP(n,j)=1, for all j.
           %
           %   See also HITM_Z for simulating actual values.
           %   Damba Lkhagvasuren, June 2005
           
           N=size(a1,1);
           znum=size(PPP,1);
           
           A=tril(ones(znum))*PPP;
           A(znum,:)=2;
           
           des=zeros(N,1);
           
           ainit=randperm(znum);
           des(1,1)=ainit(1,1);
           destemp=des(1,1);
           
           for c_ount=2:N;
               
               if a1(c_ount,1)<=A(1,destemp);
                   des(c_ount,1)=1;
               end ;
               
               for i=1:znum-1;
                   if A(i,destemp)<a1(c_ount,1)
                       if A(i+1,destemp)>=a1(c_ount,1);
                           des(c_ount,1)=i+1;
                       end
                   end
               end
               destemp=des(c_ount,1);
               
           end;
       end
       
       function [Z,Zprob] = tauchenhussey(N,mu,rho,sigma,baseSigma)
           % Function tauchenhussey
           %
           % Purpose:    Finds a Markov chain whose sample paths
           %             approximate those of the AR(1) process
           %                 z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
           %             where eps are normal with stddev sigma
           %
           % Format:     {Z, Zprob} = TauchenHussey(N,mu,rho,sigma,m)
           %
           % Input:      N         scalar, number of nodes for Z
           %             mu        scalar, unconditional mean of process
           %             rho       scalar
           %             sigma     scalar, std. dev. of epsilons
           %             baseSigma scalar, std. dev. used to calculate Gaussian
           %                       quadrature weights and nodes, i.e. to build the
           %                       grid. I recommend that you use baseSigma = w*sigma +
           %                       (1-w)*sigmaZ where sigmaZ = sigma/sqrt(1-rho^2),
           %                       and w = 0.5 + rho/4. Tauchen & Hussey recommend
           %                       baseSigma = sigma, and also mention baseSigma = sigmaZ.
           %
           % Output:     Z       N*1 vector, nodes for Z
           %             Zprob   N*N matrix, transition probabilities
           %
           %     Martin Floden, Stockholm School of Economics
           %     January 2007 (updated August 2007)
           %
           %     This procedure is an implementation of Tauchen and Hussey's
           %     algorithm, Econometrica (1991, Vol. 59(2), pp. 371-396)
           
           Z     = zeros(N,1);
           Zprob = zeros(N,N);
           
           [Z,w] = gaussnorm(N,mu,baseSigma^2);   % See note 1 below
           
           
           for i = 1:N
               for j = 1:N
                   EZprime    = (1-rho)*mu + rho*Z(i);
                   Zprob(i,j) = w(j) * norm_pdf(Z(j),EZprime,sigma^2) / norm_pdf(Z(j),mu,baseSigma^2);
               end
           end
           
           for i = 1:N
               Zprob(i,:) = Zprob(i,:) / sum(Zprob(i,:),2);
           end
           
           
           
           function c = norm_pdf(x,mu,s2)
               c = 1/sqrt(2*pi*s2) * exp(-(x-mu)^2/2/s2);
           end
           
           function [x,w] = gaussnorm(n,mu,s2)
               % Find Gaussian nodes and weights for the normal distribution
               % n  = # nodes
               % mu = mean
               % s2 = variance
               
               [x0,w0] = gausshermite(n);
               x = x0*sqrt(2*s2) + mu;
               w = w0 / sqrt(pi);
           end
           
           function [x,w] = gausshermite(n)
               % Gauss Hermite nodes and weights following "Numerical Recipes for C"
               
               MAXIT = 10;
               EPS   = 3e-14;
               PIM4  = 0.7511255444649425;
               
               x = zeros(n,1);
               w = zeros(n,1);
               
               m = floor(n+1)/2;
               for i=1:m
                   if i == 1
                       z = sqrt((2*n+1)-1.85575*(2*n+1)^(-0.16667));
                   elseif i == 2
                       z = z - 1.14*(n^0.426)/z;
                   elseif i == 3
                       z = 1.86*z - 0.86*x(1);
                   elseif i == 4
                       z = 1.91*z - 0.91*x(2);
                   else
                       z = 2*z - x(i-2);
                   end
                   
                   for iter = 1:MAXIT
                       p1 = PIM4;
                       p2 = 0;
                       for j=1:n
                           p3 = p2;
                           p2 = p1;
                           p1 = z*sqrt(2/j)*p2 - sqrt((j-1)/j)*p3;
                       end
                       pp = sqrt(2*n)*p2;
                       z1 = z;
                       z = z1 - p1/pp;
                       if abs(z-z1) <= EPS, break, end
                   end
                   if iter>MAXIT, error('too many iterations'), end
                   x(i)     = z;
                   x(n+1-i) = -z;
                   w(i)     = 2/pp/pp;
                   w(n+1-i) = w(i);
               end
               x(:) = x(end:-1:1);
           end
           
           % Note 1: If you have Miranda and Fackler's CompEcon toolbox you can use
           % their qnwnorm function to obtain quadrature nodes and weights for the
           % normal function: [Z,w] = qnwnorm(N,mu,baseSigma^2);
           % Compecon is available at http://www4.ncsu.edu/~pfackler/compecon/
           % Otherwise, use gaussnorm as here.
           
       end
   

       
       function [stlist,transmat,indlist]=makeTotalTransition(pointList,probList)
           
           N=1;
           NI=[];
           ndim=length(pointList);
           
           % total size of state space
           for i=1:ndim
               ni=length(pointList{i});
               N=N*ni;
               NI=[NI,ni];
           end
           
           indlist=grid.StateSpaceGrid.makeCombinations(NI);
           
           stlist=zeros(N,ndim);
           for d=1:ndim
               stlist(:,d)=pointList{d}(indlist(:,d));
           end
           
           transmat=zeros(N);
           for i=1:N
               dind=indlist(i,:);
               for j=1:N
                   tind=indlist(j,:);
                   tp=1;
                   for d=1:ndim
                       tp=tp*probList{d}(dind(d),tind(d));
                   end
                   transmat(i,j)=tp;
               end
           end                      
       end
       
       
        function bigstr=combineStructs(strArray)
            % strArray must be cell array of structs
            
            bigstr=struct;
            Ninstr=length(strArray);
            for i=1:Ninstr
                strtmp=strArray{i};
                fntmp=fieldnames(strtmp);
                for j=1:length(fntmp)
                    thisfield=fntmp{j};
                    bigstr.(thisfield)=strtmp.(thisfield);
                end
            end
        end
        
        function [indlist,effnames]=makeListFromNames(hashmap,names)            
            n=length(names);
            indlist=zeros(n,1);
            for i=1:n
                tmp=hashmap.get(names{i});
                if ~isempty(tmp)
                    indlist(i)=tmp;
                else
                    indlist(i)=0;
                end
            end
            effnames=names(indlist~=0);
            indlist=indlist(indlist~=0);            
        end
        
        function dum=tableExport(filename, varnames, data)
            
            dum=[];
            
            ncol=length(varnames);
            fid=fopen(filename,'w');
            for c=1:ncol-1
                fprintf(fid,'%s,',varnames{c});
            end
            fprintf(fid,'%s\n',varnames{ncol});
            fclose(fid);
            
            dlmwrite(filename,data,'-append');
            
        end
              
        function dum=scatterPoints2D(X,Y)
            dum=[];
            
            bY=max(Y);
            scale=100/bY;
            figure; hold on;
            scatter(X(:,1),X(:,2),Y*scale);            
        end
       
        
        function dum=makeHist2D(obs,freq,edges,title_str,color)
            
            dum=[];
            ncol=length(obs);
            %totfreq=zeros(length(edges),ncol);
            width=0.6./(1:ncol);
            hold on;
            for i=1:ncol
                this_obs=obs{i};
                [cnts,bins]=histc(this_obs,edges);
                if ~isempty(freq)
                    this_freq=freq{i};
                    totfreq=accumarray(bins,this_freq,[length(edges),1]);
                else
                    totfreq=cnts;
                end
                bar(edges,totfreq,width(i),color{i});
            end
            title(title_str);
            
        end
        
        
        function dum=makeImpRes(simseries,tvec,titles,colors,format,lbl,printfile, for_slides, width_ratio)
            
            % Adding font size option
            if nargin < 9
                width_ratio = 1.0;
                if nargin < 8
                    for_slides = 0;
                end
            end
            
            dum=[];
            sims=size(simseries);
            if length(sims)==2
                npp=1;
                npl=sims(2);
                nper=sims(1);
                simseries=reshape(simseries,[1,size(simseries)]);
            else
                npp=sims(1);
                npl=sims(3);
                nper=sims(2);
            end
            if npl~=prod(format)
                disp('format does not match number of plots');
                return;
            end
            fig = figure;
            for p=1:npl
                h=subplot(format(1),format(2),p);
                hold on;
                for pp=1:npp
                    plot(tvec,squeeze(simseries(pp,:,p)),colors{pp},'LineWidth',1.5);
                end
                
                if for_slides
                    set(gca, 'FontSize', 16);
                else
                    set(gca, 'FontSize', 14);
                end
                
                title(h,titles{p});
                set(h,'XLim',[min(tvec),max(tvec)],'XTick',min(tvec):round((tvec(end)-tvec(1))/5):max(tvec),'XGrid','on','GridLineStyle','-');
            end
            if ~isempty(lbl)
                l=legend(lbl{:},'Location','best');
                l.FontSize=14;
            end
            
            set(gcf,'PaperUnits','inches');
            set(gcf,'PaperSize',[width_ratio*4*format(2) 4*format(1)]);
            set(gcf,'PaperPosition',[0 0 width_ratio*4*format(2) 4*format(1)]);
            set(gcf,'PaperPositionMode','manual');
            %set(gcf,'PaperSize',[10*format(1), 3*format(2)]);
            if ~isempty(printfile)
                print('-dpdf','-r0',[printfile,'.pdf']);
                print('-depsc',[printfile,'.eps']);
            end
            
            
        end
        
       
   end
%--------------------------------------------------------------------------
% End static methods 
%--------------------------------------------------------------------------
    
    
    
end