function [c,V,r]=Spherelets(X,d)
% Find the best  d-dimensional sphere to fit data X
% input: X=data matrix
% output: c=center of the spherelet
%              r=radius of the spherelet
%              V=the subspace where the sphere lies in the affine space c+V
%              SS=sum of square error
[n,m]=size(X); % n=sample size

if n>d+1   % if there are enough samples, fit the data by a sphere or a hyperplane
    
    % do d+1 dimensional PCA first
    [V,lambda,mu] = pca_new(X,d+1);
    Y=ones(n,1)*mu+(X-ones(n,1)*mu)*V*V.'; % projection of X onto the d+1 dimensional affine space mu+V
    l=zeros(n,1);

    for i=1:n
        l(i)=norm(Y(i,:))^2;
    end
    lbar=mean(l);
    H=zeros(m,m);
    f=zeros(m,1);
    for i=1:n
        H=H+(mu-Y(i,:)).'*(mu-Y(i,:));
        f=f+(l(i)-lbar)*((mu-Y(i,:)).');
    end
    c=mu.'+V*V.'*(-0.5*pinv(H)*f-mu.');  % center of the sphere
    
    Riemd=zeros(n,1); %distance between sample and center
    for i=1:n
        Riemd(i)=sqrt((c.'-Y(i,:))*(c.'-Y(i,:)).');
    end
    r=mean(Riemd); 

end
return