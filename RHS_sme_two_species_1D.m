%the right hand sides of the discrete mean equation (DMEs) to be solved using
%ode45.

function [dCdt]=RHS_sme_two_species(C,r_a,r_b,p,q,r,rho,K)


dCdt=zeros(2,K);

C=reshape(C,[],K);

for i=1:K

    if i>=3 && i<=K-2
        
        dCdt(:,i)=[r_a/2*C(1,i-1)*(1-C(1,i)-C(2,i))*(1-p)^(C(1,i-2))*(1-p)^(C(1,i))*(1-r)^(C(2,i-2))*(1-r)^(C(2,i))+...
            r_a/2*C(1,i+1)*(1-C(1,i)-C(2,i))*(1-p)^(C(1,i+2))*(1-p)^(C(1,i))*(1-r)^(C(2,i+2))*(1-r)^(C(2,i))-...
            r_a/2*C(1,i)*(1-C(1,i-1)-C(2,i-1))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))-...
            r_a/2*C(1,i)*(1-C(1,i+1)-C(2,i+1))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))+...
            (r_a+r_b)/2*rho*C(1,i-1)*C(2,i)*(1-p)^(C(1,i-2))*(1-p)^(C(1,i))*(1-r)^(C(2,i-2))*(1-r)^(C(2,i))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))+...
            (r_a+r_b)/2*rho*C(1,i+1)*C(2,i)*(1-p)^(C(1,i+2))*(1-p)^(C(1,i))*(1-r)^(C(2,i+2))*(1-r)^(C(2,i))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))-...
            (r_a+r_b)/2*rho*C(1,i)*C(2,i-1)*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))*(1-q)^(C(2,i-2))*(1-q)^(C(2,i))*(1-r)^(C(1,i-2))*(1-r)^(C(1,i))-...
            (r_a+r_b)/2*rho*C(1,i)*C(2,i+1)*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))*(1-q)^(C(2,i+2))*(1-q)^(C(2,i))*(1-r)^(C(1,i+2))*(1-r)^(C(1,i))

            
            r_b/2*C(2,i-1)*(1-C(2,i)-C(1,i))*(1-q)^(C(2,i-2))*(1-q)^(C(2,i))*(1-r)^(C(1,i-2))*(1-r)^(C(1,i))+...
            r_b/2*C(2,i+1)*(1-C(2,i)-C(1,i))*(1-q)^(C(2,i+2))*(1-q)^(C(2,i))*(1-r)^(C(1,i+2))*(1-r)^(C(1,i))-...
            r_b/2*C(2,i)*(1-C(2,i-1)-C(1,i-1))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))-...
            r_b/2*C(2,i)*(1-C(2,i+1)-C(1,i+1))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))+...
            (r_a+r_b)/2*rho*C(2,i-1)*C(1,i)*(1-q)^(C(2,i-2))*(1-q)^(C(2,i))*(1-r)^(C(1,i-2))*(1-r)^(C(1,i))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))+...
            (r_a+r_b)/2*rho*C(2,i+1)*C(1,i)*(1-q)^(C(2,i+2))*(1-q)^(C(2,i))*(1-r)^(C(1,i+2))*(1-r)^(C(1,i))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))-...
            (r_a+r_b)/2*rho*C(2,i)*C(1,i-1)*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))*(1-p)^(C(1,i-2))*(1-p)^(C(1,i))*(1-r)^(C(2,i-2))*(1-r)^(C(2,i))-...
            (r_a+r_b)/2*rho*C(2,i)*C(1,i+1)*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))*(1-p)^(C(1,i+2))*(1-p)^(C(1,i))*(1-r)^(C(2,i+2))*(1-r)^(C(2,i))...
            ];

    
    elseif i==2
        
        dCdt(:,i)=[r_a/2*C(1,i-1)*(1-C(1,i)-C(2,i))*(1-p)^(0)*(1-p)^(C(1,i))*(1-r)^(0)*(1-r)^(C(2,i))+...
            r_a/2*C(1,i+1)*(1-C(1,i)-C(2,i))*(1-p)^(C(1,i+2))*(1-p)^(C(1,i))*(1-r)^(C(2,i+2))*(1-r)^(C(2,i))-...
            r_a/2*C(1,i)*(1-C(1,i-1)-C(2,i-1))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))-...
            r_a/2*C(1,i)*(1-C(1,i+1)-C(2,i+1))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))+...
            (r_a+r_b)/2*rho*C(1,i-1)*C(2,i)*(1-p)^(0)*(1-p)^(C(1,i))*(1-r)^(0)*(1-r)^(C(2,i))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))+...
            (r_a+r_b)/2*rho*C(1,i+1)*C(2,i)*(1-p)^(C(1,i+2))*(1-p)^(C(1,i))*(1-r)^(C(2,i+2))*(1-r)^(C(2,i))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))-...
            (r_a+r_b)/2*rho*C(1,i)*C(2,i-1)*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))*(1-q)^(0)*(1-q)^(C(2,i))*(1-r)^(0)*(1-r)^(C(1,i))-...
            (r_a+r_b)/2*rho*C(1,i)*C(2,i+1)*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))*(1-q)^(C(2,i+2))*(1-q)^(C(2,i))*(1-r)^(C(1,i+2))*(1-r)^(C(1,i))

        
            r_b/2*C(2,i-1)*(1-C(2,i)-C(1,i))*(1-q)^(0)*(1-q)^(C(2,i))*(1-r)^(0)*(1-r)^(C(1,i))+...
            r_b/2*C(2,i+1)*(1-C(2,i)-C(1,i))*(1-q)^(C(2,i+2))*(1-q)^(C(2,i))*(1-r)^(C(1,i+2))*(1-r)^(C(1,i))-...
            r_b/2*C(2,i)*(1-C(2,i-1)-C(1,i-1))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))-...
            r_b/2*C(2,i)*(1-C(2,i+1)-C(1,i+1))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))+...
            (r_a+r_b)/2*rho*C(2,i-1)*C(1,i)*(1-q)^(0)*(1-q)^(C(2,i))*(1-r)^(0)*(1-r)^(C(1,i))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))+...
            (r_a+r_b)/2*rho*C(2,i+1)*C(1,i)*(1-q)^(C(2,i+2))*(1-q)^(C(2,i))*(1-r)^(C(1,i+2))*(1-r)^(C(1,i))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))-...
            (r_a+r_b)/2*rho*C(2,i)*C(1,i-1)*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))*(1-p)^(0)*(1-p)^(C(1,i))*(1-r)^(0)*(1-r)^(C(2,i))-...
            (r_a+r_b)/2*rho*C(2,i)*C(1,i+1)*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))*(1-p)^(C(1,i+2))*(1-p)^(C(1,i))*(1-r)^(C(2,i+2))*(1-r)^(C(2,i))
            ];
        


        

    elseif i==1
        
        dCdt(1,i)=r_a/2*C(1,i+1)*(1-C(1,i)-C(2,i))*(1-p)^(C(1,i+2))*(1-p)^(C(1,i))*(1-r)^(C(2,i+2))*(1-r)^(C(2,i))-...
            r_a/2*C(1,i)*(1-C(1,i+1)-C(2,i+1))*(1-p)^(0)*(1-p)^(C(1,i+1))*(1-r)^(0)*(1-r)^(C(2,i+1))+...
            (r_a+r_b)/2*rho*C(1,i+1)*C(2,i)*(1-p)^(C(1,i+2))*(1-p)^(C(1,i))*(1-r)^(C(2,i+2))*(1-r)^(C(2,i))*(1-q)^(0)*(1-q)^(C(2,i+1))*(1-r)^(0)*(1-r)^(C(1,i+1))-...
            (r_a+r_b)/2*rho*C(1,i)*C(2,i+1)*(1-p)^(0)*(1-p)^(C(1,i+1))*(1-r)^(0)*(1-r)^(C(2,i+1))*(1-q)^(C(2,i+2))*(1-q)^(C(2,i))*(1-r)^(C(1,i+2))*(1-r)^(C(1,i));
        
        dCdt(2,i)=r_b/2*C(2,i+1)*(1-C(2,i)-C(1,i))*(1-q)^(C(2,i+2))*(1-q)^(C(2,i))*(1-r)^(C(1,i+2))*(1-r)^(C(1,i))-...
            r_b/2*C(2,i)*(1-C(2,i+1)-C(1,i+1))*(1-q)^(0)*(1-q)^(C(2,i+1))*(1-r)^(0)*(1-r)^(C(1,i+1))+...
            (r_a+r_b)/2*rho*C(2,i+1)*C(1,i)*(1-q)^(C(2,i+2))*(1-q)^(C(2,i))*(1-r)^(C(1,i+2))*(1-r)^(C(1,i))*(1-p)^(0)*(1-p)^(C(1,i+1))*(1-r)^(0)*(1-r)^(C(2,i+1))-...
            (r_a+r_b)/2*rho*C(2,i)*C(1,i+1)*(1-q)^(0)*(1-q)^(C(2,i+1))*(1-r)^(0)*(1-r)^(C(1,i+1))*(1-p)^(C(1,i+2))*(1-p)^(C(1,i))*(1-r)^(C(2,i+2))*(1-r)^(C(2,i));

        
        
        

    elseif i==K-1
        
        dCdt(:,i)=[r_a/2*C(1,i-1)*(1-C(1,i)-C(2,i))*(1-p)^(C(1,i-2))*(1-p)^(C(1,i))*(1-r)^(C(2,i-2))*(1-r)^(C(2,i))+...
            r_a/2*C(1,i+1)*(1-C(1,i)-C(2,i))*(1-p)^(0)*(1-p)^(C(1,i))*(1-r)^(0)*(1-r)^(C(2,i))-...
            r_a/2*C(1,i)*(1-C(1,i-1)-C(2,i-1))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))-...
            r_a/2*C(1,i)*(1-C(1,i+1)-C(2,i+1))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))+...
            (r_a+r_b)/2*rho*C(1,i-1)*C(2,i)*(1-p)^(C(1,i-2))*(1-p)^(C(1,i))*(1-r)^(C(2,i-2))*(1-r)^(C(2,i))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))+...
            (r_a+r_b)/2*rho*C(1,i+1)*C(2,i)*(1-p)^(0)*(1-p)^(C(1,i))*(1-r)^(0)*(1-r)^(C(2,i))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))-...
            (r_a+r_b)/2*rho*C(1,i)*C(2,i-1)*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))*(1-q)^(C(2,i-2))*(1-q)^(C(2,i))*(1-r)^(C(1,i-2))*(1-r)^(C(1,i))-...
            (r_a+r_b)/2*rho*C(1,i)*C(2,i+1)*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))*(1-q)^(0)*(1-q)^(C(2,i))*(1-r)^(0)*(1-r)^(C(1,i))
            
            r_b/2*C(2,i-1)*(1-C(2,i)-C(1,i))*(1-q)^(C(2,i-2))*(1-q)^(C(2,i))*(1-r)^(C(1,i-2))*(1-r)^(C(1,i))+...
            r_b/2*C(2,i+1)*(1-C(2,i)-C(1,i))*(1-q)^(0)*(1-q)^(C(2,i))*(1-r)^(0)*(1-r)^(C(1,i))-...
            r_b/2*C(2,i)*(1-C(2,i-1)-C(1,i-1))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))-...
            r_b/2*C(2,i)*(1-C(2,i+1)-C(1,i+1))*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))+...
            (r_a+r_b)/2*rho*C(2,i-1)*C(1,i)*(1-q)^(C(2,i-2))*(1-q)^(C(2,i))*(1-r)^(C(1,i-2))*(1-r)^(C(1,i))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))+...
            (r_a+r_b)/2*rho*C(2,i+1)*C(1,i)*(1-q)^(0)*(1-q)^(C(2,i))*(1-r)^(0)*(1-r)^(C(1,i))*(1-p)^(C(1,i-1))*(1-p)^(C(1,i+1))*(1-r)^(C(2,i-1))*(1-r)^(C(2,i+1))-...
            (r_a+r_b)/2*rho*C(2,i)*C(1,i-1)*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))*(1-p)^(C(1,i-2))*(1-p)^(C(1,i))*(1-r)^(C(2,i-2))*(1-r)^(C(2,i))-...
            (r_a+r_b)/2*rho*C(2,i)*C(1,i+1)*(1-q)^(C(2,i-1))*(1-q)^(C(2,i+1))*(1-r)^(C(1,i-1))*(1-r)^(C(1,i+1))*(1-p)^(0)*(1-p)^(C(1,i))*(1-r)^(0)*(1-r)^(C(2,i))
            ];
        
        

        
    elseif i==K
        
        dCdt(:,i)=[r_a/2*C(1,i-1)*(1-C(1,i)-C(2,i))*(1-p)^(C(1,i-2))*(1-p)^(C(1,i))*(1-r)^(C(2,i-2))*(1-r)^(C(2,i))-...
            r_a/2*C(1,i)*(1-C(1,i-1)-C(2,i-1))*(1-p)^(C(1,i-1))*(1-p)^(0)*(1-r)^(C(2,i-1))*(1-r)^(0)+...
            (r_a+r_b)/2*rho*C(1,i-1)*C(2,i)*(1-p)^(C(1,i-2))*(1-p)^(C(1,i))*(1-r)^(C(2,i-2))*(1-r)^(C(2,i))*(1-q)^(C(2,i-1))*(1-q)^(0)*(1-r)^(C(1,i-1))*(1-r)^(0)-...
            (r_a+r_b)/2*rho*C(1,i)*C(2,i-1)*(1-p)^(C(1,i-1))*(1-p)^(0)*(1-r)^(C(2,i-1))*(1-r)^(0)*(1-q)^(C(2,i-2))*(1-q)^(C(2,i))*(1-r)^(C(1,i-2))*(1-r)^(C(1,i))
            
            r_b/2*C(2,i-1)*(1-C(2,i)-C(1,i))*(1-q)^(C(2,i-2))*(1-q)^(C(2,i))*(1-r)^(C(1,i-2))*(1-r)^(C(1,i))-...
            r_b/2*C(2,i)*(1-C(2,i-1)-C(1,i-1))*(1-q)^(C(2,i-1))*(1-q)^(0)*(1-r)^(C(1,i-1))*(1-r)^(0)+...
            (r_a+r_b)/2*rho*C(2,i-1)*C(1,i)*(1-q)^(C(2,i-2))*(1-q)^(C(2,i))*(1-r)^(C(1,i-2))*(1-r)^(C(1,i))*(1-p)^(C(1,i-1))*(1-p)^(0)*(1-r)^(C(2,i-1))*(1-r)^(0)-...
            (r_a+r_b)/2*rho*C(2,i)*C(1,i-1)*(1-q)^(C(2,i-1))*(1-q)^(0)*(1-r)^(C(1,i-1))*(1-r)^(0)*(1-p)^(C(1,i-2))*(1-p)^(C(1,i))*(1-r)^(C(2,i-2))*(1-r)^(C(2,i))
            ];

        

    end
end


dCdt=dCdt(:);
